import sys
sys.path.append('../../commonfiles/python')
import data_collector_plugin as my_plugin

import logging.config
from datetime import datetime
from pytz import timezone
if sys.version_info[0] < 3:
  import ConfigParser
else:
  import configparser as ConfigParser
import traceback
import time
from yapsy.IPlugin import IPlugin

#import multiprocessing
#multiprocessing.set_start_method('fork')
from multiprocessing import Process


from wqXMRGProcessing import wqXMRGProcessingGP



class nexrad_collector_plugin(my_plugin.data_collector_plugin):

  def initialize_plugin(self, **kwargs):
    try:
      Process.__init__(self)
      IPlugin.__init__(self)

      logger = logging.getLogger(self.__class__.__name__)
      self.plugin_details = kwargs['details']
      self.ini_file = self.plugin_details.get('Settings', 'ini_file')
      self.log_config = self.plugin_details.get("Settings", "log_config")
      self.xmrg_workers_logfile = self.plugin_details.get("Settings", "xmrg_log_file")
      self._logger_name = 'nexrad_mp_logging'


      return True
    except Exception as e:
      logger.exception(e)
    return False

  def run(self):
    logger = None
    try:
      start_time = time.time()
      logging.config.fileConfig(self.log_config)
      logger = logging.getLogger()
      '''
      mp_logging = MainLogConfig(log_filename=self.xmrg_workers_logfile,
                                 logname=self._logger_name,
                                 level=logging.DEBUG,
                                 disable_existing_loggers=True)
      mp_logging.setup_logging()
      '''
      #logger = mp_logging.getLogger()
      logger.debug("run started.")

      config_file = ConfigParser.RawConfigParser()
      config_file.read(self.ini_file)
      backfill_hours = config_file.getint('nexrad_database', 'backfill_hours')
      fill_gaps = config_file.getboolean('nexrad_database', 'fill_gaps')
      logger.debug("Backfill hours: %d Fill Gaps: %s" % (backfill_hours, fill_gaps))

    except (ConfigParser.Error, Exception) as e:
      traceback.print_exc(e)
      if logger is not None:
        logger.exception(e)
    else:
      try:
        logging_config = {
            'version': 1,
            'disable_existing_loggers': False,
            'formatters': {
                'f': {
                    'format': "%(asctime)s,%(levelname)s,%(funcName)s,%(lineno)d,%(message)s",
                    'datefmt': '%Y-%m-%d %H:%M:%S'
                }
            },
            'handlers': {
                'stream': {
                    'class': 'logging.StreamHandler',
                    'formatter': 'f',
                    'level': logging.DEBUG
                },
                'file_handler': {
                    'class': 'logging.handlers.RotatingFileHandler',
                    'filename': self.xmrg_workers_logfile,
                    'formatter': 'f',
                    'level': logging.DEBUG
                }
            },
            'root': {
                'handlers': ['file_handler', 'stream'],
                'level': logging.NOTSET,
                'propagate': False
            }
        }

        geopandas_xmrg = wqXMRGProcessingGP(logger=True, logger_name=self._logger_name, logger_config=logging_config)
        geopandas_xmrg.load_config_settings(config_file = self.ini_file)

        start_date_time = timezone('US/Eastern').localize(datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)).astimezone(timezone('UTC'))
        #start_date_time = timezone('US/Eastern').localize(datetime.strptime('2021-07-08', "%Y-%m-%d").replace(hour=0, minute=0, second=0, microsecond=0)).astimezone(timezone('UTC'))

        if fill_gaps:
          logger.info("Fill gaps Start time: %s Prev Hours: %d" % (start_date_time, backfill_hours))
          geopandas_xmrg.fill_gaps(start_date_time, backfill_hours)
        else:
          logger.info("Backfill N Hours Start time: %s Prev Hours: %d" % (start_date_time, backfill_hours))
          file_list = geopandas_xmrg.download_range(start_date_time, backfill_hours)
          geopandas_xmrg.import_files(file_list)

      except Exception as e:
        logger.exception(e)
      logger.debug("run finished in %f seconds" % (time.time()-start_time))

      #mp_logging.shutdown_logging()
    return

  def finalize(self):
    return
