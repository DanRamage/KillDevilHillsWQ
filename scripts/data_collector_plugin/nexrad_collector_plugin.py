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
from multiprocessing import Process

from wqXMRGProcessing import wqXMRGProcessing

class nexrad_collector_plugin(my_plugin.data_collector_plugin):

  def initialize_plugin(self, **kwargs):
    try:
      self._start_plugin_time = time.time()

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
                'handlers': ['file_handler'],
                'level': logging.NOTSET,
                'propagate': False
            }
        }

        xmrg_proc = wqXMRGProcessing(logger=True, logger_name=self._logger_name, logger_config=logging_config)
        xmrg_proc.load_config_settings(config_file = self.ini_file)

        start_date_time = timezone('US/Eastern').localize(datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)).astimezone(timezone('UTC'))
        if fill_gaps:
          logger.info("Fill gaps Start time: %s Prev Hours: %d" % (start_date_time, backfill_hours))
          xmrg_proc.fill_gaps(start_date_time, backfill_hours)
        else:
          logger.info("Backfill N Hours Start time: %s Prev Hours: %d" % (start_date_time, backfill_hours))
          file_list = xmrg_proc.download_range(start_date_time, backfill_hours)
          xmrg_proc.import_files(file_list)

      except Exception as e:
        logger.exception(e)
      logger.debug("run finished in %f seconds" % (time.time()-start_time))

    return

  def finalize(self):
    logger = logging.getLogger(self.__class__.__name__)
    logger.debug("finalize called, total time: %f" % (time.time() - self._start_plugin_time))

    return
