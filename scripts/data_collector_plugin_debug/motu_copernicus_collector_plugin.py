import sys
sys.path.append('../../commonfiles/python')
import logging.config
from data_collector_plugin import data_collector_plugin
from datetime import datetime, timedelta
from pytz import timezone
import ConfigParser
import traceback
import time
from yapsy.IPlugin import IPlugin
from multiprocessing import Process

from motuclient import motu_api

class MotoParameters:
  def __init__(self, **kwargs):
    self.__dict__.update(kwargs)
class motu_copernicus_collector_plugin(data_collector_plugin):
  def __init__(self):
    Process.__init__(self)
    IPlugin.__init__(self)

  def initialize_plugin(self, **kwargs):
    data_collector_plugin.initialize_plugin(self, **kwargs)
    try:
      logger = logging.getLogger(self.__class__.__name__)
      self.plugin_details = kwargs['details']
      self.ini_file = self.plugin_details.get('Settings', 'ini_file')
      self.begin_date = kwargs['begin_date']
      self.output_directory = self.plugin_details.get('Settings', 'temp_directory')
      self.log_config = self.plugin_details.get("Settings", "log_config")
      return True
    except Exception as e:
      logger.exception(e)
    return False

  def run(self):
    try:
      start_time = time.time()
      #self.logging_client_cfg['disable_existing_loggers'] = True
      #logging.config.dictConfig(self.logging_client_cfg)
      logging.config.fileConfig(self.log_config)
      logger = logging.getLogger(self.__class__.__name__)
      logger.debug("run started.")
      utc_tz = timezone('UTC')
      #Give ourselves  a day ahead since we can't specify time periods, just dates.
      utc_end_date = self.begin_date.astimezone(utc_tz) + timedelta(hours=24)
      start_date = self.begin_date.astimezone(utc_tz) - timedelta(hours=216)

      motu_parms = MotoParameters(
        auth_mode=motu_api.AUTHENTICATION_MODE_CAS,
        user= 'dramage',
        pwd= 'DanCMEMS2018',
        motu= 'http://nrtcmems.mercator-ocean.fr/motu-web/Motu',
        service_id= 'GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS',
        product_id= 'global-analysis-forecast-phy-001-024',
        longitude_min=-80.0,
        longitude_max=-70.0,
        latitude_min=30.0,
        latitude_max=37.0,
        date_min= start_date.strftime('%Y-%m-%d'),
        date_max= utc_end_date.strftime('%Y-%m-%d'),
        depth_min= 0.493,
        depth_max= 0.4942,
        variable= ['so'],
        out_dir= self.output_directory,
        out_name= 'motu_data',
        block_size=65536,
        log_level=logging.DEBUG,
        outputWritten=None,
        proxy_server=None,
        proxy_user=None,
        proxy_pwd=None,
        user_agent=None,
        describe=None,
        size=None,
        sync=None,
        console_mode=None,
        socket_timeout=None
      )
      motu_log = logging.getLogger("motu_api")
      handlers = logger.handlers
      for hndlr in handlers:
        motu_log.addHandler(hndlr)
      motu_log.setLevel(logging.DEBUG)
      motu_api.execute_request(motu_parms)

    except (ConfigParser.Error, Exception) as e:
      if logger is not None:
        logger.exception(e)
      else:
        traceback.print_exc(e)
    finally:
      logger.debug("run finished in %f seconds" % (time.time()-start_time))
    return
