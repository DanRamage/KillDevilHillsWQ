import sys
sys.path.append('../commonfiles/python')
sys.path.append('./data_collector_plugins_debug')
import os

import logging.config
from datetime import datetime, timedelta
from pytz import timezone
import traceback
import time
import optparse
import ConfigParser

from yapsy.PluginManager import PluginManager
from multiprocessing import Queue

from wq_prediction_tests import predictionTest,predictionLevels

from output_plugin import output_plugin
from data_collector_plugin import data_collector_plugin

from wq_prediction_engine import wq_prediction_engine
from wq_sites import wq_sample_sites
#from wq_sample_data_collector_plugin import wq_sample_data_collector_plugin
from data_result_types import data_result_types


class bacteria_sample_test(predictionTest):
  def __init__(self, name):
    self.predictionLevel = predictionLevels(predictionLevels.NO_TEST)
    self.name = name
    self.test_time = None
    self.enabled = True
    self.sample_value = None

  """
  @property
  def predictionLevel(self):
    return self.predictionLevel
  @property
  def name(self):
    return self.name
  @property
  def test_time(self):
    return self.test_time
  @property
  def sample_value(self):
    return self.sample_value
  """
  def set_category_limits(self, low_limit, high_limit):
    self.low_limit = low_limit
    self.high_limit = high_limit

  def runTest(self, data, test_time):
    self.sample_value = data
    self.test_time = test_time
    self.categorize_result()
  """
  Function: categorize_result
  Purpose: For the bacteria sample value, this catergorizes the value.
  Parameters:
    None
  Return:
    A predictionLevels value.
  """
  def categorize_result(self):
    if self.enabled:
      self.predictionLevel.value = predictionLevels.NO_TEST
      if self.sample_value is not None:
        if self.sample_value < self.low_limit:
          self.predictionLevel.value = predictionLevels.LOW
        elif self.sample_value >= self.high_limit:
          self.predictionLevel.value = predictionLevels.HIGH
    else:
      self.predictionLevel.value = predictionLevels.DISABLED


class kdh_prediction_engine(wq_prediction_engine):
  def __init__(self):
    #wq_prediction_engine.__init__(self)
    self.logger = logging.getLogger(type(self).__name__)

    self.bacteria_sample_data = None
  '''
  Function: build_test_objects
  Purpose: Builds the models used for doing the predictions.
  Parameters:
    config_file - ConfigParser object
    site_name - The name of the site whose models we are building.
    use_logging - Flag to specify if we are to use logging.
  Return:
    A list of models constructed.
  '''

  def build_test_objects(self, config_file, wq_sites):
    logger = logging.getLogger(__name__)

    test_list = []
    # Get the sites test configuration ini, then build the test objects.
    try:
      entero_lo_limit = config_file.getint('limits', 'limit_lo')
      entero_hi_limit = config_file.getint('limits', 'limit_hi')
    except ConfigParser.Error, e:
      if logger:
        logger.exception(e)
    else:
      for site in wq_sites:
        if site.name in self.bacteria_sample_data:
          results = self.bacteria_sample_data[site.name]
          for result in results:
            test_obj = bacteria_sample_test(site.name)
            test_obj.set_category_limits(entero_lo_limit, entero_hi_limit)
            test_obj.runTest(result.value, result.date_time)
            test_list.append(test_obj)

    return test_list

  def run_wq_models(self, **kwargs):
    today_date = datetime.now()
    try:
      config_file = ConfigParser.RawConfigParser()
      config_file.read(kwargs['config_file_name'])

      boundaries_location_file = config_file.get('boundaries_settings', 'boundaries_file')
      sites_location_file = config_file.get('boundaries_settings', 'sample_sites')
      wq_sites = wq_sample_sites()
      wq_sites.load_sites(file_name=sites_location_file, boundary_file=boundaries_location_file)

      enable_data_collector_plugins = config_file.getboolean('data_collector_plugins', 'enable_plugins')
      data_collector_plugin_directories=config_file.get('data_collector_plugins', 'plugin_directories').split(',')
      if enable_data_collector_plugins:
        self.collect_data(data_collector_plugin_directories=data_collector_plugin_directories)

      test_list = self.build_test_objects(config_file, wq_sites)

      enable_output_plugins = config_file.getboolean('output_plugins', 'enable_plugins')
      output_plugin_dirs=config_file.get('output_plugins', 'plugin_directories').split(',')

      #email_settings_ini = config_file.get('email_settings', 'settings_ini')
      #email_ini_cfg = ConfigParser.RawConfigParser()
      #email_ini_cfg.read(email_settings_ini)

    except (ConfigParser.Error, Exception) as e:
      self.logger.exception(e)
    else:
      try:
        if len(test_list):
          sample_date = None
          failed_sites = []
          sampling_date = today_date - timedelta(hours=24)
          for test in test_list:
            if sample_date is None:
              sample_date = test.test_time
            if test.test_time.date() == sampling_date.date():
              if test.predictionLevel.value == predictionLevels.HIGH:
                for site in wq_sites:
                  if site.name == test.name:
                    failed_sites.append({
                      'test_result':  test,
                      'wq_site': site
                    })
                    break
          if enable_output_plugins:
            self.output_results(output_plugin_directories=output_plugin_dirs,
                                failed_sites = failed_sites,
                                feedback_email=feedback_email,
                                sample_date=sample_date)
        else:
          self.logger.debug("No sites/data found to create test objects.")
      except Exception as e:
        self.logger.exception(e)

  def collect_data(self, **kwargs):
    self.logger.info("Begin collect_data")
    try:
      simplePluginManager = PluginManager()
      logging.getLogger('yapsy').setLevel(logging.DEBUG)
      simplePluginManager.setCategoriesFilter({
         "DataCollector": data_collector_plugin
         })

      # Tell it the default place(s) where to find plugins
      self.logger.debug("Plugin directories: %s" % (kwargs['data_collector_plugin_directories']))
      logging.getLogger('yapsy').setLevel(logging.DEBUG)
      simplePluginManager.setPluginPlaces(kwargs['data_collector_plugin_directories'])

      simplePluginManager.collectPlugins()

      output_queue = Queue()
      plugin_cnt = 0
      plugin_start_time = time.time()
      for plugin in simplePluginManager.getAllPlugins():
        self.logger.info("Starting plugin: %s" % (plugin.name))
        if plugin.plugin_object.initialize_plugin(details=plugin.details,
                                                  queue=output_queue):
          plugin.plugin_object.start()
        else:
          self.logger.error("Failed to initialize plugin: %s" % (plugin.name))
        plugin_cnt += 1

      #Wait for the plugings to finish up.
      self.logger.info("Waiting for %d plugins to complete." % (plugin_cnt))
      for plugin in simplePluginManager.getAllPlugins():
        plugin.plugin_object.join()

      while not output_queue.empty():
        results = output_queue.get()
        if results[0] == data_result_types.SAMPLING_DATA_TYPE:
          self.bacteria_sample_data = results[1]

      self.logger.info("%d Plugins completed in %f seconds" % (plugin_cnt, time.time() - plugin_start_time))
    except Exception as e:
      self.logger.exception(e)

  def output_results(self, **kwargs):
    self.logger.info("Begin run_output_plugins")

    simplePluginManager = PluginManager()
    logging.getLogger('yapsy').setLevel(logging.DEBUG)
    simplePluginManager.setCategoriesFilter({
       "OutputResults": output_plugin
       })

    # Tell it the default place(s) where to find plugins
    self.logger.debug("Plugin directories: %s" % (kwargs['output_plugin_directories']))
    simplePluginManager.setPluginPlaces(kwargs['output_plugin_directories'])

    simplePluginManager.collectPlugins()

    plugin_cnt = 0
    plugin_start_time = time.time()
    for plugin in simplePluginManager.getAllPlugins():
      self.logger.info("Starting plugin: %s" % (plugin.name))
      if plugin.plugin_object.initialize_plugin(details=plugin.details):
        plugin.plugin_object.emit(sampling_date=kwargs['sample_date'],
                                  failed_sites=kwargs['failed_sites'],
                                  feedback_email=kwargs['feedback_email'])
        plugin_cnt += 1
      else:
        self.logger.error("Failed to initialize plugin: %s" % (plugin.details))
    self.logger.debug("%d output plugins run in %f seconds" % (plugin_cnt, time.time() - plugin_start_time))
    self.logger.info("Finished collect_data")

    return

def main():
  parser = optparse.OptionParser()
  parser.add_option("-c", "--ConfigFile", dest="config_file",
                    help="INI Configuration file." )
  parser.add_option("-s", "--StartDateTime", dest="start_date_time",
                    help="A date to re-run the predictions for, if not provided, the default is the current day. Format is YYYY-MM-DD HH:MM:SS." )

  (options, args) = parser.parse_args()

  if(options.config_file is None):
    parser.print_help()
    sys.exit(-1)

  try:
    config_file = ConfigParser.RawConfigParser()
    config_file.read(options.config_file)

    logger = None
    use_logging = False
    logConfFile = config_file.get('logging', 'prediction_engine')
    if logConfFile:
      logging.config.fileConfig(logConfFile)
      logger = logging.getLogger(__name__)
      logger.info("Log file opened.")
      use_logging = True

  except ConfigParser.Error, e:
    traceback.print_exc(e)
    sys.exit(-1)
  else:
    dates_to_process = []
    if options.start_date_time is not None:
      #Can be multiple dates, so let's split on ','
      collection_date_list = options.start_date_time.split(',')
      #We are going to process the previous day, so we get the current date, set the time to midnight, then convert
      #to UTC.
      eastern = timezone('US/Eastern')
      try:
        for collection_date in collection_date_list:
          est = eastern.localize(datetime.strptime(collection_date, "%Y-%m-%dT%H:%M:%S"))
          #Convert to UTC
          begin_date = est.astimezone(timezone('UTC'))
          dates_to_process.append(begin_date)
      except Exception,e:
        if logger:
          logger.exception(e)
    else:
      #We are going to process the previous day, so we get the current date, set the time to midnight, then convert
      #to UTC.
      est = datetime.now(timezone('US/Eastern'))
      est = est.replace(hour=0, minute=0, second=0,microsecond=0)
      #Convert to UTC
      begin_date = est.astimezone(timezone('UTC'))
      dates_to_process.append(begin_date)

    try:
      for process_date in dates_to_process:
        pred_engine = kdh_prediction_engine()
        pred_engine.run_wq_models(begin_date=process_date,
                        config_file_name=options.config_file)
    except Exception, e:
      logger.exception(e)

  if logger:
    logger.info("Log file closed.")

  return

if __name__ == "__main__":
  main()
