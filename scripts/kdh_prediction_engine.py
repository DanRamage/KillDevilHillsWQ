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
import sys
if sys.version_info[0] < 3:
  import ConfigParser
else:
  import configparser as ConfigParser
from collections import OrderedDict

from yapsy.PluginManager import PluginManager
from multiprocessing import Queue

from wq_prediction_tests import predictionTest,predictionLevels,wqEquations

from output_plugin import output_plugin
from data_collector_plugin import data_collector_plugin

from wq_prediction_engine import wq_prediction_engine
from wq_sites import wq_sample_sites
from data_result_types import data_result_types
from kdh_wq_data import kdh_wq_data
from stats import stats
from enterococcus_wq_test import EnterococcusPredictionTest,EnterococcusPredictionTestEx


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

  def build_test_objects(self, **kwargs):
    config_file = kwargs['config_file']
    site_name = kwargs['site_name']

    model_list = []
    #Get the sites test configuration ini, then build the test objects.
    try:
      test_config_file = config_file.get(site_name, 'prediction_config')
      entero_lo_limit = config_file.getint('entero_limits', 'limit_lo')
      entero_hi_limit = config_file.getint('entero_limits', 'limit_hi')
    except ConfigParser.Error as e:
        self.logger.exception(e)
    else:
      self.logger.debug("Site: %s Model Config File: %s" % (site_name, test_config_file))

      model_config_file = ConfigParser.RawConfigParser()
      model_config_file.read(test_config_file)
      #Get the number of prediction models we use for the site.
      model_count = model_config_file.getint("settings", "model_count")
      self.logger.debug("Site: %s Model count: %d" % (site_name, model_count))

      for cnt in range(model_count):
        model_name = model_config_file.get("model_%d" % (cnt+1), "name")
        model_equation = model_config_file.get("model_%d" % (cnt+1), "formula")
        self.logger.debug("Site: %s Model name: %s equation: %s" % (site_name, model_name, model_equation))

        test_obj = EnterococcusPredictionTestEx(formula=model_equation,
                                                site_name=site_name,
                                                model_name=model_name)
        test_obj.set_category_limits(entero_lo_limit, entero_hi_limit)
        model_list.append(test_obj)

    return model_list

  def run_wq_models(self, **kwargs):
    prediction_testrun_date = datetime.now()
    try:

      begin_date = kwargs['begin_date']
      config_file = ConfigParser.RawConfigParser()
      config_file.read(kwargs['config_file_name'])

      boundaries_location_file = config_file.get('boundaries_settings', 'boundaries_file')
      sites_location_file = config_file.get('boundaries_settings', 'sample_sites')
      wq_sites = wq_sample_sites()
      wq_sites.load_sites(file_name=sites_location_file, boundary_file=boundaries_location_file)

      enable_output_plugins = config_file.getboolean('output_plugins', 'enable_plugins')
      output_plugin_dirs=config_file.get('output_plugins', 'plugin_directories').split(',')

      enable_data_collector_plugins = config_file.getboolean('data_collector_plugins', 'enable_plugins')
      data_collector_plugin_directories = config_file.get('data_collector_plugins', 'plugin_directories').split(',')

    except (ConfigParser.Error, Exception) as e:
      self.logger.exception(e)
    else:
      try:

        #Run any data collector plugins we have.
        if enable_data_collector_plugins:
          self.collect_data(data_collector_plugin_directories=data_collector_plugin_directories,
                            begin_date=begin_date)

        site_data = OrderedDict()
        total_time = 0
        tide_offsets = []
        reset_site_specific_data_only = True
        wq_data = kdh_wq_data(config_file=kwargs['config_file_name'])

        site_model_ensemble = []
        tide_offsets = []

        for site in wq_sites:
          try:
            # Get all the models used for the particular sample site.
            model_list = self.build_test_objects(config_file=config_file, site_name=site.name)
            if len(model_list):
              # Create the container for all the models.
              site_equations = wqEquations(site.name, model_list, True)
            else:
              self.logger.error("No models found for site: %s" % (site.name))
          except (ConfigParser.Error, Exception) as e:
            self.logger.exception(e)
          else:
            try:
              #Get site specific settings.
              hycom_data_prefix = config_file.get(site.name, 'hycom_prefix')
              copernicus_data_prefix = config_file.get(site.name, 'copernicus_prefix')
              rutgers_data_prefix = config_file.get(site.name, 'rutgers_prefix')
              rutgers_cell_point = tuple(
                float(pt) for pt in config_file.get(site.name, 'rutgers_cell_loc').split(','))
              # Get the platforms the site will use
              platforms = config_file.get(site.name, 'platforms').split(',')
              platform_nfo = []

              for platform in platforms:
                obs_uoms = config_file.get(platform, 'observation').split(';')
                obs_uom_nfo = []
                for nfo in obs_uoms:
                  obs, uom = nfo.split(',')
                  obs_uom_nfo.append({'observation': obs,
                                      'uom': uom})
                platform_nfo.append({'platform_handle': config_file.get(platform, 'handle'),
                                     'observations': obs_uom_nfo})

              offset_tide_station = config_file.get(site.name, 'offset_tide_station')
              offset_param = "%s_tide_data" % (offset_tide_station)
              # We use the virtual tide sites as there no stations near the sites.
              tide_station = config_file.get(site.name, 'tide_station')
              tide_station_settings = {
                'tide_station': tide_station,
                'offset_tide_station': config_file.get(offset_param, 'station_id'),
                'hi_tide_time_offset': config_file.getint(offset_param, 'hi_tide_time_offset'),
                'lo_tide_time_offset': config_file.getint(offset_param, 'lo_tide_time_offset'),
                'hi_tide_height_offset': config_file.getfloat(offset_param, 'hi_tide_height_offset'),
                'lo_tide_height_offset': config_file.getfloat(offset_param, 'lo_tide_height_offset')
              }
              tide_offsets.append(tide_station_settings)

            except ConfigParser.Error as e:
              self.logger.exception(e)
            else:
              wq_data.reset(site=site,
                             tide_station=tide_station,
                             tide_offset_params=tide_offsets,
                             hycom_prefix=hycom_data_prefix,
                             copernicus_prefix=copernicus_data_prefix,
                             start_date=begin_date,
                             rutgers_cell_point=rutgers_cell_point,
                             rutgers_prefix=rutgers_data_prefix,
                             platform_info=platform_nfo)

              site_data['station_name'] = site.name
              wq_data.query_data(begin_date,
                                 begin_date,
                                 site_data)

              reset_site_specific_data_only = True
              site_equations.runTests(site_data)
              total_test_time = sum(testObj.test_time for testObj in site_equations.tests)
              self.logger.debug("Site: %s total time to execute models: %f ms" % (site.name, total_test_time * 1000))
              total_time += total_test_time

              #Calculate some statistics on the entero results. This is making an assumption
              #that all the tests we are running are calculating the same value, the entero
              #amount.
              entero_stats = None
              if len(site_equations.tests):
                entero_stats = stats()
                for test in site_equations.tests:
                  if test.mlrResult is not None:
                    entero_stats.addValue(test.mlrResult)
                entero_stats.doCalculations()

              site_model_ensemble.append({'metadata': site,
                                          'models': site_equations,
                                          'statistics': entero_stats,
                                          'entero_value': None})

        if enable_output_plugins:
          self.output_results(output_plugin_directories=output_plugin_dirs,
                              site_model_ensemble=site_model_ensemble,
                              prediction_date=kwargs['begin_date'],
                              prediction_run_date=prediction_testrun_date)
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
      yapsy_logger = logging.getLogger('yapsy')
      yapsy_logger.setLevel(logging.DEBUG)
      #yapsy_logger.parent.level = logging.DEBUG
      yapsy_logger.disabled = False

      simplePluginManager.setPluginPlaces(kwargs['data_collector_plugin_directories'])

      simplePluginManager.collectPlugins()

      output_queue = Queue()
      plugin_cnt = 0
      plugin_start_time = time.time()
      for plugin in simplePluginManager.getAllPlugins():
        plugin_start_time = time.time()
        self.logger.info("Starting plugin: %s" % (plugin.name))
        if plugin.plugin_object.initialize_plugin(details=plugin.details,
                                                  queue=output_queue,
                                                  begin_date=kwargs['begin_date']):
          plugin.plugin_object.start()
          self.logger.info("Waiting for %s plugin to complete." % (plugin.name))
          plugin.plugin_object.join()
          self.logger.info("%s plugin to completed in %f seconds." % (plugin.name, time.time()-plugin_start_time))
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
  """
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
        plugin.plugin_object.emit(sampling_date=kwargs['prediction_date'],
                                  failed_sites=kwargs['failed_sites'],
                                  feedback_email=kwargs['feedback_email'])
        plugin_cnt += 1
      else:
        self.logger.error("Failed to initialize plugin: %s" % (plugin.details))
    self.logger.debug("%d output plugins run in %f seconds" % (plugin_cnt, time.time() - plugin_start_time))
    self.logger.info("Finished collect_data")

    return
  """
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

  except ConfigParser.Error as e:
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
      except Exception as e:
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
    except Exception as e:
      logger.exception(e)

  if logger:
    logger.info("Log file closed.")

  return

if __name__ == "__main__":
  main()
