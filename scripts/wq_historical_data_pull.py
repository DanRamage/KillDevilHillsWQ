import os
import sys
sys.path.append('../commonfiles/python')
from datetime import datetime, timedelta
import time
from pytz import timezone
import logging.config
import ConfigParser
import optparse
import csv
import netCDF4 as nc
import numpy as np
from bisect import bisect_left,bisect_right
from pyoos.collectors.coops.coops_sos import CoopsSos
from pyoos.collectors.ndbc.ndbc_sos import NdbcSos
import requests
from unitsConversion import uomconversionFunctions
from build_tide_file import create_tide_data_file_mp
from wqXMRGProcessing import wqXMRGProcessing
from wq_sites import wq_sample_sites

from wqDatabase import wqDB


"""
Stations:
  NWS:
    -KFFA
      Wind Speed/Dir
      Air Temp

  NDBC:
    -41062
      Water Temp
      Salinity
      WInd Speed/Dir
    -44095
      Water Temp
      Waves
  NOS/Tide:
    -8651370(primary)
      Tide
    -8652226(subordinate)
      Offset:
        High: -5 Min 1.04ft
        Low: 1 Min 1.43ft

    -8651605(subordinate from 8651370)
      Offset:
        High: -1 Min 1.01ft
        Low: 2 Min 1.43ft
"""
nos_sites = ['8651370']
nos_obs = {
  "sea_water_temperature": {
                             "units": "celsius",
                             "xenia_name": "water_temperature",
                             "xenia_units": "celsius"
                           },
  "wind_speed": {
                  "units": "m_s-1",
                  "xenia_name": "wind_speed",
                  "xenia_units": "m_s-1"

                },
  "wind_from_direction": {
    "units": "degrees_true",
    "xenia_name": "wind_from_direction",
    "xenia_units": "degrees_true"

  }
}
#ndbc_sites = ['44095']
ndbc_sites = ['41062']
"""
  {
    'sos_obs_query': 'wind_speed',
    'sites': ['41062'],
    'xenia_obs': [
      {
        "sos_obs_name": "wind_speed",
        "units": "m_s-1",
        "xenia_name": "wind_speed",
        "xenia_units": "m_s-1"
      }
    ]
  },
  {
    'sos_obs_query': 'wind_from_direction',
    'sites': ['41062'],
    'xenia_obs': [
      {
        "sos_obs_name": "wind_from_direction",
        "units": "degrees_true",
        "xenia_name": "wind_from_direction",
        "xenia_units": "degrees_true"
      }
    ]
  },

"""
ndbc_obs = [
  {
    'sos_obs_query': 'sea_water_temperature',
    'sites': ['44095'],
    'xenia_obs': [
      {
        "sos_obs_name": "sea_water_temperature",
        "units": "celsius",
        "xenia_name": "water_temperature",
        "xenia_units": "celsius"
      }
    ]
  },
  {
    'sos_obs_query': 'sea_water_salinity',
    'sites': ['41062'],
    'xenia_obs': [
      {
        "sos_obs_name": "sea_water_salinity",
        "units": "psu",
        "xenia_name": "salinity",
        "xenia_units": "psu"
      }
    ]
  },
  {
    'sos_obs_query': 'waves',
    'sites': ['44095'],
    'xenia_obs': [
      {
        "sos_obs_name": "sea_surface_wave_to_direction",
        "units": "degree",
        "xenia_name": "sea_surface_wave_to_direction",
        "xenia_units": "degree"
      }
    ]
  }
]

nws_site=['KFFA']
nws_obs = {
  "p01i": {
    "units": "inch",
    "xenia_name": "precipitation",
    "xenia_units": "millimeter"
  },
  "tmpf": {
             "units": "fahrenheit",
             "xenia_name": "air_temperature",
             "xenia_units": "celsius"
           },
  "sknt": {
                  "units": "knots",
                  "xenia_name": "wind_speed",
                  "xenia_units": "m_s-1"

                },
  "drct": {
    "units": "degrees_true",
    "xenia_name": "wind_from_direction",
    "xenia_units": "degrees_true"

  }
}
"""
"""


def get_nos_data(site, dates, units_coverter, db_obj):
  start_time = time.time()
  logger = logging.getLogger(__name__)
  logger.debug("Starting get_nos_data")

  row_entry_date = datetime.now()
  utc_tz = timezone('UTC')
  eastern_tz= timezone('US/Eastern')

  platform_handle = 'nos.%s.met' % (site)
  if db_obj.platformExists(platform_handle) == -1:
    obs_list = []
    for single_obs in nos_obs:
      obs_list.append({'obs_name': nos_obs[single_obs]['xenia_name'],
                       'uom_name': nos_obs[single_obs]['xenia_units'],
                       's_order': 1})
    db_obj.buildMinimalPlatform(platform_handle, obs_list)

  nos_query = CoopsSos()
  #dates.sort(reverse=True)
  for rec_date in dates:
    logger.debug("Query site: %s for date: %s" % (site, rec_date))
    nos_query.clear()
    utc_end_date = rec_date.astimezone(utc_tz) + timedelta(hours=24)
    start_date = rec_date.astimezone(utc_tz) - timedelta(hours=24)

    for single_obs in nos_obs:
      obs_type = nos_obs[single_obs]['xenia_name']
      uom_type = nos_obs[single_obs]['xenia_units']
      s_order = 1

      nos_query.filter(features=[site], start=start_date, end=utc_end_date, variables=[single_obs])
      try:
        #results = nos_query.collect()
        response = nos_query.raw(responseFormat="text/csv")
      except Exception as e:
        logger.exception(e)
      else:
        csv_reader = csv.reader(response.split('\n'), delimiter=',')
        line_cnt = 0
        for row in csv_reader:
          if line_cnt > 0 and len(row):
            obs_date = datetime.strptime(row[4], '%Y-%m-%dT%H:%M:%SZ')
            obs_val = float(row[5])
            logger.debug("Adding obs: %s(%s) Date: %s Value: %s S_Order: %d" %\
                         (obs_type, uom_type, obs_date, obs_val, s_order))

            if not db_obj.addMeasurement(obs_type,
                                    uom_type,
                                    platform_handle,
                                    obs_date.strftime('%Y-%m-%dT%H:%M:%S'),
                                    float(row[2]),
                                    float(row[3]),
                                    0,
                                    [obs_val],
                                    sOrder=s_order,
                                    autoCommit=True,
                                    rowEntryDate=row_entry_date ):
              logger.error(db_obj.lastErrorMsg)


          line_cnt += 1

  logger.debug("Finished get_nos_data in %f seconds" % (time.time() - start_time))

  return
def process_nos_data(**kwargs):
  out_directory = kwargs['output_directory']
  all_dates = kwargs['all_dates']
  db_obj = kwargs['db_obj']
  units_converter = kwargs['units_converter']

  if kwargs.get('query_nos', False):
    for site in nos_sites:
      get_nos_data(site, all_dates, units_converter, db_obj)
  return

def process_ndbc_data(**kwargs):
  out_directory = kwargs['output_directory']
  all_dates = kwargs['all_dates']
  db_obj = kwargs['db_obj']
  units_converter = kwargs['units_converter']

  if kwargs.get('query_ndbc', False):
    for site in ndbc_sites:
      get_ndbc_data(site, all_dates, units_converter, db_obj)
  return

def get_ndbc_data(site, dates, units_coverter, db_obj):
  start_time = time.time()
  logger = logging.getLogger(__name__)
  logger.debug("Starting get_ndbc_data")

  row_entry_date = datetime.now()
  utc_tz = timezone('UTC')
  eastern_tz= timezone('US/Eastern')

  platform_handle = 'ndbc.%s.met' % (site)
  if db_obj.platformExists(platform_handle) == -1:
    obs_list = []
    for obs_setup in ndbc_obs:
      if site in obs_setup['sites']:
        for xenia_obs in obs_setup['xenia_obs']:
          obs_list.append({'obs_name': xenia_obs['xenia_name'],
                           'uom_name': xenia_obs['xenia_units'],
                           's_order': 1})
    db_obj.buildMinimalPlatform(platform_handle, obs_list)

  sos_query = NdbcSos()
  #dates.sort(reverse=True)
  dates.sort(reverse=True)
  for rec_date in dates:
    logger.debug("Query site: %s for date: %s" % (site, rec_date))
    sos_query.clear()
    utc_end_date = rec_date.astimezone(utc_tz) + timedelta(hours=24)
    start_date = rec_date.astimezone(utc_tz) - timedelta(hours=24)

    for obs_setup in ndbc_obs:
      if site in obs_setup['sites']:
        date_ndx = None
        value_ndx = None
        lat_ndx = None
        lon_ndx = None
        depth_ndx = None

        sos_query.filter(features=[site], start=start_date, end=utc_end_date, variables=[obs_setup['sos_obs_query']])
        try:
          #results = nos_query.collect()
          response = sos_query.raw(responseFormat="text/csv")
        except Exception as e:
          logger.exception(e)
        else:
          csv_reader = csv.reader(response.split('\n'), delimiter=',')
          line_cnt = 0

          for row in csv_reader:
            for xenia_obs_setup in obs_setup['xenia_obs']:
              obs_type = xenia_obs_setup['xenia_name']
              uom_type = xenia_obs_setup['xenia_units']
              s_order = 1

              if line_cnt > 0 and len(row):
                obs_date = datetime.strptime(row[date_ndx], '%Y-%m-%dT%H:%M:%SZ')
                try:
                  obs_val = float(row[value_ndx])
                except ValueError as e:
                  logger.exception(e)
                  obs_val = 0.0
                logger.debug("Adding obs: %s(%s) Date: %s Value: %s S_Order: %d" %\
                             (obs_type, uom_type, obs_date, obs_val, s_order))
                depth = 0
                if depth_ndx is not None:
                  depth = float(row[depth_ndx])
                if not db_obj.addMeasurement(obs_type,
                                        uom_type,
                                        platform_handle,
                                        obs_date.strftime('%Y-%m-%dT%H:%M:%S'),
                                        float(row[lat_ndx]),
                                        float(row[lon_ndx]),
                                        depth,
                                        [obs_val],
                                        sOrder=s_order,
                                        autoCommit=True,
                                        rowEntryDate=row_entry_date ):
                  logger.error(db_obj.lastErrorMsg)
              else:
                if value_ndx is None:
                  for ndx,val in enumerate(row):
                    if val.lower().find(xenia_obs_setup['sos_obs_name']) != -1:
                      value_ndx = ndx
                    if val.lower().find('date_time') != -1:
                      date_ndx = ndx
                    if val.lower().find('latitude') != -1:
                      lat_ndx = ndx
                    if val.lower().find('longitude') != -1:
                      lon_ndx = ndx
                    if val.lower().find('depth') != -1:
                      depth_ndx = ndx
            line_cnt += 1

  logger.debug("Finished get_ndbc_data in %f seconds" % (time.time() - start_time))

def process_nws_data(**kwargs):
  out_directory = kwargs['output_directory']
  all_dates = kwargs['all_dates']
  db_obj = kwargs['db_obj']
  units_converter = kwargs['units_converter']
  file_list = kwargs['file_list']
  for site in nws_site:
    for file in file_list:
      if file.find(site) != -1:
        get_nws_data(site, all_dates, units_converter, db_obj, out_directory, file)
  return
def get_nws_data(site, dates, units_coverter, db_obj, out_directory, file):
  logger = logging.getLogger(__name__)
  header =[
    "station","valid","lon","lat","tmpf","dwpf","relh","drct","sknt","p01i","alti","mslp","vsby","gust","skyc1","skyc2","skyc3","skyc4","skyl1","skyl2","skyl3","skyl4","wxcodes","metar"
  ]
  platform_handle = 'nws.%s.met' % (site)
  if db_obj.platformExists(platform_handle) == -1:
    obs_list = []
    for obs_name in nws_obs:
      xenia_obs = nws_obs[obs_name]
      obs_list.append({'obs_name': xenia_obs['xenia_name'],
                       'uom_name': xenia_obs['xenia_units'],
                       's_order': 1})
    db_obj.buildMinimalPlatform(platform_handle, obs_list)

  utc_tz = timezone('UTC')
  row_entry_date = datetime.now()
  with open(file, "r") as nws_data_file:
    nws_csv_rdr = csv.DictReader(nws_data_file, header)
    for row_ndx,row in enumerate(nws_csv_rdr):
      if row_ndx != 0:
        try:
          date_time = utc_tz.localize(datetime.strptime(row['valid'], "%Y-%m-%d %H:%M"))
        except (ValueError, Exception) as e:
          logger.exception(e)
        else:
          for obs_name in nws_obs:
            try:
              xenia_obs = nws_obs[obs_name]
              obs_val = units_coverter.measurementConvert(float(row[obs_name]), xenia_obs['units'], xenia_obs['xenia_units'])
              if obs_val is not None:
                if not db_obj.addMeasurement(xenia_obs['xenia_name'],
                                             xenia_obs['xenia_units'],
                                              platform_handle,
                                              date_time.strftime('%Y-%m-%dT%H:%M:%S'),
                                              float(row['lat']),
                                              float(row['lon']),
                                              0,
                                              [obs_val],
                                              sOrder=1,
                                              autoCommit=True,
                                              rowEntryDate=row_entry_date):
                  logger.error(db_obj.lastErrorMsg)
                else:
                  logger.debug("Platform: %s Date: %s %s: %f(%s)" % (platform_handle, date_time.strftime('%Y-%m-%dT%H:%M:%S'), xenia_obs['xenia_name'], obs_val, xenia_obs['xenia_units']))
              else:
                logger.error("Unable to convert obs: %s: %f" % (obs_name, float(row[obs_name])))

            except ValueError as e:
              logger.exception(e)

tide_stations = ['8651370']
def process_tide_data(**kwargs):
  start_time = time.time()
  logger = logging.getLogger(__name__)
  logger.debug("Starting process_tide_data")

  out_directory = kwargs['output_directory']
  all_dates = kwargs['all_dates']
  #db_obj = kwargs['db_obj']
  #units_converter = kwargs['units_converter']
  log_conf_file = kwargs['log_config_file']
  #eastern_tz = timezone('US/Eastern')
  tide_dates = []
  for tide_station in tide_stations:
    for date_rec in all_dates:
      #Add 24 hours since we want to make sure we have +/- 24 hours around our date. This
      #way we can have enough data to use if we want the sample times starting at midnight
      #or we want to use the actual sample time. Instead of getting the data for each
      #time for a given sample day, just do a more coarse approach.
      #tide_date = date_rec + timedelta(hours=24)
      #tide_date = tide_date.replace(hour=0, minute=0, second=0)
      tide_dates.append(date_rec)

  tide_output_file = os.path.join(out_directory, "%s.csv" % (tide_station))
  create_tide_data_file_mp(tide_station,
                           tide_dates,
                           tide_output_file,
                           1,
                           log_conf_file,
                           False)

  logger.debug("Finished process_tide_data in %f seconds" % (time.time()-start_time))


def process_ctd_data(**kwargs):
  logger = logging.getLogger(__name__)
  out_directory = kwargs['output_directory']
  all_dates = kwargs['all_dates']
  db_obj = kwargs['db_obj']
  units_converter = kwargs['units_converter']
  file_name = os.path.join(out_directory, 'ctd_data.csv')

  if kwargs.get('query_ctd', False):
    get_ctd_data(file_name)
  platform_handle = 'usace.frf_pier.ctd'
  if db_obj.platformExists(platform_handle) == -1:
    obs_list = []
    for obs_name in nws_obs:
      xenia_obs = nws_obs[obs_name]
      obs_list.append({'obs_name': 'salinity',
                       'uom_name': 'psu',
                       's_order': 1})
    db_obj.buildMinimalPlatform(platform_handle, obs_list)

  utc_tz = timezone('UTC')
  row_entry_date = datetime.now()

  with open(file_name, "r") as ctd_data:
    header = ["date_time","Lon", "Lat", "Salinity @ 0.50 m","Salinity @ 1.00 m","Salinity @ 1.50 m","Salinity @ 2.00 m","Salinity @ 2.50 m","Salinity @ 3.00 m","Salinity @ 3.50 m","Salinity @ 4.00 m","Salinity @ 4.50 m","Salinity @ 5.00 m","Salinity @ 5.50 m","Salinity @ 6.00 m","Salinity @ 6.50 m","Salinity @ 7.00 m","Salinity @ 7.50 m","Salinity @ 8.00 m"]
    ctd_csv_reader = csv.DictReader(ctd_data, fieldnames=header)
    for row_ndx, row in enumerate(ctd_csv_reader):
      if row_ndx > 0:
        if len(row['Salinity @ 0.50 m']):
          if not db_obj.addMeasurement('salinity',
                                       'psu',
                                       platform_handle,
                                       row['date_time'],
                                       float(row['Lat']),
                                       float(row['Lon']),
                                       0.5,
                                       [float(row['Salinity @ 0.50 m'])],
                                       sOrder=1,
                                       autoCommit=True,
                                       rowEntryDate=row_entry_date):
            logger.error(db_obj.lastErrorMsg)
          else:
            logger.debug("Platform: %s Date: %s %s: %f(%s)" % (platform_handle, row['date_time'], 'salinity', float(row['Salinity @ 0.50 m']),'psu'))


def get_ctd_data(output_filename):
  full_url = "https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/ctd/eop-ctd/eop-ctd.ncml"
  logger = logging.getLogger(__name__)
  logger.debug("Connecting to thredds endpoint for ctd data: %s" % (full_url))
  #Connect to netcdf file for retrieving data from c10 buoy. To speed up retrieval, we connect
  #only once and retrieve the times.
  ncObj = nc.Dataset(full_url)
  ctd_times = ncObj.variables['time'][:]
  ctd_salinity = ncObj.variables['salinity'][:]
  ctd_depth = ncObj.variables['depth'][:]
  lon = ncObj.variables['lon'][:]
  lat = ncObj.variables['lat'][:]
  with open(output_filename, "w") as ctd_csv:
    header = ["date_time", "Lon", "Lat"]
    for depth in ctd_depth:
      header.append('\"Salinity @ %.2f m\"' % (depth))
    ctd_csv.write(",".join(header))
    ctd_csv.write("\n")
    for ndx,time in enumerate(ctd_times):
      row = []
      date_time = datetime.utcfromtimestamp(time)
      logger.debug("Processing date: %s" % (date_time))
      row.append("\"%s\"" % (date_time.strftime("%Y-%m-%d %H:%M:%S")))
      row.append(str(lon))
      row.append(str(lat))
      salinity_data = ctd_salinity[ndx]
      for salinity_val in salinity_data:
        try:
          np.testing.assert_equal(salinity_val, np.nan)
        except AssertionError:
          #Value is non NaN
          val = "%.2f" % (salinity_val)
        else:
          val = ""
        row.append(val)
      ctd_csv.write(",".join(row))
      ctd_csv.write("\n")


def get_dates(**kwargs):
  start_time = time.time()
  est_tz = timezone('US/Eastern')
  utc_tz = timezone('UTC')
  logger = logging.getLogger(__name__)
  logger.debug("Starting get_dates")
  dates = []
  sample_data_files = os.listdir(kwargs['data_directory'])
  start_date = kwargs.get('start_date', None)
  for file in sample_data_files:
    name,ext = os.path.splitext(file)
    if ext == '.csv':
      full_path = os.path.join(kwargs['data_directory'], file)
      logger.debug("Getting dates from: %s" % (file))
      header = [
        'site_id',
        'site_name',
        'date',
        'entero ssm',
        'entero gm'
      ]
      with open(full_path, "r") as data_file:
        csv_reader = csv.DictReader(data_file, header)
        row_ndx = 0
        for row in csv_reader:
          if row_ndx > 0:
            date_obj = (utc_tz.localize(datetime.strptime(row['date'], '%Y-%m-%d %H:%M:%S'))).astimezone(est_tz)
            if start_date is not None and date_obj >= start_date:
              date_in_list = [date for date in dates if date_obj == date]
              if len(date_in_list) == 0:
                logger.debug("Adding date: %s to list" % (date_obj))
                dates.append(date_obj)
              #else:
                #logger.debug("Date: %s aready in list" % (date_obj))
          row_ndx += 1
  dates.sort()
  logger.debug("Finished get_dates in %f seconds." % (time.time()-start_time))
  return dates

def main():
  parser = optparse.OptionParser()
  parser.add_option("--ConfigFile", dest="config_file", default=None,
                    help="INI Configuration file." )
  parser.add_option("--OutputDirectory", dest="out_dir", default='./')
  parser.add_option("--ProcessNWSData", dest="process_nws_data", action="store_true",
                    help="If set this will process NWS data, either by querying the server or processing the json files.")
  parser.add_option("--NWSDataFile", dest="nws_file",
                    help="Historical metar file to process.")

  parser.add_option("--ProcessNOSData", dest="process_nos_data", action="store_true",
                    help="If set this will process NOS data, either by querying the server or processing the json files.")
  parser.add_option("--QueryNOS", dest="query_nos", action="store_true",
                    help="If set this will query the NDBC web service to get the historical data.")
  parser.add_option("--ProcessNDBCData", dest="process_ndbc_data", action="store_true",
                    help="If set this will process NDBC data, either by querying the server or processing the json files.")
  parser.add_option("--QueryNDBC", dest="query_ndbc", action="store_true",
                    help="If set this will query the NOS web service to get the historical data.")
  parser.add_option("--ProcessTideData", dest="process_tide_data", action="store_true",
                    help="If set this will process Tide data")
  parser.add_option("--StartDate", dest="starting_date", default=None)
  parser.add_option("--QueryCTD", dest="query_ctd", action="store_true",
                    help="If set this will query the web service to get the historical data.")
  parser.add_option("--ProcessCTDData", dest="process_ctd_data", action="store_true",
                    help="If set this will process CTD data, either by querying the server or processing the json files.")

  (options, args) = parser.parse_args()

  if options.config_file is None:
    parser.print_help()
    sys.exit(-1)

  try:
    config_file = ConfigParser.RawConfigParser()
    config_file.read(options.config_file)
  except Exception as e:
    raise
  else:
    try:
      logConfFile = config_file.get('logging', 'config_file')
      logging.config.fileConfig(logConfFile)
      logger = logging.getLogger(__name__)
      logger.info("Log file opened.")
      boundaries_location_file = config_file.get('boundaries_settings', 'boundaries_file')
      sites_location_file = config_file.get('boundaries_settings', 'sample_sites')
      sample_data_dir = config_file.get('historic_sample_data', 'directory')
      historical_db_name = config_file.get("database", "name")
      units_file = config_file.get("units_conversion", "config_file")
      if options.starting_date is not None:
        start_date = timezone('US/Eastern').localize(datetime.strptime(options.starting_date, "%Y-%m-%d"))
    except ConfigParser.Error as e:
      if logger:
        logger.exception(e)
    else:
      units_conversion = uomconversionFunctions(units_file)
      historic_db = wqDB(historical_db_name, __name__)
      wq_sites = wq_sample_sites()
      wq_sites.load_sites(file_name=sites_location_file, boundary_file=boundaries_location_file)

    #Get all dates.
    all_dates = get_dates(data_directory=sample_data_dir, station_data=wq_sites, start_date=start_date)
    if options.process_ndbc_data:
      process_ndbc_data(output_directory=options.out_dir,
                        query_ndbc=options.query_ndbc,
                        all_dates=all_dates,
                        db_obj=historic_db,
                        units_converter=units_conversion)
    if options.process_nos_data:
      process_nos_data(output_directory=options.out_dir,
                        query_nos=options.query_nos,
                        all_dates=all_dates,
                        db_obj=historic_db,
                        units_converter=units_conversion)
    if options.process_tide_data:
      process_tide_data(output_directory=options.out_dir, all_dates=all_dates, log_config_file=logConfFile)

    if options.process_nws_data:
      process_nws_data(output_directory=options.out_dir,
                        all_dates=all_dates,
                        db_obj=historic_db,
                        units_converter=units_conversion,
                        file_list=options.nws_file.split(","))
    if options.process_ctd_data:
      process_ctd_data(output_directory=options.out_dir,
                        query_ctd=options.query_ctd,
                        all_dates=all_dates,
                        db_obj=historic_db,
                        units_converter=units_conversion)

  logger.info("Closing log.")
  #historical_wq_data = kdh_historical_wq_data()
  #wq_data = {}
  #historical_wq_data.get_ctd_data(datetime.now(),wq_data )
  return

if __name__ == "__main__":
  main()