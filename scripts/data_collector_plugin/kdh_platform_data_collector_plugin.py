import sys
sys.path.append('../../commonfiles/python')
import os
import logging.config
import time
from datetime import datetime, timedelta
from pytz import timezone
import csv
import traceback
from yapsy.IPlugin import IPlugin
from multiprocessing import Process

import data_collector_plugin as my_plugin

if sys.version_info[0] < 3:
  import ConfigParser
else:
  import configparser as ConfigParser

from pyoos.collectors.ndbc.ndbc_sos import NdbcSos
from pyoos.collectors.coops.coops_sos import CoopsSos
from erddapy import ERDDAP
import pandas as pd
import requests
from unitsConversion import uomconversionFunctions
from wq_sites import wq_sample_sites

#from wqDatabase import wqDB
from xeniaSQLiteAlchemy import xeniaAlchemy as sqliteAlchemy
from xeniaSQLiteAlchemy import multi_obs


# FROM pyoos SOS handling
# Convenience function to build record style time series representation
def flatten_element(p):
    rd = {'time':p.time}
    for m in p.members:
      if m['standard'] is not None:
        rd[m['standard']] = m['value']
      else:
        rd[m['name']] = m['value']
      unit = "%s_unit" % (m['name'])
      rd[unit] = m['unit']

    return rd

erddap_sites = ['44095']
erddap_obs = [
  {
    'erdap_obs_query': 'sea_water_temperature',
    'sites': ['44095', '44100'],
    'urls': ['https://erddap.secoora.org/erddap', 'https://erddap.secoora.org/erddap'],
    'xenia_obs': [
      {
        "erdap_obs_name": "sea_water_temperature",
        "units": "celsius",
        "xenia_name": "water_temperature",
        "xenia_units": "celsius",
        "dataset_ids": ['192-oregon-inlet-nc-44095', 'edu_ucsd_cdip_44100']
      }
    ]
  },
  {
    'erdap_obs_query': 'waves',
    'sites': ['44095', '44100'],
    'urls': ['https://erddap.secoora.org/erddap', 'https://erddap.secoora.org/erddap'],
    'xenia_obs': [
      {
        "erdap_obs_name": "sea_surface_wave_to_direction",
        "units": "degree",
        "xenia_name": "sea_surface_wave_to_direction",
        "xenia_units": "degree",
        "dataset_ids": ['192-oregon-inlet-nc-44095', 'edu_ucsd_cdip_44100']
      }
    ]
  }
]
ndbc_sites=['44095', '44100']
ndbc_obs = [
  {
    'sos_obs_query': 'sea_water_temperature',
    'sites': ['44095','44100'],
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
    'sos_obs_query': 'waves',
    'sites': ['44095','44100'],
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
nos_sites = ['8651370']
nos_obs = [
  {
    'sos_obs_query': 'sea_water_temperature',
    'sites': ['8651370'],
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
    'sos_obs_query': 'wind_speed',
    'sites': ['8651370'],
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
    'sites': ['8651370'],
    'xenia_obs': [
      {
        "sos_obs_name": "wind_from_direction",
        "units": "degrees_true",
        "xenia_name": "wind_from_direction",
        "xenia_units": "degrees_true"
      }
    ]
  }
]

nws_site=['KFFA']
nws_obs = [
  {
    'sites': ['KFFA'],
    'xenia_obs': [
      {
        "obs_name": "precip_in",
        "units": "inch",
        "xenia_name": "precipitation",
        "xenia_units": "millimeter"
      }
    ]
  },
  {
    'sites': ['KFFA'],
    'xenia_obs': [
      {
        "obs_name": "temp_c",
        "units": "celsius",
        "xenia_name": "air_temperature",
        "xenia_units": "celsius"
      }
    ]
  },
  {
    'sites': ['KFFA'],
    'xenia_obs': [
      {
        "obs_name": "wind_speed_kt",
        "units": "knots",
        "xenia_name": "wind_speed",
        "xenia_units": "m_s-1"
      }
    ]
  },
  {
    'sites': ['KFFA'],
    'xenia_obs': [
      {
        "obs_name": "wind_dir_degrees",
        "units": "degrees_true",
        "xenia_name": "wind_from_direction",
        "xenia_units": "degrees_true"
      }
    ]
  }

]

class kdh_platform_data_collector_plugin(my_plugin.data_collector_plugin):
  def __init__(self):
    Process.__init__(self)
    IPlugin.__init__(self)

  def initialize_plugin(self, **kwargs):
    #data_collector_plugin.initialize_plugin(self, **kwargs)
    try:
      logger = logging.getLogger(self.__class__.__name__)
      self.plugin_details = kwargs['details']
      self.begin_date = kwargs['begin_date']
      self.temp_directory = self.plugin_details.get("Settings", "temp_directory")
      self.log_config = self.plugin_details.get("Settings", "log_config")
      return True
    except Exception as e:
      logger.exception(e)
    return False

  def run(self):
    start_time = time.time()
    logger = None
    try:
      xenia_db = None
      #self.logging_client_cfg['disable_existing_loggers'] = True
      #logging.config.dictConfig(self.logging_client_cfg)
      logging.config.fileConfig(self.log_config)
      logger = logging.getLogger(self.__class__.__name__)

      logger.setLevel(logging.DEBUG)
      logger.debug("run started.")

      config_file = ConfigParser.RawConfigParser()
      config_file.read(self.plugin_details.get('Settings', 'ini_file'))

      xenia_obs_db_name = config_file.get('database', 'name')

      #xenia_db = wqDB(xenia_obs_db_name, self.__class__.__name__)
      xenia_db = sqliteAlchemy()
      xenia_db.connectDB(databaseType='sqlite',
                        dbHost=xenia_obs_db_name,
                         dbUser=None,
                         dbPwd=None,
                         dbName=None,
                         printSQL=False)

      units_file = config_file.get("units_conversion", "config_file")
      units_conversion = uomconversionFunctions(units_file)

      boundaries_location_file = config_file.get('boundaries_settings', 'boundaries_file')
      sites_location_file = config_file.get('boundaries_settings', 'sample_sites')
      wq_sites = wq_sample_sites()
      wq_sites.load_sites(file_name=sites_location_file, boundary_file=boundaries_location_file)

      for site in erddap_sites:
        try:
          self.get_erdap_data(site, erddap_obs, self.begin_date, units_conversion, xenia_db)
        except Exception as e:
          logger.exception(e)
      '''
      for site in ndbc_sites:
        try:
          self.get_ndbc_data(site, ndbc_obs, self.begin_date, units_conversion, xenia_db)
        except Exception as e:
          logger.exception(e)
      '''
      for site in nws_site:
        try:
          self.get_nws_data(site, nws_obs, self.begin_date, units_conversion, xenia_db)
        except Exception as e:
          logger.exception(e)
      for site in nos_sites:
        try:
          self.get_nos_data(site, nos_obs, self.begin_date, units_conversion, xenia_db)
        except Exception as e:
          logger.exception(e)

    except Exception as e:
      if logger is not None:
        logger.exception(e)
      else:
        traceback.print_exc(e)
    finally:
      if xenia_db is not None:
        xenia_db.disconnect()
      if logger is not None:
        logger.debug("run finished in %f seconds" % (time.time()-start_time))

    return

  def get_erdap_data(self, site, observations, begin_date, units_coverter, db_obj):
    start_time = time.time()
    logger = logging.getLogger(self.__class__.__name__)
    logger.debug("Starting get_erdap_data")

    row_entry_date = datetime.now()
    utc_tz = timezone('UTC')
    eastern_tz = timezone('US/Eastern')

    utc_end_date = begin_date.astimezone(utc_tz)
    start_date = begin_date.astimezone(utc_tz) - timedelta(hours=24)

    platform_handle = 'ndbc.%s.met' % (site)
    if db_obj.platformExists(platform_handle) is None:
      obs_list = []
      for obs_setup in observations:
        if site in obs_setup['sites']:
          for xenia_obs in obs_setup['xenia_obs']:
            obs_list.append({'obs_name': xenia_obs['xenia_name'],
                             'uom_name': xenia_obs['xenia_units'],
                             's_order': 1})
      db_obj.buildMinimalPlatform(platform_handle, obs_list)
    #Build sensor_id and m_type_id list.
    for obs_setup in observations:
      if site in obs_setup['sites']:
        ndx = obs_setup['sites'].index(site)
        url = obs_setup['urls'][ndx]
        for xenia_obs in obs_setup['xenia_obs']:
          m_type_id = db_obj.mTypeExists(xenia_obs['xenia_name'],xenia_obs['xenia_units'])
          sensor_id = db_obj.sensorExists(xenia_obs['xenia_name'],
                                          xenia_obs['xenia_units'],
                                          platform_handle,
                                          1)
          dataset_id = xenia_obs['dataset_ids'][ndx]
          erddap_query = ERDDAP(
            server= url,
            protocol="tabledap",
            response="csv")
          erddap_query.dataset_id = dataset_id
          erddap_query.variables = [obs_setup['erdap_obs_query'], 'time', 'longitude', 'latitude', 'z']
          erddap_query.constraints = {
            "time>=": start_date.strftime('%Y-%m-%dT%H:%M:%SZ'),
            "time<": utc_end_date.strftime('%Y-%m-%dT%H:%M:%SZ')
          }
          df = erddap_query.to_pandas(parse_dates=True).dropna()
          for index, row in df.iterrows():
            date_time = datetime.strptime(row[1], '%Y-%m-%dT%H:%M:%SZ')
            obs_rec = multi_obs(row_entry_date=row_entry_date,
                                platform_handle=platform_handle,
                                sensor_id=sensor_id,
                                m_type_id=m_type_id,
                                m_date=date_time.strftime('%Y-%m-%dT%H:%M:%S'),
                                m_lon=row[2],
                                m_lat=row[3],
                                m_z=row[4],
                                m_value=row[0],
                                )

            rec_id = db_obj.addRec(obs_rec, True)
            if rec_id is not None:
              logger.debug("Adding obs: %s(%s) Date: %s Value: %s S_Order: %d" % \
                           (obs_setup['erdap_obs_query'], xenia_obs['xenia_units'], date_time, row[0], row[4]))
            else:
              logger.error("Failed adding obs: %s(%s) Date: %s Value: %s S_Order: %d" % \
                           (obs_setup['erdap_obs_query'], xenia_obs['xenia_units'], date_time, row[0], row[4]))


    logger.debug("Finished get_erddap_data in %f seconds" % (time.time() - start_time))

    return
  def get_ndbc_data(self, site, observations, begin_date, units_coverter, db_obj):
    start_time = time.time()
    logger = logging.getLogger(self.__class__.__name__)
    logger.debug("Starting get_ndbc_data")

    row_entry_date = datetime.now()
    utc_tz = timezone('UTC')
    eastern_tz = timezone('US/Eastern')

    platform_handle = 'ndbc.%s.met' % (site)
    #if db_obj.platformExists(platform_handle) == -1:
    if db_obj.platformExists(platform_handle) is None:
      obs_list = []
      for obs_setup in observations:
        if site in obs_setup['sites']:
          for xenia_obs in obs_setup['xenia_obs']:
            obs_list.append({'obs_name': xenia_obs['xenia_name'],
                             'uom_name': xenia_obs['xenia_units'],
                             's_order': 1})
      db_obj.buildMinimalPlatform(platform_handle, obs_list)
    #Build sensor_id and m_type_id list.
    for obs_setup in observations:
      if site in obs_setup['sites']:
        for xenia_obs in obs_setup['xenia_obs']:
          m_type_id = db_obj.mTypeExists(xenia_obs['xenia_name'],xenia_obs['xenia_units'])
          sensor_id = db_obj.sensorExists(xenia_obs['xenia_name'],
                                          xenia_obs['xenia_units'],
                                          platform_handle,
                                          1)
          xenia_obs['m_type_id'] = m_type_id
          xenia_obs['sensor_id'] = sensor_id
    sos_query = NdbcSos()
    # dates.sort(reverse=True)

    logger.debug("Query site: %s for date: %s" % (site, begin_date))
    sos_query.clear()
    #utc_end_date = begin_date.astimezone(utc_tz) + timedelta(hours=24)
    utc_end_date = begin_date.astimezone(utc_tz)
    start_date = begin_date.astimezone(utc_tz) - timedelta(hours=24)

    for obs_setup in ndbc_obs:
      if site in obs_setup['sites']:
        date_ndx = None
        value_ndx = None
        lat_ndx = None
        lon_ndx = None
        depth_ndx = None

        sos_query.filter(features=[site], start=start_date, end=utc_end_date, variables=[obs_setup['sos_obs_query']])
        try:
          # results = nos_query.collect()
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
                try:
                  if depth_ndx is not None:
                    depth = float(row[depth_ndx])
                except ValueError as e:
                  logger.exception(e)
                  depth = 0
                try:
                  if lat_ndx is not None:
                    latitude = float(row[lat_ndx])
                  if lon_ndx is not None:
                    longitude = float(row[lon_ndx])
                except ValueError as e:
                  logger.exception(e)
                  latitude=0.0
                  longitude=0.0

                obs_rec = multi_obs(row_entry_date=row_entry_date,
                                    platform_handle=platform_handle,
                                    sensor_id=xenia_obs_setup['sensor_id'],
                                    m_type_id=xenia_obs_setup['m_type_id'],
                                    m_date=obs_date.strftime('%Y-%m-%dT%H:%M:%S'),
                                    m_lon=longitude,
                                    m_lat=latitude,
                                    m_z=depth,
                                    m_value=obs_val,
                                    )

                rec_id = db_obj.addRec(obs_rec, True)
                if rec_id is not None:
                  logger.debug("Adding obs: %s(%s) Date: %s Value: %s S_Order: %d" % \
                               (obs_type, uom_type, obs_date, obs_val, s_order))
                else:
                  logger.error("Failed adding obs: %s(%s) Date: %s Value: %s S_Order: %d" % \
                               (obs_type, uom_type, obs_date, obs_val, s_order))
              else:
                if value_ndx is None:
                  for ndx, val in enumerate(row):
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

  def get_nos_data(self, site, observations, begin_date, units_coverter, db_obj):
    start_time = time.time()
    logger = logging.getLogger(self.__class__.__name__)
    logger.debug("Starting get_nos_data")

    row_entry_date = datetime.now()
    utc_tz = timezone('UTC')
    eastern_tz= timezone('US/Eastern')

    platform_handle = 'nos.%s.met' % (site)
    if db_obj.platformExists(platform_handle) is None:
      obs_list = []
      for obs_setup in observations:
        if site in obs_setup['sites']:
          for xenia_obs in obs_setup['xenia_obs']:
            obs_list.append({'obs_name': xenia_obs['xenia_name'],
                             'uom_name': xenia_obs['xenia_units'],
                             's_order': 1})
      db_obj.buildMinimalPlatform(platform_handle, obs_list)
    #Build sensor_id and m_type_id list.
    for obs_setup in observations:
      if site in obs_setup['sites']:
        for xenia_obs in obs_setup['xenia_obs']:
          m_type_id = db_obj.mTypeExists(xenia_obs['xenia_name'],xenia_obs['xenia_units'])
          sensor_id = db_obj.sensorExists(xenia_obs['xenia_name'],
                                          xenia_obs['xenia_units'],
                                          platform_handle,
                                          1)
          xenia_obs['m_type_id'] = m_type_id
          xenia_obs['sensor_id'] = sensor_id

    sos_query = CoopsSos()
    #dates.sort(reverse=True)

    logger.debug("Query site: %s for date: %s" % (site, begin_date))
    sos_query.clear()
    #utc_end_date = begin_date.astimezone(utc_tz) + timedelta(hours=24)
    utc_end_date = begin_date.astimezone(utc_tz)
    start_date = begin_date.astimezone(utc_tz) - timedelta(hours=24)

    for obs_setup in observations:
      if site in obs_setup['sites']:
        date_ndx = None
        value_ndx = None
        lat_ndx = None
        lon_ndx = None
        depth_ndx = None

        sos_query.filter(features=[site], start=start_date, end=utc_end_date, variables=[obs_setup['sos_obs_query']])
        try:
          # results = nos_query.collect()
          response = sos_query.raw(responseFormat="text/csv")
        except Exception as e:
          logger.exception(e)
        else:
          csv_reader = csv.reader(response.decode('utf-8').split('\n'), delimiter=',')
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
                try:
                  if depth_ndx is not None:
                    depth = float(row[depth_ndx])
                  else:
                    depth = 0
                except ValueError as e:
                  logger.exception(e)
                  depth = 0
                try:
                  if lat_ndx is not None:
                    latitude = float(row[lat_ndx])
                  if lon_ndx is not None:
                    longitude = float(row[lon_ndx])
                except ValueError as e:
                  logger.exception(e)
                  latitude=0.0
                  longitude=0.0

                obs_rec = multi_obs(row_entry_date=row_entry_date,
                                    platform_handle=platform_handle,
                                    sensor_id=xenia_obs_setup['sensor_id'],
                                    m_type_id=xenia_obs_setup['m_type_id'],
                                    m_date=obs_date.strftime('%Y-%m-%dT%H:%M:%S'),
                                    m_lon=longitude,
                                    m_lat=latitude,
                                    m_z=depth,
                                    m_value=obs_val,
                                    )

                rec_id = db_obj.addRec(obs_rec, True)
                if rec_id is not None:
                  logger.debug("Adding obs: %s(%s) Date: %s Value: %s S_Order: %d" % \
                               (obs_type, uom_type, obs_date, obs_val, s_order))
                else:
                  logger.error("Failed adding obs: %s(%s) Date: %s Value: %s S_Order: %d" % \
                               (obs_type, uom_type, obs_date, obs_val, s_order))
              else:
                if value_ndx is None:
                  for ndx, val in enumerate(row):
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


    """
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



          line_cnt += 1
    """
    logger.debug("Finished get_nos_data in %f seconds" % (time.time() - start_time))

    return

  def get_nws_data(self, site, observations, begin_date, units_coverter, db_obj):
    start_time = time.time()
    logger = logging.getLogger(self.__class__.__name__)
    logger.debug("Starting get_nws_data")

    row_entry_date = datetime.now()
    utc_tz = timezone('UTC')
    eastern_tz = timezone('US/Eastern')

    header = ["raw_text","station_id","observation_time","latitude","longitude","temp_c","dewpoint_c","wind_dir_degrees","wind_speed_kt","wind_gust_kt","visibility_statute_mi","altim_in_hg","sea_level_pressure_mb","corrected","auto","auto_station","maintenance_indicator_on","no_signal","lightning_sensor_off","freezing_rain_sensor_off","present_weather_sensor_off","wx_string","sky_cover","cloud_base_ft_agl","sky_cover","cloud_base_ft_agl","sky_cover","cloud_base_ft_agl","sky_cover","cloud_base_ft_agl","flight_category","three_hr_pressure_tendency_mb","maxT_c","minT_c","maxT24hr_c","minT24hr_c","precip_in","pcp3hr_in","pcp6hr_in","pcp24hr_in","snow_in","vert_vis_ft","metar_type","elevation_m"]

    platform_handle = 'nws.%s.met' % (site)
    if db_obj.platformExists(platform_handle) is None:
      obs_list = []
      for obs_setup in observations:
        if site in obs_setup['sites']:
          for xenia_obs in obs_setup['xenia_obs']:
            obs_list.append({'obs_name': xenia_obs['xenia_name'],
                             'uom_name': xenia_obs['xenia_units'],
                             's_order': 1})
      db_obj.buildMinimalPlatform(platform_handle, obs_list)
    # Build sensor_id and m_type_id list.
    for obs_setup in observations:
      if site in obs_setup['sites']:
        for xenia_obs in obs_setup['xenia_obs']:
          m_type_id = db_obj.mTypeExists(xenia_obs['xenia_name'], xenia_obs['xenia_units'])
          sensor_id = db_obj.sensorExists(xenia_obs['xenia_name'],
                                          xenia_obs['xenia_units'],
                                          platform_handle,
                                          1)
          xenia_obs['m_type_id'] = m_type_id
          xenia_obs['sensor_id'] = sensor_id

    #awc_query = AwcRest()
    # dates.sort(reverse=True)

    logger.debug("Query site: %s for date: %s" % (site, begin_date))
    #awc_query.clear()
    # utc_end_date = begin_date.astimezone(utc_tz) + timedelta(hours=24)
    utc_end_date = begin_date.astimezone(utc_tz)
    start_date = begin_date.astimezone(utc_tz) - timedelta(hours=24)

    url = "https://aviationweather.gov/cgi-bin/data/metar.php?ids={station}&hours={hours}&format=csv"\
      .format(station=site,hours=48)
    '''
    url = "https://www.aviationweather.gov/adds/dataserver_current/httpparam?" \
          "dataSource=metars&requestType=retrieve&format=csv&stationString={station}&hoursBeforeNow={hours}".format(
      station=site,
      hours=48
    )
    '''
    '''
    #awc_query.filter(bbox=(-77,34,-74,36),stationString=site, start=start_date, end=utc_end_date)
    awc_query.filter(bbox=(0,0,0,0), start=start_date, end=utc_end_date)
    try:
      # response = awc_query.collect()
      temp_file = os.path.join(self.temp_directory, "nws_data.csv")
      with open(temp_file, "w") as nws_file_obj:
        response = awc_query.raw(format="csv", stationString=site)
        nws_file_obj.write(response[0])
    except Exception as e:
      logger.exception(e)
    '''

    try:
      temp_file = os.path.join(self.temp_directory, "nws_data.csv")
      with open(temp_file, "w") as nws_file_obj:
        req = requests.get(url)
        nws_file_obj.write(req.text)
    except Exception as e:
      logger.exception(e)
    else:
      with open(temp_file, "rU") as nws_data_file:
        line_cnt = 0
        header_row = 5
        nws_csv_rdr = csv.DictReader(nws_data_file, header)
        for row_ndx, row in enumerate(nws_csv_rdr):
          if row_ndx > header_row:
            try:
              obs_date = utc_tz.localize(datetime.strptime(row['observation_time'], "%Y-%m-%dT%H:%M:%SZ"))
            except (ValueError, Exception) as e:
              logger.exception(e)
            else:
              for obs in nws_obs:
                try:
                  for xenia_obs in obs["xenia_obs"]:
                    src_obs_name = xenia_obs['obs_name']
                    src_obs_uom = xenia_obs['units']
                    xenia_obs_name = xenia_obs['xenia_name']
                    xenia_uom = xenia_obs['xenia_units']
                    sensor_id = xenia_obs['sensor_id']
                    m_type_id = xenia_obs['m_type_id']
                    try:
                      obs_value = float(row[src_obs_name])
                    except ValueError as e:
                      obs_value = 0.0
                    try:
                      latitude = float(row["latitude"])
                      longitude = float(row["longitude"])
                    except ValueError as e:
                      latitude = 0.0
                      longitude = 0.0
                    try:
                      elevation = float(row['elevation_m'])
                    except ValueError as e:
                      elevation = 0

                    obs_val = units_coverter.measurementConvert(obs_value,
                                                                src_obs_uom,
                                                                xenia_uom)
                    if obs_val is not None:
                      obs_rec = multi_obs(row_entry_date=row_entry_date,
                                          platform_handle=platform_handle,
                                          sensor_id=sensor_id,
                                          m_type_id=m_type_id,
                                          m_date=obs_date.strftime('%Y-%m-%dT%H:%M:%S'),
                                          m_lon=longitude,
                                          m_lat=latitude,
                                          m_z=elevation,
                                          m_value=obs_val
                                          )

                      rec_id = db_obj.addRec(obs_rec, True)
                      if rec_id is not None:
                        logger.debug("Adding obs: %s(%s) Date: %s Value: %s S_Order: 1" % \
                                     (xenia_obs_name, xenia_uom, obs_date, obs_val))
                      else:
                        logger.error("Failed adding obs: %s(%s) Date: %s Value: %s S_Order: 1" % \
                                     (xenia_obs_name, xenia_uom, obs_date, obs_val))

                except (Exception,ValueError) as e:
                  logger.exception(e)

    logger.debug("Finished get_nws_data in %f seconds" % (time.time() - start_time))

    return

