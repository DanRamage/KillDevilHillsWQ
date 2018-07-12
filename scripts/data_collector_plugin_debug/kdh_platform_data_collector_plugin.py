import sys
sys.path.append('../../commonfiles/python')
import os
import logging.config
import requests
import time
import xlrd
from datetime import datetime, timedelta
from pytz import timezone
import csv


from wq_sites import wq_sample_sites
from data_collector_plugin import data_collector_plugin
import ConfigParser
from pyoos.collectors.ndbc.ndbc_sos import NdbcSos
from pyoos.collectors.coops.coops_sos import CoopsSos

import requests
from unitsConversion import uomconversionFunctions
from build_tide_file import create_tide_data_file_mp
from wqXMRGProcessing import wqXMRGProcessing
from wq_sites import wq_sample_sites

#from wqDatabase import wqDB
from xeniaSQLiteAlchemy import xeniaAlchemy as sqliteAlchemy
from xeniaSQLiteAlchemy import multi_obs

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



class kdh_platform_data_collector_plugin(data_collector_plugin):

  def initialize_plugin(self, **kwargs):
    data_collector_plugin.initialize_plugin(self, **kwargs)
    try:
      logger = logging.getLogger(self.__class__.__name__)
      self.plugin_details = kwargs['details']
      self.begin_date = kwargs['begin_date']
      return True
    except Exception as e:
      logger.exception(e)
    return False

  def run(self):
    start_time = time.time()
    try:

      self.logging_client_cfg['disable_existing_loggers'] = True
      logging.config.dictConfig(self.logging_client_cfg)
      logger = logging.getLogger(self.__class__.__name__)
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

      for site in ndbc_sites:
        self.get_ndbc_data(site, ndbc_obs, self.begin_date, units_conversion, xenia_db)

      logger.debug("run finished in %f seconds" % (time.time()-start_time))
    except ConfigParser.Error, e:
      print("No log configuration file given, logging disabled.")
    except Exception,e:
      import traceback
      traceback.print_exc(e)
      sys.exit(-1)
    finally:
      xenia_db.disconnect()
    return

  def get_ndbc_data(self, site, observations, begin_date, units_coverter, db_obj):
    start_time = time.time()
    logger = logging.getLogger(__name__)
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
                                    m_date=obs_date,
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

  def get_nos_data(self, site, dates, units_coverter, db_obj):
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
              """
              row_id           = Column(Integer,primary_key=True)                     
              row_entry_date   = Column(String)    
              row_update_date  = Column(String)    
              platform_handle  = Column(String(100))    
              sensor_id        = Column(Integer,ForeignKey(sensor.row_id))                     
              m_type_id        = Column(Integer,ForeignKey(m_type.row_id))                     
              m_date           = Column(String) 
              m_lon            = Column(Float)            
              m_lat            = Column(Float)             
              m_z              = Column(Float)             
              m_value          = Column(Float)             
              m_value_2        = Column(Float)             
              m_value_3        = Column(Float)             
              m_value_4        = Column(Float)             
              m_value_5        = Column(Float)             
              m_value_6        = Column(Float)             
              m_value_7        = Column(Float)             
              m_value_8        = Column(Float)             
              qc_metadata_id   = Column(Integer)                     
              qc_level         = Column(Integer)                     
              qc_flag          = Column(String(100))         
              qc_metadata_id_2 = Column(Integer)                     
              qc_level_2       = Column(Integer)                     
              qc_flag_2        = Column(String(100))         
              metadata_id      = Column(Integer)                     
              d_label_theta    = Column(Integer)                     
              d_top_of_hour    = Column(Integer)                     
              d_report_hour    = Column(String) 
              """
              """
              obs_rec = multi_obs(row_entry_date=row_entry_date,
                                platform_handle=platform_handle,
                                sensor_id=
                                m_type_id
                                m_date
                                m_lon
                                m_lat
                                m_z
                                m_value
                                m_value_2
                                m_value_3
                                m_value_4
                                m_value_5
                                m_value_6
                                m_value_7
                                m_value_8
                                qc_metadata_id
                                qc_level
                                qc_flag
                                qc_metadata_id_2
                                qc_level_2
                                qc_flag_2
                                metadata_id
                                d_label_theta
                                d_top_of_hour
                                d_report_hour

              )
              """
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
