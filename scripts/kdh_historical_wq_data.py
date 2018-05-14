import sys
sys.path.append('../commonfiles/python')


import os
import logging.config
import optparse
import ConfigParser
import csv
from datetime import datetime, timedelta
import time
from pytz import timezone
from shapely.geometry import Polygon
import netCDF4 as nc
import numpy as np
from bisect import bisect_left,bisect_right
import csv
from collections import OrderedDict
import math

from NOAATideData import noaaTideData
from wqHistoricalData import tide_data_file
from wqHistoricalData import tide_data_file_ex,station_geometry,sampling_sites, wq_defines, geometry_list
from wq_output_results import wq_sample_data,wq_samples_collection,wq_advisories_file,wq_station_advisories_file
from wq_sites import wq_sample_sites

from wqHistoricalData import wq_data
from romsTools import closestCellFromPtInPolygon,bbox2ij,closestCellFromPt,closestLonLatFromPt
from nc_sample_data import nc_wq_sample_data
from xeniaSQLiteAlchemy import xeniaAlchemy as sl_xeniaAlchemy, multi_obs as sl_multi_obs, func as sl_func
from sqlalchemy import or_
from xenia import qaqcTestFlags
from stats import calcAvgSpeedAndDir

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

from wqDatabase import wqDB

from date_time_utils import *


def find_le(a, x):
  'Find rightmost ndx less than or equal to x'
  i = bisect_right(a, x)
  if i:
    return i - 1
  raise ValueError


def find_ge(a, x):
  'Find leftmost ndx greater than or equal to x'
  i = bisect_left(a, x)
  if i != len(a):
    return i
  raise ValueError


class kdh_historical_wq_data(wq_data):
  def __init__(self, **kwargs):
    wq_data.__init__(self, **kwargs)
    try:
      config_file = ConfigParser.RawConfigParser()
      config_file.read(kwargs['config_file'])
      xenia_database_name = config_file.get('database', 'name')
      hycom_model_bbox = config_file.get('hycom_model_data', 'bbox').split(',')
      model_bbox = [float(hycom_model_bbox[0]), float(hycom_model_bbox[2]),
                    float(hycom_model_bbox[1]), float(hycom_model_bbox[3])]

      poly_string = config_file.get('hycom_model_data', 'within_polygon').split(',')
      hycom_within_poly = [(float(lon_lat.split(' ')[0]), float(lon_lat.split(' ')[1])) for lon_lat in poly_string]
      hycom_endpoint_dates = config_file.get('hycom_model_data', 'endpoint_dates').split(';')
      #hycom_start_end = [(hycom.split(',')[0],hycom.split(',')[1]) for hycom in hycom_endpoint_dates]
      self.hycom_start_end = []
      self.hycomm_endpoints = []
      self.hycom_current_date_ndx = -1
      self.hycom_180_longitude = True
      utc_tz = timezone('UTC')
      for hycom_date in hycom_endpoint_dates:
        if len(hycom_date):
          parts = hycom_date.split(',')
          start = utc_tz.localize(datetime.strptime(parts[0], "%Y-%m"))
          end = ''
          if len(parts[1]) > 0:
            end = utc_tz.localize(datetime.strptime(parts[1], "%Y-%m"))
          else:
            end = None
          self.hycom_start_end.append((start, end))
          key = "%s_%s_hycom" % (parts[0], parts[1])
          self.hycomm_endpoints.append(config_file.get(key, 'thredds_url'))

      rutgers_endpoint_dates = config_file.get('rutgers_roms_model_data', 'endpoint_dates').split(';')
      self.rutgers_start_end = []
      self.rutgers_endpoints = []
      self.rutgers_current_date_ndx = -1
      utc_tz = timezone('UTC')
      for rutgers_date in rutgers_endpoint_dates:
        if len(rutgers_date):
          parts = rutgers_date.split(',')
          start = utc_tz.localize(datetime.strptime(parts[0], "%Y"))
          end = None
          if len(parts[1]) > 0:
            end = utc_tz.localize(datetime.strptime(parts[1], "%Y"))
          self.rutgers_start_end.append((start, end))
          key = "%s_%s_rutgers" % (parts[0], parts[1])
          self.rutgers_endpoints.append(config_file.get(key, 'thredds_url'))

      copernicus_file = config_file.get('copernicus_model_data', 'datafile')

      tide_file_name = config_file.get('tide_station', 'tide_file')

    except (ConfigParser.Error, Exception) as e:
      self.logger.exception(e)
      raise
    else:
      self.site = None
      #The main station we retrieve the values from.
      self.tide_station =  None
      #These are the settings to correct the tide for the subordinate station.
      self.tide_offset_settings = None
      self.tide_data_obj = None

      #Course bounding box to clip data from.
      self.model_bbox = hycom_model_bbox
      self.model_within_polygon = Polygon(hycom_within_poly)

      self.logger.debug("Connecting to copernicus data: %s" % (copernicus_file))
      self.copernicus_model = nc.Dataset(copernicus_file)
      self.copernicus_model_time = self.copernicus_model.variables['time'][:]
      #Determine the bounding box indexes.
      lons = self.copernicus_model.variables['longitude'][:]
      lats = self.copernicus_model.variables['latitude'][:]
      # latitude lower and upper index
      self.copernicus_latli = np.argmin(np.abs(lats - model_bbox[2]))
      self.copernicus_latui = np.argmin(np.abs(lats - model_bbox[3]))
      # longitude lower and upper index
      self.copernicus_lonli = np.argmin(np.abs(lons - model_bbox[0]))
      self.copernicus_lonui = np.argmin(np.abs(lons - model_bbox[1]))
      self.copernicus_lon_array = self.copernicus_model.variables['longitude'][self.copernicus_lonli:self.copernicus_lonui]
      self.copernicus_lat_array = self.copernicus_model.variables['latitude'][self.copernicus_latli:self.copernicus_latui]
      self.copernicus_data_prefix = None

      #self.connect_to_hycom(self.hycomm_endpoints[self.hycom_current_date_ndx])
      self.hycom_data_prefix = None

      self.rutgers_data_prefix = None
      #We define the cell lat/lon a site uses since the rutgers model also models in the
      #back side of the island and those cells can be closer than the ocean side so we
      #make sure we get the ocean cell.
      self.rutgers_cell_pt = None

      try:
        self.tide_data_obj = tide_data_file_ex()
        self.tide_data_obj.open(tide_file_name)
      except (IOError, Exception) as e:
        self.logger.exception(e)
        raise
      if self.logger:
        self.logger.debug("Connection to xenia db: %s" % (xenia_database_name))
      self.nexrad_db = wqDB(xenia_database_name, type(self).__name__)
      try:
        #Connect to the xenia database we use for observations aggregation.
        self.xenia_obs_db = sl_xeniaAlchemy()
        if self.xenia_obs_db.connectDB('sqlite', None, None, xenia_database_name, None, False):
          self.logger.info("Succesfully connect to DB: %s" %(xenia_database_name))
        else:
          self.logger.error("Unable to connect to DB: %s." %(xenia_database_name))


      except Exception,e:
        self.logger.exception(e)
        raise

  def connect_to_hycom(self, endpoint):
    self.logger.debug(
      "Connecting to thredds endpoint for hycom data: %s" % (endpoint))
    self.hycom_model = nc.Dataset(self.hycomm_endpoints[self.hycom_current_date_ndx])

    self.hycom_model_time = self.hycom_model.variables['time'][:]

    model_bbox = [float(self.model_bbox[0]), float(self.model_bbox[2]),
                  float(self.model_bbox[1]), float(self.model_bbox[3])]

    # Determine the bounding box indexes.
    lons = self.hycom_model.variables['lon'][:]
    lats = self.hycom_model.variables['lat'][:]
    # latitude lower and upper index
    self.hycom_latli = np.argmin(np.abs(lats - model_bbox[2]))
    self.hycom_latui = np.argmin(np.abs(lats - model_bbox[3]))
    # longitude lower and upper index
    #Check how the longitude is laid out.
    #Longitudes are -180 to 180
    if lons[0] == -180.0:
      self.hycom_lonli = np.argmin(np.abs(lons - model_bbox[0]))
      self.hycom_lonui = np.argmin(np.abs(lons - model_bbox[1]))
      self.hycom_180_longitude = True
    #Longitudes are 0-360.
    else:
      lon_li = model_bbox[0] % 360
      lon_ui = model_bbox[1] % 360
      self.hycom_lonli = np.argmin(np.abs(lons - lon_li))
      self.hycom_lonui = np.argmin(np.abs(lons - lon_ui))
      self.hycom_180_longitude = False


    self.hycom_lon_array = self.hycom_model.variables['lon'][self.hycom_lonli:self.hycom_lonui]
    self.hycom_lat_array = self.hycom_model.variables['lat'][self.hycom_latli:self.hycom_latui]
    if not self.hycom_180_longitude:
      #Convert our longitudes into -180 to 180.
      self.hycom_lon_array = np.array([(((lon+180) % 360)-180) for lon in self.hycom_lon_array])
    return

  def connect_to_rutgers(self, endpoint):
    self.logger.debug(
      "Connecting to thredds endpoint for rutgers data: %s" % (endpoint))
    self.rutgers_model = nc.Dataset(self.rutgers_endpoints[self.rutgers_current_date_ndx])
    #Seconds since 2006-01-01 00:00:00
    self.rutgers_ocean_time = self.rutgers_model.variables['ocean_time'][:]

    #In the new(2013-present) endpoint, the time parameter is time and not ocean_time.
    self.rutgers_time = None
    if 'time'in self.rutgers_model.variables:
      self.rutgers_time = self.rutgers_model.variables['time'][:]

    model_bbox = [float(self.model_bbox[0]), float(self.model_bbox[2]),
                  float(self.model_bbox[1]), float(self.model_bbox[3])]

    # Determine the bounding box indexes.
    lons = self.rutgers_model.variables['lon_rho'][:]
    lats = self.rutgers_model.variables['lat_rho'][:]

    self.rutgers_lonli, self.rutgers_lonui, self.rutgers_latli, self.rutgers_latui = bbox2ij(lons, lats, model_bbox)

    self.rutgers_lon_array = self.rutgers_model.variables['lon_rho'][self.rutgers_latli:self.rutgers_latui, self.rutgers_lonli:self.rutgers_lonui]
    self.rutgers_lat_array = self.rutgers_model.variables['lat_rho'][self.rutgers_latli:self.rutgers_latui, self.rutgers_lonli:self.rutgers_lonui]

    return

  def reset(self, **kwargs):
    if self.site is None or self.site != kwargs['site']:
      self.site = kwargs['site']
      #The main station we retrieve the values from.
      self.tide_station = kwargs['tide_station']
      #These are the settings to correct the tide for the subordinate station.
      self.tide_offset_settings = kwargs['tide_offset_params']

      #self.tide_data_obj = None
      #if 'tide_data_obj' in kwargs and kwargs['tide_data_obj'] is not None:
      #  self.tide_data_obj = kwargs['tide_data_obj']

      self.platforms_info = kwargs['platform_info']

      self.hycom_data_prefix = kwargs['hycom_prefix']

      self.copernicus_data_prefix = kwargs['copernicus_prefix']

      self.rutgers_data_prefix = kwargs['rutgers_prefix']
      self.rutgers_cell_pt = kwargs['rutgers_cell_point']

    start_date = kwargs['start_date']
    #Check the date and reconnect endpoints as needed

    for ndx, hycom_date in enumerate(self.hycom_start_end):
      if hycom_date[1] != None:
        if start_date > hycom_date[0] and start_date <= hycom_date[1]:
          hycom_endpoint_ndx = ndx
          break
      else:
        if start_date > hycom_date[0]:
          hycom_endpoint_ndx = ndx
          break
    if self.hycom_current_date_ndx != hycom_endpoint_ndx:
      self.hycom_current_date_ndx = hycom_endpoint_ndx
      self.connect_to_hycom(self.hycomm_endpoints[self.hycom_current_date_ndx])

    #Check the date and reconnect endpoints as needed
    in_between = False
    for ndx, rutgers_date in enumerate(self.rutgers_start_end):
      if start_date > rutgers_date[0]:
        if rutgers_date[1] is not None:
          if start_date <= rutgers_date[1]:
            in_between = True
        else:
          in_between = True
        if in_between:
          rutgers_endpoint_ndx = ndx
          break
    if in_between and self.rutgers_current_date_ndx != ndx:
      self.rutgers_current_date_ndx = ndx
      self.connect_to_rutgers(self.rutgers_endpoints[self.rutgers_current_date_ndx])

  """
  Function: initialize_return_data
  Purpose: INitialize our ordered dict with the data variables and assign a NO_DATA
    initial value.
  Parameters:
    wq_tests_data - An OrderedDict that is initialized.
  Return:
    None
  """
  def initialize_return_data(self, wq_tests_data):
    if self.logger:
      self.logger.debug("Creating and initializing data dict.")
    #Build variables for the base tide station.
    var_name = 'tide_range_%s' % (self.tide_station)
    wq_tests_data[var_name] = wq_defines.NO_DATA
    var_name = 'tide_hi_%s' % (self.tide_station)
    wq_tests_data[var_name] = wq_defines.NO_DATA
    var_name = 'tide_lo_%s' % (self.tide_station)
    wq_tests_data[var_name] = wq_defines.NO_DATA
    var_name = 'tide_stage_%s' % (self.tide_station)
    wq_tests_data[var_name] = wq_defines.NO_DATA

    #Build variables for the subordinate tide station.
    var_name = 'tide_range_%s' % (self.tide_offset_settings['tide_station'])
    wq_tests_data[var_name] = wq_defines.NO_DATA
    var_name = 'tide_hi_%s' % (self.tide_offset_settings['tide_station'])
    wq_tests_data[var_name] = wq_defines.NO_DATA
    var_name = 'tide_lo_%s' % (self.tide_offset_settings['tide_station'])
    wq_tests_data[var_name] = wq_defines.NO_DATA

    for platform_nfo in self.platforms_info:
      handle = platform_nfo['platform_handle']
      for obs_nfo in platform_nfo['observations']:
        var_name = '%s_avg_%s' % (handle.replace('.', '_'), obs_nfo['observation'])
        wq_tests_data[var_name] = wq_defines.NO_DATA

    for boundary in self.site.contained_by:
      if len(boundary.name):
        for prev_hours in range(24, 192, 24):
          clean_var_boundary_name = boundary.name.lower().replace(' ', '_')
          var_name = '%s_nexrad_summary_%d' % (clean_var_boundary_name, prev_hours)
          wq_tests_data[var_name] = wq_defines.NO_DATA

        var_name = '%s_nexrad_dry_days_count' % (clean_var_boundary_name)
        wq_tests_data[var_name] = wq_defines.NO_DATA

        var_name = '%s_nexrad_rainfall_intensity' % (clean_var_boundary_name)
        wq_tests_data[var_name] = wq_defines.NO_DATA

        var_name = '%s_nexrad_total_1_day_delay' % (clean_var_boundary_name)
        wq_tests_data[var_name] = wq_defines.NO_DATA
        var_name = '%s_nexrad_total_2_day_delay' % (clean_var_boundary_name)
        wq_tests_data[var_name] = wq_defines.NO_DATA
        var_name = '%s_nexrad_total_3_day_delay' % (clean_var_boundary_name)
        wq_tests_data[var_name] = wq_defines.NO_DATA


    for hour in range(24,48,24):
      #wq_tests_data['%s_time_%d' % (self.hycom_data_prefix,hour)] = wq_defines.NO_DATA
      wq_tests_data['%s_avg_salinity_%d' % (self.hycom_data_prefix,hour)] = wq_defines.NO_DATA
      #wq_tests_data['%s_min_salinity_%d' % (self.hycom_data_prefix,hour)] = wq_defines.NO_DATA
      #wq_tests_data['%s_max_salinity_%d' % (self.hycom_data_prefix,hour)] = wq_defines.NO_DATA

      #wq_tests_data['%s_time_%d' % (self.copernicus_data_prefix,hour)] = wq_defines.NO_DATA
      wq_tests_data['%s_avg_salinity_%d' % (self.copernicus_data_prefix,hour)] = wq_defines.NO_DATA
      #wq_tests_data['%s_min_salinity_%d' % (self.copernicus_data_prefix,hour)] = wq_defines.NO_DATA
      #wq_tests_data['%s_max_salinity_%d' % (self.copernicus_data_prefix,hour)] = wq_defines.NO_DATA

      #wq_tests_data['%s_time_%d' % (self.rutgers_data_prefix,hour)] = wq_defines.NO_DATA
      wq_tests_data['%s_avg_salinity_%d' % (self.rutgers_data_prefix,hour)] = wq_defines.NO_DATA
      #wq_tests_data['%s_min_salinity_%d' % (self.rutgers_data_prefix,hour)] = wq_defines.NO_DATA
      #wq_tests_data['%s_max_salinity_%d' % (self.rutgers_data_prefix,hour)] = wq_defines.NO_DATA


    #wq_tests_data['%s_avg_water_temp_24' % (self.hycom_data_prefix)] = wq_defines.NO_DATA
    #wq_tests_data['%s_min_water_temp' % (self.hycom_data_prefix)] = wq_defines.NO_DATA
    #wq_tests_data['%s_max_water_temp' % (self.hycom_data_prefix)] = wq_defines.NO_DATA

    if self.logger:
      self.logger.debug("Finished creating and initializing data dict.")

    return

  def query_data(self, start_date, end_date, wq_tests_data):
    if self.logger:
      self.logger.debug("Site: %s start query data for datetime: %s" % (self.site.name, start_date))

    self.initialize_return_data(wq_tests_data)

    for platform in self.platforms_info:
      for obs_nfo in platform['observations']:
        self.get_platform_data(platform['platform_handle'],
                               obs_nfo['observation'], obs_nfo['uom'],
                               start_date,
                               wq_tests_data)

    self.get_nexrad_data(start_date, wq_tests_data)
    self.get_tide_data(start_date, wq_tests_data)
    self.get_rutgers_model_data(start_date, wq_tests_data)
    self.get_copernicus_model_data(start_date, wq_tests_data)
    self.get_hycom_model_data(start_date, wq_tests_data)

    if self.logger:
      self.logger.debug("Site: %s Finished query data for datetime: %s" % (self.site.name, start_date))


  def get_rutgers_model_data(self, start_date, wq_tests_data):
    if self.rutgers_current_date_ndx != -1:
      #start_date = timezone('UTC').localize(datetime.strptime('2009-12-25 00:00:00', "%Y-%m-%d %H:%M:%S"))
      # Ocean Time is seconds since 2006-01-01 00:00:00
      model_time = timezone('UTC').localize(datetime.strptime('2006-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'))
      #Run time is hours since 2013-05-18T00:00:00Z
      if self.rutgers_time is not None:
        model_time = timezone('UTC').localize(datetime.strptime('2013-05-18T00:00:00', '%Y-%m-%dT%H:%M:%S'))

      begin_date = start_date - timedelta(hours=192)
      end_date = start_date
      closest_start_ndx = \
        closest_end_ndx = -1
      if self.rutgers_time is None:
        if begin_date >= (model_time + timedelta(seconds=int(self.rutgers_ocean_time[0]))):

          start_time_delta = begin_date - model_time
          end_time_delta = end_date - model_time

          # The time dimension in the model is hours offset since the beginning_time above.
          offset_start = start_time_delta.total_seconds()
          offset_end = end_time_delta.total_seconds()
          closest_start_ndx = bisect_left(self.rutgers_ocean_time, offset_start)
          closest_end_ndx = bisect_left(self.rutgers_ocean_time, offset_end)

      else:
        if begin_date >= (model_time + timedelta(hours=int(self.rutgers_time[0]))):
          start_time_delta = begin_date - model_time
          end_time_delta = end_date - model_time

          # The time dimension in the model is hours offset since the beginning_time above.
          offset_start = start_time_delta.total_seconds() / (60.0 * 60.0)
          offset_end = end_time_delta.total_seconds() / (60.0 * 60.0)
          closest_start_ndx = bisect_left(self.rutgers_time, offset_start)
          closest_end_ndx = bisect_left(self.rutgers_time, offset_end)

      if closest_start_ndx != -1 and closest_end_ndx != -1:
        salinity_data = self.rutgers_model.variables['salt'][closest_start_ndx:closest_end_ndx + 1, 0,
                        self.rutgers_latli:self.rutgers_latui, self.rutgers_lonli:self.rutgers_lonui]
        """
        with open("/Users/danramage/tmp/kdh_salinity_rutgers_model.csv", "w") as out_file:
          #i0,i1,j0,j1 = bbox2ij(self.rutgers_model.variables['lon_rho'][:], self.rutgers_model.variables['lat_rho'][:], [-75.728228,-74.985642,35.152562,35.746631])
          lon_array = self.rutgers_model.variables['lon_rho'][self.rutgers_latli:self.rutgers_latui, self.rutgers_lonli:self.rutgers_lonui]
          lat_array = self.rutgers_model.variables['lat_rho'][self.rutgers_latli:self.rutgers_latui, self.rutgers_lonli:self.rutgers_lonui]

          #flat_lon_array = lon_array.flatten()
          #flat_lat_array = lat_array.flatten()
          lons_i,lons_j = np.shape(lon_array)
          lats_i, lats_j = np.shape(lat_array)
          out_file.write("Longitude,Latitude\n")
          out_file.write("%f,%f\n" % (lon_array[0][0], lat_array[0][0]))
          out_file.write("%f,%f\n" % (lon_array[0][-1], lat_array[0][-1]))
          out_file.write("%f,%f\n" % (lon_array[-1][0], lat_array[-1][0]))
          out_file.write("%f,%f\n" % (lon_array[-1][-1], lat_array[-1][-1]))

        """
        try:
          fill_value =  self.rutgers_model.variables['salt']._FillValue
        except Exception as e:
          fill_value = None

        pt_cell = closestCellFromPt(self.rutgers_cell_pt[0],self.rutgers_cell_pt[1],
                               self.rutgers_lon_array, self.rutgers_lat_array,
                               salinity_data[0],
                               fill_value,
                               self.rutgers_model.variables['mask_rho'][self.rutgers_latli:self.rutgers_latui, self.rutgers_lonli:self.rutgers_lonui])

        adjusted_end_ndx = None
        for hour in range(24, 48, 24):
          # Calculate the date we will start at referenced to the beginning_time of the data.
          begin_delta = (start_date - timedelta(hours=hour)) - model_time
          # Get the hours count.
          if self.rutgers_time is None:
            begin_cnt = begin_delta.total_seconds()
            # Find the index in the data our days count is closest to.
            begin_ndx = bisect_left(self.rutgers_ocean_time, begin_cnt)
          else:
            begin_cnt = (begin_delta.total_seconds() / (60.0 * 60.0))
            # Find the index in the data our days count is closest to.
            begin_ndx = bisect_left(self.rutgers_time, begin_cnt)
          # Offset indexes from the larger dataset to our subset.
          adjusted_begin_ndx = begin_ndx - closest_start_ndx
          if adjusted_end_ndx is None:
            adjusted_end_ndx = closest_end_ndx - closest_start_ndx

          # Get the last 24 hour average salinity data
          avg_salinity_pts = salinity_data[adjusted_begin_ndx:adjusted_end_ndx, int(pt_cell.x), int(pt_cell.y)]
          """
          if self.rutgers_time is None:
            wq_tests_data['%s_time_%d' % (self.rutgers_data_prefix, hour)] = model_time + timedelta(
              seconds=int(self.rutgers_ocean_time[begin_ndx]))
          else:
            wq_tests_data['%s_time_%d' % (self.rutgers_data_prefix, hour)] = model_time + timedelta(
            hours=int(self.rutgers_time[begin_ndx]))
          """
          wq_tests_data['%s_avg_salinity_%d' % (self.rutgers_data_prefix, hour)] = float(
            np.average(avg_salinity_pts))
          #wq_tests_data['%s_min_salinity_%d' % (self.rutgers_data_prefix, hour)] = float(avg_salinity_pts.min())
          #wq_tests_data['%s_max_salinity_%d' % (self.rutgers_data_prefix, hour)] = float(avg_salinity_pts.max())

          adjusted_end_ndx = adjusted_begin_ndx
        """
        with open("/Users/danramage/tmp/kdh_sites_salinity_rutgers.csv", "wa") as out_file:
          out_file.write("DateTime,Longitude,Latitude,Salinity\n")
          if self.rutgers_time is None:
            date_time = model_time + timedelta(seconds=self.rutgers_ocean_time[closest_start_ndx])
          else:
            date_time = model_time + timedelta(hours=self.rutgers_time[closest_start_ndx])
          out_file.write("%s,%f,%f,%f\n" % (date_time, self.rutgers_lon_array[int(pt_cell.x)][int(pt_cell.y)], self.rutgers_lat_array[int(pt_cell.x)][int(pt_cell.y)], salinity_data[0][int(pt_cell.x)][int(pt_cell.y)]))
        """

    return

  def get_copernicus_model_data(self, start_date, wq_tests_data):
    self.logger.debug("Start retrieving copernicus model data: %s" % (start_date))
    #Time is hours since 1950-01-01 00:00:00
    beginning_time = timezone('UTC').localize(datetime.strptime('1950-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'))
    begin_date = start_date - timedelta(hours=192)
    end_date = start_date
    if begin_date >= (beginning_time + timedelta(hours=int(self.copernicus_model_time[0]))):
      start_time_delta = begin_date - beginning_time
      end_time_delta = end_date - beginning_time

      #The time dimension in the model is hours offset since the beginning_time above.
      offset_start = start_time_delta.total_seconds() / (60.0 * 60.0)
      #offset_start = (start_time_delta.days) + (start_time_delta.seconds / (60.0 * 60.0))
      offset_end = end_time_delta.total_seconds() / (60.0 * 60.0)
      #offset_end = (end_time_delta.days) + (end_time_delta.seconds / (60.0 * 60.0))
      closest_start_ndx = bisect_left(self.copernicus_model_time, offset_start)
      closest_end_ndx = bisect_left(self.copernicus_model_time, offset_end)
      if closest_start_ndx != -1 and closest_end_ndx != -1:
        """
        with open("/Users/danramage/tmp/kdh_salinity_copernicus_model.csv", "w") as out_file:
          out_file.write("Longitude,Latitude,Lat Ndx,Lon Ndx\n")
          for lon_ndx in range(0,len(self.copernicus_lon_array)):
            for lat_ndx in range(0,len(self.copernicus_lat_array)):
              out_file.write("%f,%f,%d,%d\n" % (self.copernicus_lon_array[lon_ndx], self.copernicus_lat_array[lat_ndx], lon_ndx, lat_ndx))
        """
        try:
          if self.logger:
            self.logger.debug("Retrieving copernicus salinity data.")
          salinity_data = self.copernicus_model.variables['so'][closest_start_ndx:closest_end_ndx+1,0,self.copernicus_latli:self.copernicus_latui,self.copernicus_lonli:self.copernicus_lonui]
          pt = closestCellFromPtInPolygon(self.site.object_geometry,
                                          self.copernicus_lon_array, self.copernicus_lat_array,
                                          salinity_data[0],
                                          self.copernicus_model.variables['so']._FillValue,
                                          self.model_within_polygon)
          adjusted_end_ndx = None

          for hour in range(24, 48, 24):
            #Calculate the date we will start at referenced to the beginning_time of the data.
            begin_delta = (start_date - timedelta(hours=hour)) - beginning_time
            #Get the hours count.
            begin_hours_cnt = (begin_delta.total_seconds() / (60.0 * 60.0))
            #begin_days_cnt = (begin_delta.days) + (begin_delta.seconds / (60.0 * 60.0 * 24.0))
            #Find the index in the data our days count is closest to.
            begin_ndx = bisect_left(self.copernicus_model_time, begin_hours_cnt)
            #Offset indexes from the larger dataset to our subset.
            adjusted_begin_ndx = begin_ndx-closest_start_ndx
            if adjusted_end_ndx is None:
              adjusted_end_ndx = closest_end_ndx-closest_start_ndx
            #Validate we're getting the same dates.
            #if self.hycom_model_time[begin_ndx] != times_192[adjusted_begin_ndx]:
            #  if self.logger:
            #    self.logger.error("Times do not match")
            #Get the last 24 hour average salinity data
            avg_salinity_pts = salinity_data[adjusted_begin_ndx:adjusted_end_ndx,int(pt.y),int(pt.x)]

            #wq_tests_data['%s_time_%d' % (self.copernicus_data_prefix, hour)] = beginning_time + timedelta(hours=int(self.copernicus_model_time[begin_ndx]))
            wq_tests_data['%s_avg_salinity_%d' % (self.copernicus_data_prefix, hour)] = float(np.average(avg_salinity_pts))
            #wq_tests_data['%s_min_salinity_%d' % (self.copernicus_data_prefix, hour)] = float(avg_salinity_pts.min())
            #wq_tests_data['%s_max_salinity_%d' % (self.copernicus_data_prefix, hour)] = float(avg_salinity_pts.max())

            adjusted_end_ndx = adjusted_begin_ndx
          """
          with open("/Users/danramage/tmp/copernicus_salinity.csv", "w") as out_file:
            out_file.write("Date,Salinity\n")
            row = ["StartDateTime", "Lon", "Lat"]
            for hour in range(24, 192, 24):
              row.append("Model_Time_%d" % (hour))
              row.append("Salinity_%d" % (hour))
            out_file.write(",".join(row))
            out_file.write("\n")
            del row[:]

            row.append(str(start_date))
            row.append(str(self.copernicus_lon_array[int(pt.x)]))
            row.append(str(self.copernicus_lat_array[int(pt.y)]))
            for hour in range(24, 192, 24):
              row.append(str(wq_tests_data['%s_time_%d' % (self.copernicus_data_prefix, hour)]))
              row.append(str(wq_tests_data['%s_avg_salinity_%d' % (self.copernicus_data_prefix, hour)]))
            out_file.write(",".join(row))
            out_file.write("\n")
          """
        except Exception as e:
          self.logger.exception(e)
    return

  def get_hycom_model_data(self, start_date, wq_tests_data):
    if self.logger:
      self.logger.debug("Start retrieving hycom model data: %s" % (start_date))

    #if start_date >= self.
    #Hycom time is referenced in hours since 2000-01-01.
    beginning_time = timezone('UTC').localize(datetime.strptime('2000-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'))
    begin_date = start_date - timedelta(hours=192)
    end_date = start_date
    #Verify that the date we are interested should be in the hycom model data.
    if begin_date >= (beginning_time + timedelta(hours=self.hycom_model_time[0])):
      start_time_delta = begin_date - beginning_time
      end_time_delta = end_date - beginning_time

      #The time dimension in the model is hours offset since the beginning_time above.
      offset_start = start_time_delta.total_seconds() / (60.0 * 60.0)
      #offset_start = (start_time_delta.days) + (start_time_delta.seconds / (60.0 * 60.0))
      offset_end = end_time_delta.total_seconds() / (60.0 * 60.0)
      #offset_end = (end_time_delta.days) + (end_time_delta.seconds / (60.0 * 60.0))
      closest_start_ndx = bisect_left(self.hycom_model_time, offset_start)
      closest_end_ndx = bisect_left(self.hycom_model_time, offset_end)
      if closest_start_ndx != -1 and closest_end_ndx != -1:
        """
        with open("/Users/danramage/tmp/kdh_salinity_hycom_model.csv", "w") as out_file:
          out_file.write("Longitude,Latitude,Lat Ndx,Lon Ndx\n")
          for lon_ndx in range(0,len(self.hycom_lon_array)):
            for lat_ndx in range(0,len(self.hycom_lat_array)):
              out_file.write("%f,%f,%d,%d\n" % (self.hycom_lon_array[lon_ndx], self.hycom_lat_array[lat_ndx], lon_ndx, lat_ndx))
        """
        try:
          if self.logger:
            self.logger.debug("Retrieving hycom salinity data.")
          salinity_data = self.hycom_model.variables['salinity'][closest_start_ndx:closest_end_ndx+1,0,self.hycom_latli:self.hycom_latui,self.hycom_lonli:self.hycom_lonui]
          pt = closestCellFromPtInPolygon(self.site.object_geometry,
                                          self.hycom_lon_array, self.hycom_lat_array,
                                          salinity_data[0],
                                          self.hycom_model.variables['salinity']._FillValue,
                                          self.model_within_polygon)
          adjusted_end_ndx = None

          for hour in range(24, 48, 24):
            #Calculate the date we will start at referenced to the beginning_time of the data.
            begin_delta = (start_date - timedelta(hours=hour)) - beginning_time
            #Get the hours count.
            begin_hours_cnt = (begin_delta.total_seconds() / (60.0 * 60.0))
            #begin_days_cnt = (begin_delta.days) + (begin_delta.seconds / (60.0 * 60.0 * 24.0))
            #Find the index in the data our days count is closest to.
            begin_ndx = bisect_left(self.hycom_model_time, begin_hours_cnt)
            #Offset indexes from the larger dataset to our subset.
            adjusted_begin_ndx = begin_ndx-closest_start_ndx
            if adjusted_end_ndx is None:
              adjusted_end_ndx = closest_end_ndx-closest_start_ndx
            #Validate we're getting the same dates.
            #if self.hycom_model_time[begin_ndx] != times_192[adjusted_begin_ndx]:
            #  if self.logger:
            #    self.logger.error("Times do not match")
            #Get the last 24 hour average salinity data
            avg_salinity_pts = salinity_data[adjusted_begin_ndx:adjusted_end_ndx,int(pt.y),int(pt.x)]

            #wq_tests_data['%s_time_%d' % (self.hycom_data_prefix, hour)] = beginning_time + timedelta(hours=self.hycom_model_time[begin_ndx])
            wq_tests_data['%s_avg_salinity_%d' % (self.hycom_data_prefix, hour)] = float(np.average(avg_salinity_pts))
            #wq_tests_data['%s_min_salinity_%d' % (self.hycom_data_prefix, hour)] = float(avg_salinity_pts.min())
            #wq_tests_data['%s_max_salinity_%d' % (self.hycom_data_prefix, hour)] = float(avg_salinity_pts.max())

            adjusted_end_ndx = adjusted_begin_ndx
          """
          with open("/Users/danramage/tmp/hycom_salinity.csv", "w") as out_file:
            out_file.write("Date,Salinity\n")
            row = ["StartDateTime", "Lon", "Lat"]
            for hour in range(24, 192, 24):
              row.append("Hycom_Time_%d" % (hour))
              row.append("Salinity_%d" % (hour))
            out_file.write(",".join(row))
            out_file.write("\n")
            del row[:]

            row.append(str(start_date))
            row.append(str(self.hycom_lon_array[int(pt.x)]))
            row.append(str(self.hycom_lat_array[int(pt.y)]))
            for hour in range(24, 192, 24):
              row.append(str(wq_tests_data['%s_time_%d' % (self.hycom_data_prefix, hour)]))
              row.append(str(wq_tests_data['%s_avg_salinity_%d' % (self.hycom_data_prefix, hour)]))
            out_file.write(",".join(row))
            out_file.write("\n")
          """
          #Now get the rate of change from 24-192 hours back
          """
          ndx_list = []
          for time_ndx in range(24, 192, 24):
            date_ndx = start_date - timedelta(hours=time_ndx)
            time_delta = date_ndx - beginning_time
            #The time dimension in the model is hours offset since the beginning_time above.
            #offset_hours = (time_delta.days * 24) + (time_delta.seconds / (60 * 60 * 24))
            offset_ndx = (time_delta.days) + (time_delta.seconds / (60.0 * 60.0 * 24.0))
            ndx_list.append(bisect_left(self.hycom_model_time, offset_ndx))
          """
          """   
          if self.logger:
            self.logger.debug("Retrieving hycom water temp data.")

          last_24_delta = (start_date - timedelta(hours=24)) - beginning_time
          last_24_start = (last_24_delta.days) + (last_24_delta.seconds / (60.0 * 60.0 * 24.0))
          last_24_ndx = bisect_left(self.hycom_model_time, last_24_start)

          water_temp_data = self.hycom_model.variables['temperature'][last_24_ndx:closest_end_ndx,0,self.hycom_latli:self.hycom_latui,self.hycom_lonli:self.hycom_lonui]
          avg_water_temp_pts = water_temp_data[0:last_24_ndx-closest_start_ndx,pt.y,pt.x]
          #Get the water temperature data
          wq_tests_data['%s_avg_water_temp_24' % (self.hycom_data_prefix)] = float(np.average(avg_water_temp_pts))
          wq_tests_data['%s_min_water_temp' % (self.hycom_data_prefix)] = float(avg_water_temp_pts.min())
          wq_tests_data['%s_max_water_temp' % (self.hycom_data_prefix)] = float(avg_water_temp_pts.max())
          if self.logger:
            cell_lat = self.hycom_lat_array[pt.y]
            cell_lon = self.hycom_lon_array[pt.x]
            begin_dt = beginning_time + timedelta(days=(self.hycom_model_time[closest_start_ndx]))
            end_dt = beginning_time + timedelta(days=(self.hycom_model_time[closest_end_ndx]))
            self.logger.debug("Site: %s Dates: %s to %s closest cell@ Lat: %f(%d) Lon: %f(%d) Salinity 24 hrAvg: %f Water Temp Avg: %f"\
                              % (self.site.name,\
                                 begin_dt.strftime('%Y-%m-%dT%H:%M:%S'), end_dt.strftime('%Y-%m-%dT%H:%M:%S'),\
                                 cell_lat, pt.x, cell_lon, pt.y,\
                                 wq_tests_data['%s_avg_salinity_24' % (self.hycom_data_prefix)],wq_tests_data['%s_avg_water_temp_24'  % (self.hycom_data_prefix)]))
          """
        except Exception, e:
          if self.logger:
            self.logger.exception(e)

      else:
        if self.logger:
          self.logger.error("Cannot find start: %s or end: %s date range." % (offset_start, offset_end))
    else:
      if self.logger:
        self.logger.error("Date: %s out of range of data source." % (start_date))

    if self.logger:
      self.logger.debug("Finished retrieving hycom model data: %s" % (start_date))
    return

  def get_nexrad_data(self, start_date, wq_tests_data):
    start_time = time.time()
    if self.logger:
      self.logger.debug("Start retrieving nexrad data datetime: %s" % (start_date.strftime('%Y-%m-%d %H:%M:%S')))

    # Collect the radar data for the boundaries.
    for boundary in self.site.contained_by:
      clean_var_bndry_name = boundary.name.lower().replace(' ', '_')

      platform_handle = 'nws.%s.radarcoverage' % (boundary.name)
      if self.logger:
        self.logger.debug("Start retrieving nexrad platfrom: %s" % (platform_handle))
      # Get the radar data for previous 8 days in 24 hour intervals
      for prev_hours in range(24, 192, 24):
        var_name = '%s_nexrad_summary_%d' % (clean_var_bndry_name, prev_hours)
        radar_val = self.nexrad_db.getLastNHoursSummaryFromRadarPrecip(platform_handle,
                                                                       start_date,
                                                                       prev_hours,
                                                                       'precipitation_radar_weighted_average',
                                                                       'mm')
        if radar_val != None:
          # Convert mm to inches
          wq_tests_data[var_name] = radar_val * 0.0393701
        else:
          if self.logger:
            self.logger.error("No data available for boundary: %s Date: %s. Error: %s" % (
            var_name, start_date, self.nexrad_db.getErrorInfo()))

      # calculate the X day delay totals
      if wq_tests_data['%s_nexrad_summary_48' % (clean_var_bndry_name)] != wq_defines.NO_DATA and \
          wq_tests_data['%s_nexrad_summary_24' % (clean_var_bndry_name)] != wq_defines.NO_DATA:
        wq_tests_data['%s_nexrad_total_1_day_delay' % (clean_var_bndry_name)] = wq_tests_data['%s_nexrad_summary_48' % (
        clean_var_bndry_name)] - wq_tests_data['%s_nexrad_summary_24' % (clean_var_bndry_name)]

      if wq_tests_data['%s_nexrad_summary_72' % (clean_var_bndry_name)] != wq_defines.NO_DATA and \
          wq_tests_data['%s_nexrad_summary_48' % (clean_var_bndry_name)] != wq_defines.NO_DATA:
        wq_tests_data['%s_nexrad_total_2_day_delay' % (clean_var_bndry_name)] = wq_tests_data['%s_nexrad_summary_72' % (
        clean_var_bndry_name)] - wq_tests_data['%s_nexrad_summary_48' % (clean_var_bndry_name)]

      if wq_tests_data['%s_nexrad_summary_96' % (clean_var_bndry_name)] != wq_defines.NO_DATA and \
          wq_tests_data['%s_nexrad_summary_72' % (clean_var_bndry_name)] != wq_defines.NO_DATA:
        wq_tests_data['%s_nexrad_total_3_day_delay' % (clean_var_bndry_name)] = wq_tests_data['%s_nexrad_summary_96' % (
        clean_var_bndry_name)] - wq_tests_data['%s_nexrad_summary_72' % (clean_var_bndry_name)]

      prev_dry_days = self.nexrad_db.getPrecedingRadarDryDaysCount(platform_handle,
                                                                   start_date,
                                                                   'precipitation_radar_weighted_average',
                                                                   'mm')
      if prev_dry_days is not None:
        var_name = '%s_nexrad_dry_days_count' % (clean_var_bndry_name)
        wq_tests_data[var_name] = prev_dry_days

      rainfall_intensity = self.nexrad_db.calcRadarRainfallIntensity(platform_handle,
                                                                     start_date,
                                                                     60,
                                                                     'precipitation_radar_weighted_average',
                                                                     'mm')
      if rainfall_intensity is not None:
        var_name = '%s_nexrad_rainfall_intensity' % (clean_var_bndry_name)
        wq_tests_data[var_name] = rainfall_intensity

      if self.logger:
        self.logger.debug("Finished retrieving nexrad platfrom: %s" % (platform_handle))

    if self.logger:
      self.logger.debug("Finished retrieving nexrad data datetime: %s in %f seconds" % (start_date.strftime('%Y-%m-%d %H:%M:%S'),
                                                                                        time.time() - start_time))

  def get_tide_data(self, start_date, wq_tests_data):
    if self.logger:
      self.logger.debug("Start retrieving tide data for station: %s date: %s" % (self.tide_station, start_date))

    use_web_service = True
    if self.tide_data_obj is not None:
      use_web_service = False
      date_key = start_date.strftime('%Y-%m-%dT%H:%M:%S')
      if date_key in self.tide_data_obj:
        tide_rec = self.tide_data_obj[date_key]
        if tide_rec['range'] is not None:
          wq_tests_data['tide_range_%s' % (self.tide_station)] = tide_rec['range']

        if tide_rec['hh'] is not None:
          wq_tests_data['tide_hi_%s' % (self.tide_station)] = tide_rec['hh']

        if tide_rec['ll'] is not None:
          wq_tests_data['tide_lo_%s' % (self.tide_station)] = tide_rec['ll']

        if tide_rec['tide_stage'] is not None:
          wq_tests_data['tide_stage_%s' % (self.tide_station)] = tide_rec['tide_stage']

    #Save subordinate station values
    if wq_tests_data['tide_hi_%s'%(self.tide_station)] != wq_defines.NO_DATA:
      offset_hi = wq_tests_data['tide_hi_%s'%(self.tide_station)] * self.tide_offset_settings['hi_tide_height_offset']
      offset_lo = wq_tests_data['tide_lo_%s'%(self.tide_station)] * self.tide_offset_settings['lo_tide_height_offset']

      tide_station = self.tide_offset_settings['tide_station']
      wq_tests_data['tide_range_%s' % (tide_station)] = offset_hi - offset_lo
      wq_tests_data['tide_hi_%s' % (tide_station)] = offset_hi
      wq_tests_data['tide_lo_%s' % (tide_station)] = offset_lo

    if self.logger:
      self.logger.debug("Finished retrieving tide data for station: %s date: %s" % (self.tide_station, start_date))

    return

  def get_platform_data(self, platform_handle, variable, uom, start_date, wq_tests_data):
    start_time = time.time()
    try:
      self.logger.debug("Platform: %s Obs: %s(%s) Date: %s query" % (platform_handle, variable, uom, start_date))

      station = platform_handle.replace('.', '_')
      var_name = '%s_avg_%s' % (station, variable)
      end_date = start_date
      begin_date = start_date - timedelta(hours=24)
      dir_id = None
      sensor_id = self.xenia_obs_db.sensorExists(variable, uom, platform_handle, 1)
      if variable == 'wind_speed':
        dir_id = self.xenia_obs_db.sensorExists('wind_from_direction', 'degrees_true', platform_handle, 1)

      if sensor_id is not -1 and sensor_id is not None:
        recs = self.xenia_obs_db.session.query(sl_multi_obs) \
          .filter(sl_multi_obs.m_date >= begin_date.strftime('%Y-%m-%dT%H:%M:%S')) \
          .filter(sl_multi_obs.m_date < end_date.strftime('%Y-%m-%dT%H:%M:%S')) \
          .filter(sl_multi_obs.sensor_id == sensor_id) \
          .filter(or_(sl_multi_obs.qc_level == qaqcTestFlags.DATA_QUAL_GOOD, sl_multi_obs.qc_level == None)) \
          .order_by(sl_multi_obs.m_date).all()
        if dir_id is not None:
          dir_recs = self.xenia_obs_db.session.query(sl_multi_obs) \
            .filter(sl_multi_obs.m_date >= begin_date.strftime('%Y-%m-%dT%H:%M:%S')) \
            .filter(sl_multi_obs.m_date < end_date.strftime('%Y-%m-%dT%H:%M:%S')) \
            .filter(sl_multi_obs.sensor_id == dir_id) \
            .filter(or_(sl_multi_obs.qc_level == qaqcTestFlags.DATA_QUAL_GOOD, sl_multi_obs.qc_level == None)) \
            .order_by(sl_multi_obs.m_date).all()

        if len(recs):
          if variable == 'wind_speed':
            if sensor_id is not None and dir_id is not None:
              wind_dir_tuples = []
              direction_tuples = []
              scalar_speed_avg = None
              speed_count = 0
              for wind_speed_row in recs:
                for wind_dir_row in dir_recs:
                  if wind_speed_row.m_date == wind_dir_row.m_date:
                    # self.logger.debug("Building tuple for Speed(%s): %f Dir(%s): %f" % (
                    # wind_speed_row.m_date, wind_speed_row.m_value, wind_dir_row.m_date, wind_dir_row.m_value))
                    if scalar_speed_avg is None:
                      scalar_speed_avg = 0
                    scalar_speed_avg += wind_speed_row.m_value
                    speed_count += 1
                    # Vector using both speed and direction.
                    wind_dir_tuples.append((wind_speed_row.m_value, wind_dir_row.m_value))
                    # Vector with speed as constant(1), and direction.
                    direction_tuples.append((1, wind_dir_row.m_value))
                    break

              if len(wind_dir_tuples):
                avg_speed_dir_components = calcAvgSpeedAndDir(wind_dir_tuples)
                self.logger.debug("Platform: %s Avg Wind Speed: %f(m_s-1) %f(mph) Direction: %f" % (platform_handle,
                                                                                                    avg_speed_dir_components[
                                                                                                      0],
                                                                                                    avg_speed_dir_components[
                                                                                                      0],
                                                                                                    avg_speed_dir_components[
                                                                                                      1]))

                # Unity components, just direction with speeds all 1.
                avg_dir_components = calcAvgSpeedAndDir(direction_tuples)
                scalar_speed_avg = scalar_speed_avg / speed_count
                wq_tests_data[var_name] = scalar_speed_avg
                wind_dir_var_name = '%s_avg_%s' % (station, 'wind_from_direction')
                wq_tests_data[wind_dir_var_name] = avg_dir_components[1]
                self.logger.debug(
                  "Platform: %s Avg Scalar Wind Speed: %f(m_s-1) %f(mph) Direction: %f" % (platform_handle,
                                                                                           scalar_speed_avg,
                                                                                           scalar_speed_avg,
                                                                                           avg_dir_components[1]))
          #Calculate vector direction.
          elif variable == 'sea_surface_wave_to_direction':
            direction_tuples = []
            for dir_row in recs:
              # Vector with speed as constant(1), and direction.
              direction_tuples.append((1, dir_row.m_value))

            if len(direction_tuples):
              # Unity components, just direction with speeds all 1.
              avg_dir_components = calcAvgSpeedAndDir(direction_tuples)
              wq_tests_data[var_name] = avg_dir_components[1]
              self.logger.debug(
                "Platform: %s Avg Scalar Direction: %f" % (platform_handle,
                                                           avg_dir_components[1]))

          else:
            wq_tests_data[var_name] = sum(rec.m_value for rec in recs) / len(recs)
            self.logger.debug("Platform: %s Avg %s: %f Records used: %d" % (
              platform_handle, variable, wq_tests_data[var_name], len(recs)))

            if variable == 'water_conductivity':
              water_con = wq_tests_data[var_name]
              #if uom == 'uS_cm-1':
              water_con = water_con / 1000.0
              salinity_var = '%s_avg_%s' % (station, 'salinity')
              wq_tests_data[salinity_var] = 0.47413 / (math.pow((1 / water_con), 1.07) - 0.7464 * math.pow(10, -3))
              self.logger.debug("Platform: %s Avg %s: %f Records used: %d" % (
                platform_handle, 'salinity', wq_tests_data[salinity_var], len(recs)))
        else:
          self.logger.error(
            "Platform: %s sensor: %s(%s) Date: %s had no data" % (platform_handle, variable, uom, start_date))
      else:
        self.logger.error("Platform: %s sensor: %s(%s) does not exist" % (platform_handle, variable, uom))
      self.logger.debug("Platform: %s query finished in %f seconds" % (platform_handle, time.time()-start_time))
    except Exception as e:
      self.logger.exception(e)
      return False

    return True

def parse_file(**kwargs):
  start_time = time.time()
  #est_tz = timezone('US/Eastern')
  utc_tz = timezone('UTC')
  logger = logging.getLogger(__name__)
  logger.debug("Starting parse_file")
  dates = []
  start_date = kwargs.get('start_date', None)
  sample_collection = kwargs['samples_collection']
  sample_data_file = kwargs['data_file']

  logger.debug("Getting data from: %s" % (sample_data_file))
  header = [
    'site_id',
    'site_name',
    'date',
    'entero ssm',
    'entero gm'
  ]
  with open(sample_data_file, "r") as data_file:
    csv_reader = csv.DictReader(data_file, header)
    for row_ndx,row in enumerate(csv_reader):
      if row_ndx > 0:
        add_rec = False
        sample_data = nc_wq_sample_data()
        date_obj = (utc_tz.localize(datetime.strptime(row['date'], '%Y-%m-%d %H:%M:%S')))
        if start_date is not None:
          if date_obj >= start_date:
            add_rec = True
        else:
          add_rec = False
        if add_rec:
          sample_data.date_time = date_obj
          sample_data.station = row['site_id']
          sample_data.entero_gm = row['entero gm']
          sample_data.entero_ssm = row['entero ssm']
          sample_collection.append(sample_data)

  logger.debug("Finished processing file in %f seconds." % (time.time()-start_time))
  return

def main():
  parser = optparse.OptionParser()
  parser.add_option("--ConfigFile", dest="config_file", default=None,
                    help="INI Configuration file." )
  parser.add_option("--OutputDirectory", dest="output_dir", default=None,
                    help="Directory to save the historical data site files." )
  (options, args) = parser.parse_args()


  try:
    config_file = ConfigParser.RawConfigParser()
    config_file.read(options.config_file)

    logConfFile = config_file.get('logging', 'config_file')

    logging.config.fileConfig(logConfFile)
    logger = logging.getLogger('build_historical_logger')
    logger.info("Log file opened.")


    boundaries_location_file = config_file.get('boundaries_settings', 'boundaries_file')
    sites_location_file = config_file.get('boundaries_settings', 'sample_sites')
    wq_historical_db = config_file.get('database', 'name')

  except ConfigParser.Error, e:
    import traceback
    traceback.print_exc(e)
    sys.exit(-1)

  else:
    #Load the sample site information. Has name, location and the boundaries that contain the site.
    wq_sites = wq_sample_sites()
    wq_sites.load_sites(file_name=sites_location_file, boundary_file=boundaries_location_file)

    wq_historical_data = kdh_historical_wq_data(config_file=options.config_file)

    sample_data_directory = '/Users/danramage/Documents/workspace/WaterQuality/NorthCarolina-OuterBanks/data/historical/sample_data'
    historical_sample_files = os.listdir(sample_data_directory)
    #start_date = timezone('UTC').localize(datetime.strptime('2005-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'))
    utc_tz = timezone('UTC')
    est_tz = timezone('US/Eastern')
    data_start_date = utc_tz.localize(datetime.strptime('2005-01-01 00:00:00', '%Y-%m-%d %H:%M:%S'))
    for site in wq_sites:
      out_file = os.path.join(options.output_dir, "%s_historical_data.csv" % (site.name))
      write_header = True
      with open(out_file, 'w') as site_data_file:
        try:
          hycom_data_prefix = config_file.get(site.description, 'hycom_prefix')
          copernicus_data_prefix = config_file.get(site.description, 'copernicus_prefix')
          rutgers_data_prefix = config_file.get(site.description, 'rutgers_prefix')
          rutgers_cell_point = tuple(float(pt) for pt in config_file.get(site.description, 'rutgers_cell_loc').split(','))

          # Get the station specific tide stations
          tide_station = config_file.get(site.description, 'tide_station')
          offset_tide_station = config_file.get(site.description, 'offset_tide_station')
          offset_key = "%s_tide_data" % (offset_tide_station)
          tide_offset_settings = {
            'tide_station': config_file.get(offset_key, 'station_id'),
            'hi_tide_time_offset': config_file.getint(offset_key, 'hi_tide_time_offset'),
            'lo_tide_time_offset': config_file.getint(offset_key, 'lo_tide_time_offset'),
            'hi_tide_height_offset': config_file.getfloat(offset_key, 'hi_tide_height_offset'),
            'lo_tide_height_offset': config_file.getfloat(offset_key, 'lo_tide_height_offset')
          }
          #Get the platforms the site will use
          platforms = config_file.get(site.description, 'platforms').split(',')
          platform_nfo = []
          for platform in platforms:
            obs_uoms = config_file.get(platform,'observation').split(';')
            obs_uom_nfo = []
            for nfo in obs_uoms:
              obs,uom = nfo.split(',')
              obs_uom_nfo.append({'observation': obs,
                                  'uom': uom})
            platform_nfo.append({'platform_handle': config_file.get(platform,'handle'),
                                 'observations': obs_uom_nfo})

        except ConfigParser.Error, e:
          if logger:
            logger.exception(e)

        file_name = site.name
        for file in historical_sample_files:
          if file.find(file_name) != -1:
            samples_collection = wq_samples_collection()
            full_path = os.path.join(sample_data_directory, file)
            parse_file(data_file=full_path,
                       samples_collection=samples_collection,
                       start_date=data_start_date)

            try:
              sample_recs = samples_collection[site.name]
            except (KeyError, Exception) as e:
              logger.exception(e)
            else:
              sample_recs.sort(key=lambda x: x.date_time, reverse=False)
              auto_num = 1
              for sample_data in sample_recs:
                start_date = sample_data.date_time
                try:
                  wq_date_time_local = sample_data.date_time.astimezone(est_tz)
                  site_data = OrderedDict([
                    ('autonumber', auto_num),
                    ('station_name', site.name),
                    ('station_desc', site.description),
                    ('sample_datetime', wq_date_time_local),
                    ('sample_datetime_utc', sample_data.date_time),
                    ('enterococcus_value', sample_data.entero_ssm),
                  ])
                  wq_historical_data.reset(site=site,
                                           tide_station=tide_station,
                                           tide_offset_params=tide_offset_settings,
                                           hycom_prefix=hycom_data_prefix,
                                           copernicus_prefix=copernicus_data_prefix,
                                           start_date=sample_data.date_time,
                                           rutgers_cell_point=rutgers_cell_point,
                                           rutgers_prefix=rutgers_data_prefix,
                                           platform_info=platform_nfo)
                  wq_historical_data.query_data(sample_data.date_time, sample_data.date_time, site_data)

                  header_buf = []
                  data = []
                  for key in site_data:
                    if write_header:
                      header_buf.append(key)
                    if site_data[key] != wq_defines.NO_DATA:
                      data.append(str(site_data[key]))
                    else:
                      data.append("")
                  if write_header:
                    site_data_file.write(",".join(header_buf))
                    site_data_file.write('\n')
                    header_buf[:]
                    write_header = False

                  site_data_file.write(",".join(data))
                  site_data_file.write('\n')
                  site_data_file.flush()
                  data[:]

                  auto_num += 1
                except Exception, e:
                  if logger:
                    logger.exception(e)
                  sys.exit(-1)
    """
    site_data = OrderedDict([('autonumber', 1),
                             ('station_name', row['SPLocation']),
                             ('sample_datetime', wq_date.strftime("%Y-%m-%d %H:%M:%S")),
                             ('sample_datetime_utc', wq_utc_date.strftime("%Y-%m-%d %H:%M:%S")),
                             ('County', row['County']),
                             ('enterococcus_value', row['enterococcus']),
                             ('enterococcus_code', row['enterococcus_code'])])
    """
  return

if __name__ == "__main__":
  main()