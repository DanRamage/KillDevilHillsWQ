import sys
sys.path.append('../../commonfiles/python')
import os
import logging.config
import requests
import time
import xlrd
from datetime import datetime
from pytz import timezone
from yapsy.IPlugin import IPlugin
from multiprocessing import Process
import traceback

from wq_output_results import wq_sample_data,wq_samples_collection,wq_advisories_file,wq_station_advisories_file
from nc_sample_data import nc_wq_sample_data

from wq_sites import wq_sample_sites
from data_collector_plugin import data_collector_plugin
import ConfigParser

"""
https://reports.ncdenr.org/BOE/OpenDocument/opendoc/CrystalReports/viewrpt.cwr?
id=1887482&
apsuser=REC_U&
apspassword=mf_pr0_b0&
apsauthtype=secEnterprise&
cmd=export&
export_fmt=U2FXLS:4&
prompt0=Date(2018,1,1)&
prompt1=Date(2018,12,31)&
prompt2=*&
prompt3=*&
prompt4=*&
prompt5=Public
"""
def download_historical_sample_data(**kwargs):
  logger=logging.getLogger(kwargs['logger_name'])

  func_start_time = time.time()
  start_year = kwargs['start_year']
  end_year = kwargs['end_year']
  base_url = kwargs['url']
  logger.info("download_historical_sample_data started")
  output_directory = kwargs['output_directory']
  params = {'id': '1887482',
    'apsuser': 'REC_U',
    'apspassword': 'mf_pr0_b0',
    'apsauthtype': 'secEnterprise',
    'cmd': 'export',
    'export_fmt': 'U2FXLS:4',
    'prompt0': '',
    'prompt1': '',
    'prompt2': '*',
    'prompt3': '*',
    'prompt4': '*',
    'prompt5': 'Public'
  }

  for year in range(start_year, end_year, -1):
    start_time = time.time()
    params['prompt0'] = 'Date(%d,1,1)' % (year)
    params['prompt1'] = 'Date(%d,12,31)' % (year)
    logger.debug("Preparing to retrieve data year %d" % (year))
    try:
      req = requests.get(base_url, params=params)
      if req.status_code == 200:
        logger.debug("Successfully retrieved data year %d in %f seconds. url: %s" % (year, time.time()-start_time, req.url))
        out_filename = os.path.join(output_directory, 'sample_data_%d.xls' % (year))
        logger.debug("Saving sample data to: %s" % (out_filename))
        start_time = time.time()
        try:
          with open(out_filename, 'wb') as out_data_file:
            for chunk in req.iter_content(chunk_size=1024):
              out_data_file.write(chunk)
            logger.debug("Saved sample data in %f seconds" % (time.time()-start_time))
        except (IOError, Exception) as e:
          logger.exception(e)
      else:
        logger.error("Failed to retrieve data, status code: %d" % (req.status_code))
    except Exception as e:
      logger.exception(e)

    logger.info("download_historical_sample_data finished in: %f seconds" % (time.time()-func_start_time))

  return

def parse_files(**kwargs):
  logger=logging.getLogger(kwargs['logger_name'])
  wq_sites = kwargs['sample_sites']
  src_data_directory = kwargs['src_data_directory']
  output_directory = kwargs['output_directory']
  monitoring_sites = [
    {
      'name': "drain pipe at oregon st",
      'site_id': '14A'
    },
    {
      'name': "drain pipe at lake dr beach access",
      'site_id': '15',
    },
    {
      'name': "drain pipe at curlew street",
      'site_id': '16A',
    },
    {
      'name': "drain pipe at mile post 10.5",
      'site_id': '16',
    },
    {
      'name': "drain pipe at mile post 12.5",
      'site_id': '17A'
    },
    {
      'name': "drain pipe at conch street",
      'site_id': '17'
    },
    {
      'name': "drain pipe at s nags head/federal park border",
      'site_id': '22'
    },
    {
      'name':"drain pipe at martin street",
      'site_id': '85A'
    },
    {
      'name': "drain pipe at mp 8 3/4",
      'site_id': '85'
    }
  ]

  header = [
    'site_id',
    'site_name',
    'date',
    'entero ssm',
    'entero gm'
  ]
  historical_files = os.listdir(src_data_directory)
  historical_files.reverse()

  for file in historical_files:
    try:
      path,ext = os.path.splitext(file)
      if ext == '.xls':
        wq_data_collection = wq_samples_collection()
        parse_excel_data(os.path.join(src_data_directory,file), monitoring_sites, wq_data_collection, kwargs['logger_name'])

        json_results_file = os.path.join(output_directory, 'kdh_beach_advisories.json')
        logger.debug("Creating beach advisories file: %s" % (json_results_file))
        try:
          current_advisories = wq_advisories_file(wq_sites)
          current_advisories.create_file(json_results_file, wq_data_collection)
        except Exception as e:
          logger.exception(e)
        try:
          for site in wq_sites:
            logger.debug("Creating site: %s advisories file" % (site.name))
            site_advisories = wq_station_advisories_file(site)
            site_advisories.create_file(output_directory, wq_data_collection)
        except Exception as e:
          logger.exception(e)
        """
        for wq_data_key in wq_data_collection:
          wq_data = wq_data_collection[wq_data_key]
          wq_data.sort(key=lambda x: x.date_time, reverse=True)
          clean_site_name = wq_data[0].station.replace(' ', '_').replace('/', '-')
          site_data_file = os.path.join(output_directory, "%s.csv" % (wq_data[0].site_id))
          write_header = False
          if os.path.isfile(site_data_file) is False:
            write_header = True
          with open(site_data_file, "a") as site_file:
            if write_header:
              site_file.write(",".join(header))
              site_file.write("\n")
            for data_rec in wq_data:
              site_file.write('%s,%s,%s,%s,%s\n' % (data_rec.site_id, data_rec.station, data_rec.date_time.strftime('%Y-%m-%d %H:%M:%S'), str(data_rec.entero_ssm), str(data_rec.entero_gm)))
        """
    except Exception as e:
      logger.exception(e)
  return

def parse_excel_data(file, monitoring_sites, wq_data_collection, logger_name):
  logger=logging.getLogger(logger_name)

  logger.debug("Parsing file: %s" % (file))
  wb = xlrd.open_workbook(filename = file)
  sheet = wb.sheet_by_name('Sheet1')
  row_headers = []
  results_ndx = \
    station_ndx = \
    county_ndx = \
    site_id_ndx = \
    entero_gm_ndx = \
    entero_ssm_ndx = \
    date_ndx = None

  est_tz = timezone('US/Eastern')
  utc_tz = timezone('UTC')
  current_excel_station = None
  current_site_station_ndx = None
  for row_ndx, data_row in enumerate(sheet.get_rows()):
    if row_ndx != 0:
      if data_row[county_ndx].value.strip() != "Dare":
        continue
      station_name = data_row[station_ndx].value.strip().lower()
      site_id = data_row[site_id_ndx].value.strip()
      logger.debug("Processing row: %d Site: %s-%s" % (row_ndx, station_name, site_id))
      try:
        matched = False
        if current_excel_station is None or current_excel_station != station_name:
          for ndx, site_nfo in enumerate(monitoring_sites):
            site_name = site_nfo['name']
            id = site_nfo['site_id']
            if site_id == id:
              matched = True
              station_name = id
              current_excel_station = station_name
              current_site_station_ndx = ndx
        else:
          matched = True
      except ValueError as e:
        pass
      else:
        if matched:
          try:
            wq_sample_rec = nc_wq_sample_data()
            wq_sample_rec.station = station_name
            wq_sample_rec.site_id = site_id
            try:
              date_val = xlrd.xldate.xldate_as_datetime(data_row[date_ndx].value, wb.datemode)
            except Exception as e:
                logger.error("Date format error on line: %d" % (row_ndx))
                logger.exception(e)
                break
            wq_sample_rec.date_time = (est_tz.localize(date_val)).astimezone(utc_tz)
            wq_sample_rec.entero_gm = data_row[entero_gm_ndx].value
            wq_sample_rec.entero_ssm = data_row[entero_ssm_ndx].value
            wq_sample_rec.value = wq_sample_rec.entero_ssm
            logger.debug("Site: %s Date: %s SSM: %s GM: %s" % (wq_sample_rec.station,
                                                          wq_sample_rec.date_time,
                                                          wq_sample_rec.entero_ssm,
                                                          wq_sample_rec.entero_gm))
            wq_data_collection.append(wq_sample_rec)
          except Exception as e:
            logger.exception(e)
    else:
      """
      Date	Area	Site	County	Description	Run	Location	Tier	24-hr Precipitation	Salinity	Water Temp	Tide	Current	Wind	Entero MPN1	Entero MPN2	Entero MPN3	Entero SSM	Entero GM	"Entero
      CFU1"	Entero CFU2	Entero CFU3	Entero SSM CFU	Entero_GM2	# of Days	E.Coli MPN1	E.Coli MPN2	E.Coli MPN3	E.Coli SSM	E.Coli GM	# of Days	Fecal MPN1	Fecal MPN2	Fecal MPN3	Fecal SSM	Fecal GM	# of Days	E.Coli CFU1	E.Coli CFU2	E.Coli CFU3	E.Coli SSF	E.Coli GM2	Fecal CFU1	Fecal CFU2	Fecal CFU3	Fecal SSC	Fecal GM2      
      """
      for cell in data_row:
        row_headers.append(cell.value)
      #Save the indexes for quicker access
      station_ndx = row_headers.index('Description')
      date_ndx = row_headers.index('Date')
      entero_gm_ndx = row_headers.index('Entero GM')
      entero_ssm_ndx = row_headers.index('Entero SSM')
      county_ndx = row_headers.index('County')
      site_id_ndx = row_headers.index('Site')
  return

class kdh_sample_data_collector_plugin(data_collector_plugin):
  def __init__(self):
    Process.__init__(self)
    IPlugin.__init__(self)

  def initialize_plugin(self, **kwargs):
    #data_collector_plugin.initialize_plugin(self, **kwargs)
    try:
      logger = logging.getLogger(self.__class__.__name__)
      plugin_details = kwargs['details']
      self.base_url = plugin_details.get('Settings', 'url')
      self.source_directory = plugin_details.get('Settings', 'download_directory')
      self.sample_site_directory = plugin_details.get('Settings', 'sample_site_directory')
      self.boundaries_file = plugin_details.get('Settings', 'boundaries_file')
      self.sample_sites_file = plugin_details.get('Settings', 'sample_sites')
      self.log_config = plugin_details.get("Settings", "log_config")

      return True
    except Exception as e:
      logger.exception(e)
    return False

  def run(self):
    logger = None
    start_time = time.time()
    try:

      #self.logging_client_cfg['disable_existing_loggers'] = True
      #logging.config.dictConfig(self.logging_client_cfg)
      logging.config.fileConfig(self.log_config)
      logger = logging.getLogger(self.__class__.__name__)
      logger.debug("run started.")

      try:
        wq_sites = wq_sample_sites()
        if wq_sites.load_sites(file_name=self.sample_sites_file, boundary_file=self.boundaries_file):

          start_year = datetime.now().date().year

          download_historical_sample_data(output_directory=self.source_directory,
                                          url=self.base_url,
                                          start_year=start_year,
                                          end_year=start_year - 1,
                                          logger_name=self.__class__.__name__)

          parse_files(sample_sites=wq_sites,
                      src_data_directory=self.source_directory,
                      output_directory=self.sample_site_directory,
                      logger_name=self.__class__.__name__)
        else:
          logger.error("Failed to load sites file: %s %s" % (self.sample_sites_file, self.boundaries_file))

      except (IOError, Exception) as e:
        if logger is not None:
          logger.exception(e)

    except Exception as e:
      if logger is not None:
        logger.exception(e)
      else:
        traceback.print_exc(e)
    finally:
      if logger is not None:
        logger.debug("run finished in %f seconds" % (time.time()-start_time))

    return
