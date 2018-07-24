import sys
sys.path.append('../../commonfiles/python')
import logging.config
import json

from output_plugin import output_plugin

class json_output_plugin(output_plugin):
  def __init__(self):
    output_plugin.__init__(self)
    self.logger = logging.getLogger(__name__)

  def initialize_plugin(self, **kwargs):
    try:
      details = kwargs['details']

      self.json_outfile = details.get("Settings", "json_outfile")
      return True
    except Exception as e:
      self.logger.exception(e)
    return False

  def emit(self, **kwargs):
    if self.logger:
      self.logger.debug("Starting emit for json output.")

    ensemble_data = kwargs['ensemble_tests']
    try:
      with open(self.json_outfile, 'w') as json_output_file:
        station_data = {'features' : [],
                        'type': 'FeatureCollection'}
        features = []
        for rec in ensemble_data:
          site_metadata = rec['metadata']
          test_results = rec['models']
          if 'statistics' in rec:
            stats = rec['statistics']
          test_data = []
          for test in test_results.tests:
            test_data.append({
              'name': test.model_name,
              'p_level': test.predictionLevel.__str__(),
              'p_value': test.mlrResult,
              'data': test.data_used
            })
          features.append({
            'type': 'Feature',
            'geometry' : {
              'type': 'Point',
              'coordinates': [site_metadata.object_geometry.x, site_metadata.object_geometry.y]
            },
            'properties': {
              'desc': site_metadata.name,
              'ensemble': str(test_results.ensemblePrediction),
              'station': site_metadata.name,
              'tests': test_data
            }
          })
        station_data['features'] = features
        json_data = {
          'status': {'http_code': 200},
          'contents': {
            'run_date': kwargs['execution_date'],
            'testDate': kwargs['prediction_date'],
            'stationData': station_data
          }
        }
        try:
          json_output_file.write(json.dumps(json_data, sort_keys=True))
        except Exception,e:
          if self.logger:
            self.logger.exception(e)
    except IOError,e:
      if self.logger:
        self.logger.exception(e)
    if self.logger:
      self.logger.debug("Finished emit for json output.")
    return
