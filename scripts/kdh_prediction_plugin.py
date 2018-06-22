import sys
from os.path import dirname, realpath
sys.path.append('../commonfiles/python')
sys.path.append(dirname(realpath(__file__)))
import time
from wq_prediction_plugin import wq_prediction_engine_plugin
from kdh_prediction_engine import kdh_prediction_engine

class kdh_prediction_plugin(wq_prediction_engine_plugin):

  def do_processing(self, **kwargs):
    start_do_processing_time = time.time()
    self.logger.debug("Starting do_processing")
    dates_to_process = kwargs.get('processing_dates', [])
    for process_date in dates_to_process:
      chs_engine = kdh_prediction_engine()
      chs_engine.run_wq_models(begin_date=process_date,
                      config_file_name=self.config_file)
    self.logger.debug("do_processing took: %f seconds" % (time.time() - start_do_processing_time))
    self.logger.debug("Finished do_processing")

