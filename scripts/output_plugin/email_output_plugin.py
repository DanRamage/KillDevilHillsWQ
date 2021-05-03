import sys
sys.path.append('../../commonfiles/python')

from mako.template import Template
from mako import exceptions as makoExceptions
from smtp_utils import smtpClass
import os
import logging.config
from output_plugin import output_plugin

class email_output_plugin(output_plugin):
  def __init__(self):
    output_plugin.__init__(self)
    self.logger = logging.getLogger(__name__)


  def initialize_plugin(self, **kwargs):
    try:
      details = kwargs['details']
      self.mailhost = details.get("Settings", "mailhost")
      self.mailport = None
      self.fromaddr = details.get("Settings", "fromaddr")
      self.toaddrs = details.get("Settings", "toaddrs").split(',')
      self.subject = details.get("Settings", "subject")
      self.user = details.get("Settings", "user")
      self.password = details.get("Settings", "password")
      self.result_outfile = details.get("Settings", "results_outfile")
      self.results_template = details.get("Settings", "results_template")
      self.report_url = details.get("Settings", "report_url")
      return True
    except Exception as e:
      self.logger.exception(e)
    return False

  def emit(self, **kwargs):
    if self.logger:
      self.logger.debug("Starting emit for email output.")
    try:
      mytemplate = Template(filename=self.results_template)
      file_ext = os.path.splitext(self.result_outfile)
      file_parts = os.path.split(file_ext[0])
      #Add the prediction date into the filename
      file_name = "%s-%s%s" % (file_parts[1], kwargs['prediction_date'].replace(':', '_').replace(' ', '-'), file_ext[1])
      out_filename = os.path.join(file_parts[0], file_name)
      with open(out_filename, 'w') as report_out_file:
        report_url = '%s/%s' % (self.report_url, file_name)
        results_report = mytemplate.render(ensemble_tests=kwargs['ensemble_tests'],
                                                prediction_date=kwargs['prediction_date'],
                                                execution_date=kwargs['execution_date'],
                                                report_url=report_url)
        report_out_file.write(results_report)
    except TypeError as e:
      if self.logger:
        self.logger.exception(makoExceptions.text_error_template().render())
    except (IOError,AttributeError,Exception) as e:
      if self.logger:
        self.logger.exception(e)
    else:
      try:
        subject = self.subject % (kwargs['prediction_date'])
        #Now send the email.
        smtp = smtpClass(host=self.mailhost, user=self.user, password=self.password)
        smtp.rcpt_to(self.toaddrs)
        smtp.from_addr(self.fromaddr)
        smtp.subject(subject)
        smtp.message(results_report)
        smtp.send(content_type="html")
      except Exception as e:
        if self.logger:
          self.logger.exception(e)
    if self.logger:
      self.logger.debug("Finished emit for email output.")

