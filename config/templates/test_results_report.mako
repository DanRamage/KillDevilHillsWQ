<!DOCTYPE html>

<html lang="en">
    <head>
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <link href="http://howsthebeach.org/static/css/bootstrap/css/bootstrap.min.css" rel="stylesheet">
      <link href=http://howsthebeach.org//static/css/bootstrap/css/css/bootstrap-theme.min.css" rel="stylesheet">

      <title>Kill Devil Hills Water Quality Prediction Results</title>
    </head>
    <body>
        <style>
          .high_bacteria {
            background-color: #ff3633;
          }
          .medium_bacteria {
            background-color: #fff45c;
          }

        </style>
        <div class="container">
            <div class="row">
              <div class="col-xs-12">
                <h1>Kill Devil Hills Quality Prediction Results</h1>
                <h2>Prediction for: ${prediction_date}</h2>
                <h3>Prediction executed: ${execution_date}</h3>
              </div>
            </div>
            </br>
            <div class = "row">
              <div class="col-xs-12">
                <h4><a href='${report_url}'>Report File</a> Click here if the report below displays incorrectly.</h4>
              </div>
            </div>
            % for site_data in ensemble_tests:
                <div class="row">
                    <div class="col-xs-6">
                      <h2>Site: ${site_data['metadata'].name}</h2>
                    </div>
                </div>
                <div class="row">
                    <div class="col-xs-6">
                      <h3>
                        % if str(site_data['models'].ensemblePrediction) == 'LOW':
                          <span>
                        % elif str(site_data['models'].ensemblePrediction) == 'MEDIUM':
                          <span class="medium_bacteria">
                        % else:
                          <span class="high_bacteria">
                        % endif
                          Overall Prediction: ${str(site_data['models'].ensemblePrediction)}
                        </span>
                      </h3>
                    </div>
                </div>
                % if site_data['entero_value'] is not None:
                    <div class="row">
                        <table class="table table-bordered">
                            <tr>
                                <th> Site Sampled Entero Value </th>
                            </tr>
                            <tr>
                                <td>
                                    ${"%.2f" % (site_data['entero_value'])}
                                </td>
                            </tr>
                        </table>
                    </div>
                % endif
                % if site_data['statistics'] is not None:
                    <div class="row">
                        <table class="table table-bordered">
                            <tr>
                                <th>Minimum Entero</th>
                                <th>Maximum Entero</th>
                                <th>Average Entero</th>
                                <th>Median Entero</th>
                                <th>Geometric Mean</th>
                                <th>StdDev Entero</th>
                            </tr>
                            <tr>
                              % if site_data['statistics'].minVal is not None:
                                <td>${"%.2f" % (site_data['statistics'].minVal)}</td>
                              % else:
                                <td></td>
                              % endif
                              % if site_data['statistics'].maxVal is not None:
                                <td>${"%.2f" % (site_data['statistics'].maxVal)}</td>
                              % else:
                                <td></td>
                              % endif
                              % if site_data['statistics'].average is not None:
                                  <td>${"%.2f" % (site_data['statistics'].average)}</td>
                              % else:
                                <td></td>
                              % endif
                              % if site_data['statistics'].median is not None:
                                  <td>${"%.2f" % (site_data['statistics'].median)}</td>
                              % else:
                                <td></td>
                              % endif
                              % if site_data['statistics'].geometric_mean is not None:
                                  <td>${"%.2f" % (site_data['statistics'].geometric_mean)}</td>
                              % else:
                                <td></td>
                              % endif
                              % if site_data['statistics'].stdDev is not None:
                                  <td>${"%.2f" % (site_data['statistics'].stdDev)}</td>
                              % else:
                                <td></td>
                              % endif
                            </tr>
                        </table>
                    </div>
                % endif
                <div class="row">
                    <table class="table table-bordered">
                        <tr>
                            <th>Model Name</th>
                            <th>Prediction Level</th>
                            <th>Prediction Value</th>
                            <th>Data Used</th>
                        </tr>
                        % for test_obj in site_data['models'].tests:
                            % if test_obj is not None:
                                %if test_obj.mlrResults is not None:
                                    % if test_obj.mlrResult < 36:
                                      <tr>
                                    % elif test_obj.mlrResult >= 36 and test_obj.mlrResult < 104:
                                      <tr class="medium_bacteria">
                                    % else:
                                      <tr class="high_bacteria">
                                    % endif
                                %else:
                                  <td>
                                      NO TEST
                                  </td>
                                %endif

                                  <td>
                                    ${test_obj.model_name}
                                  </td>
                                  <td>
                                    ${test_obj.predictionLevel.__str__()}
                                  </td>
                                  <td>
                                    % if test_obj.mlrResult is not None:
                                      ${"%.2f" % (test_obj.mlrResult)}
                                    % else:
                                      NO TEST
                                    % endif
                                  </td>
                                  <td>
                                    % for key in test_obj.data_used:
                                      % if test_obj.data_used[key] != -9999:
                                        ${key}: ${test_obj.data_used[key]}
                                      % else:
                                        ${key}: Data unavailable
                                      % endif
                                        </br>
                                    % endfor
                                  </td>
                                </tr>
                            %else:
                                <tr>
                                    <td>
                                        NO TEST
                                    </td>
                                </tr>
                            % endif
                        % endfor
                    </table>
                </div>
            % endfor
        </div>
    </body>
</html>
