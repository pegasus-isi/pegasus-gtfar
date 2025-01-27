/*
 * #  Copyright 2007-2012 University Of Southern California
 * #
 * #  Licensed under the Apache License, Version 2.0 (the "License");
 * #  you may not use this file except in compliance with the License.
 * #  You may obtain a copy of the License at
 * #
 * #  http://www.apache.org/licenses/LICENSE-2.0
 * #
 * #  Unless required by applicable law or agreed to in writing,
 * #  software distributed under the License is distributed on an "AS IS" BASIS,
 * #  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * #  See the License for the specific language governing permissions and
 * #  limitations under the License.
 * # Revision : $Revision: 2012 $
 */

/**
 * Created by dcbriggs on 7/28/14.
 */

define(["angular"],
    function (angular) {
        "use strict"
        var runsCreationController = function ($scope, $window, $http, $state) {

            var file;
            var FIVE_GIGABYTES = {{ config['MAX_UPLOAD_SIZE'] }};
            $scope.alerts = [];
            $scope.getAlertIcon = function (status) {
                return (status == 'success') ? "text-success fa fa-check-circle" : "text-danger fa fa-exclamation-triangle";
            };
            $scope.closeAlert = function (index) {
                $scope.alerts.splice(index, 1);
            };
            $scope.uploadProgress = null;
            $scope.addingRun = null;
            $scope.speciesOptions = [
                {% for s, o in species.iteritems() %}
                    "{{ s }}"
                    {% if not loop.last %}
                        ,
                    {% endif %}
                {% endfor %}
            ];
            $scope.strandRuleOptions = [
                {rule: "unstranded"},
                {rule: "sense"},
                {rule: "anti-sense"}
            ];
            // Defaults

            // We have to define this out of the run object because we need to link it to .rule in the end
            $scope.species = $scope.speciesOptions[0];
            $scope.strandRule = $scope.strandRuleOptions[0];

            $scope.run = {
                status: -1, // default status to running
                userName: "genericUser",
                readLength: 100,
                mismatches: 3,
                trimUnmapped: false,
                mapFiltered: false,
                genSplice: false
            };

            $scope.sample = $window.apiLinks.sample;

            $scope.uploadFile = function () {
                /*
                 * Upload the file
                 * Angular does not natively work with a formData object so we have to override
                 * the default behavior.  Refer to the link below
                 * http://uncorkedstudios.com/blog/multipartformdata-file-upload-with-angularjs
                 */
                var inputFile = document.getElementById("inputFile");

                if (inputFile && !inputFile.value) {
                    return;
                }

                if (inputFile.files[0].size > FIVE_GIGABYTES) {
                    $scope.$apply(function () {
                        $scope.alerts.push({
                            "type": "danger",
                            "message": "That file is too large, please select one smaller than five gigabytes"
                        });
                    });

                    return;
                }

                $scope.alerts = [];

                var fileData = new FormData(document.getElementById("runForm"));
                var xhr = new XMLHttpRequest();
                // TODO: disable the browse button
                xhr.upload.onprogress = function (e) {
                    $scope.$apply(function () {
                        if (e.lengthComputable) {
                            $scope.uploadProgress = Math.round(e.loaded / e.total * 100);
                            $scope.uploadProgressStyle = {width: $scope.uploadProgress + '%'};
                        }
                    });
                };

                xhr.upload.onload = function (e) {
                    // Event listener for when the file completed uploading
                    $scope.$apply(function () {
                        $scope.uploadProgress = null;
                        $scope.alerts.push({
                            "type": "success",
                            "message": "File uploaded successfully"
                        });
                    });
                };

                xhr.onerror = function (e) {
                    $scope.$apply(function () {
                        $scope.uploadProgress = null;
                        $scope.alerts.push({
                            "type": "danger",
                            "message": "File failed to upload"
                        });
                    });
                };

                xhr.open('POST', $window.apiLinks.upload + "/" + $window.tempFolder);
                xhr.send(fileData);

            };

            $scope.cancel = function () {
                $state.go('runs');
            };

            $scope.addRun = function () {
                $scope.alerts = [];
                var invalid = false;
                if ($scope.form.$invalid) {
                    if ($scope.form.name.$invalid) {
                        invalid = true;
                        $scope.alerts.push({'type': 'danger', 'message': 'The name field is required and must be an alphanumeric value (underscores acceptable)'});
                    }
                    if ($scope.form.mismatches.$invalid) {
                        invalid = true;
                        $scope.alerts.push({'type': 'danger', 'message': 'The mismatches field is required and must be an integer between 0 and 8 (inclusive)'});
                    }
                    if ($scope.form.readLength.$invalid) {
                        invalid = true;
                        $scope.alerts.push({'type': 'danger', 'message': 'The read length field is required and must be an integer between 50 and 128 (inclusive)'});
                    }
                    // Check for the empty string because we could have a case where the user enters and then erases an email
                    if ($scope.form.email.$invalid && $scope.form.email.$viewValue != "") {
                        if ($scope.form.email.$viewValue.indexOf(',') != -1) { // Multiple emails input
                            var emails = $scope.form.email.$viewValue.split(",");
                            for (var i = 0; i < emails.length; i++) {
                                if (emails[i].indexOf("@") == -1) {
                                    invalid = true;
                                    $scope.alerts.push({'type': 'danger', 'message': 'Emails must be properly formatted example@domain.com'});
                                    return;
                                }
                            }
                        }
                        else {
                            invalid = true;
                            $scope.alerts.push({'type': 'danger', 'message': 'Emails must be properly formatted example@domain.com'});
                        }
                    }
                    if (invalid) {
                        return;
                    }
                }

                if ($scope.run["genSplice"] && $scope.run["readLength"] < 75) {
                    $scope.alerts.push({'type': 'danger', 'message': 'When the Generate New Splice Candidates option is true, the read length must be >= 75'});
                    return;
                }

                $scope.run["strandRule"] = $scope.strandRule.rule;
                $scope.run["species"] = $scope.species;
                $scope.run["uploadFolder"] = $window.tempFolder;
                $scope.run.email = $scope.form.email.$viewValue;
                try {
                    $scope.run["filename"] = document.getElementById("inputFile").files[0].name;
                }
                catch (error) {
                    return; // An error message should show to the user so we just need to return
                }
                $scope.addingRun = true;
                // Add the run to the database, the runId will be used to create the folder
                $http.post($window.apiLinks.runs, $scope.run).success(function (data) {
                    $scope.addingRun = null;
                    $state.go('runDetails', {name: data.name});
                }).error(function (data, status) {
                    $scope.addingRun = null;
                    if (data.code && data.message) {
                        $scope.alerts.push({
                            'message': data.message,
                            'type': 'danger'
                        });
                    }
                    else if (status == 400) {
                        for (var i = 0; i < data.errors.length; i++) {
                            $scope.alerts.push({'type': 'danger', 'message': data.errors[i].message});
                        }
                    }
                    else if (data.code && data.code == 503) {
                        $scope.alerts.push({'type': 'danger', 'message': 'Could not connect to server, please check your connection.'});
                    }
                    else {
                        $scope.alerts.push({'type': 'danger', 'message': 'Unknown error occurred while retrieving error analysis.  Please contact the developers'});
                    }

                });


            }
        };

        function getFullConstructor() {
            return ["$scope", "$window", "$http", "$state", runsCreationController];
        }

        // For testing so we can mock the injected variables
        function getMinimalConstructor() {
            return runsCreationController;
        }

        return {
            getFullConstructor: getFullConstructor,
            getName: getMinimalConstructor
        };

    });
