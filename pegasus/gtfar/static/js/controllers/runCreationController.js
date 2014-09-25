/*
 * Copyright 2007-2014 University Of Southern California
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/**
 * Created by dcbriggs on 7/28/14.
 */

define(["angular"],
function(angular) {
    "use strict"
    var runsCreationController = function($scope, $window, $http, $state) {

        var file;
        var FIVE_GIGABYTES = 5000000000;
        $scope.alerts = [];
        $scope.getAlertIcon = function(status) {
            return (status == 'success') ? "text-success fa fa-check-circle" : "text-danger fa fa-exclamation-triangle";
        };
        $scope.closeAlert = function(index) {
            $scope.alerts.splice(index, 1);
        };
        $scope.uploadProgress = null;
        $scope.addingRun = null;
        $scope.strandRuleOptions = [{rule: "unstranded"},{rule: "sense"},{rule: "anti-sense"}];
        // Defaults

        // We have to define this out of the run object becaues we need to link it to .rule in the end
        $scope.strandRule = $scope.strandRuleOptions[0];

        $scope.run = {
            genome : "Genome", // TODO: add this in a bit $scope.genome,
            gtf : "gencode.v19.annotation.gtf", // TODO: Add this in a bit $scope.gtf,
            status : -1, // default status to running
            userName : "genericUser",
            readLength : 100,
            mismatches : 3,
            trimUnmapped : false,
            mapFiltered : false,
            genSplice : false
        };

        $scope.uploadFile = function() {
            /*
             * Upload the file
             * Angular does not natively work with a formData object so we have to override
             * the default behavior.  Refer to the link below
             * http://uncorkedstudios.com/blog/multipartformdata-file-upload-with-angularjs
             */
            var  inputFile = document.getElementById("inputFile");

            if (inputFile && !inputFile.value) {
                return;
            }

            if( inputFile.files[0].size > FIVE_GIGABYTES) {
                $scope.alerts.push({
                    "type" : "danger",
                    "message" : "That file is too large, please select one smaller than five gigabytes"
                });
                return;
            }

            $scope.alerts = [];

            var fileData = new FormData(document.getElementById("runForm"));
            var xhr = new XMLHttpRequest();
            // TODO: disable the browse button
            xhr.upload.onprogress = function(e) {
                $scope.$apply(function() {
                    if (e.lengthComputable) {
                        $scope.uploadProgress = Math.round(e.loaded / e.total * 100);
                    }
                });
            };

            xhr.upload.onload = function(e) {
                // Event listener for when the file completed uploading
                $scope.$apply(function() {
                    $scope.uploadProgress = null;
                    $scope.alerts.push({
                        "type" : "success",
                        "message" : "File uploaded successfully"
                    });
                    // TODO: make the browse button clickable again
                });
            };

            xhr.onerror = function(e) {
                $scope.$apply(function() {
                    $scope.uploadProgress = null;
                    $scope.alerts.push({
                        "type" : "danger",
                        "message" : "File failed to upload"
                    });
                    // TODO: make the browse button clickable again
                });
            };

            xhr.open('POST', $window.apiLinks.upload);
            xhr.send(fileData);

        };

        $scope.cancel = function() {
            $state.go('runs');
        };

        $scope.addRun = function() {
            $scope.alerts = [];
            var invalid = false;
            if($scope.form.$invalid) {
                if($scope.form.name.$invalid) {
                    invalid = true;
                    $scope.alerts.push({'type' : 'danger', 'message' : 'The name field is required and must be an alphanumeric value'});
                }
                if($scope.form.mismatches.$invalid) {
                    invalid = true;
                    $scope.alerts.push({'type' : 'danger', 'message' : 'The mismatches field is required and must be an integer between 0 and 8 (inclusive)'});
                }
                if($scope.form.readLength.$invalid) {
                    invalid = true;
                    $scope.alerts.push({'type' : 'danger', 'message' : 'The read length field is required and must be an integer between 50 and 128 (inclusive)'});
                }
                if($scope.form.email.$invalid) {
                    if($scope.form.email.$viewValue.indexOf(',') != -1) { // Multiple emails input
                        var emails = $scope.form.email.$viewValue.split(",");
                        for(var i = 0; i < emails.length; i++) {
                            if(emails[i].indexOf("@") == -1) {
                                invalid = true;
                                $scope.alerts.push({'type' : 'danger', 'message' : 'Emails must be properly formatted user@location.domain'});
                                return;
                            }
                        }
                    }
                    else {
                        invalid = true;
                        $scope.alerts.push({'type' : 'danger', 'message' : 'Emails must be properly formatted user@location.domain'});
                    }
                }
                if(invalid) {
                    return;
                }
            }

            $scope.run["strandRule"] = $scope.strandRule.rule;
            $scope.run.email = $scope.form.email.$viewValue;
            try {
                $scope.run["filename"] = document.getElementById("inputFile").files[0].name;
            }
            catch(error)
            {
                return; // An error message should show to the user so we just ned to return
            }
            $scope.addingRun = true;
            // Add the run to the database, the runId will be used to create the folder
            $http.post($window.apiLinks.runs, $scope.run).success(function(data) {
                $scope.addingRun = null;
                $scope.alerts.push({'type' : 'success', 'message' : 'Run successfully created and starting now, please wait'});
                $state.go('runDetails', {id : data.id});
             }).error(function(data) {
                $scope.addingRun = null;
                if(data.message && data.message == "400: Bad Request") {
                    $http.get($window.apiLinks.errors).success(function(data) {
                        if(data.errors) {
                            for (var i = 0; i < data.errors.length; i++) {
                                $scope.alerts.push({'type': 'danger', 'message': data.errors[i].message});
                            }
                        }
                    }).error(function(data){});
                }
                else if(data.code && data.code == 503) {
                    $scope.alerts.push({'type' : 'danger', 'message' : 'Could not connect to server, please check your connection.'});
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
        getFullConstructor : getFullConstructor,
        getName : getMinimalConstructor
    };

});
