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
 * Created by dcbriggs on 7/14/14.
 */
define(["moment"],
    function (moment) {
        "use strict"

        var runDetailsController = function ($scope, $state, $stateParams, $http, $window) {
            $scope.alerts = [];
            $scope.getAlertIcon = function (status) {
                return (status == 'success') ? "text-success fa fa-check-circle" : "text-danger fa fa-exclamation-triangle";
            };
            $scope.closeAlert = function (index) {
                $scope.alerts.splice(index, 1);
            };

            $scope.inputDownload = $window.apiLinks.runs + "/" + $stateParams.name + "/input/";
            $scope.outputDownload = $window.apiLinks.runs + "/" + $stateParams.name + "/output/";


            $scope.getFormattedDate = function (dateString) {
                var createdOn = moment.utc(dateString).local();
                var s = createdOn.format('dddd MMMM DD, YYYY hh:mm:ss A');
                s += ' (' + createdOn.fromNow() + ')';
                return s;
            };

            $scope.GTFAR_S3_BUCKET = $window.GTFAR_S3_BUCKET;

            $scope.analyzeURL = $window.apiLinks.runs + "/" + $stateParams.name + $window.apiLinks.analyze;

            function getRun() {
                $http.get($window.apiLinks.runs + "/" + $stateParams.name).success(function (data) {
                    $scope.run = data;
                    $scope.run.emails = null;
                    if ($scope.run.email) {
                        $scope.run.emails = $scope.run.email.split(',');
                    }
                    // If the run was stopped then we don't need to get the status
                    if ($scope.run.status != 256) {
                        setTimeout(getStatus, 1500);

                    }
                    else {
                        getOutputFiles();
                        getErrorReport();
                    }
                }).error(function (data) {
                    console.error(data);
                });
            }

            function getStatus() {
                $http.get($window.apiLinks.runs + "/" + $stateParams.name + $window.apiLinks.status).success(function (data) {
                    if (statusChanged(data)) {
                        $scope.status = data;
                        $scope.progressStyle = {width: $scope.status.progress + '%'};

                        getOutputFiles();
                        if($scope.status.failed) {
                            getErrorReport();
                        }
                        else {
                            clearErrorReport();
                        }
                    }
                    // We don't want the calls to continue when we leave this state
                    if ($state.includes("runDetails")) {
                        if ($scope.status.state <= 3) {
                            setTimeout(getStatus, 30000);
                        }
                    }

                }).error(function (data) {
                    if(data.code && data.message) {
                        $scope.alerts.push({
                            'message': data.message + ".  Error Code: " + data.code,
                            'type': 'danger'
                        });
                    }
                    else {
                        $scope.alerts.push({'type': 'danger', 'message': 'Unknown error occurred while retrieving run status.  Please contact the developers'});
                    }
                });
            }

            function statusChanged(newStatus) {
                if (!$scope.status)
                    return true;
                for (var key in $scope.status) {
                    if ($scope.status[key] != newStatus[key])
                        return true;
                }
                return false;
            }

            function getOutputFiles() {
                $http.get($window.apiLinks.runs + "/" + $stateParams.name + $window.apiLinks.outputs).success(function (data) {
                    $scope.outputFiles = data.objects;
                    // TODO: need to present links to download the output files
                }).error(function (data) {
                    $scope.alerts.push({'type': 'danger', 'message': 'Error retrieving output files. Please contact the developers'});
                });
            }

            function getErrorReport() {
                $http.get($window.apiLinks.runs + "/" + $stateParams.name + $window.apiLinks.analyze).success(function(data) {
                    $scope.errorReport = [];

                    var summarySplit = data.match(/\*{4,}([^\*]+)\*{4,}([^\*]+)\*{4,}([^\*]+)\*{4,}([\s\S]+)/);
                    $scope.errorReport.push({'section' : summarySplit[1], 'details' : [], 'subsections' : []});
                    var summaryDetails = summarySplit[2].trim().split('\n');
                    for(var i = 0; i < summaryDetails.length; i++) {
                        $scope.errorReport[0].details.push({'text' : summaryDetails[i]});
                    }

                    // This will give us all the names of the files and the lines to split on in order to isolate the details
                    var fileDetails = summarySplit[4];
                    var fileSplits = fileDetails.trim().match(/={4,}[^=]+={4,}/g);

                    for(var i = 0; i < fileSplits.length; i++) {
                        $scope.errorReport.push({'section' : fileSplits[i].match(/[^=]+/)[0], 'details' : [], 'subsections' : []});
                        var headerRemoved = fileDetails.split(fileSplits[i]);
                        var previousFileDetails = headerRemoved[0].trim();
                        // This will be true for all splits except for the first one, and because split separates the parts before
                        // and after the splitting string on the ith iteration we will populate the details for the i-1th file
                        // but because we have the summary section as the 0th index of the errorReport
                        if(previousFileDetails != "") {
                            // To make sure we add data to the proper section of the report if some parts are missing
                            var previousSectionIndex = -1;
                            var summaryHeader = previousFileDetails.match(/-{4,}(.+?Summary)-{4,}/);
                            var splitByHeader;
                            if(summaryHeader) {
                                splitByHeader = previousFileDetails.split(summaryHeader[0]);

                                // The previous section
                                fillErrorDetails($scope.errorReport[i], previousSectionIndex, splitByHeader[0].trim().split('\n'));
                                $scope.errorReport[i].subsections.push({'section' : summaryHeader[1], 'details' : []});
                                previousFileDetails = splitByHeader[1];
                                previousSectionIndex += 1;
                            }

                            var stdoutHeader = previousFileDetails.match(/-{4,}(.+?stdout)-{4,}/);
                            if(stdoutHeader) {
                                splitByHeader = previousFileDetails.split(stdoutHeader[0]);

                                // Add the previous section
                                fillErrorDetails($scope.errorReport[i], previousSectionIndex, splitByHeader[0].trim().split('\n'));
                                $scope.errorReport[i].subsections.push({'section' : stdoutHeader[1], 'details' : []});
                                previousFileDetails = splitByHeader[1];
                                previousSectionIndex += 1;
                            }

                            var stderrHeader = previousFileDetails.match(/-{4,}(.+?stderr)-{4,}/);
                            if(stderrHeader) {
                                splitByHeader = previousFileDetails.split(stderrHeader[0]);

                                // Add the previous section
                                fillErrorDetails($scope.errorReport[i], previousSectionIndex, splitByHeader[0].trim().split('\n'));
                                $scope.errorReport[i].subsections.push({'section': stderrHeader[1], 'details' : []});
                                previousFileDetails = splitByHeader[1];
                                previousSectionIndex += 1;
                            }
                            // Add the final section
                            fillErrorDetails($scope.errorReport[i], previousSectionIndex, previousFileDetails.trim().split('\n'));
                        }

                        // Special case for the last file
                        if(i == fileSplits.length - 1) {
                            var currentFileDetails = headerRemoved[1].trim();
                            var previousSectionIndex = -1;

                            var summaryHeader = currentFileDetails.match(/-{4,}(.+?Summary)-{4,}/);
                            var splitByHeader;
                            if(summaryHeader) {
                                splitByHeader = currentFileDetails.split(summaryHeader[0]);
                                // The first section will give any general details on the file
                                fillErrorDetails($scope.errorReport[i+1], previousSectionIndex, splitByHeader[0].trim().split('\n'));
                                $scope.errorReport[i+1].subsections.push({'section' : summaryHeader[1], 'details' : []});
                                currentFileDetails = splitByHeader[1];
                                previousSectionIndex += 1;
                            }

                            var stdoutHeader = currentFileDetails.match(/-{4,}(.+?stdout)-{4,}/);
                            if(stdoutHeader) {
                                splitByHeader = currentFileDetails.split(stdoutHeader[0]);

                                // Add the previous section
                                fillErrorDetails($scope.errorReport[i+1], previousSectionIndex, splitByHeader[0].trim().split('\n'));
                                $scope.errorReport[i+1].subsections.push({'section' : stdoutHeader[1], 'details' : []});
                                currentFileDetails = splitByHeader[1];
                                previousSectionIndex += 1;
                            }

                            var stderrHeader = currentFileDetails.match(/-{4,}(.+?stderr)-{4,}/);
                            if(stderrHeader) {
                                splitByHeader = currentFileDetails.split(stderrHeader[0]);

                                // Add the previous section
                                fillErrorDetails($scope.errorReport[i+1], previousSectionIndex, splitByHeader[0].trim().split('\n'));
                                $scope.errorReport[i+1].subsections.push({'section': stderrHeader[1], 'details' : []});
                                currentFileDetails = splitByHeader[1];
                                previousSectionIndex += 1;
                            }
                            // Add the final section
                            fillErrorDetails($scope.errorReport[i+1], previousSectionIndex, currentFileDetails.trim().split('\n'));
                        }

                        // This should be everything after the last header we split by
                        fileDetails = headerRemoved[1].trim();
                    }
                }).error(function(data) {
                    if(data.code && data.message) {
                        $scope.alerts.push({
                            'message': data.message + ".  Error Code: " + data.code,
                            'type': 'danger'
                        });
                    }
                    else {
                        $scope.alerts.push({'type': 'danger', 'message': 'Unknown error occurred while retrieving error analysis.  Please contact the developers'});
                    }
                });
            }

            function fillErrorDetails(report, sectionIndex, details) {
                var reportSection = (sectionIndex == -1) ? report.details : report.subsections[sectionIndex].details;
                for(var i = 0; i < details.length; i++){
                    reportSection.push({'text' : details[i]});
                }
            }

            function clearErrorReport() {
                $scope.errorReport = [];
            }

            $scope.getStatusColor = function (prefix) {
                if ($scope.run && $scope.run.status > 0) {
                    if ($scope.run.status == 256) {
                        return prefix + "-warning";
                    }
                    return prefix + "-danger";
                }
                if ($scope.status) { // We have to make sure we've gotten data back from the server
                    if ($scope.status.state == 5) { // Failed
                        return prefix + "-danger";
                    }
                    if ($scope.status.state == 4) { // Success
                        return prefix + "-success";
                    }
                    if ($scope.status.failed > 0) { // Failing
                        return prefix + "-warning";
                    }
                    // Running
                    return prefix + "-info";
                }
                return prefix + "-danger";

            };

            // Pass in either the option 'bar' or 'alert'
            $scope.showBarOrAlert = function (option) {
                if ($scope.run && $scope.run.status > 0) {
                    return option == "alert";
                }
                if ($scope.status) {
                    if ($scope.status.state > 3) {
                        return option == "alert";
                    }
                    return option == "bar";
                }
                return false; // Hides everything until we have data
            };

            $scope.getStatusText = function (short) {
                if ($scope.run && $scope.run.status == 256) {
                    if (short) {
                        return "Aborted"
                    }
                    return "Run aborted before completion!";
                }
                if ($scope.status) { // We have to make sure we've gotten data back from the server
                    if ($scope.status.state == 5) { // Failed
                        if (short) {
                            return "Failed";
                        }
                        return "Run Failed!";
                    }
                    if ($scope.status.state == 4) { // Success
                        if (short) {
                            return "Successful";
                        }
                        return "Run Succeeded!";
                    }
                    if ($scope.status.failed > 0) { // Failing
                        if (short) {
                            return "Failing";
                        }
                        return "Failures found, but run still progressing.";
                    }
                    // Running
                    if (short) {
                        return "Running";
                    }
                    return "Still Running.";
                }
            };

            $scope.stopRun = function () {
                var doubleCheck = confirm("Are you sure you want to stop " + $scope.run.name + "?  You will not be able to restart it.");
                if(!doubleCheck){
                    return;
                }
                $http.get($window.apiLinks.runs + "/" + $stateParams.name + $window.apiLinks.stop).success(function (data) {
                    $scope.run.status = 256;

                }).error(function (data) {
                    if(data.code && data.message){
                        $scope.alerts.push({
                            'message': data.message,
                            'type': 'danger'
                        });
                    }
                    else {
                        $scope.alerts.push({
                            'message': "Unknown error occurred while attempting to stop the run.  Please contact the development team",
                            'type': 'danger'
                        });
                    }
                    $scope.run.status = 0;
                });
            };

            $scope.gotoList = function () {
                $state.go("runs");
            };

            getRun();
        };

        // Return the array so the dependecy injections are accounted for
        return ["$scope", "$state", "$stateParams", "$http", "$window", runDetailsController];
    });
