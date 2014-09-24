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
function(moment) {
    "use strict"

    var runDetailsController = function($scope, $state, $stateParams, $http, $window) {

        $scope.alerts = [];
        $scope.getAlertIcon = function(status) {
            return (status == 'success') ? "text-success fa fa-check-circle" : "text-danger fa fa-exclamation-triangle";
        };
        $scope.closeAlert = function(index) {
            $scope.alerts.splice(index, 1);
        };

        $scope.inputDownload = $window.apiLinks.download + "/" + $stateParams.id + "/input/";
        $scope.outputDownload = $window.apiLinks.download + "/" + $stateParams.id + "/output/";


        $scope.getFormattedDate = function(dateString) {
            var createdOn = moment.utc(dateString).local();
            var s = createdOn.format('dddd MMMM DD, YYYY hh:mm:ss A');
            s += ' (' + createdOn.fromNow() + ')';
            return s;
        };

        function getRun() {
            $http.get($window.apiLinks.runs + "/" + $stateParams.id).success(function(data) {
                $scope.run = data;
                $scope.run.emails = null;
                if($scope.run.email) {
                    $scope.run.emails = $scope.run.email.split(',');
                }
                // If the run was stopped then we don't need to get the status
                if($scope.run.status <= 0) {
                    setTimeout(getStatus, 500);

                }
                else {
                    getOutputFiles();
                }
            }).error(function(data) {
                console.error(data);
            });
        }

        function getStatus() {
            $http.get($window.apiLinks.runs + "/" + $stateParams.id + $window.apiLinks.status).success(function(data) {
                if(statusChanged(data)) {
                    $scope.status = data;
                    getOutputFiles();
                }
                // We don't want the calls to continue when we leave this state
                if($state.includes("runDetails")) {
                    if($scope.status.state <= 3) {
                        setTimeout(getStatus, 2000);
                    }
                }

            }).error(function(data) {
                console.error(data);
            });
        }

        function statusChanged(newStatus) {
            if(!$scope.status)
                return true;
            for(var key in $scope.status) {
                if($scope.status[key] != newStatus[key])
                    return true;
            }
            return false;
        }

        function getOutputFiles() {
            $http.get($window.apiLinks.runs + "/" + $stateParams.id + $window.apiLinks.outputs).success(function(data) {
                $scope.outputFiles = data.objects;
                // TODO: need to present links to download the output files
            }).error(function(data) {
                $scope.alerts.push({'type' : 'danger', 'message' : 'Error retrieving output files. Please contact the developers'});
            });
        }

        $scope.getStatusColor = function(prefix) {
            if($scope.run && $scope.run.status > 0) {
                return prefix + "-danger";
            }
            if($scope.status) { // We have to make sure we've gotten data back from the server
                if($scope.status.state == 5) { // Failed
                    return prefix + "-danger";
                }
                if($scope.status.state == 4) { // Success
                    return prefix + "-success";
                }
                if($scope.status.failed > 0) { // Failing
                    return prefix + "-warning";
                }
                // Running
                return prefix + "-info";
            }
            return prefix + "-danger";

        };

        // Pass in either the option 'bar' or 'alert'
        $scope.showBarOrAlert = function(option) {
            if($scope.run && $scope.run.status > 0) {
                return option == "alert";
            }
            if($scope.status) {
                if($scope.status.state > 3) {
                    return option == "alert";
                }
                return option == "bar";
            }
            return false; // Hides everything until we have data
        };

        $scope.getStatusText = function(short) {
            if($scope.run && $scope.run.status > 0) {
                if(short) {
                    return "stopped"
                }
                return "Run stopped before completion!";
            }
            if($scope.status) { // We have to make sure we've gotten data back from the server
                if($scope.status.state == 5) { // Failed
                    if(short) {
                        return "fail";
                    }
                    return "Run Failed!";
                }
                if($scope.status.state == 4) { // Success
                    if(short) {
                        return "success";
                    }
                    return "Run Succeeded!";
                }
                if($scope.status.failed > 0) { // Failing
                    if(short) {
                        return "failing";
                    }
                    return "Failures found, but run still progressing.";
                }
                // Running
                if(short) {
                    return "running";
                }
                return "Still Running.";
            }
        };

        $scope.stopRun = function() {
            $http.get($window.apiLinks.runs + "/" + $stateParams.id + $window.apiLinks.stop).success(function(data) {
                $scope.run.status = 1;
                $http.put($window.apiLinks.runs + "/" + $stateParams.id, {"status" : 1}).success(function(data) {
                    // No need to do anything
                }).error(function(data){
                    $scope.run.status = 0;
                });
            }).error(function(data) {
                $scope.alerts.push({
                    'message' : data.reason + "Error Code: " + data.status,
                    'type' : 'danger'
                });
                $scope.run.status = 0;
            });
        };

        $scope.gotoList = function() {
            $state.go("runs");
        };

        getRun();




    };

    // Return the array so the dependecy injections are accounted for
    return ["$scope", "$state", "$stateParams", "$http", "$window", runDetailsController];
});
