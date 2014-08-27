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
define([],
function() {
    "use strict"

    var runDetailsController = function($scope, $state, $stateParams, $http, $window) {


        function getRun() {
            $http.get($window.apiLinks.runs + "/" + $stateParams.id).success(function(data) {
                $scope.run = data;
            }).error(function(data) {
                console.error(data);
            });
        }


        function getStatus() {
            $http.get($window.apiLinks.runs + "/" + $stateParams.id + $window.apiLinks.status).success(function(data) {
                if(statusChanged(data)) {
                    $scope.status = data;
                    $scope.statusText =
                    getOutputFiles();
                }
                // We don't want the calls to continue when we leave this state
                if($state.includes('runDetails')) {
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
            console.log("WE NEED TO CHECK THE OUTPUT FILES!");
            $http.get($window.apiLinks.runs + "/" + $stateParams.id + $window.apiLinks.outputs).success(function(data) {
                $scope.outputFiles = data.objects;
                // TODO: need to present links to download the output files
            }).error(function(data) {
                console.error(data);
            });
        }

        $scope.getStatusColor = function(prefix) {
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

        };

        /*function downloadFile(section, filename) {
            $http.get($window.apiLinks.download + "/" + $stateParams.id + "/" + section + "/" + filename).success(function(data) {

            })
        }*/

        getRun();
        getStatus();

        $scope.inputDownload = $window.apiLinks.download + "/" + $stateParams.id + "/input/";
        $scope.outputDownload = $window.apiLinks.download + "/" + $stateParams.id + "/outputs/";


    };

    // Return the array so the dependecy injections are accounted for
    return ["$scope", "$state", "$stateParams", "$http", "$window", runDetailsController];
});