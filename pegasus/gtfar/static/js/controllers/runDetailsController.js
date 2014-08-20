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

    var runDetailsController = function($scope, $state, $stateParams, $window) {


        function getRun() {
            $httop.get($window.apiLInks.runs + "")
        }
        function getStatus() {
            $http.get($window.apiLinks.runs + "/" + $stateParams.id + $window.apiLinks.status).success(function(data) {
                if(statusChanged(data.object)) {
                    $scope.status = data.object;
                    getOutputs();
                }
                setTimeout(getStatus, 2000);
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

        function getOutputs() {
            $http.get($window.apiLinks.runs + "/" + $stateParams.id + $window.apiLinks.outputs).success(function(data) {
                $scope.outputs = data.objects;
                // TODO: need to present links to download the output files
            }).error(function(data) {
                console.error(data);
            });
        }

        //getStatus(); uncomment when we have the backend hooked up properly


    };

    // Return the array so the dependecy injections are accounted for
    return ["$scope", "$state", "$stateParams", "$window", runDetailsController];
});