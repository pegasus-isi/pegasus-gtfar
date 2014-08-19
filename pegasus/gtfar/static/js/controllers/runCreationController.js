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

        $scope.strandRuleOptions = [{rule: "unstranded"},{rule: "same"},{rule: "opposite"}];
        // Defaults

        // We have to define this out of the run object becaues we need to link it to .rule in the end
        $scope.strandRule = $scope.strandRuleOptions[0];

        $scope.run = {
            genome : "testGenome", // TODO: add this in a bit $scope.genome,
            gtf : "testGTF", // TODO: Add this in a bit $scope.gtf,
            status : 0,
            userName : "genericUser",
            readLength : 100,
            mismatches : 3,
            trimUnmapped : false,
            mapFiltered : false
        };

        $scope.uploadFile = function() {
            /*
             * Upload the file
             * Angular does not natively work with a formData object so we have to override
             * the default behavior.  Refer to the link below
             * http://uncorkedstudios.com/blog/multipartformdata-file-upload-with-angularjs
             */
            var fileData = new FormData(document.getElementById("runForm"));
            $http.post($window.apiLinks.upload, fileData, {
                transformRequest : angular.identity,
                headers : {
                    "Content-Type" : undefined
                }
            }).success(function(data) {
                console.log("file uploaded");
            }).error(function(data) {
                console.error("upload failed!");
                console.error(data);
            });
        };

        $scope.cancel = function() {
            $state.go('runs');
        };

        $scope.addRun = function() {
            if($scope.form.$invalid) {
                return;
            }

            try {
                $scope.run["filename"] = document.getElementById("inputFile").files[0].name;
            }
            catch(error)
            {
                return; // An error message should show to the user so we just ned to return
            }
            // Add the run to the database, the runId will be used to create the folder
            $http.post($window.apiLinks.runs, $scope.run).success(function(data) {
                console.log(data);
                $state.go('runs');
             }).error(function(data) {
                console.error(data);
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