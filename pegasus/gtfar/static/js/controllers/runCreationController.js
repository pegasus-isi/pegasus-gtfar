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
        var reader = new FileReader();

        $scope.strandRuleOptions = [{rule: "unstranded"},{rule: "same"},{rule: "opposite"}];
        // Defaults
        $scope.readLength = 100;
        $scope.mismatches = 3;
        $scope.trimUnmapped = false;
        $scope.mapFiltered = false;
        $scope.strandRule = $scope.strandRuleOptions[0];

        reader.onload = function(event) {
            // TODO: We could update a progress bar with this event
        };
        reader.onloadend = function(event) {
            $scope.fileData = event.target.result;
            //$scope.inputFile = event.target.result;
        };

        $scope.uploadFile = function() {
            // Angular does not support file upload so we have to go by the element name
            file = document.getElementById("inputFile").files[0];
            // readAsArrayBuffer would never give any contents so I'm using readAsBinaryString
            reader.readAsBinaryString(file);
        };
        $scope.addRun = function() {

            // TODO: Figure out where to add validation checks for all the inputs

            var runData = {
                "name" : $scope.runName,
                "userName" : "genericUser",
                "status" : 0,  // For now we just set the status to running
                "filename" : document.getElementById("inputFile").files[0].name,
                "readLength" : $scope.readLength,
                "genome" : "testGenome", // TODO: add this in a bit $scope.genome,
                "gtf" : "testGTF", // TODO: Add this in a bit $scope.gtf,
                "mismatches" : $scope.mismatches,
                "trimUnmapped" : $scope.trimUnmapped,
                "mapFiltered" : $scope.mapFiltered,
                "strandRule" : $scope.strandRule.rule,
                "email" : $scope.email
            };
            //console.log($scope.fileData);
            var fileData = new FormData(document.getElementById("runForm"));

            // Add the run to the database, the runId will be used to create the folder
            $http.post($window.apiLinks.runs, runData).success(function(data) {
             console.log(data);
             $state.go('runs');
             }).error(function(data) {
             console.error(data);
             });

            /*
             * Upload the file
             * Angular does not natively work with a formData object so we have to override
             * the default behavior.  Refer to the link below
             * http://uncorkedstudios.com/blog/multipartformdata-file-upload-with-angularjs
             */
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
        getMinimalConstructor : getMinimalConstructor
    };

});