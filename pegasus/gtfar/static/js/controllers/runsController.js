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
    var runsController = function($scope, $window, $http, $state) {

        function getRuns() {
            $http.get($window.apiLinks.runs).success(function(data) {
                $scope.runs = data.objects;
            }).error(function(data) {
                console.error();
                console.error(data);
            });
        }
        getRuns();

        $scope.refresh = getRuns;

        // We have to define this with the CSS hardcoded in here because the system does not seem
        // to like inlining the -'s that bootstrap classes use
        $scope.getbgColor = function(exitCode) {
            return (exitCode == 0) ? "bg-success" : (exitCode == -1)  ? "bg-info" : "bg-danger";
        };

        $scope.createRun = function() {
            $state.go('createRun');
        }

        $scope.runsGrid = {
            data : "runs",
            enablePaging : true,
            showGroupPanel : true,
            showSelectionCheckbox : true,
            selectWithCheckboxOnly : true,
            showFilter : true,
            columnDefs : [
                {field : "name", displayName : "Name"},
                {field : "filename", displayName : "File"},
                {field : "created", displayName : "Created On"}
            ],
            rowTemplate : '<div style="height: 100%" ng-class="getbgColor(row.getProperty(\'exitcode\'))"><div ng-style="{ \'cursor\': row.cursor }" ng-repeat="col in renderedColumns" ng-class="col.colIndex()" class="ngCell ">' +
                '<div class="ngVerticalBar" ng-style="{height: rowHeight}" ng-class="{ ngVerticalBarVisible: !$last }"> </div>' +
                '<div ng-cell></div>' +
                '</div></div>'
        }

    };

    function getFullConstructor() {
        return ["$scope", "$window", "$http", "$state", runsController];
    }

    // For testing so we can mock the injected variables
    function getMinimalConstructor() {
        return runsController;
    }

    return {
        getFullConstructor : getFullConstructor,
        getMinimalConstructor : getMinimalConstructor
    };

});