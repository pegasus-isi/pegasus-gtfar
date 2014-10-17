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
        var runsController = function ($scope, $window, $http, $state) {

            // Adjust this to change the number of pages you have to request from the server
            var resultsPerPage = 100;

            $scope.alerts = [];
            $scope.getAlertIcon = function (status) {
                return (status == 'success') ? "text-success fa fa-check-circle" : "text-danger fa fa-exclamation-triangle";
            };

            $scope.closeAlert = function (index) {
                $scope.alerts.splice(index, 1);
            };

            function getRuns() {
                $http.get($window.apiLinks.runs, {
                    params: {
                        results_per_page: resultsPerPage
                    }
                }).success(function (data) {
                    $scope.runs = data.objects;
                    // Get all the results after the first page
                    for (var currentPage = 2; currentPage <= data.total_pages; currentPage++) {
                        $http.get($window.apiLinks.runs, {
                            params: {
                                results_per_page: resultsPerPage,
                                page: currentPage
                            }
                        }).success(function (data) {
                            $scope.runs = $scope.runs.concat(data.objects);
                        }).error(function (data) {
                            $scope.alerts.push({
                                "type": "danger",
                                "message": "Unable to get runs, please contact the Pegasus team if this issue persists"
                            });
                        })
                    }
                }).error(function (data) {
                    console.error(data);
                });
            }

            getRuns();

            $scope.refresh = getRuns;

            // We have to define this with the CSS hardcoded in here because the system does not seem
            // to like inlining the -'s that bootstrap classes use
            $scope.getbgColor = function (exitCode) {
                return (exitCode == 0) ? "bg-success" : (exitCode == -1) ? "bg-info" : "bg-danger";
            };

            $scope.getExitCodeIcon = function (status) {
                switch (status) {
                    case -1:
                        return "text-primary fa fa-check-circle";
                    case 0:
                        return "text-success fa fa-check-circle";
                    case 256:
                        return "text-warning fa fa-warning";
                    default:
                        return "text-danger fa fa-exclamation-triangle";
                }
            };

            $scope.createRun = function () {
                $state.go('createRun');
            };

            $scope.deleteSelected = function () {

                var selected = $scope.runsGrid.selectedItems;
                var deleteTotal = $scope.runsGrid.selectedItems.length;
                var runName = selected[0].name;
                if (deleteTotal > 1) {
                    runName = "these " + deleteTotal + " runs"
                }
                var doubleCheck = confirm("Are you sure you want to delete " + runName + "?  This will delete all output files as well.");
                if (!doubleCheck) {
                    return;
                }

                var deleteCount = 0;
                for (var i = 0; i < deleteTotal; i++) {
                    $http.delete($window.apiLinks.runs + "/" + selected[i].name)
                        .success(function () {
                            // No specifications aside from the then clause
                        }).error(function () {
                            // TODO: display an error popup
                            deleteCount++;

                        })
                        .then(function () {
                            deleteCount++;
                            if (deleteCount == deleteTotal) {
                                $scope.runsGrid.selectedItems.length = 0;
                                getRuns();
                            }
                        })
                }
            };

            $scope.getFormattedDate = function (dateString) {
                var createdOn = moment.utc(dateString).local();
                var s = createdOn.format('ddd MMM DD, YYYY hh:mm:ss A');
                s += ' (' + createdOn.fromNow() + ')';
                return s;
            };


            $scope.runsGrid = {
                data: "runs",
                enablePaging: true,
                showGroupPanel: true,
                showSelectionCheckbox: true,
                selectWithCheckboxOnly: true,
                showFilter: true,
                columnDefs: [
                    {field: "name", displayName: "Name", cellTemplate: "<a class='gridCenter' title='{{row.getProperty(col.field)}}' ui-sref='runDetails({name : row.getProperty(\"name\")})'><div class='ngCellText'>{{row.getProperty(col.field)}}</div></a>" },
                    {field: "status", displayName: "Status", width: 60, cellTemplate: "<i class='ngCellText' ng-class='getExitCodeIcon(row.getProperty(col.field))'></i>" },
                    {field: "readLength", displayName: "Read Length", width: '10%', cellTemplate: "<div class='ngCellText'>{{row.getProperty(col.field)}}</div>" },
                    {field: "mismatches", displayName: "Mismatches", width: '10%', cellTemplate: "<div class='ngCellText'>{{row.getProperty(col.field)}}</div>" },
                    {field: "filename", displayName: "File", cellTemplate: "<a class='gridCenter' ng-href='{{getDownloadLink(row.getProperty(\"name\"), row.getProperty(col.field))}}' title='{{row.getProperty(col.field)}}'><div class='ngCellText'><i class='fa fa-file-text-o'></i> {{ row.getProperty(col.field) }}</div></a>" },
                    {field: "created", displayName: "Started On", width: '30%', cellTemplate: "<div class='ngCellText'>{{getFormattedDate(row.getProperty(col.field))}}</div>" }
                ],
                sortInfo: {fields: ['created'], directions: ['desc']},
                selectedItems: []
            };


            $scope.getDownloadLink = function (id, filename) {
                return $window.apiLinks.runs + "/" + id + "/input/" + filename;
            };

        };

        function getFullConstructor() {
            return ["$scope", "$window", "$http", "$state", runsController];
        }

        // For testing so we can mock the injected variables
        function getMinimalConstructor() {
            return runsController;
        }

        return {
            getFullConstructor: getFullConstructor,
            getMinimalConstructor: getMinimalConstructor
        };

    });
