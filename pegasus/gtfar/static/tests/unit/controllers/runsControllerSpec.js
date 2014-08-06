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
 * Created by dcbriggs on 7/23/14.
 */

define(["angular", "angularMocks", "js/controllers/runsController.js"],
    function(angular, mocks, runsController) {
    'use strict'
    describe('runsController Unit Tests', function() {
        var scope, httpBackend, objects, stateProvider;

        objects = [
            {rid : 1, filename : "myFile.txt", exitCode : 0},
            {rid : 2, filename : "test.zip", exitCode : 0},
            {rid : 3, filename : "test2.tar", exitCode : 0}
        ];

        beforeEach(function() {
            var mockApp = angular.module("gtfarApp", ['stateMock']);
            mockApp.controller("runsController", runsController.getMinimalConstructor());
        });

        beforeEach(inject(function ($rootScope, $controller, $httpBackend, $http, $state) {
            scope = $rootScope.$new();
            window.apiLinks = {"runs" : "/mock/api/runs"};
            httpBackend = $httpBackend;
            httpBackend.when("GET", "/mock/api/runs").respond({
                "objects"  : objects
            });
            $controller(runsController.getMinimalConstructor(), {
                $scope : scope,
                $window : window,
                $http : $http,
                $state : $state
            });
            httpBackend.flush();
        }));

        it("Test that controller automatically gets the runs", function() {
            // Because of the way our controller files are setup, we have to specify the
            // last element in the array instead of just
            expect(scope.runs.length).toBe(3);
        });

        it("Test that controller can refresh the runs", function() {
            expect(scope.refresh).toBeTruthy();
        });

        it("Test that when refresh is called but the data is the same that nothing changes", function() {

            scope.refresh();
            httpBackend.flush();
            expect(scope.runs.length).toBe(3);
        });

        it("Test that when refresh is called with new data that the changes are noticed", function () {
            objects.shift();
            objects.unshift({rid : 4, filename : "hello.tar.bz", exitCode : 1});
            scope.refresh();
            httpBackend.flush();
            expect(scope.runs.length).toBe(3);
            expect(scope.runs[0].rid).toBe(4);
            expect(scope.runs[1].rid).toBe(2);
            expect(scope.runs[2].rid).toBe(3);
        });


    });
});