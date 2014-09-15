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
 * This module handles hooking up all the controllers to the proper route within the app.
 *
 * Created by dcbriggs on 7/14/14.
 */
define(["./controllers/runsController", "./controllers/runDetailsController", "./controllers/runCreationController"],
function(runsController, runDetailsController, runCreationController){

    /*var routes = function($routeProvider) {
        $routeProvider.when('/workflows', {
            templateUrl : "partials/runs.html",
            controller : "runsController"

        }).when('/workflows/:Id', {
            templateUrl : "partials/runDetails.html",
            controller : "runDetailsController"
        }).when('/createRun', {
         templateUrl : "partials/runCreator.html",
         controller : "runCreationController"
         }).otherwise({
            redirectTo : '/workflows'
        });
    };*/

    var routes = function($stateProvider, $urlRouterProvider) {
        //
        // For any unmatched url, redirect to /state1
        $urlRouterProvider.otherwise("/runs");
        //
        // Now set up the states
        $stateProvider
            .state('runs', {
                url: "/runs",
                templateUrl: "partials/runs.html",
                controller: "runsController"

            })
            .state('runDetails', {
                url: "/runs/details/{id:[^/]*}",
                templateUrl: "partials/runDetails.html",
                controller: "runDetailsController"
            })
            .state('createRun', {
                url: "/createRun",
                templateUrl: "partials/runCreator.html",
                controller : "runCreationController"
            }).state('help', {
                url: "/help",
                templateUrl : "partials/help.html"
            });
    };

    /*
     * This will be used in the .config of the app.
     * Return an array with the variables we want
     * to inject so the caller does not need to worry about that.
     */
    function getFullConstructor() {
        return ["$stateProvider", "$urlRouterProvider", routes];
    }

    // For testing so that we can mock the variables you normally inject
    function getMinimalConstructor() {
        return routes;
    }

    return {
        getFullConstructor : getFullConstructor,
        getName : getMinimalConstructor
    }




});