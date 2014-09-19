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
 * Instantiates the main app and attaches all services and controllers to it.
 *
 * Created by dcbriggs on 7/14/14.
 */
// Keep angular files at the end because they do not export anything
define(["angular", "js/controllers/runsController", "js/controllers/runDetailsController", "js/controllers/runCreationController",
        "js/router", "uiRouter", "angularGrid", "bootstrapUI", "moment"],
function(angular, runsController, runDetailsController, runCreationController,
         router) {
    'use strict'

    var appName = "gtfarApp";
    var app = angular.module(appName, ["ui.router", "ngGrid", "ui.bootstrap"]);

    // For form validation
    var INTEGER_REGEX = /^\-?\d+$/;
    app.directive('integer', function() {
        return {
            require: 'ngModel',
            link: function(scope, elm, attrs, ctrl) {
                ctrl.$parsers.unshift(function(viewValue) {
                    if (INTEGER_REGEX.test(viewValue)) {
                        // it is valid
                        ctrl.$setValidity('integer', true);
                        return viewValue;
                    } else {
                        // it is invalid, return undefined (no model update)
                        ctrl.$setValidity('integer', false);
                        return undefined;
                    }
                });
            }
        };
    });

    var ALPHANUM_REGEX = /^\w+$/i;
    app.directive('alphanumeric', function() {
        return {
            require: 'ngModel',
            link: function(scope, elm, attrs, ctrl) {
                ctrl.$parsers.unshift(function(viewValue) {
                    if (ALPHANUM_REGEX.test(viewValue)) {
                        // it is valid
                        ctrl.$setValidity('alphanumeric', true);
                        return viewValue;
                    } else {
                        // it is invalid, return undefined (no model update)
                        ctrl.$setValidity('alphanumeric', false);
                        return undefined;
                    }
                });
            }
        };
    });

    app.config(router.getFullConstructor());

    /*app.run(["$window", function($window) {
        console.log($window.apiLinks.runs);
        $window.apiLinksJSON = JSON.parse($window.apiLinks);
    }]);*/

    app.controller("runsController", runsController.getFullConstructor());
    app.controller("runDetailsController", runDetailsController);
    app.controller("runCreationController", runCreationController.getFullConstructor());


    function getName() {
        return appName;
    }

    return {
        app : app,
        getName : getName
    };

});
