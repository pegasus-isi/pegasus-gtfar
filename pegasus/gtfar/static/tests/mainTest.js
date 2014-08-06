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
(function() {
    'use strict'
    require.config({

        //"baseURL": "js",
        "paths": {
            js : '/js',
            angular: "//cdnjs.cloudflare.com/ajax/libs/angular.js/1.2.18/angular",
            uiRouter : "//cdnjs.cloudflare.com/ajax/libs/angular-ui-router/0.2.10/angular-ui-router",
            bootstrapUI : "//cdnjs.cloudflare.com/ajax/libs/angular-ui-bootstrap/0.11.0/ui-bootstrap-tpls",
            angularGrid : '//cdnjs.cloudflare.com/ajax/libs/ng-grid/2.0.11/ng-grid.min',
            jquery : "//cdnjs.cloudflare.com/ajax/libs/jquery/2.1.1/jquery",
            angularMocks: "//cdnjs.cloudflare.com/ajax/libs/angular.js/1.2.18/angular-mocks",
            jasmine: "//cdnjs.cloudflare.com/ajax/libs/jasmine/2.0.0/jasmine",
            jasmineHtml: "//cdnjs.cloudflare.com/ajax/libs/jasmine/2.0.0/jasmine-html",
            jasmineBoot: "//cdnjs.cloudflare.com/ajax/libs/jasmine/2.0.0/boot"
        },
        /*
         * Some libraries do not use requirejs, which means that any libraries they require are not guaranteed to be loaded
         * beforehand, it also means they do not have a pointer to call in our files, to fix this add to the shim for any
         * library that gets added that is not AMD so that it can be properly used in this project.
         */
        shim: {
            angularMocks: {
                deps: ["angular"],
                exports : "mocks"
            },
            uiRouter : {
                deps : ["angular"]
                // No export because this just adds fields into angular
            },
            bootstrapUI : {
                deps : ["angular"]
                // No export because this just adds fields into angular
            },
            angularGrid : {
                deps : ["jquery","angular"]
                // No export because this just adds fields into angular
            },
            angular: {
                exports: "angular"
            },
            jasmine: {
                exports: "window.jasmineRequire"
            },
            jasmineHtml: {
                deps: ["jasmine"],
                exports: "window.jasmineRequire"
            },
            jasmineBoot: {
                deps: ["jasmine", "jasmineHtml"],
                exports: "window.jasmineRequire"
            },
            jquery : {
                exports : "$"
            }
        }
    });

    // Place all test files in this array
    var specs = [
        "./unit/controllers/runsControllerSpec"
    ];

    /*
     * Hack that calls window.onload after require is setup even though window.onload will have
     * already been called when the page first loaded.  This is so we can still have Jasmine working
     * in our require system.  This is only a problem with version 2.0 of Jasmine, I imagine (hope)
     * that this issue gets resolved in a later version.
     */
    require(["jasmineBoot"], function () {

        // Load the specs
        require(specs, function () {

            // Initialize the HTML Reporter and execute the environment (setup by `boot.js`)
            window.onload();
        });
    });
})();