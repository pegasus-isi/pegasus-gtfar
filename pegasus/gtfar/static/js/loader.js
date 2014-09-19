/**
 * Copyright 2007-2014 University Of Southern California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Created by dcbriggs on 7/9/14.
 */


/*
 * Uses require js to load the libraries as modules.  This will make it so that our Javascript has to explicitly
 * require other modules in order to use them and cannot just access any function in any file.
 */
requirejs.config({

    //"baseURL" : "js",
    "paths" : {
        js : "/js",
        angular : "//cdnjs.cloudflare.com/ajax/libs/angular.js/1.2.18/angular",
        uiRouter : "//cdnjs.cloudflare.com/ajax/libs/angular-ui-router/0.2.10/angular-ui-router",
        bootstrapUI : "//cdnjs.cloudflare.com/ajax/libs/angular-ui-bootstrap/0.11.0/ui-bootstrap-tpls",
        angularGrid : '//cdnjs.cloudflare.com/ajax/libs/ng-grid/2.0.11/ng-grid.min',
        jquery : "//cdnjs.cloudflare.com/ajax/libs/jquery/2.1.1/jquery",
        "moment" : "//cdnjs.cloudflare.com/ajax/libs/moment.js/2.8.3/moment.min"
    },
    /*
     * Some libraries do not use requirejs, which means that any libraries they require are not guaranteed to be loaded
     * beforehand, it also means they do not have a pointer to call in our files, to fix this add to the shim for any
     * library that gets added that is not AMD so that it can be properly used in this project.
     */
    shim : {
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
        angular : {
            exports : "angular"
        },
        jquery : {
            exports : "$"
        }
    }
});

/*
 * Load the main module
 */

requirejs(["angular", "app"], function(angular, app) {
    angular.element(document).ready(function () {
        angular.bootstrap(document, [app.getName()]);
    });
});
