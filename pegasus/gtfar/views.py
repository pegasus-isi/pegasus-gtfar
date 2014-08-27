# Copyright 2007-2014 University Of Southern California
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__author__ = 'dcbriggs'

import os
from flask import render_template, request, redirect, url_for, jsonify, send_from_directory
from werkzeug import secure_filename


from pegasus.workflow import wrapper
from pegasus.gtfar import app
from pegasus.gtfar.models import isValidFile

@app.route("/")
def index():
    """
    Loads up the main page
    :return the template for the main page:
    """
    apiLinks = '{"runs" : "/api/runs", "upload" : "/api/upload", "status" : "/status", "outputs" : "/outputs", "logs" : "/logs", "download" : "/api/download"}'
    return render_template("mainView.html", apiLinks=apiLinks)


@app.route("/api/upload", methods=["POST"])
def upload():
    file = request.files["file"]
    if file and isValidFile(file.filename):
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config["UPLOAD_FOLDER"], filename))
        return redirect(url_for("index"))

@app.route("/api/runs/<int:id>/status", methods=["GET"])
def getStatus(id):
    workflow = wrapper.PegasusWorkflow(app.config["PEGASUS_HOME"], app.config["STORAGE_DIR"] + os.sep + str(id) + os.sep + "submit")
    status = workflow.monitor(['-l'])
    # we have to change the state to be a basic data type
    return jsonify(status)

@app.route("/api/download/<path:file>")
def download(file):
    return send_from_directory(app.config["STORAGE_DIR"], file)

@app.route("/api/runs/<int:id>/outputs", methods=["GET"])
def getOutputFiles(id):
    files = {"objects" : []}
    for filename in os.listdir(app.config["STORAGE_DIR"] + os.sep + str(id) + os.sep + "outputs"):
        files["objects"].append({"name" : filename})
    return jsonify(files)

@app.route("/api/runs/<int:id>/logs", methods=["GET"])
def getLogs(id):
    pass

@app.route("/tests")
def tests():
    """
    Loads up the testing environment
    :return: the template for the test page
    """
    return render_template("testRunner.html")

if __name__ == "__main__":
    app.run(debug=True)
