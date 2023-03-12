import json
import re
import os

import pandas as pd
import numpy as np
from IPython.display import SVG
import matplotlib.pyplot as plt

html = '''<!DOCTYPE html>
<html>
  <head>
    <meta http-equiv="content-type" content="text/html; charset=UTF-8">
    <title>NPMINE</title>
    <script src="https://code.jquery.com/jquery-3.3.1.js" integrity="sha256-2Kok7MbOyxpgUVvAk/HJ2jigOSYS2auK4Pfzbm7uH60=" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.6.1/css/all.css" integrity="sha384-gfdkjb5BdAXd+lj+gudLWI+BXq4IuLW5IT+brZEZsLFm++aCMlF1V92rMkPaX4PP" crossorigin="anonymous">
</head>

  <body>
    <!-- A grey horizontal navbar that becomes vertical on small screens -->
    <nav class="navbar navbar-expand-sm bg-light navbar-light">
        <a class="navbar-brand" href="">
        </a>

      <!-- Links-->
      <ul class="navbar-nav">
        <li class="nav-item">
          <a class="nav-link" href="http://ccbl.fcfrp.usp.br">CCBL</a>
        </li>
      </ul>

    </nav>
    <div id="content">

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.6.1/css/all.css" integrity="sha384-gfdkjb5BdAXd+lj+gudLWI+BXq4IuLW5IT+brZEZsLFm++aCMlF1V92rMkPaX4PP" crossorigin="anonymous">

    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>

    <script src="https://cdn.datatables.net/1.10.4/js/jquery.dataTables.min.js"></script>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.10.4/css/jquery.dataTables.min.css">
    <script>
        var dataSet = REPLACE;
        console.log(dataSet);
        $(document).ready(function() {
        $('#resTable').DataTable( {
            data : dataSet,
            // add column definitions to map your json to the table
            "columns": [
            {title: "row ID"},
            {title: "row m/z"},
            {title: "row retention time"},
            {title: "XIC"}
            ]
        } );
        });
    </script>
    <div class="m-5">
        <table id="resTable" class="table table-striped" style="width:100%" >
		<thead>
        <tr>
		<th>row ID</th>
		<th>row m/z</th>
		<th>row retention time</th>
		<th>XIC</th>
        </tr>
        </thead>
        </table>
    </div>

    </div>
  </body>
</html>'''


def create_report(report_print, out_file='regfilter_report.html'):
    """Creates an html report from NPMINE's results
    Parameters
    ----------
    report_print: pd.DataFrame
        DataFrame containing columns 'doi', 'pubchem', 'ExactMolWt', 'smiles','source'.
    dois: str or list
        Directory of doi link files or link list.
    out_file: str
        HTML output file name.
    useSVG: bool
        If svg format should be used. Default is png.
    Returns
        Report html file.
    -------
    """

    for i in report_print.index:
        fig = "./figs/%s.png" % report_print.loc[i, 'row ID']
        report_print.loc[i, 'XIC'] = '<img src="%s" width="120" height="120">' % fig

    html_local = re.sub('REPLACE', json.dumps(report_print.apply(lambda a: a.tolist(), axis=1).tolist()), html)
    with open(out_file, 'w+') as f:
        f.write(html_local)

