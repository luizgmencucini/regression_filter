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
    <link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.3.6/css/buttons.dataTables.min.css">
    <script type="text/javascript" src="https://code.jquery.com/jquery-3.5.1.js"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/buttons/2.3.6/js/dataTables.buttons.min.js"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.html5.min.js"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/buttons/2.3.6/js/buttons.print.min.js"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/select/1.6.2/js/dataTables.select.min.js"></script>

    <script>
        var dataSet = REPLACE;
        console.log(dataSet);
        $(document).ready(function() {
        $('#resTable').DataTable( {
            dom: 'Bfrtip',
            buttons: [
                {
                    extend: 'csv',
                    text: 'Export Selected',
                    exportOptions: {
                        modifier: {
                            selected: true
                        }
                    }
                }
            ],
            data : dataSet,
            // add column definitions to map your json to the table
            "columns": [
            {title: "row ID"},
            {title: "row m/z"},
            {title: "row retention time"},
            {title: "XIC"}
            ],
            select: true
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
        fig = "./images/%s.png" % report_print.loc[i, 'row ID']
        report_print.loc[i, 'XIC'] = '<img src="%s" width="190" height="190">' % fig

    html_local = re.sub('REPLACE', json.dumps(report_print.apply(lambda a: a.tolist(), axis=1).tolist()), html)
    with open(out_file, 'w+') as f:
        f.write(html_local)

#Exclusion list -> It's a .csv file generated by html report
def create_SPL(exclusion_list, peaks_info):
    df = pd.read_csv(exclusion_list)
    lista = df.iloc[:, 0].tolist()
    
    dfx = peaks_info.reset_index()
    dfx.drop(['index'], axis=1, inplace=True)
    dfx.set_index('row ID', inplace = True)
    dfx.index = dfx.index.astype(int)
    dfx.drop(lista, 0, inplace = True)

    name = []
    time = []
    time_tol = []
    mass_begin = []
    mass_end = []

    for item in dfx.index:
        name.append('')
        time.append(dfx['row retention time'][item]*60)
        time_tol.append(30)
        mass_begin.append(dfx['row m/z'][item])
        mass_end.append('')
    
    col_names = ["Compound Name", "Time", "Time Tolerance", "Mass Begin", "Mass End"]
    SPL = pd.DataFrame(zip(name,time,time_tol,mass_begin,mass_end), columns = col_names)
    SPL.to_csv("SPL_list.csv")