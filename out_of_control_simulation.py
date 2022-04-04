import numpy as np
import pandas as pd
import ssDMEWMAC.functions as spc
import os
from os import path

# import the environment variable for the method to be simulated
method = os.environ['METH']
print(method)

# import the environment variable for the dimension to be simulated
dim = os.environ['DIME']

# import the environment variable for the iterations to be simulated
ite = int(os.environ['ITERATIONS'])

# import the environment variable for the smoothing parameter lambda to be simulated
lam = os.environ['LAMBDA']

# import the environment variable for the smoothing parameter lambda to be targeted
a = os.environ['ARL']

# import the environment variable for the magnitude of the shift to be induced
delta = os.environ['DELTA']

# import the environment variable for the wait time before perturbing
wait = os.environ['WAIT_TIME']

# build file paths to store the simulation results
f = "arl1_" + str(method) + "_" + str(dim) + "_" + str(lam).replace('.', '_', 1) + "_" + str(a) + "_" + str(delta).replace('.', '_', 1) + "_" + str(wait) + '.json'

pa = "/out/arl1_" + str(method) + "_" + str(dim) + "_" + str(lam).replace('.', '_', 1) + "_" + str(a) + "_" + str(delta).replace('.', '_', 1) + "_" + str(wait) + '.json'
print(pa)

ex = path.exists(pa)
# print(ex)

# Read in the completed h_lookup worksheet
h_lookup = pd.read_csv('/data/h_lookup.csv')

# setting up an escape clause to not run again if the output file already exists
if ex == False:
    # set the h_value to use in the out of control simulation
    h = h_lookup[(h_lookup['Method'] == method) & (h_lookup['Dimension'] == int(dim)) & (
                h_lookup['Lambda'] == float(lam)) & (h_lookup['ARL'] == int(a))]['h'].iloc[0]

    arl1 = spc.find_out_of_control_run_length(method=method, iterations=ite, dimension=dim,
                                             smoothing=float(lam), h_val=h, delta=float(delta), change_after=int(wait))

    out = {'Method': [method], 'Dimension': [dim], 'Lambda': [lam], 'ARL': [a], 'h': [h], 'Delta': [delta],  'Wait': [wait], 'ARL_observed': [arl1]}
    out = pd.DataFrame.from_dict(out)

    out.to_json(path_or_buf=pa)

if os.environ['LOCAL_OR_GCP'] == 'GCP':
    spc.upload_blob(bucket_name="fresch_ssdmewmac", source_file_name=pa, destination_blob_name=f)
