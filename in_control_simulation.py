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

# set the list of dimensions to simulate
dims = list((2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 50))
lambdas = list((0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95))
target_arls = list((100, 200, 500, 1000, 2000))

# Read in the completed h_init worksheet
h_init_lookup = pd.read_csv('/data/h_init_lookup.csv')


# Create nested loops to go over these combos and pipe the results to json
for l in lambdas:
    for a in target_arls:
        f = "h_"+str(method)+"_"+str(dim)+"_"+str(l).replace('.', '_', 1)+"_"+str(a)+'.json'

        pa = "/out/h_"+str(method)+"_"+str(dim)+"_"+str(l).replace('.', '_', 1)+"_"+str(a)+'.json'
        print(pa)

        ex = path.exists(pa)
        #print(ex)


        if ex == False:
            h_init = h_init_lookup[(h_init_lookup['Method'] == method) & (h_init_lookup['Dimension'] == int(dim)) & (h_init_lookup['Lambda'] == l) & (h_init_lookup['ARL'] == a)]['h'].iloc[0]
            #print(h_init)
            h, arl0 = spc.find_h(method=method, dimension=int(dim), smoothing=l, target_arl=a, h_init=h_init, tolerance=.001, iterations=ite)

            out = {'Method' : [method] , 'Dimension' : [dim], 'Lambda' : [l], 'ARL' : [a], 'h' : [h], 'ARL_observed' : [arl0]}
            out = pd.DataFrame.from_dict(out)

            out.to_json(path_or_buf=pa)
                
        if os.environ['LOCAL_OR_GCP'] == 'GCP':
            spc.upload_blob(bucket_name="fresch_ssdmewmac", source_file_name=pa, destination_blob_name=f)
