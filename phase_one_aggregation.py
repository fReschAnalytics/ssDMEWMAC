"""
The following script is a simple script to scan the output directory for phase I output files.
It then concatenates the results into a single lookup file for use in phase II.
"""

import os
import pandas as pd

h_base = {'Method': ["(ss)(D)MEWM(A)(C)"], 'Dimension': [0], 'Lambda': [0], 'ARL': [0], 'h': [0], 'ARL_observed': [0]}
h_base = pd.DataFrame.from_dict(h_base)

for filename in os.listdir("data/"):
    if filename.endswith(".json") and filename.startswith("h_"):
        # read the json file
        print("Accessing file: %s" % filename)
        h = pd.read_json("data/"+filename)

        # concatenate with existing h_base dataframe
        h_base = pd.concat((h_base, h))

        continue
    else:
        continue

h_base = h_base.sort_values(by=["Method", "Dimension", "Lambda", "ARL"])
h_base.to_csv("data/h_lookup.csv", header=True, index=None)