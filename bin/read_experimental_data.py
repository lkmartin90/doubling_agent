import numpy as np
import pandas as pd
import argparse
import os

# This code takes a data file takes a text file as an input, the text file should be the output from
# Qupath containing details of the location of each cell.
# eg '/Users/martinl/Documents/XDF/Code/doubling_agent/experimental_data/GBM_test_3/GBM_test3_cells.txt'

# input parameters
parser = argparse.ArgumentParser(description='Process input parameters for plotting')
parser.add_argument('--datafile', type=str, help='File containing the data', required=True)

args = parser.parse_args()


data_file = args.datafile
data_df = pd.read_csv(data_file, sep="\t")
coords_df = data_df[['Centroid X px', 'Centroid Y px', 'Cell: Max caliper', 'Cell: Min caliper']]

# Then want to combine data into the same form as the simulation data, i.e a column with the coordinates
# followed by the sate. The state will be the same for all the cells.

# At some point will need to normalise for average cell size in terms of pixels - will need metadata for this.
# or maybe not... use the "cell capliper" data? take average ?
diameter = (np.mean(coords_df['Cell: Max caliper'].values) + np.mean(coords_df['Cell: Min caliper'].values))/2
#print(diameter)

corrected_df = coords_df.copy()/diameter
corrected_df = corrected_df[['Centroid X px', 'Centroid Y px']]
corrected_df = corrected_df.rename(columns={'Centroid X px': "x", 'Centroid Y px': "y"})
corrected_df['state'] = 0
corrected_df['count'] = 1
corrected_df = corrected_df[['count', 'x', 'y', 'state']]

corrected_df.to_csv(os.path.dirname(os.path.abspath(args.datafile)) + '/repeat_0',
                    header=None, index=None, sep=',', mode='a')

# Then should be able to use the same code to analyse the data as in the simulation case, specifying --exp in the input

# But for this to work will also have to write the correct metadata so that it can be processed.
params = [0, 0, 0, 0]
dim = 2

with open(os.path.dirname(os.path.abspath(args.datafile)) + '/metadata.txt', 'a+') as f:
    f.write('parameters = ' + str(params) + '\nmotility parameters = '
            + 'None' + '\ntimestep = ' + '0.1' +
            '\nlength in days = ' + '1' + '\nnumber of repeats = ' + '1' +
            '\ndimensions = ' + str(dim))
