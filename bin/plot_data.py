from doubling_agent.image_analysis_functions import *
from doubling_agent.common_functions import *
import argparse


class StoreAsArray(argparse._StoreAction):
    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)

# sate min lacuarity box size
min_lac_box = 4

# Parse all of the arguments given by the user

parser = argparse.ArgumentParser(description='Process input parameters for plotting')
parser.add_argument('--folder', type=str, help='Folder containing the data you wish to plot.', required=True)
parser.add_argument('--repeat', action=StoreAsArray, type=int, nargs='+', help="Which repeats to plot"
                                            "can take multiple values.",
                    default=np.array([0]))
parser.add_argument('--subset', action=StoreAsArray, type=float, nargs='+', help="Do you wish to look at a subset of "
                                            "the data? If so enter the min and max value, subset must be square.")
parser.add_argument('--exp', help='Is this experimental data include flag if yes', action = 'store_true')
parser.add_argument('--mapping', help='Include flag to map state 0 onto state 1 for analysis', action = 'store_true')
args = parser.parse_args()

# open the file with the metadata describing the length of the simulation, number of dimensions, etc.

with open(args.folder + '/metadata.txt') as f:
    meta_data = f.readlines()

if int(meta_data[5].strip().split()[2]) == 2:
    switch_3d = False
else:
    switch_3d = True

# determine if cells are motile

if meta_data[1].split()[3] == 'None' or meta_data[1].split()[3] == '[0.0]' :
    from doubling_agent.basic_functions import *
    mot = False
else:
    from doubling_agent.motility_functions import *
    mot = True

# find the number of repeats
repeats = int(meta_data[4].strip().split()[4])

# mapping

mapping = args.mapping

# set subrange of data to look at from input parameters
subrange = args.subset
if subrange is not None:
    if len(subrange) != 2:
        print('Invalid subset length')
        exit()
time_step = float(meta_data[2].strip().split()[2])

for r in args.repeat:
    if args.exp:
        print('experimental data')
    else:
        data = read_from_file(args.folder + '/repeat_' + str(r), switch_3d)

    if subrange is None:
        folder_name = args.folder + '/plots_repeat_' + str(r)
    else:
        folder_name = args.folder + '/plots_repeat_' + str(r) + '_subset_' + str(subrange[0]) + '_' + str(subrange[1])

    # don't waste computational effort by doing the same thing twice if the plots already exist
    if os.path.isdir(folder_name):
        print('Plots already exist for repeat ', str(r))
    else:
        ensure_dir(folder_name)
        # data comes back as a list of touples. count, x,y, state

        if args.exp:
            data_df = pd.read_csv(args.folder + '/repeat_' + str(r), sep = ',', names=['count', 'x', 'y', 'state'])
        else:
            if mot and switch_3d:
                data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'z', 'state', 'mot'])
            elif mot:
                data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'state', 'mot'])
            elif switch_3d:
                data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'z', 'state'])
            else:
                data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'state'])

        # optional mapping of state 0 onto state 1 for analysis
        if mapping:
            mapping = {0: 1}
            data_df = data_df.replace({'state': mapping})

        counts = data_df.drop_duplicates(subset=['count'], keep='last')
        # want the option here to only look at a subset of the data, user specify subset range, then can do analysis
        # on the full tumour followed by the analysis of a densly populated central slice, as that will be easier
        # to find/ study in real data.
        #print(data_df)
        if subrange is not None and switch_3d:
            data_df = data_df.loc[(data_df['x'] < subrange[1]) & (data_df['x'] > subrange[0])]
            data_df = data_df.loc[(data_df['y'] < subrange[1]) & (data_df['y'] > subrange[0])]
            data_df = data_df.loc[(data_df['z'] < subrange[1]) & (data_df['z'] > subrange[0])]
        elif subrange is not None:
            data_df = data_df.loc[(data_df['x'] < subrange[1]) & (data_df['x'] > subrange[0])]
            data_df = data_df.loc[(data_df['y'] < subrange[1]) & (data_df['y'] > subrange[0])]

        # here each "value" will be a number of timesteps after which data was recorded.
        for value in counts['count'].values:

            # If the simulation is 3D and the last repeat to be analysed and the last time point at which
            # data was recorded, only then look at the 2D slices, otherwise will take too long.
            if switch_3d and r == np.max(args.repeat) and value == np.max(counts['count'].values):
                print('Plotting...')
                plot_cells_3d(data_df, value, folder_name, r, time_step)
                print('Plotting 2D slices of repeat ' + str(r) + ' and day ' + str(int(value * time_step)))
                plot_2d_slice(folder_name, data_df, value, time_step, r)
                print('2D slice analysis...')
                density_analysis_2d_slice(folder_name, data_df, value, 8, time_step, r, subrange, min_lac_box)
            elif switch_3d:
                print('Plotting...')
                plot_cells_3d(data_df, value, folder_name, r, time_step)
            else:
                plot_cells(data_df, value, folder_name, r, time_step)

            if len(data_df.loc[data_df['count'] == value]) > 30:
                print('Analysis of repeat ' + str(r) + ' and day ' + str(int(value * time_step)))
                density_analysis(folder_name, data_df, value, 8, switch_3d, time_step, r)

                distance_from_centre(folder_name, data_df, value, switch_3d, time_step, r)
                # At the momennt FFT just projects onto 2d is you give it a 3d dataset
                fft_analysis(data_df, value, folder_name, r, time_step, 'full')

            if subrange is not None:
                # At the moment if you send a 3d dataset to the lacunarity function it simply projects it onto 2d
                # it doesn't do 2d slices.
                print('Lacunarity for repeat ' + str(r) + ' and day ' + str(int(value * time_step)))
                lacunarity(folder_name, data_df, value, time_step, r, subrange, 'full', min_lac_box)

max_r = np.max(args.repeat)
min_r = np.min(args.repeat)

#fit_analysis(args.folder, min_r, max_r, time_step)

# where a 3d tumour has been split into 2d segments, want to compare data across the tumour.
print('Section analysis...')
cell_dict = {0: "Stem cell", 1: "Progenitor cell", 2: "Differentiated cell", 3: "Quiescent cell"}
section_analysis(args.folder, 'sect_*_cell_numbers.txt', 'number of cells', cell_dict, min_r, max_r)
fft_dict = {1: 'Mean freq', 0: 'Std of freq'}
section_analysis_fft(args.folder, 'fft_data_day_*_sect_*.txt', 'FFT frequencies', cell_dict, fft_dict, min_r, max_r)

# average the lacunarity across the 3D slices if the simulation was in 3D
if args.subset is not None and switch_3d is True:
    lac_analysis(args.folder, 'lac_tot_day_*_sect_*_state_0.txt',
                'Stem cell average lacunarity', min_r, max_r, min_lac_box)
    lac_analysis(args.folder, 'lac_tot_day_*_sect_*_state_1.txt',
                 'Progenitor cell average lacunarity', min_r, max_r, min_lac_box)
    lac_analysis(args.folder, 'lac_tot_day_*_sect_*_state_2.txt',
                 'Differentiated cell average lacunarity', min_r, max_r, min_lac_box)
    lac_analysis(args.folder, 'lac_tot_day_*_sect_*_state_3.txt',
                 'Quiescent cell average lacunarity', min_r, max_r, min_lac_box)
