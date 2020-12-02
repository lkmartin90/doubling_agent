from doubling_agent.image_steve_functions import *
from doubling_agent.steve_functions import *
import argparse


class StoreAsArray(argparse._StoreAction):
    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)


# sate minimum lacunarity box size
min_lac_box = 4

# deal with the user specified parameters
parser = argparse.ArgumentParser(description='Process input parameters for plotting')
parser.add_argument('--folder', type=str, help='Folder containing the data you wish to plot.', required=True)
parser.add_argument('--repeat', action=StoreAsArray, type=int, nargs='+', help="Which repeats to plot"
                                            "can take multiple values.",
                    default=np.array([0]))
parser.add_argument('--subset', action=StoreAsArray, type=float, nargs='+', help="Do you wish to look at a subset of "
                                            "the data? If so enter the min and max value, subset must be square.")
args = parser.parse_args()

with open(args.folder + '/metadata.txt') as f:
    meta_data = f.readlines()

if int(meta_data[5].strip().split()[2]) == 2:
    switch_3d = False
else:
    switch_3d = True

repeats = int(meta_data[4].strip().split()[4])

time_step = float(meta_data[2].strip().split()[2])

subrange = args.subset
if subrange is not None:
    if len(subrange) != 2:
        print('Invalid subset length')
        exit()

print('Analysing data as continuous...')
for r in args.repeat:
    data = read_from_file(args.folder + '/repeat_' + str(r), switch_3d)

    if subrange is None:
        folder_name_cont = args.folder + '/plots_repeat_cont' + str(r)
    else:
        folder_name_cont = args.folder + '/plots_repeat_cont' + str(r) + '_subset_' \
                           + str(subrange[0]) + '_' + str(subrange[1])

    if os.path.isdir(folder_name_cont):
        print('Continuous plots already exist for repeat ', str(r))
    else:
        ensure_dir(folder_name_cont)
        # data comes back as a list of touples. count, x,y, state

        if switch_3d:
            data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'z', 'state'])
        else:
            data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'state'])

        counts = data_df.drop_duplicates(subset=['count'], keep='last')

        for value in counts['count'].values:

            if switch_3d and r == np.max(args.repeat) and value == np.max(counts['count'].values):
                print('Plotting...')
                plot_cells_3d_cont(data_df, value, folder_name_cont, r, time_step)
                plot_2d_slice_cont(folder_name_cont, data_df, value, time_step, r)
            elif switch_3d:
                print('Plotting...')
                plot_cells_3d_cont(data_df, value, folder_name_cont, r, time_step)
            else:
                plot_cells_cont(data_df, value, folder_name_cont, r, time_step)

# Then want to preform the analysis by splitting the sliding scale into 3: Dormant, primed and dividing.
# Should be able to do this by mapping the values from the recorded data onto one of three states, 0, 1, and 2
# Then can use all the plotting tools from the previous analysis.
# Allow the map to be user defined so that it is simple to change between the different states.

print('Analysing data as discrete...')
for r in args.repeat:
    data = read_from_file(args.folder + '/repeat_' + str(r), switch_3d)

    if subrange is None:
        folder_name = args.folder + '/plots_repeat_' + str(r)
    else:
        folder_name = args.folder + '/plots_repeat_' + str(r) + '_subset_' + str(subrange[0]) + '_' + str(subrange[1])

    if os.path.isdir(folder_name):
        print('Plots already exist for repeat ', str(r))
    else:
        ensure_dir(folder_name)
        # data comes back as a list of touples. count, x,y, state

        if switch_3d:
            data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'z', 'state'])
        else:
            data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'state'])

        # now map continuous state column to 0, 1 or 2
        mapping = {range(0, 30): 0, range(30, 60): 1, range(60, 100): 2}

        data_df['state'] = data_df['state'].apply(lambda x: next((v for k, v in mapping.items() if x in k), 0))

        counts = data_df.drop_duplicates(subset=['count'], keep='last')

        if subrange is not None and switch_3d:
            data_df = data_df.loc[(data_df['x'] < subrange[1]) & (data_df['x'] > subrange[0])]
            data_df = data_df.loc[(data_df['y'] < subrange[1]) & (data_df['y'] > subrange[0])]
            data_df = data_df.loc[(data_df['z'] < subrange[1]) & (data_df['z'] > subrange[0])]
        elif subrange is not None:
            data_df = data_df.loc[(data_df['x'] < subrange[1]) & (data_df['x'] > subrange[0])]
            data_df = data_df.loc[(data_df['y'] < subrange[1]) & (data_df['y'] > subrange[0])]

        for value in counts['count'].values:

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

# fit_analysis(args.folder, min_r, max_r, time_step)

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
