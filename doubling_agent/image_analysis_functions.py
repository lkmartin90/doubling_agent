import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.spatial.distance import euclidean
import matplotlib.patches as mpatches
from scipy import optimize
import pandas as pd
import os
import fnmatch
plt.style.use('ggplot')


def plot_cells(data_df, value, folder_name, r, time_step):
    # basic function to plot the cells at a given time snapshot
    df_to_plot = data_df.loc[data_df['count'] == value]
    col = df_to_plot.state.map({0: 'b', 1: 'r', 2: 'g', 3: 'k'})
    df_to_plot.plot.scatter(x='x', y='y', c=col, s=8)
    blue_patch = mpatches.Patch(color='blue', label='Stem cell')
    red_patch = mpatches.Patch(color='red', label='Progentior cell')
    green_patch = mpatches.Patch(color='green', label='Differentiated cell')
    black_patch = mpatches.Patch(color='black', label='Quiescent cell')
    plt.legend(handles=[red_patch, blue_patch, green_patch, black_patch])
    plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                str(int(value*time_step)) + '.png')
    plt.cla()
    plt.close('all')


def plot_cells_3d(data_df, value, folder_name, r, time_step):
    # basic function to plot the cells at a given time snapshot (for 3D data)
    df_to_plot_3d = data_df.loc[data_df['count'] == value]
    col = df_to_plot_3d.state.map({0: 'b', 1: 'r', 2: 'g', 3: 'k'})
    fig = plt.figure()
    threedee = fig.gca(projection='3d')
    #print(df_to_plot.state.values)
    threedee.scatter(df_to_plot_3d.x.values, df_to_plot_3d.y.values, df_to_plot_3d.z.values, c=col)
    threedee.set_xlabel('x')
    threedee.set_ylabel('y')
    threedee.set_zlabel('z')
    blue_patch = mpatches.Patch(color='blue', label='Stem cell')
    red_patch = mpatches.Patch(color='red', label='Progentior cell')
    green_patch = mpatches.Patch(color='green', label='Differentiated cell')
    black_patch = mpatches.Patch(color='black', label='Quiescent cell')
    plt.legend(handles=[red_patch, blue_patch, green_patch, black_patch], loc='upper left')
    plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                str(int(value*time_step)) + '.png')
    plt.cla()
    plt.close('all')


def plot_2d_slice(folder_name, data_df, value, time_step, r):
    # Plot 2D slices of 3D data
    df_2d_slice = data_df.loc[data_df['count'] == value].copy()
    # will take slices in x to get 2d analysis
    x_values = df_2d_slice['z'].values
    unique, counts = np.unique(x_values, return_counts=True)
    tot_dict = dict(zip(unique, counts))

    for sect in unique:
        if tot_dict.get(sect) > 10:
            print(sect)
            df_for_image = df_2d_slice.loc[(data_df['z'] == sect)].copy()
            col = df_for_image.state.map({0: 'b', 1: 'r', 2: 'g', 3: 'k'})
            df_for_image.plot.scatter(x='x', y='y', c=col, s=8)
            blue_patch = mpatches.Patch(color='blue', label='Stem cell')
            red_patch = mpatches.Patch(color='red', label='Progentior cell')
            green_patch = mpatches.Patch(color='green', label='Differentiated cell')
            black_patch = mpatches.Patch(color='black', label='Quiescent cell')
            plt.legend(handles=[red_patch, blue_patch, green_patch, black_patch])
            plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                        str(int(value * time_step)) + '_sect_' + str(sect) + '.png')
            plt.cla()
            plt.close('all')


def fft_analysis(data_df, value,  folder_name, r, time_step, sect):
    # wanted to look at the FFT of the cell data to determine if there is a difference in frequencies present in
    # different simulations. For a comparison to experiment would perhapse have to use the ratio between the
    # frequencies for the different cell types
    df_fft = data_df.loc[data_df['count'] == value]
    label_dict = {0: "Stem", 1: "Progenitor", 2: "Differentiated", 3: "Quiescent"}
    # only look at the snapshot that match the value parameter

    # will plot figure with data from all 4 cell types
    fig, ax = plt.subplots(nrows=3, ncols=4)

    # Faff to set the scale of the plots
    if abs(df_fft['x'].max()) > abs(df_fft['x'].min()):
        x_len = df_fft['x'].max() + 1
    else:
        x_len = -1 * df_fft['x'].min() + 1

    if abs(df_fft['y'].max()) > abs(df_fft['y'].min()):
        y_len = df_fft['y'].max() + 2
    else:
        y_len = -1 * df_fft['y'].min() + 2

    if x_len > y_len:
        length = int(x_len)
    else:
        length = int(y_len)

    fft_stats = []

    # If there are no Quiescent cells just produces empty plot in 4th column
    for state_type in range(0, 4):
        df_for_image = df_fft.loc[df_fft['state'] == state_type]

        image = np.zeros((2 * length, 2 * length))
        for i in range(len(df_for_image)):
            image[length + int(np.round(df_for_image['x'].iloc[i]))][length + int(np.round(df_for_image['y'].iloc[i]))] = 1

        # Do the FFT
        ftimage = np.fft.fft2(image)
        ftimage = np.fft.fftshift(ftimage)
        freq = np.abs(ftimage)

        ax[0, state_type].hist(freq.ravel(), bins=100, range=(0, 30))
        ax[0, state_type].set_title(str(label_dict.get(state_type)))
        im1 = ax[1, state_type].imshow(freq, interpolation="none", cmap='jet', aspect="auto")
        fig.colorbar(im1, ax=ax[1, state_type])
        im2 = ax[2, state_type].imshow(image, interpolation="none", cmap='jet', aspect="auto")
        fig.colorbar(im2, ax=ax[2, state_type])
        fig.tight_layout()

        fft_stats.append(str(state_type))
        fft_stats.append(str(np.mean(freq)))
        fft_stats.append(str(np.std(freq)))

    plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                str(int(value * time_step)) + '_sect_' + str(sect) + 'FFT.png')
    plt.cla()
    plt.close('all')

    with open(folder_name + '/fft_data_day_' +
                        str(int(value * time_step)) + '_sect_' + str(sect) + '.txt', 'w') as f:
        f.write("\n".join(fft_stats))


def density_analysis(folder_name, data_df, value, k, switch_3d, time_step, r):
    # For the whole tumour look at the density of each cell type through the distance to a number of its
    # nearest neighbours, specified by k. In an attempt to quantify this further the data is then fitted
    # to with a simple function.

    df_dens = data_df.loc[data_df['count'] == value].copy()
    label_dict = {0: "Stem cell", 1: "Progenitor cell", 2: "Differentiated cell", 3: "Quiescent cell"}
    tot_cells = df_dens.shape[0]

    # If there are more cells than this the code will take a prohibatively long time to run
    if tot_cells > 30000:
        return

    # Distance between the array and itself
    plt.figure()

    sta = df_dens.drop_duplicates(subset=['state'], keep='last')
    for state_type in sta['state'].values:
        df_dens_state = df_dens.loc[data_df['state'] == state_type].copy()

        # BE WARNED THESE DATAFRAMES ARE COPIES, CHANGING THEM WILL CHANGE ORIGINAL

        if switch_3d:
            df_dens_state['coords'] = df_dens_state.loc[:, ['x', 'y', 'z']].values.tolist()
            df_dens_state['coords'] = df_dens_state.loc[:, 'coords'].apply(np.array)
        else:
            df_dens_state['coords'] = df_dens_state.loc[:, ['x', 'y']].values.tolist()
            df_dens_state['coords'] = df_dens_state.loc[:, 'coords'].apply(np.array)

        data = np.stack(df_dens_state['coords'].values, axis=0)
        dists = cdist(data, data)
        # Sort by distances
        k_nearest = np.sort(dists)[:, 1:k + 1]
        mean_k_nearest = np.mean(k_nearest, axis=1)

        distances = np.sort(mean_k_nearest)

        with open(folder_name + '/cell_numbers.txt', 'a+') as f:
            f.write(str(len(distances)) + '\n')

        plt.scatter(y=distances, x=np.arange(len(distances))/tot_cells, label=label_dict.get(state_type),
                    marker='+', alpha=0.5)

        if len(distances) > 2:
            try:
                fit_params, fit_params_covariance = optimize.curve_fit(fit_func,
                                    np.arange(1, len(distances)+1)/tot_cells,
                                    distances, p0=[1, 0.2, (len(distances)+1)/tot_cells])
                # print(fit_params)
                plt.plot(np.arange(len(distances))/tot_cells,
                                    fit_func(np.arange(len(distances))/tot_cells, fit_params[0],
                                    fit_params[1], fit_params[2]), label='Fitted function')

                with open(folder_name + '/fitting.txt', 'a+') as f:
                    # write the data to file in a noce, human readable way
                    f.write('value =' + str(value) + '\n' + str(label_dict.get(state_type))
                            + ' = ' + str(fit_params[0])
                            + ' ' + str(fit_params[1]) + ' ' + str(fit_params[2]) + '\n')

                with open(folder_name + '/fitting_dat.txt', 'a+') as f:
                    # write the data to file in a way which makes it easier to process later
                    f.write(str(value) + ' ' + str(state_type) + ' ' + str(fit_params[0])
                            + ' ' + str(fit_params[1]) + ' ' + str(fit_params[2]) + '\n')

            except RuntimeError or ValueError:
                print('Didnt find a good fit')
                with open(folder_name + '/fitting_dat.txt', 'a+') as f:
                    f.write(str(value) + ' ' + str(state_type) + ' ' + str(0)
                            + ' ' + str(0) + ' ' + str(0) + '\n')

    plt.xlabel('Cells, ordered from smallest to largest mean distance')
    plt.ylabel('Mean distance to 8 nearest neighbours')
    plt.legend(loc="upper right")
    plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                str(int(value*time_step)) + 'density.png')
    plt.close('all')


def distance_from_centre(folder_name, data_df, value, switch_3d, time_step, r):
    # find the distance of each type of cell from the tumour centre and plot

    df_dist = data_df.loc[data_df['count'] == value].copy()
    label_dict = {0: "Stem cell", 1: "Progenitor cell", 2: "Differentiated cell", 3: "Quiescent cell"}

    x_mean = np.mean(df_dist.loc[:, ['x']].values)
    y_mean = np.mean(df_dist.loc[:, ['y']].values)

    # find the centre of the tumour
    if switch_3d:
        z_mean = np.mean(df_dist.loc[:, ['z']].values)
        cent = [x_mean, y_mean, z_mean]
    else:
        cent = [x_mean, y_mean]

    print('Tumour center is ' + str(cent))
    plt.figure()

    # loop through the different cell states
    for state_type in range(0, df_dist.state.max()+1):
        df_for_image = df_dist.loc[data_df['state'] == state_type].copy()

        # BE WARNED THESE DATAFRAMES ARE COPIES, CHANGING THEM WILL CHANGE ORIGINAL

        if switch_3d:
            df_for_image['coords'] = df_for_image.loc[:, ['x', 'y', 'z']].values.tolist()
            df_for_image['coords'] = df_for_image.loc[:, 'coords'].apply(np.array)
        else:
            df_for_image['coords'] = df_for_image.loc[:, ['x', 'y']].values.tolist()
            df_for_image['coords'] = df_for_image.loc[:, 'coords'].apply(np.array)

        dist_list = []
        for i in range(len(df_for_image['coords'].values)):
            dist_list.append(euclidean(df_for_image['coords'].values[i], cent))

        plt.hist(dist_list, label=label_dict.get(state_type), alpha=0.5)

    plt.xlabel('Distance from tumour centre')
    plt.ylabel('Number of cells')
    plt.legend(loc="upper right")
    plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' + str(
        int(value * time_step)) + 'distance.png')
    plt.close('all')


def fit_func(x, a, b, c):
    # define the fitting function for density analysis (currently used)
    return a + b/(c - x)


def fit_func2(x, a, b, c, d):
    # define the more complex fitting function for density analysis (not currently used)
    return a + b/(c - x) - d/x


def fit_analysis(main_folder, min_r, max_r, time_step):
    # Compares the fitting parameters over a number of different repeats of the simulation
    # Compares the whole tumour data, not subsets or slices

    label_dict = {0: "Stem cell", 1: "Progenitor cell", 2: "Differentiated cell", 3: "Quiescent cell"}

    # loop over the repeats that we're dealing with
    for r in range(min_r, max_r+1):
        with open(main_folder + '/plots_repeat_' + str(r) + '/fitting_dat.txt', 'r') as f:
            fit_data = f.readlines()

        data_mat = []
        for i in range(len(fit_data)):
            dat = fit_data[i].strip().split()
            dat.append(str(r))
            dat = np.array(dat).astype(float)
            data_mat.append(dat)

        if r == min_r:
            data_df = pd.DataFrame(data_mat)
        else:
            data_df = data_df.append(pd.DataFrame(data_mat))
        # extract the values

    data_df = data_df.rename(columns={0: "value", 1: "state"})

    # find the points at which data was saved
    values = data_df.drop_duplicates(subset=['value'])['value'].values

    tot_len = len(data_df.columns)

    # find all different cell states present
    states = data_df.drop_duplicates(subset=['state'])['state'].values

    # loop over the possible "values", which are the time steps at which data was taken.
    for val in values:
        print(val)
        plt.figure()

        # loop over the diffent cell states present
        for sta in states:
            to_plot = data_df.loc[(data_df['value'] == val) & (data_df['state'] == sta)][tot_len-2].values
            plt.plot(to_plot, label=label_dict.get(sta))

        plt.legend(loc="upper left")
        plt.xlabel('repeat')
        plt.ylabel('fraction of cells in each state')
        plt.savefig(main_folder + '/day_' + str(int(val * time_step)) + 'fraction.png')
        plt.close('all')

    for val in values:
        print(val)
        plt.figure()

        for sta in states:
            to_plot = data_df.loc[(data_df['value'] == val) & (data_df['state'] == sta)][tot_len-3].values
            plt.plot(to_plot, label=label_dict.get(sta))
        plt.legend(loc="upper left")
        plt.xlabel('repeat')
        plt.ylabel('second fitting parameter')
        plt.savefig(main_folder + '/day_' + str(int(val * time_step)) + 'dens_change.png')
        plt.close('all')


def density_analysis_2d_slice(folder_name, data_df, value, k, time_step, r, subset, min_box_size):
    # performs the density analysis for 2D slices of the 3D data

    df_dens_2d = data_df.loc[data_df['count'] == value].copy()
    # will take slices in x to get 2d analysis
    z_values = df_dens_2d['z'].values

    # unique is the unique z values, and counts is the number of cells at this z value
    unique, counts = np.unique(z_values, return_counts=True)
    tot_dict = dict(zip(unique, counts))

    label_dict = {0: "Stem cell", 1: "Progenitor cell", 2: "Differentiated cell", 3: "Quiescent cell"}

    # here, sect is a z value identified in unique
    for sect in unique:

        # if there are more then 30 cell in this 3 axis slice...
        if tot_dict.get(sect) > 30:
            print(sect)

            # intitialise figure
            plt.figure()

            tot_cells = df_dens_2d.loc[df_dens_2d['z'] == sect].shape[0]
            # for each cell type at this z axis value...
            for state_type in np.unique(df_dens_2d.loc[df_dens_2d['z'] == sect]['state'].values):
                # takes a copy of the data in which cells are in this state
                df_for_image = df_dens_2d.loc[(data_df['state'] == state_type) & (data_df['z'] == sect)].copy()

                # BE WARNED THESE DATAFRAMES ARE COPIES, CHANGING THEM WILL CHANGE ORIGINAL

                df_for_image['coords'] = df_for_image.loc[:, ['x', 'y', 'z']].values.tolist()
                df_for_image['coords'] = df_for_image.loc[:, 'coords'].apply(np.array)

                data = np.stack(df_for_image['coords'].values, axis=0)
                # Distance between the array and itself
                dists = cdist(data, data)
                # Sort by distances
                k_nearest = np.sort(dists)[:, 1:k + 1]
                mean_k_nearest = np.mean(k_nearest, axis=1)
                # print(mean_k_nearest)
                distances = np.sort(mean_k_nearest)

                with open(folder_name + '/sect_' + str(sect) + '_cell_numbers.txt', 'a+') as f:
                    f.write(str(len(distances)) + '\n')

                plt.scatter(y=distances, x=np.arange(len(distances))/tot_cells, label=label_dict.get(state_type),
                            marker='+', alpha=0.5)

                # if there are more than 3 of this cell type fit to the data
                if np.isnan(distances[0]) == False and len(distances) > 3:

                    try:
                        fit_params, fit_params_covariance = optimize.curve_fit(fit_func,
                                                                np.arange(1, len(distances) + 1)/tot_cells,
                                                                distances, p0=[1, 0.2, (len(distances) + 1)/tot_cells])

                        plt.plot(np.arange(len(distances))/tot_cells, fit_func(np.arange(len(distances))/tot_cells,
                                                                     fit_params[0], fit_params[1],
                                                                     fit_params[2]), label='Fitted function')

                        # write the dat to fine in a way which makes it easier to process later on
                        with open(folder_name + '/sect_' + str(sect) +
                                  '_fitting_dat.txt', 'a+') as f:
                            f.write(str(value) + ' ' + str(state_type) + ' ' + str(fit_params[0])
                                    + ' ' + str(fit_params[1]) + ' ' + str(fit_params[2]) + '\n')

                    except RuntimeError or ValueError or TypeError:
                        with open(folder_name + '/sect_' + str(sect) + '_fitting_dat.txt', 'a+') as f:
                            f.write(str(value) + ' ' + str(state_type) + ' ' + str(0)
                                    + ' ' + str(0) + ' ' + str(0) + '\n')
                        print('Didnt find a good fit')

            plt.xlabel('Cells, ordered from smallest to largest mean distance')
            plt.ylabel('Mean distance to 8 nearest neighbours')
            plt.legend(loc="upper right")
            plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                        str(int(value * time_step)) + 'sect_' + str(sect) + 'density.png')
            plt.close('all')

            # Here as the data is already in the correct format, just pass to fft analysis.
            df_for_analysis = df_dens_2d.loc[(data_df['z'] == sect)].copy()
            if subset is not None:
                lacunarity(folder_name, df_for_analysis, value, time_step, r, subset, sect, min_box_size)
            fft_analysis(df_for_analysis, value, folder_name, r, time_step, sect)


def lacunarity(folder_name, data_df, value, time_step, r, subset, sect, min_box_size):
    # The idea here is to split the image into boxes of many sizes and to count the number of cells in each box,
    # the standard deviation in relation to the mean tells you about the clumping, and by changing the scale of
    # the box this gives you the scale of the variation
    # want this function to be able to take 2D or 3D input, if input is 3D then o the same thing for each slice.
    # print(value)
    label_dict = {0: "Stem cell", 1: "Progenitor cell", 2: "Differentiated cell", 3: "Quiescent cell"}

    # need to decide on box size, Feel like this won't work for small tumours, need a large area that is evenly populated.
    min_data_size = 20
    min_grid_size = min_box_size
    max_grid_size = (subset[1] - subset[0])/4

    df_lac = data_df.loc[data_df['count'] == value].copy()
    # only use for tumours that have been subsetted
    if subset[1] - subset[0] < min_data_size:
        print('Subset specified too small for lacunarity analysis')
    else:
        plt.figure()
        for state_type in np.unique(df_lac['state'].values):
            # no idea how computationally intensive this will be or what is a sensible number of boxes to use.
            # can simply bin data into boxes of the correct size with different starting points?
            # need to specify the start points of the bins
            df_process = df_lac.loc[(df_lac['state'] == state_type)].copy()
            lac_r = []
            r_size = []

            for bin_size in range(int(min_grid_size), int(max_grid_size)+1):

                lac_r_data = np.array([])
                r_size.append(bin_size)
                for i in range(bin_size):
                    x_edge = np.arange(subset[0] + i, subset[1] - bin_size, bin_size)
                    # can then set the range of the histogram based on these values
                    x = df_process['x'].values
                    y = df_process['y'].values
                    hist = np.histogram2d(x, y, bins=x_edge,
                                         range=[[np.min(x_edge), np.max(x_edge)+bin_size],
                                                  [np.min(x_edge), np.max(x_edge)+bin_size]])

                    lac_r_data = np.append(lac_r_data, np.ndarray.flatten(hist[0]))

                lac_r.append((np.std(lac_r_data)/np.mean(lac_r_data))**2)

            plt.plot(r_size, lac_r, label=label_dict.get(state_type))

            with open(folder_name + '/lac_tot_day_' + str(int(value * time_step)) + '_sect_' + str(sect) +
                      '_state_' + str(state_type) + '.txt','a+') as f:
                f.write(str(r_size).strip('[]') + '\n')
                f.write(str(lac_r).strip('[]'))

        plt.legend()
        plt.xlabel('box size (r)')
        plt.ylabel('lacunarity')
        plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                        str(int(value * time_step)) + '_subset_' + str(subset[0]) + '_' + str(subset[1])
                    + '_sect_' + str(sect) + '_lac.png')
        plt.close('all')



def section_analysis(folder, file_pattern, quantity, label_type, minr, maxr):
    # analyses data from 2D slices of 3D data, can be used to plot these quantities as a function of position along the
    # z axis. Different quanities can be plotted with this function by changing "file pattern" and "quantity"
    for r in range(minr, maxr+1):

        for filename in os.listdir(str(folder)):
            if fnmatch.fnmatch(filename, 'plots_repeat_' + str(r) + '*'):
                cell_nos = []
                sect = []
                for filename2 in os.listdir(str(folder) + '/' + filename):
                    if fnmatch.fnmatch(filename2, file_pattern):
                        sect.append(int(filename2.split('_')[1]))

                        with open(str(folder) + '/' + filename + '/' + filename2) as f:
                            lines = f.read().splitlines()
                        cell_nos.append(lines)

                cell_nos_df = pd.DataFrame(cell_nos)

                if len(sect) > 2:
                    plt.figure()
                    for i in range(len(cell_nos_df.columns)):
                        print(i)
                        plt.scatter(x=sect, y=cell_nos_df[i].astype(float).values, linestyle='None',
                                    label=label_type.get(i))
                    plt.xlabel('Distance of sect along z axis')
                    plt.ylabel(quantity)
                    plt.legend()
                    plt.savefig(str(folder) + '/' + filename + '/' + str(quantity))


def section_analysis_fft(folder, file_pattern, quantity, cell_tye, label_type, minr, maxr):
    # Analyses the FFT data for a number of 2D slices over 3D data. FFt has already been computed for these
    #slices and saved to file

    # loop over the number of repeats to be analysed
    for r in range(minr, maxr+1):
        # find the data from the file name in which the data is stored
        for filename in os.listdir(str(folder)):
            if fnmatch.fnmatch(filename, 'plots_repeat_' + str(r) + '*'):

                fft_stats = []
                sect = []
                for filename2 in os.listdir(str(folder) + '/' + filename):
                    if fnmatch.fnmatch(filename2, file_pattern) and filename2.split('_')[5].split('.')[0] != 'full':
                        sect.append(int(filename2.split('_')[5].split('.')[0]))

                        with open(str(folder) + '/' + filename + '/' + filename2) as f:
                            lines = f.read().splitlines()
                        fft_stats.append(lines)

                fft_stats_df = pd.DataFrame(fft_stats)

                if len(sect) > 2:
                    plt.figure()
                    # The indices containing the data we wish to plot, omittied indices denote cell type
                    for i in [1, 2, 4, 5, 7, 8, 10, 11]:

                        plt.scatter(x=sect, y=fft_stats_df[i].astype(float).values, linestyle='None',
                                    label=label_type.get(i % 2) + cell_tye.get(np.floor(i/3)))
                    plt.xlabel('Distance of slice along z axis')
                    plt.ylabel(quantity)
                    plt.legend()
                    plt.savefig(str(folder) + '/' + filename + '/' + str(quantity))


def lac_analysis(folder, file_pattern, quantity, minr, maxr, min_box_size):
    # Analyses the lacunarity data for many slices of a 3D data ser, taking the mean for each box size
    # and plotting this for each cell type. The lacunarity data for these slices has already been computed
    # and saved to file.

    print('lacunarity analysis')
    # loops over the repeats to analyse
    for r in range(minr, maxr+1):
        # finds the data files with the correct names for the lacunarity data of these repeats
        for filename in os.listdir(str(folder)):
            if fnmatch.fnmatch(filename, 'plots_repeat_' + str(r) + '*'):
                lac_nos = []
                for filename2 in os.listdir(str(folder) + '/' + filename):
                    if fnmatch.fnmatch(filename2, file_pattern) and filename2.split('_')[5].split('.')[0] != 'full':
                        with open(str(folder) + '/' + filename + '/' + filename2) as f:
                            lines = f.read().splitlines()[0].split(',')
                        lac_nos.append([float(i) for i in lines])

                # take mean of data
                lac_mean = np.mean(lac_nos, axis=0)

                # find the correct box size for plotting
                box_size = np.arange(min_box_size, len(lac_mean)+min_box_size)

                plt.figure()
                plt.scatter(box_size, lac_mean)
                plt.xlabel('Box size')
                plt.ylabel('Mean lacunarity')
                plt.savefig(str(folder) + '/' + filename + '/' + str(quantity))

                np.savetxt(str(folder) + '/' + filename + '/' + str(quantity) + '.txt',
                           np.concatenate((box_size, lac_mean), axis=0), fmt="%s")
