from doubling_agent.image_analysis_functions import *
plt.style.use('ggplot')

def plot_cells_cont(data_df, value, folder_name, r, time_step):
    # basic function to plot the cells at a given time snapshot, when cell states are continuous
    df_to_plot = data_df.loc[data_df['count'] == value]
    df_to_plot.plot.scatter(x='x', y='y', c='state', s=8, cmap='viridis', vmin=0, vmax=100)
    plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                str(int(value*time_step)) + '_cont.png')
    plt.cla()
    plt.close('all')


def plot_cells_3d_cont(data_df, value, folder_name, r, time_step):
    # basic function to plot the cells at a given time snapshot, when the cell states are continuous
    df_to_plot_3d = data_df.loc[data_df['count'] == value]
    fig = plt.figure()
    threedee = fig.gca(projection='3d')
    #print(df_to_plot.state.values)
    p = threedee.scatter(df_to_plot_3d.x.values, df_to_plot_3d.y.values, df_to_plot_3d.z.values, c=df_to_plot_3d.state,
                     cmap='viridis', vmin=0, vmax=100)
    threedee.set_xlabel('x')
    threedee.set_ylabel('y')
    threedee.set_zlabel('z')
    fig.colorbar(p)
    plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                str(int(value*time_step)) + '_cont.png')
    plt.cla()
    plt.close('all')


def plot_2d_slice_cont(folder_name, data_df, value, time_step, r):
    # plot 2d slices of the data when the data is continuous

    df_2d_slice = data_df.loc[data_df['count'] == value].copy()
    # will take slices in x to get 2d analysis
    x_values = df_2d_slice['z'].values
    unique, counts = np.unique(x_values, return_counts=True)
    tot_dict = dict(zip(unique, counts))

    for sect in unique:
        # if there are more than 10 cells at this value of z, then plot
        if tot_dict.get(sect) > 10:
            print(sect)
            df_for_image = df_2d_slice.loc[(data_df['z'] == sect)].copy()
            plt.figure()
            df_for_image.plot.scatter(x='x', y='y', c='state', s=8, cmap='viridis', vmin=0, vmax=100)
            plt.savefig(folder_name + '/repeat_' + str(r) + '_day_' +
                        str(int(value * time_step)) + '_sect_' + str(sect) + '_cont.png')
            plt.close('all')
