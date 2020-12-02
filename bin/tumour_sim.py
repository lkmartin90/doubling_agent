#!/usr/bin/env python
import numpy as np
import argparse
from datetime import datetime
from doubling_agent.common_functions import *
from matplotlib import pyplot as plt
import pandas as pd

########################################
# including the user defined parameters
########################################

class StoreAsArray(argparse._StoreAction):
    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)


parser = argparse.ArgumentParser(description='Process input parameters for simulation')
parser.add_argument('--params', action=StoreAsArray, type=float, nargs='+', help="Parameters for the simulation, "
                                            "please ensure they are entered in the correct order.",
                    default=np.array([0.15, 0.15, 1.0, 0.48]))
parser.add_argument('--dim', type=int, nargs='?', help='Dimension for the simulation, 2 or 3.', const=2, default=2)
parser.add_argument('--motparams', action=StoreAsArray, type=float, nargs='+',
                    help=" Motility parameters for the simulation, "
                    "please ensure they are entered in the correct order, if no movement desired leave blank.")
parser.add_argument('--length', type=int, nargs='?', help='Length of the simulation in days',
                    const=200, default=200)
parser.add_argument('--repeats', type=int, nargs='?', help='Number of repeats for the simulation',
                    const=1, default=1)
parser.add_argument('--saving', help='Do you require data saving? include flag if no', action = 'store_false')
parser.add_argument('--killcells', action=StoreAsArray, type=float, nargs='+',
                    help=" Include if you want to kill a population of cells"
                    "requires 2 parameters, cell type to kill and day of death")
args = parser.parse_args()
print(args)

########################################
# interpreting user defined parameters
########################################

FORMAT = '%y%m%d%H%M'
date_and_time = str(datetime.now().strftime(FORMAT))

# if no motility import the correct function
if args.motparams is None:
    from doubling_agent.basic_functions import *
else:
    from doubling_agent.motility_functions import *

if len(args.params) == 4:
    sim_type = "Basic"
elif len(args.params) == 8:
    sim_type = "Quiescent"
else:
    raise ValueError("Unrecognised number of parameters!")

params = []
if sim_type == "Basic":
    params = ParametersBasic(args.params[0], args.params[1], args.params[2], args.params[3])
elif sim_type == "Quiescent":
    params = ParametersQuiescent(args.params[0], args.params[1], args.params[2], args.params[3],
                                 args.params[4], args.params[5], args.params[6], args.params[7])

mot_params = []
if args.motparams is not None:
    mot_params = MotilityParameters(args.motparams[0], args.motparams[1], args.motparams[2])

time_step = calculate_timestep(params, mot_params)

if args.dim == 3:
    switch_3d = True
elif args.dim == 2:
    switch_3d = False
else:
    raise ValueError('Not a recognised number of dimensions')

repeats = args.repeats
sim_steps = int(round(args.length/time_step))
sim_length = args.length
killcells = args.killcells
if killcells is None:
    print('No treatment')
    killcells = [-1, -1]
elif len(killcells) != 2:
    print('Incorrect number of parameters in killcells, try again')
    exit()

# For each repeat want to create a file with the data stored in it. Each run should have a different
# folder in output_data
folder_name = 'output_data/' + sim_type + "_" + date_and_time + '_' + str(args.dim) + 'D_' + str(repeats) \
              + 'r_' + str(sim_length) + 'l'
if args.saving:
    ensure_dir(folder_name)

# How often to save data based on timestep size and length

if sim_length <= 100:
    record_after = 1 / time_step
#elif sim_length < 100:
#    record_after = 10 / time_step
elif sim_length < 10000:
    record_after = 100 / time_step
else:
    record_after = 1000 / time_step

########################################
# write description of data
########################################
if args.saving:
    with open(folder_name + '/metadata.txt', 'w') as f:
        f.write('parameters = ' + str(args.params) + '\nmotility parameters = '
                + str(args.motparams) + '\ntimestep = ' + str(time_step) +
                '\nlength in days = ' + str(sim_length) + '\nnumber of repeats = ' + str(repeats) +
                '\ndimensions = ' + str(args.dim))

######################################
# perform simulation
########################################
tot_cells = []
s0 = []
s1 = []
s2 = []
if sim_type == "Quiescent":
    s3 = []

plt.figure()
for r in range(repeats):
    print('repeat ' + str(r + 1))
    cells = {}
    cancer_seed_single(cells, switch_3d)
    #print(cells)
    # "For bus arrivals consistent with a Poisson process,
    # the expected wait time for a passenger is identical to the average interval between arrivals."
    count = 0
    timing_update_all(cells, params, mot_params)

    # Now need to deal with motility. There are two states, moving (1) and proliferating (0). There is a
    # rate of switching between the two. In the (0) state the cell progression is as it was in the basic
    # simulation. If the cell is in the (1) state it can either move or die, equally there is a rate for
    # cell movement.

    # There are now more things that can happen to each cell type, this must be reflected in the timing
    # update, does the cell divide, or does it become mobile first? Want to use the same trick as in the
    # quiescent case, combining the two rates and then deciding which event has happened.

    # each cell is now described by a set of coordinates and the vector [cell type, time, motility state]

    no_of_cells = []
    while count < sim_steps + 1:
        print(count)
        if count % record_after == 0 and args.saving or count == sim_steps and args.saving:
            #print(count)
            for cell in cells:
                # want data in the form: time step, x, y, state
                if args.motparams is None:
                    if switch_3d:
                        data = [count, cell[0], cell[1], cell[2], cells.get(cell)[0]]
                    else:
                        data = [count, cell[0], cell[1], cells.get(cell)[0]]
                else:
                    if switch_3d:
                        data = [count, cell[0], cell[1], cell[2], cells.get(cell)[0], cells.get(cell)[2]]
                    else:
                        data = [count, cell[0], cell[1], cells.get(cell)[0], cells.get(cell)[2]]

                write_to_file(data, folder_name + '/repeat_' + str(r))

        if count == killcells[1]/time_step:
            print(count)
            kill_cells(cells, int(killcells[0]))

        timing_update_all(cells, params, mot_params)
        # then want to use a map to select cells relevant in this time step, as all cells will not be updated
        sorted_cells = sorted(cells.items(), key=lambda cells: cells[1][1])

        for cell in sorted_cells:
            #print(cell[0])
            t = cells.get(cell[0])[1]

            if t < time_step:
                if sim_type == "Quiescent":
                    update_cell_quiescent(cells, cell[0], params, switch_3d, mot_params)
                else:
                    update_cell_basic(cells, cell[0], params, switch_3d, mot_params)
            else:
                break

        count = count + 1
        no_of_cells.append(len(cells))

    tot_cells.append([no_of_cells[-1]])
    cells_df = pd.DataFrame(cells.values())
    cells_df.rename(columns={0: 'state'}, inplace=True)
    #print(cells_df)
    s0.append(cells_df.loc[cells_df['state'] == 0].shape[0])
    s1.append(cells_df.loc[cells_df['state'] == 1].shape[0])
    s2.append(cells_df.loc[cells_df['state'] == 2].shape[0])
    if sim_type == "Quiescent":
        s3.append(cells_df.loc[cells_df['state'] == 3].shape[0])

    plt.plot(no_of_cells)
    plt.xlabel('simulation step (length/timestep)')
    plt.ylabel('Number of cells')
    if args.saving:
        plt.savefig(folder_name + '/repeat_' + str(r) + '_tumour_growth.png')
    plt.cla()

if args.saving:
    with open(folder_name + '/metadata.txt', 'a+') as f:
        f.write("\nmean number of cells = " + str(np.mean(tot_cells)) +
                "\nstd of cell numbers = " + str(np.std(tot_cells)) +
                "\nnumber of S0 = " + str(np.mean(s0)) + "\nstd of S0 = " + str(np.std(s0)) +
                "\nnumber of S1 = " + str(np.mean(s1)) + "\nstd of S1 = " + str(np.std(s1)) +
                "\nnumber of S2 = " + str(np.mean(s2)) + "\nstd of S2 = " + str(np.std(s2)))
        if sim_type == "Quiescent":
            f.write("\nnumber of S3 = " + str(np.mean(s3)) + "\nstd of S3 = " + str(np.std(s3)))

