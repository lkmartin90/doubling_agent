#!/usr/bin/env python
import numpy as np
import argparse
from datetime import datetime
from doubling_agent.state_functions import *
from matplotlib import pyplot as plt

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
                    default=np.array([0., 0., 0.5, 0.5, 0.3, 1., 0., 0.]))
parser.add_argument('--dim', type=int, nargs='?', help='Dimension for the simulation, 2 or 3.', const=2, default=2)
parser.add_argument('--motparams', action=StoreAsArray, type=float, nargs='+',
                    help=" Motility parameters for the simulation, if no movement desired leave blank.",
                    default=[0.])
parser.add_argument('--length', type=int, nargs='?', help='Length of the simulation in days',
                    const=200, default=200)
parser.add_argument('--repeats', type=int, nargs='?', help='Number of repeats for the simulation',
                    const=1, default=1)
parser.add_argument('--saving', help='Do you require data saving? include flag if no', action='store_false')

args = parser.parse_args()
print(args)

########################################
# interpreting user defined parameters
########################################

FORMAT = '%y%m%d%H%M'
date_and_time = str(datetime.now().strftime(FORMAT))

if len(args.params) == 8:
    sim_type = "State"
else:
    raise ValueError("Unrecognised number of parameters!")

params = []
if sim_type == "State":
    params = ParametersState(args.params[0], args.params[1], args.params[2], args.params[3],
                             args.params[4], args.params[5], args.params[6], args.params[7])

mot_params = []
if args.motparams is not None:
    mot_params = MotilityParameters(args.motparams[0])

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

# For each repeat want to create a file with the data stored in it. Each run should have a different
# folder in output_data
folder_name = 'output_data/' + sim_type + "_" + date_and_time + '_' + str(args.dim) + 'D_' + str(repeats) \
              + 'r_' + str(sim_length) + 'l'
if args.saving:
    ensure_dir(folder_name)

# How often to save data based on timestep size and length

if sim_length <= 10:
    record_after = 1 / time_step
elif sim_length < 400:
    record_after = 10 / time_step
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

# In this simulation will include reaction-diffusion and will slightly re-write. This time the timing of cells
# will only be updated once something has happened to them.

tot_cells = []
plt.figure()
for r in range(repeats):
    print('repeat ' + str(r + 1))
    cells = {}
    cancer_seed_single(cells, switch_3d)
    #print(cells)
    # "For bus arrivals consistent with a Poisson process,
    # the expected wait time for a passenger is identical to the average interval between arrivals."
    count = 0
    timing_update_all(cells, params, mot_params, time_step)

    # Things are adjusted slightly for this model. There are 3 columns per cell: time to next event, cell state (0,1,2)
    # and "division time". Once a cell has comitted
    # to dividing it can't move for a certain amount of time, this 'division time' will be updated on every time step.
    # If it contains a none 0 value the cell can neither divide or move again. Should make model slightly more realistic.

    no_of_cells = []
    while count < sim_steps + 1:
        print(count)
        #print(cells)
        if count % record_after == 0 or count == sim_steps:
            print(count)
            for cell in cells:
                # want data in the form: time step, x, y, state
                if switch_3d:
                    data = [count, cell[0], cell[1], cell[2], cells.get(cell)[0]]
                else:
                    data = [count, cell[0], cell[1], cells.get(cell)[0]]

                if args.saving:
                    write_to_file(data, folder_name + '/repeat_' + str(r))

        timing_update_all(cells, params, mot_params, time_step)
        # then want to use a map to select cells relevant in this time step, as all cells will not be updated
        sorted_cells = sorted(cells.items(), key=lambda cells: cells[1][1])

        for cell in sorted_cells:
            # print(cell[0])
            t = cells.get(cell[0])[1]

            if t < time_step:
                update_cell_diff(cells, cell[0], params, switch_3d, mot_params)
            else:
                break

        count = count + 1
        no_of_cells.append(len(cells))

    tot_cells.append([no_of_cells[-1]])

    plt.plot(no_of_cells)
    plt.xlabel('simulation step (length/timestep)')
    plt.ylabel('Number of cells')
    if args.saving:
        plt.savefig(folder_name + '/repeat_' + str(r) + '_tumour_growth.png')
    plt.cla()

if args.saving:
    with open(folder_name + '/metadata.txt', 'a+') as f:
        f.write("\nmean number of cells = " + str(np.mean(tot_cells)))