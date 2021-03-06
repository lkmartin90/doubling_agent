from random import random
from random import choice
import numpy as np
import plotly.express as px
import struct
import os

###
# Broadly the same as "basic_functions.py" but updated to include motility
# intentionally trying to keep them separate so as not to slow down the basic version
###


class MotilityParameters:
    def __init__(self, motility_rate):
        self.rate = motility_rate


class ParametersDens:
    def __init__(self, g1, g2, d1, d2, b, h, m, div_time):
        # - g2: [0] -> [1]
        # - g1: [1] -> [0]
        # - d1: [2] -> [1]
        # - d2: [1] -> [2]
        # - b: apoptosis rate (same for all cells)
        # - h: division rate of [2] cells
        # - m: division rate of [1] cells
        # - div_time: The time for which cells cannot move or divide after they have just divided.
        self.g1 = g1
        self.g2 = g2
        self.d1 = d1
        self.d2 = d2
        self.b = b
        self.h = h
        self.m = m
        self.t = div_time

        # want a dictionary for each cell type to access parameters
        # for the middle state the rate of transition down is first
        self.state_trans = {0: [g2], 1: [g1, d2], 2: [d1]}
        self.div_rate = {0: 0, 1: m, 2: h}


def cancer_seed_single(cells, switch_3d):
    # created initial cancer stem cell at [0,0]
    if switch_3d:
        cells.update({(0, 0, 0): [2, 0, 0, 0]})
    else:
        cells.update({(0,0): [2, 0, 0, 0]})


def timing_update_all(cells, params, mot_params, time_step, transition):
    # update second entry in dict to give a timing based on the first entry, the state
    # time is log(1/rand_no)/rate

    # here include the effect of local density on transition rates.
    for k in cells.keys():

        state = cells.get(k)[0]
        div_time = cells.get(k)[2]
        density = cells.get(k)[3]
        transition_dict = return_transition_dict(transition, density)
        move = mot_params.rate  # move rate
        death = params.b  # death rate

        #print(np.array(params.state_trans.get(state)))
        #print(transition_dict.get(state))
        # np.sum(params.state_trans.get(state)) is the sum of the transition rates for each cell. Want also to include
        # the density effects here.
        transition_sum = np.sum(np.array(params.state_trans.get(state))*np.array(transition_dict.get(state)))

        if div_time > 0:
            div_time_new = div_time - time_step
            rate = transition_sum + death

        elif div_time < 0:
            div_time_new = 0.
            rate = transition_sum + death + params.div_rate.get(state)

        else:
            div_time_new = div_time
            rate = transition_sum + death + params.div_rate.get(state)

        rate = rate + move

        cells.update({k: [state, np.log(1/random())/rate, div_time_new, density]})


def choose_new_pos(pos, cells):
    # Identifies a free position for a cell to divide or move into. In this function a 2d square grid is used
    # space is searched for in the surrounding area, by random number generator, if there is already a cell
    # occupying the space then that space is excluded from possible locations and a new random number is generated.
    i = pos[0]
    j = pos[1]

    neighbours = [(i+1, j), (i-1, j), (i, j-1), (i, j+1)]
    options = [0, 1, 2, 3]
    cont = 0
    new_pos = 0
    while cont == 0 and len(options) > 0:
        pick = choice(options)
        check = neighbours[pick]
        if check in cells:
            options.remove(pick)
        else:
            cont = 1
            new_pos = check

    return new_pos


def choose_new_pos_eq(pos, cells):
    # chooses a new position by identifying all the free spaces first and then assigning them all equal probability
    i = pos[0]
    j = pos[1]

    neighbours = [(i+1, j), (i-1, j), (i, j-1), (i, j+1)]
    options = [0, 1, 2, 3]

    for n in range(len(neighbours)):
        if neighbours[n] in cells:
            options.remove(n)

    if len(options) > 0:
        new_pos = neighbours[choice(options)]
    else:
        new_pos = 0

    return new_pos


def choose_new_pos_3d(pos, cells):
    # 3d version of "choose_new_pos", the same method is used
    i = pos[0]
    j = pos[1]
    k = pos[2]
    # this currently assumes only square transitions on the cubic grid, may want to alter
    neighbours = [(i + 1, j, k), (i - 1, j, k), (i, j + 1, k), (i, j - 1, k), (i, j, k + 1), (i, j, k - 1)]
    options = [0, 1, 2, 3, 4, 5]
    cont = 0
    new_pos = 0
    while cont == 0 and len(options) > 0:
        pick = choice(options)
        check = neighbours[pick]
        if check in cells:
            options.remove(pick)
        else:
            cont = 1
            new_pos = check

    return new_pos


def choose_new_pos_3d_eq(pos, cells):
    # 3d version of "choose_new_pos", the same method is used
    i = pos[0]
    j = pos[1]
    k = pos[2]
    # this currently assumes only square transitions on the cubic grid, may want to alter
    neighbours = [(i + 1, j, k), (i - 1, j, k), (i, j + 1, k), (i, j - 1, k), (i, j, k + 1), (i, j, k - 1)]
    options = [0, 1, 2, 3, 4, 5]
    cont = 0
    new_pos = 0
    while cont == 0 and len(options) > 0:
        pick = choice(options)
        check = neighbours[pick]
        if check in cells:
            options.remove(pick)
        else:
            cont = 1
            new_pos = check

    return new_pos


def move_cell(cells, pos, state, switch_3d, div_time):
    # moves the cell
    if switch_3d:
        new_location = choose_new_pos_3d(pos, cells)
    else:
        new_location = choose_new_pos(pos, cells)

    if new_location != 0:
        del cells[pos]
        cells.update({new_location: [state, 0, div_time, 0]})


def update_cell_diff(cells, pos, params, switch_3d, mot_params, transition):
    # updates a given cell based on the current state of that cell
    # pos is string describing position
    # time is from random number generator giving time of interaction
    # cells is dict describing all cells in the tumour

    state = cells.get(pos)[0]
    div_time = cells.get(pos)[2]
    density = cells.get(pos)[3]
    transition_dict = return_transition_dict(transition, density)
    transition_sum = np.sum(np.array(params.state_trans.get(state)) * np.array(transition_dict.get(state)))

    if switch_3d:
        daughter = choose_new_pos_3d(pos, cells)
    else:
        daughter = choose_new_pos(pos, cells)

    # There are a number of options, cells can either move up or down the state scale, they can divide, die, or move

    if div_time <= 0:
        # cells can divide, die, move, and change state
        r_num = random()
        tot = transition_sum + params.b + params.div_rate.get(state) + mot_params.rate

        if tot == 0:
            print('Total time = 0, Cell stuck in state, something has gone wrong')
            print(state)
            print(density)
            print(transition_sum)
            print(transition_dict)

        elif r_num < params.div_rate.get(state)/tot:
            # divide
            if daughter != 0:
                cells.update({daughter: [state, 0, params.t, 0]})
                cells.update({pos: [state, 0, params.t, 0]})

        elif r_num < (params.div_rate.get(state) + mot_params.rate)/tot:
            # move
            move_cell(cells, pos, state, switch_3d, div_time)

        elif r_num < (params.div_rate.get(state) + mot_params.rate + transition_sum)/tot:
            # change state, this is dependent on the sate the cell is currently in
            r_num_2 = random()
            if state != 1:
                # only one possibility
                cells.update({pos: [1, 0, 0, 0]})
            elif r_num_2 < (np.array(params.state_trans.get(state)) * np.array(transition_dict.get(state)))[0]/transition_sum:
                # first transition happens
                cells.update({pos: [0, 0, 0, 0]})
            else:
                # second transition happens
                cells.update({pos: [2, 0, 0, 0]})

        else:
            # die
            del cells[pos]


    else:
        # cells can't divide or move, but can die and change state
        r_num = random()
        tot = params.b + transition_sum
        if tot == 0:
            print('Total time = 0, Cell stuck in state, something has gone wrong')

        elif r_num < transition_sum/tot:
            # change state, depends on the state the cell is currently in
            r_num_2 = random()
            if state != 1:
                # only one possibility
                cells.update({pos: [1, 0, div_time, 0]})
            elif r_num_2 < (np.array(params.state_trans.get(state)) * np.array(transition_dict.get(state)))[0] / transition_sum:
                # first transition happens
                cells.update({pos: [0, 0, div_time, 0]})
            else:
                # second transition happens
                cells.update({pos: [2, 0, div_time, 0]})
        else:
            # die
            del cells[pos]



def animate(animation_df, r, name):
    # animate the simulations using plotly and save as a .html

    animation_df['coord'] = animation_df[['x', 'y']].values.tolist()
    animation_df['coord'] = animation_df['coord'].apply(lambda x: np.array(x))
    #print(animation_df)
    if len(animation_df['coord'].values[0]) > 2:
        print("currently cannot animate for 3d")
        raise ValueError()

    mapping = {0: 'stem cell', 1: 'progenitor cell', 2: 'differentiated cell', 3: 'quiescent cell'}
    animation_df = animation_df.replace({'state': mapping})

    animation_df = animation_df.append(
        {'state': 'differentiated cell', 'count': 0, 'coord': 0, 'x': 10000, 'y': 10000},
        ignore_index=True)
    animation_df = animation_df.append(
        {'state': 'progenitor cell', 'count': 0, 'coord': 0, 'x': 10000, 'y': 10000},
        ignore_index=True)
    animation_df = animation_df.append(
        {'state': 'quiescent cell', 'count': 0, 'coord': 0, 'x': 10000, 'y': 10000},
        ignore_index=True)

    fig = px.scatter(animation_df, x="x", y="y", animation_frame="count",
            color='state', size_max=55, range_x=[-50, 50], range_y=[-50, 50])
    fig.update_traces(marker=dict(size=12))
    fig.layout.updatemenus[0].buttons[0].args[1]["frame"]["duration"] = 20
    fig.show()
    fig.write_html(name + '/ani_' + str(r) + '.html')


def read_from_file(file_name, switch_3d):
    # read data from binary file in the form: time step, x, y, state, motility
    if switch_3d:
        struct_fmt = '=iiiii'  # 6 ints
    else:
        struct_fmt = '=iiii'  # 5 ints
    struct_len = struct.calcsize(struct_fmt)
    struct_unpack = struct.Struct(struct_fmt).unpack_from

    results = []
    with open(file_name, "rb") as f:
        while True:
            data = f.read(struct_len)
            if not data: break
            s = struct_unpack(data)
            results.append(s)

    return results


def calculate_timestep(params, mot_params):
    # calculates timestep based on the probability of 2 or more events happening in a timestep (<0.01)
    max_rate = params.g1 + params.g2 + params.d1 + params.m + params.b + mot_params.rate
    # playing it safe by summing max of each
    lambert = 0.135157
    step = lambert/max_rate
    print('exact timestep from calculation', step)
    if step > 0.1:
        return step // 0.1 * 0.1
    elif step > 0.01:
        return step // 0.01 * 0.01
    else:
        return step // 0.001 * 0.001


def write_to_file(data, file_name):
    # write data to binary file in the form: time step, x, y, state
    #print(data)
    s = struct.pack('i' * len(data), *data)
    with open(file_name, 'ab') as f:
        f.write(s)


def ensure_dir(file_path):
    # if the directory doesn't already exist then make it
    directory = file_path
    if not os.path.exists(directory):
        print('...making', directory)
        os.makedirs(directory)


def calc_local_density(distance, dens_state, cells, pos, switch_3d):
    # calculates the local density of a given cell type. If there are many dormant cells this may bias the surrounding
    # cells. As we know the location of the cell in question it will be quicker to search the cells dictionary for
    # the cells with a location key within the given distance of the cell in question. If this is not good enough can
    # weight the contribution to the density based on distance.
    neighbours = []
    if switch_3d:
        for i in range(-distance, distance+1):
            for j in range(-distance, distance+1):
                for k in range(-distance, distance+1):
                    neighbours.append((pos[0] + i, pos[1] + j, pos[2] + k))

    else:
        for i in range(-distance, distance+1):
            for j in range(-distance, distance+1):
                neighbours.append((pos[0] + i, pos[1] + j))

    number_of_cells = 0
    for cell in neighbours:
        for state in dens_state:
            if cell in cells:
                if cells.get(cell)[0] == state:
                    number_of_cells = number_of_cells + 1

    #print(number_of_cells)

    if switch_3d:
        return number_of_cells/(2*distance + 1)**3
    else:
        return number_of_cells/(2*distance + 1)**2


def density_update_all(distance, cells, switch_3d, dens_state):
    # Update the density parameter of all cells based on the calculated local density of relevant states
    #print('density')
    for k in cells.keys():
        state = cells.get(k)[0]
        div_time = cells.get(k)[2]
        density = calc_local_density(distance, dens_state, cells, k, switch_3d)
        cells.update({k: [state, 0, div_time, density]})


def return_transition_dict(transition, density):
    # returns the transition dictionary. Density will be a number between 0 and 1. Multiplying transition_dict by the
    # original dictionary of transition values will change them correctly.

    if transition == 'g1':
        transition_dict = {0: [1], 1: [density, 1], 2: [1]}
    elif transition == 'g2':
        transition_dict = {0: [density], 1: [1, 1], 2: [1]}
    elif transition == 'd1':
        transition_dict = {0: [1], 1: [1, 1], 2: [density]}
    elif transition == 'd2':
        transition_dict = {0: [1], 1: [1, density], 2: [1]}

    return transition_dict