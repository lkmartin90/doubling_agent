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
        # 0 motility state is proliferating, 1 is moving
        self.rate = motility_rate


class ParametersSteve:
    def __init__(self, right_rate, left_rate, division_rate, apoptosis_rate, granularity, div_time):
        self.r = right_rate
        self.d = division_rate
        self.l = left_rate
        self.a = apoptosis_rate
        self.g = int(granularity)
        self.t = div_time


def cancer_seed_single(cells, switch_3d):
    # created initial cancer stem cell at [0,0]
    if switch_3d:
        cells.update({(0, 0, 0): [100, 0, 0]})
    else:
        cells.update({(0,0): [100, 0, 0]})


def timing_update_all(cells, params, mot_params, time_step):
    # update second entry in dict to give a timing based on the first entry, the state
    # time is log(1/rand_no)/rate
    # Now want to account for fact that cells can either move or divide, depending on where on the sliding scale they sit.

    for k in cells.keys():
        state = cells.get(k)[0]
        div_time = cells.get(k)[2]
        div = params.d  # division rate
        move = mot_params.rate  # move rate
        death = params.a  # death rate

        # start by assuming that this is a linear scale, the larger the value of state the more likely everything
        # is to happen

        if div_time > 0:
            div_time_new = div_time - time_step
            rate = params.r + params.l
        elif div_time < 0:
            div_time_new = 0
            rate = (div + move + death)*state/100 + params.r + params.l
        else:
            div_time_new = div_time
            rate = (div + move + death)*state/100 + params.r + params.l

        cells.update({k: [state, np.log(1/random())/rate, div_time_new]})


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
        cells.update({new_location: [state, 0, div_time]})


def update_cell_steve(cells, pos, params, switch_3d, mot_params):
    # updates a given cell based on the current state of that cell
    # pos is string describing position
    # time is from random number generator giving time of interaction
    # cells is dict describing all cells in the tumour

    state = cells.get(pos)[0]
    div_time = cells.get(pos)[2]

    if switch_3d:
        daughter = choose_new_pos_3d(pos, cells)
    else:
        daughter = choose_new_pos(pos, cells)

    # There are a number of options, cells can either move up or down the state scale, they can divide, die, or move

    if div_time != 0:
        # if div_time is non zero there are only  options, move left or right on the sliding scale
        r_num = random()
        tot = (params.l + params.r)
        if r_num < params.l/tot:
            # left on state scale
            if state - params.g < 0:
                cells.update({pos: [0, 0, div_time]})
            else:
                cells.update({pos: [state - params.g, 0, div_time]})

        else:
            # move right on state scale
            if state - params.g > 1:
                cells.update({pos: [0, 0, div_time]})
            else:
                cells.update({pos: [state + params.g, 0, div_time]})

    else:
        # If div_time is 0 then we're ok to both move and divide so there are more options. The timing has been updated
        # to reflect this. Timing is updaed before the cells are updates
            r_num = random()
            tot = (params.d + mot_params.rate + params.a)*state/100 + params.r + params.l
            if r_num < params.l/tot:
                # left on state scale
                if state - params.g < 0:
                    cells.update({pos: [0, 0, div_time]})
                else:
                    cells.update({pos: [state - params.g, 0, div_time]})

            elif r_num < (params.l + params.r)/tot:
                # move right on state scale
                if state + params.g > 100:
                    cells.update({pos: [100, 0, div_time]})
                else:
                    cells.update({pos: [state + params.g, 0, div_time]})

            elif r_num < (params.l + params.r + params.d*(state/100))/tot:
                # divide
                if daughter != 0:
                    cells.update({daughter: [state, 0, params.t]})
                    cells.update({pos: [state, 0, params.t]})

            elif r_num < (params.l + params.r + (params.d + mot_params.rate)*(state/100))/tot:
                # move
                move_cell(cells, pos, state, switch_3d, div_time)

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
    max_rate = params.a + params.l + params.r + params.d + mot_params.rate
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
    directory = file_path
    if not os.path.exists(directory):
        print('...making', directory)
        os.makedirs(directory)