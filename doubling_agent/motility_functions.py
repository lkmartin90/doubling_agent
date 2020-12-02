from random import random
from random import choice
import numpy as np
import plotly.express as px
import struct
import operator

###
# Broadly the same as "basic_functions.py" but updated to include motility
# intentionally trying to keep them separate so as not to slow down the basic version
###


class MotilityParameters:
    def __init__(self, switch_to_m_rate, switch_to_p_rate, motility_rate):
        # 0 motility state is proliferating, 1 is moving
        self.m = switch_to_m_rate
        self.p = switch_to_p_rate
        self.rate = motility_rate
        self.dict = {0: switch_to_m_rate, 1: switch_to_p_rate + motility_rate}


class ParametersBasic:
    def __init__(self, s_division_rate, epsilon, p_division_rate, apoptosis_rate):
        self.s = s_division_rate
        self.p = p_division_rate
        self.e = epsilon
        self.a = apoptosis_rate
        self.dict = {0: s_division_rate, 1: p_division_rate, 2: apoptosis_rate}
        self.death = {0: 0, 1: 0, 2: apoptosis_rate}


class ParametersQuiescent:
    def __init__(self, k1, k2, k3, k4, k5, k6, k7, k8):
        # s>s+s :k1, s>s+p:k2, s>dead:k3, p>p+p:k4, p>dead:k5, p>Q:k6, D>dead:k7, Q>s:k8
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5
        self.k6 = k6
        self.k7 = k7
        self.k8 = k8
        # rate of something happening for each state
        # note the slight change in notation, 0 is stem cell, 1 progenitor, 2 differentiated and 3 quiescent
        self.dict = {0: k1+k2+k3, 1: k4+k5+k6, 2: k7, 3: k8}
        self.death = {0: k3, 1: k5, 2: k7, 3: 0}


def cancer_seed_single(cells, switch_3d):
    # created initial cancer stem cell at [0,0]
    if switch_3d:
        cells.update({(0, 0, 0): [0, 0, 0]})
    else:
        cells.update({(0,0): [0, 0, 0]})


def cancer_seed_single_quiescent(cells):
    # created initial cancer cell (differentiated) at [0,0]
    cells.update({(0,0): [3, 0, 0]})


def cancer_seed_single_progen(cells):
    # created initial cancer cell (differentiated) at [0,0]
    cells.update({(0,0): [1, 0, 0]})


def timing_update_all(cells, params, mot_params):
    # update second entry in dict to give a timing based on the first entry, the state
    # time is log(1/rand_no)/rate
    # Now want to account for fact that cells can either change motility state or move or divide
    # options:
    #   motility state 1: can either move, change motility state, or die
    #   motility state 0: can either change motility state or go though all division choices
    #   including death (already written)

    for k in cells.keys():
        state = cells.get(k)[0]
        mot = cells.get(k)[2]
        div = params.dict[state]  # division or death rate for motility 0
        m_or_c = mot_params.dict[mot]  # move or change rate
        mot_death = params.death[state]  # death rate for motility state 1
        if mot == 0:
            rate = div+m_or_c
        else:
            rate = m_or_c + mot_death
        cells.update({k: [state, np.log(1/random())/rate, mot]})


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
    # choses a new position by identifying all the free spaces first and then assigning them all equal probability
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


def move_cell(cells, pos, state, switch_3d):
    # moves the cell
    if switch_3d:
        new_location = choose_new_pos_3d(pos, cells)
    else:
        new_location = choose_new_pos(pos, cells)

    if new_location != 0:
        del cells[pos]
        cells.update({new_location: [state, 0, 1]})


def update_cell_basic(cells, pos, params, switch_3d, mot_params):
    # updates a given cell based on the current state of that cell
    # pos is string describing position
    # time is from random number generator giving time of interaction
    # cells is dict describing all cells in the tumour

    state = cells.get(pos)[0]
    mot = cells.get(pos)[2]
    # Once motility is included first thing is to make a decision on whether the cell moves, divides, or switches
    # to a different motility state. Need to check that an appropriate time step is being used still.

    mot_check = random()
    if mot == 1:
        # Can move, cell can either switch motility state or move or die
        if mot_check < mot_params.p/(mot_params.dict.get(mot) + params.death.get(state)):
            # then the motilty state changes
            cells.update({pos: [state, 0, abs(mot-1)]})

        elif mot_check < mot_params.dict.get(mot)/(mot_params.dict.get(mot) + params.death.get(state)):
            # The cell moves
            move_cell(cells, pos, state, switch_3d)

        else:
            # cell death
            del cells[pos]

    else:
        # No motility, can either switch state or go to division decisions
        if mot_check < mot_params.m/(mot_params.dict.get(mot) + params.dict.get(state)):
            # then the motilty state changes
            cells.update({pos: [state, 0, abs(mot - 1)]})

        # The cell divides or dies, can ust move on to that section as we have already conditioned on the
        # probability of it happening
        else:
            if switch_3d:
                daughter = choose_new_pos_3d(pos, cells)
            else:
                daughter = choose_new_pos(pos, cells)

            if state == 0:
                # if it's a stem cell there are 2 possibilities, S > S + S, S > S + P
                # generate random number to determine fate, compare to epsilon
                r_num = random()
                if r_num < params.e:
                    # divide > S + S
                    if daughter != 0:
                        cells.update({daughter: [0, 0, 0]})
                else:
                    # divide > S + P
                    if daughter != 0:
                        cells.update({daughter: [1, 0, 0]})

            elif state == 1:
                # if it's a progentior cell there are 2 possibilities, P > P + P, P > D
                # generate random number to determine fate, start by assuming each happens with equal chance
                r_num = random()
                if r_num < 0.5:
                    # P > P + P
                    if daughter != 0:
                        cells.update({daughter: [1, 0, 0]})
                else:
                    # P > D
                    cells.update({pos: [2, 0, 0]})

            else:
                # If it is differentiated cell the only possible state change is death
                del cells[pos]


def update_cell_quiescent(cells, pos, params, switch_3d, mot_params):
    # updates a given cell based on the current state of that cell
    # pos is string describing position
    # time is from random number generator giving time of interaction
    # cells is dict describing all cells in the tumour

    state = cells.get(pos)[0]
    mot = cells.get(pos)[2]
    # Once motility is included first thing is to make a decision on whether the cell moves, divides, or switches
    # to a different motility state. Need to check that an appropriate time step is being used still.

    mot_check = random()
    if mot == 1:
        # Can move, cell can either switch motility state or move or die
        if mot_check < mot_params.p / (mot_params.dict.get(mot) + params.death.get(state)):
            # then the motilty state changes
            cells.update({pos: [state, 0, abs(mot - 1)]})

        elif mot_check < mot_params.dict.get(mot) / (mot_params.dict.get(mot) + params.death.get(state)):
            # The cell moves
            move_cell(cells, pos, state, switch_3d)

        else:
            # cell death
            del cells[pos]

    else:
        # No motility, can either switch state or go to division decisions
        if mot_check < mot_params.m / (mot_params.dict.get(mot) + params.dict.get(state)):
            # then the motilty state changes
            cells.update({pos: [state, 0, abs(mot - 1)]})

        else:

            if switch_3d:
                daughter = choose_new_pos_3d(pos, cells)
            else:
                daughter = choose_new_pos(pos, cells)

            if state == 0:
                # if it's a stem cell there are 3 possibilities, S > S + S, S > S + P and S > dead
                # generate random number to determine fate
                r_num = random()
                if r_num < params.k1/params.dict.get(0):
                    # divide > S + S
                    if daughter != 0:
                        cells.update({daughter: [0, 0, 0]})
                elif r_num < (params.k1+params.k2)/params.dict.get(0):
                    # divide > S + P
                    if daughter != 0:
                        cells.update({daughter: [1, 0, 0]})
                else:
                    # die
                    del cells[pos]

            elif state == 1:
                # if it's a progentior cell there are 3 possibilities, P > P + P, P > D, P > Q
                # generate random number to determine fate
                r_num = random()
                if r_num < params.k4/params.dict.get(1):
                    # P > P + P
                    if daughter != 0:
                        cells.update({daughter: [1, 0, 0]})
                elif r_num < (params.k4+params.k5)/params.dict.get(1):
                    # P > D
                    cells.update({pos: [2, 0, 0]})
                else:
                    # P > Q
                    cells.update({pos: [3, 0, 0]})

            elif state == 2:
                # If it is differentiated cell the only possible state change is death
                del cells[pos]

            else:
                # If its Quiescent the only possible fate is to return to a stem cell
                cells.update({pos: [0, 0, 0]})


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
        struct_fmt = '=iiiiii'  # 6 ints
    else:
        struct_fmt = '=iiiii'  # 5 ints
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
    max_rate = max(params.dict.items(), key=operator.itemgetter(1))[1] + \
               max(mot_params.dict.items(), key=operator.itemgetter(1))[1]
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
