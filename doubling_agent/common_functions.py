import struct
import operator
import os
import numpy as np

def ensure_dir(file_path):
    directory = file_path
    if not os.path.exists(directory):
        print('...making', directory)
        os.makedirs(directory)


def write_to_file(data, file_name):
    # write data to binary file in the form: time step, x, y, state
    s = struct.pack('i' * len(data), *data)
    with open(file_name, 'ab') as f:
        f.write(s)

def write_exp_to_file(data, file_name):
    # write data to binary file in the form: time step, x, y, state
    s = struct.pack('f' * len(data), *data)
    with open(file_name, 'ab') as f:
        f.write(s)


def Extract(lst):
    return [item[0] for item in lst]


def kill_cells(cells, cell_type):
    # find coordinates of all cells of the specified type
    to_kill = []
    for name, type in cells.items():
        if type[0] == cell_type:
            to_kill.append(name)

    # delete these cells
    for k in to_kill:
        cells.pop(k, None)
