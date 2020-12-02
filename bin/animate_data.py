from doubling_agent.image_analysis_functions import *
from doubling_agent.common_functions import *
import argparse


class StoreAsArray(argparse._StoreAction):
    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)


parser = argparse.ArgumentParser(description='Process input parameters for plotting')
parser.add_argument('--folder', type=str, help='Folder containing the data you wish to plot.', required=True)
parser.add_argument('--repeat', action=StoreAsArray, type=int, nargs='+', help="Which repeats to plot"
                                            "can take multiple values.",
                    default=np.array([0]))
args = parser.parse_args()

with open(args.folder + '/metadata.txt') as f:
    meta_data = f.readlines()

if int(meta_data[5].strip().split()[2]) == 2:
    switch_3d = False
else:
    switch_3d = True

if meta_data[1].split()[3] == 'None':
    from doubling_agent.basic_functions import *
    mot = False
else:
    from doubling_agent.motility_functions import *
    mot = True

repeats = int(meta_data[4].strip().split()[4])

time_step = float(meta_data[2].strip().split()[2])

for r in args.repeat:
    data = read_from_file(args.folder + '/repeat_' + str(r), switch_3d)
    if os.path.isdir(args.folder + '/animation_repeat_' + str(r)):
        print('Animation already exist for repeat ', str(r))
    else:
        ensure_dir(args.folder + '/animation_repeat_' + str(r))
        # data comes back as a list of touples. count, x,y, state

        if mot:
            if switch_3d:
                data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'z', 'state', 'mot'])
            else:
                data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'state', 'mot'])
        else:
            if switch_3d:
                data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'z', 'state'])
            else:
                data_df = pd.DataFrame(np.array(data), columns=['count', 'x', 'y', 'state'])

        name = args.folder + '/animation_repeat_' + str(r)
        animate(data_df, r, name)
