## Doubling agent
Doubling agent is a simple 2D and 3D model for the
growth of a cancerous tumour consisting of a number of 
stem like cells. Cells exist on a Von Neumann grid, and 
grow, divide, die, and move based on the model and parameters
specified. If cells are surrounded on all four sides then they 
cannot move as there is no space for then to move into, 
and they cannot divide as there is no space for then to divide into. 

The parameters of the models are Poisson rate, e.g. the
rate at which a cell differentiates into a specific cell 
type. The units of these parameters is days<sup>-1</sup>. 
 
All simulations follow the same basic structure, the simulation parameter, 
number of dimension, and length of the simulation are specified by the user.

Firstly, the required time step is calculated and the simulation is initiated
by one cancerous cell. 
Each cell is stored in a dictionary, it’s key is its location in a Von Neumann grid. 
The values that do with it are [state, time], where state is the cell state 
time is the next “time” that something will happen to the cell.

"For bus arrivals consistent with a Poisson process, the expected wait time 
for a passenger is identical to the average interval between arrivals.”

Therefore a new “time” to each cell in each time step.
If the assigned time is within the timestep something happens to the cell.
Cells are “evaluated” in time order.
What happens to the cell depends on the state of the cell, another random 
number is generated and used to determine what happens to the cell once 
it has been decided that something will happen. 

If a cell divides the new cell can be created in one of the 4 grid points 
surrounding the location of the initial cell. 
All four locations are given equal probability, if the cell tries to chose 
a location which is taken that cell is removed from the list of options and 
the cell tries again. If all surrounding cells are occupied then the cell 
cannot divide.
If the cell dies it is deleted from the dictionary.
Cells are tracked for the specified number of time steps.
Information is saved to a binary file after a given number of time steps.
The program automatically decides how often to save data based on the 
length of the simulation so as not so slow the simulation down too much.

The saved data is later read and processed by a different script.

### tumour_sim.py:
In this code it is possible to implement 2 models, either with 3 cell states:
stem cell, proliferating cell and differentiated cell (referred to as the 'simple
model'); or with 4 cell states: stem cell, proliferating cell, differentiated cell,
and quiescent cells (referred to as the quiescent model).

#### Simple model rates:

- stem cell -> stem cell + stem cell
- stem cell -> stem cell + progenitor cell
- progenitor cell -> progenitor cell + progenitor cell
- progenitor cell -> differentiated cell
- differentiated cell -> dead

Parameters are [stem cell division rate, epsilon, progenitor cell division rate
, apoptosis rate]. Where epsilon is the probability that a stem cell divides symmetrically,
producing two stem cells.

#### Quiescent model rates:

- k1: stem cell -> stem cell + stem cell
- k2: stem cell -> stem cell + progenitor cell
- k3: stem cell -> dead
- k4: progenitor cell -> progenitor cell + progenitor cell
- k5: progenitor cell -> differentiated cell
- k6: progenitor cell -> quiescent cell
- k7: differentiated cell -> dead
- k8: quiescent cell -> stem cell

Parameters = [k1, k2, k3, k4, k5, k6 k7, k8]

In these simulations each cell is given a state in the form of a number:
0 - Stem cell
1 - Progenitor cell
2 - Differentiated cell
3 - Quiescent cell

Files are saved in a folder created by the simulation. The name of the folder
contains the name of the type of simulation run, the year, month, day hour and
minute of the folder creation and a few basic simulation parameters to ensure
that the data is easy to find. 

Data should later be analysed with 'plot_data.py'

The user must specify a number of parameters:

#### --params
Enter the parameters for either of the two model types supported by this simulation. 
The code will decide which model to run based ont the number of parameters entered.
Default parameters are [0.15, 0.15, 1.0, 0.48].

#### --motparams

In this model a "go or grow", motility is assumed. Cells are either in a state in 
which they can move, or a state in which they can divide. If this flag is not included 
it is assumed that cells cannot move. If the flag is included three parameters must 
also be added:
[switch to ‘go’ rate, switch to ‘grow’ rate, rate of movement]

#### --dim

The number of dimensions that the simulation should run in, the default is 2.

#### --length

Length of the simulation in days.

#### --repeats

Number of repeats of the same simulation (It's a stochastic procees, you'll get a 
different output each time).

#### --saving

This flag is mostly for code development, if it is included the output of the
simulation will not be saved.

#### --killcells

At a given point in the simulation kill all of one cell type instantaneously.
Requires two parameters, the state of the cells you want to kill and the 
day of death.

### steve_sim.py:

In this model cells do not have a distinct state, there are a large number of 
states and the model is linear, so that cells with a higher 
value of state are more likely to both move and divide. 
The 'state' of the cell is between 0 and 100 not 0 and 1 to make saving 
to binary files easier. 
When a cell divides both of the cells it divides into are in the same state.
To better model biological cell division and movement,  each cell has a 
period while it divides where it can’t move, reflecting perhaps a slightly 
more accurate go or grow model. This time length is user specified. 
 
The user specifies the granularity of the 'stemness' scale and the 
probability of a transition to move left or right on this scale (which is constant).
Cells to the right of this scale move, divide and die with greater probability. 

The simulation can currently run with the following parameters:

#### --params
Parameters, in order, are:
- rate of transition to a higher cell state
 - rate of transition to a lower cell state
 - division rate of cell with maximum diving potential
 - apoptosis_rate
 - the granularity of the 'stemness' scale
 - The amount of time for which a cell is paused after dividing (in days)
 
Default parameters are [1., 1., 1., 0.48, 10., 0.1].

#### --motparams

A single parameter discribing the movement rate of a cell with the maximum 
ability to move (cell state 100).

#### --dims

The number of dimensions for the simulation to run in, the default is 2. 

#### --length
Length of the simulation in days.

#### --repeats
Number of repeats of the simulation.

#### --saving
For development, including this flag stops the output from being saved.

The output of this simulation should be analysed with plot_steve.py. 

### state_sim.py:

In this simulation 2 or 3 cell states can be specified. The cell states are 0, 1, 
and 2. The parameters (--params) given to the simulation are [g1, g2, d1, d2, b, 
h, m, div_time]. 

- g2: [0] -> [1]
- g1: [1] -> [0]
- d1: [2] -> [1]
- d2: [1] -> [2]
- b: apoptosis rate (same for all cells)
- h: division rate of [2] cells
- m: division rate of [1] cells
- div_time: The time for which cells cannot move or divide after they have just divided.

To specify a model in which there are only 2 states, simply set g1 and g2
to 0. The initiating cell in this simulation is a dividing cell and so 
the simulation will still run. 

Currently all cells can move if movement is specified.

### density_sim.py

A slightly more advanced version of state_sim.py, which takes all of the same 
parameters. In addition the code calculates the local density for each cell. 
There is one additional input parameter, --density, which specifies the distance
from any given cell over which the density should be calculated, the cell states
which are included in the density calculation, and the cell state transition which 
this density will impact. 

The density can impact any of the 4 transitions: g1, g2, d1, d2. If the local density
of the given cell states is high then the given reaction will take place at the 
maximum possible rate (this is the rate specified in --parameters). However, if the 
local density of the specified cells is 0 then the rate of the transition will also 
be 0. 

This is very much an over simplification, but is a nice starting point. The local density is 
calculated on each timestep, and is used to update rate at which something will 
happen to each cell. 

### read_experimental_data.py:
This code takes a data file takes a text file as an input, the text file should be the output from
QuPath containing details of the location of each cell. The code then saves the 
data in the correct format to be analysed with plot_data.py, using the --exp 
flag. A 'fake' metadata file is also created to provide the correct parameters 
to plt_data.py. 

The data from QuPath specifies the location of cells in pixels, which is not
comparable to the simulation results. Due to this data is normalised by the 
average size of one cell (also included in the data from QuPath), so that the 
location of each cell is then specified in units of cell diameter, which is
comparable to simulation data.

The code requires the parameters:

#### -- datafile
The location of the datafile outputted by QuPath. The data should be tab separated
and contain the columns 'Centroid X px', 'Centroid Y px', 'Cell: Max caliper', and 'Cell: Min caliper'.


### plot_data:

This is the main plotting code. The spatial data from the simulations can be 
analysed in a number of ways:
- distance of each cell type from the centre of the tumour.
- the fourier transform of the location of the cells on a grid.
 - the density of each cell type, through the mean distance to the 
 eight (can be changed) nearest neighbours of the same cell type.
 - the lacunarity - a measure of the clustering of cells of the same
 type over a range of scales.
 
##### Density:

The code plots the mean distance of each cell to its eight nearest neighbours, 
normalised by the number of cells in the tumour. The code then fits to this
data with _y = a + b/(c - x)_, to allow a better numerical comparison between samples. 

##### Lacunarity:

Lacunarity is a measure of spatial data often used in ecology. A box is moved 
over the dataset and the number of cells of a given state within the box are 
counted. The lacunaity is then (standard deviation/ mean)<sup>2</sup>. By changing 
the box size the variation of the data can be studied over a number of scales.

The lacunarity calculation is only performed when --subset is specified, so 
that the box can be swept over a well defined rectangular region of the data. 

The code takes the parameters:

#### --folder
Location of the folder containing the data to be plotted.

#### --repeat
Which repeat of the data should be analysed.

#### --subset 
The subset of data for the lacunarity analysis. 2 parameters are require, a
lower and upper bound to the square or rectangle of data to be considered.
For example, if "--subset -20 20" is used, the lacunarity of the data would
be calculated on the grid from -20 to 20 in the x, y, and if the data is in 3d, 
z direction.

### plot_steve:

The same a plot_data.py but for the output from steve_sim.py. The code shows
plots of the data in terms of the cell state, which in this case can take many 
values between 0 and 100 (inclusive). The simulation then maps this into three
distinct cell states. For example, cells with a state between 0 and 30 are 
described as stem cells and given the cell state 0 (to match the cell states
in tumour_sim.py). 

The code then calls plot_data.py to analyse these three cell states in the 
usual way. 

The code requires the parameters:

#### --folder
Location of the folder containing the data to be plotted.

#### --repeat
Which repeat of the data should be analysed

#### --subset 
The subset of data for the lacunarity analysis
