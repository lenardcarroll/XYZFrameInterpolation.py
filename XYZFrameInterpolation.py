#All modules to be used imported
import argparse
import pandas as pd
import numpy as np
import csv
parser = argparse.ArgumentParser()
#Cheat sheet of arguments to be used with the script
parser.add_argument("-inp", "--input", dest = "input", default = "input.xyz", help="Name of input file")
parser.add_argument("-out", "--output", dest = "output", default = "output.xyz", help="Name of output file")
parser.add_argument("-rep", "--replicas", dest = "Num_Replicas", default = "10", help="Enter the number of structures to interpolate")
parser.add_argument("-ter", "--interpolate", dest = "interpolate", default = "CSpline", help="Choose an interpolation method from: Linear, CSpline, NN")
parser.add_argument("-fr", "--frozen", dest = "frozenrange", default = "0", help="The range of atoms that are frozen. If none are frozen, specify 0 or ignore the flag. Must be written as 'X-Y'")
args = parser.parse_args()

#Here the necessary interpolation module is imported
interpolationmethod = args.interpolate
if interpolationmethod == 'CSpline':
    from scipy.interpolate import CubicSpline
else:
    from scipy.interpolate import interp1d

pd.options.mode.chained_assignment = None

#Here the number of atoms is determined by reading in the first value of the xyz file, the number of lines for the whole file is determined and the number of frames calculated using the aforementioned two variables
with open(args.input) as f:
    lines = f.read()
    first = lines.split('\n', 1)[0]
num_of_atoms = int(first)
num_of_lines = len(list(open(args.input)))
num_of_frames = int(num_of_lines/(num_of_atoms+2))
Num_Replicas = 1/int(args.Num_Replicas)

#Here the headers of the xyz file is removed
skiparray = []
for i in range(num_of_frames):
    k = num_of_atoms*i+2*i
    skiparray.append(k)
    skiparray.append(k+1)
df = pd.read_csv(args.input, skiprows=skiparray, names=['Atom', 'X', 'Y', 'Z'], sep="\s+" , engine='python')

#Here the dataframe is split into multiple dataframes, each making up a frame from the file
frames = [ df.iloc[i*num_of_atoms:(i+1)*num_of_atoms].copy() for i in range(num_of_frames+1) ]

for i in range(0,num_of_frames):
    frames[i] = frames[i].reset_index(drop=True)

frozenrange, firstval, lastval, dash = args.frozenrange, -2, -2, -2 # First value is the range of frozen atoms as specificied by the user, the other values are just to initialize the other terms

#Here the range of frozen atoms is properly stored
if '-' in frozenrange:
    dash = frozenrange.index('-')
    firstval = int(frozenrange[:dash])
    lastval = int(frozenrange[dash+1:])
elif int(frozenrange)==0:
    firstval = -2
    lastval = -1
else:
    firstval = int(frozenrange)
    lastval = firstval

Time = range(num_of_frames)
Atoms = []
X_values = []
Y_values = []
Z_values = []

#Place the x, y and z coordinates of each atom in every frame into lists
for i in range(num_of_atoms):
    X = []
    Y = []
    Z = []
    for j in range(num_of_frames):
        X.append(frames[j]['X'].iloc[i])
        Y.append(frames[j]['Y'].iloc[i])
        Z.append(frames[j]['Z'].iloc[i])
    X_values.append(X)
    Y_values.append(Y)
    Z_values.append(Z)
    Atoms.append(frames[0]['Atom'].iloc[i])

X = []
Y = []
Z = []
#The X,Y,Z_values lists are then used to interpolate the coordinates 
for i in range(num_of_atoms):
    #Here the number of frames to be interpolated is calculated
    TimeNew = np.arange(0,num_of_frames-1,Num_Replicas)
    if num_of_frames-1 not in TimeNew:
        TimeNew = np.append(TimeNew,num_of_frames-1)
    #If any atoms are fixed, it is not necessary to do an interpolation with it, since it will stay the same
    if i in range(firstval-1,lastval):
        X.append([X_values[i][0]]*len(TimeNew))
        Y.append([Y_values[i][0]]*len(TimeNew))
        Z.append([Z_values[i][0]]*len(TimeNew))
    #Here the interpolation is done on every atom in every frame based on the interpolation method chosen
    else:
        if interpolationmethod == 'CSpline':
            fX = CubicSpline(Time, X_values[i], bc_type='natural')
            fY = CubicSpline(Time, Y_values[i], bc_type='natural')
            fZ = CubicSpline(Time, Z_values[i], bc_type='natural')
            X.append(fX(TimeNew))
            Y.append(fY(TimeNew))
            Z.append(fZ(TimeNew))
        elif interpolationmethod == 'Linear':
            fX = interp1d(Time, X_values[i], kind='linear')
            fY = interp1d(Time, Y_values[i], kind='linear')
            fZ = interp1d(Time, Z_values[i], kind='linear')
            X.append(fX(TimeNew))
            Y.append(fY(TimeNew))
            Z.append(fZ(TimeNew))
        else:
            fX = interp1d(Time, X_values[i], kind='nearest')
            fY = interp1d(Time, Y_values[i], kind='nearest')
            fZ = interp1d(Time, Z_values[i], kind='nearest')
            X.append(fX(TimeNew))
            Y.append(fY(TimeNew))
            Z.append(fZ(TimeNew))

X_values = []
Y_values = []
Z_values = []

f = open(args.output,"w")
f.close()

#Here the interpolated frames are saved
for i in range(len(X[0])):
    X_values = []
    Y_values = []
    Z_values = []
    for j in range(num_of_atoms):
        X_values.append(X[j][i])
        Y_values.append(Y[j][i])
        Z_values.append(Z[j][i])
    dict = {'Atom':Atoms,
    'X': X_values,
    'Y': Y_values,
    'Z': Z_values}
    df2 = pd.DataFrame(dict)  
    with open(args.output,"a") as f:
        print(num_of_atoms, file=f)
        print("XYZ file created with XYZFrameInterpolation.py", file=f)
        for j in df2.index:
            print(df2['Atom'].iloc[j], df2['X'].iloc[j], df2['Y'].iloc[j], df2['Z'].iloc[j], file=f)
