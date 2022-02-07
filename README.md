# XYZFrameInterpolation.py

This python script follows on from my LinearXYZFrameInterpolation.py (https://github.com/lenardcarroll/LinearXYZFrameInterpolation.py) script, where frames were purely linearly interpolated. What makes this python script different is that unlike the aforementioned script where the positions of atoms were linearly interpolated, here the individual x,y and z coordinates of atoms are interpolated using either a linear, cubic spline or nearest neighbour interpolation method.

The script operates generally as follows:

```
python -i XYZFrameInterpolation.py -inp <INPUT_FILE_NAME> -out <OUTPUT_FILE_NAME> -rep <NUMER_OF_INTERPOLATED_FRAMES> -ter <INTERPOLATION_METHOD> -fr <RANGE_OF_FIXED_ATOMS>
```

where ```-rep <NUMBER_OF_INTERPOLATED_FRAMES>``` corresponds to the number of frames to be interpolated in-between original frames from the .xyz file, ```-ter <INTERPOLATION_METHOD>``` must be either ```NN``` (for nearest neighbor), ```Linear``` for linear interpolation or ```CSpline``` for Cubic Spline interpolation, and ```-fr <RANGE_OF_FIXED_ATOMS>``` is the range of atoms that are fixed in the simulation/computation, if none exists, ignore the flag. 

Make sure to pip install pandas, argparse, csv and scipy.
