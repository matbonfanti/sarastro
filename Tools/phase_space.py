#! /usr/bin/env python3

import numpy as np
import sys
import os
import matplotlib.pyplot as plt

# list of colors
Colors = ["blue", "green", 'red', 'magenta', 'cyan', 'black']

# define the degree of freedom index
if not sys.argv[1]:
    print("Missing first command line argument. Please provide DoF index!")
    sys.exit()
else:
    IndDof = int(sys.argv[1])
DofFileName = "centers_pq_" + "{0:03d}".format(IndDof) + ".dat"

# define the calculation directories
# when no command argument is given, use the working dir
if sys.argv[2:]:
    DirList = sys.argv[2:]
else:
    DirList = [os.getcwd()]

# check whether the directories contains the appropriate files
FileList = [DofFileName]
for d in DirList:
    for f in FileList:
        if not os.path.isfile(os.path.join(d, f)):
            print("file " + os.path.join(d, f) + " is missing!")
            sys.exit()

graphA = plt.subplot("111")

# cycle over directories
n = 0
for d in DirList:

    dofdata = np.loadtxt(os.path.join(d, DofFileName))
    dofdata = list(zip(*dofdata))

    time = dofdata[0]
    n = 0
    for colq, colp in zip(dofdata[1::2], dofdata[2::2]):
        n += 1
        graphA.plot(colq, colp, label="x " + str(n))

plt.legend()
plt.show()
