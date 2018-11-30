#! /usr/bin/env python

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib

# list of colors
Colors = [ "#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4",
           "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff",
           "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
           "#000075", "#a9a9a9", "#ffffff", "#000000" ]

# define the calculation directories
# when no command argument is given, use the working dir
if sys.argv[1:]:
   DirList = sys.argv[1:]
else:
   DirList = [ os.getcwd() ]
   
# check whether the directories contains the appropriate files
FileList = ["integrator.log", "expectations.dat"]
for d in DirList:
  for f in FileList:
     if not os.path.isfile(os.path.join(d,f)): 
        print( "file "+os.path.join(d,f)+" is missing!" )
        sys.exit()

# define graphs
graphA = plt.subplot("211")
graphB = plt.subplot("212")
graphC = graphB.twinx()

# cycle over directories
n = 0
for d in DirList:

   # read integration step
   integdata = np.loadtxt(os.path.join(d,"integrator.log"), dtype="S", usecols=(3,5), skiprows=2)[::10]
   integdata = zip(*integdata)
   stepsize_y = [ float(t) for t in  integdata[0] ]
   stepsize_x = [ float(t) for t in  integdata[1] ]
   average = sum(stepsize_y)/len(stepsize_y) 

   graphA.plot( stepsize_x, stepsize_y, ls="--", c=Colors[n] )
   graphA.plot( (stepsize_x[0],stepsize_x[-1]), (average,average), lw=2, c=Colors[n], label = d )

   # read energy over time
   expectdata = np.loadtxt(os.path.join(d,"expectations.dat"), dtype="float", usecols=(0,1,2)  ) 
   if isinstance(expectdata[0], np.ndarray):
      expectdata = zip(*expectdata)
      norm_x = [ float(t) for t in  expectdata[0] ]
      norm_y = [ (float(t)-1.0)*1000. for t in  expectdata[1] ]
      energy_x = [ float(t) for t in  expectdata[0] ]
      energy_y = [ (float(t)-float(expectdata[2][0]))*1000. for t in  expectdata[2] ]
   else:
      norm_x = [ float(t) for t in  [ expectdata[0] ] ]
      norm_y = [ (float(t)-1.0)*1000. for t in [ expectdata[1]] ]
      energy_x = [ float(t) for t in  [expectdata[0]] ]
      energy_y = [ (float(t)-float(expectdata[2]))*1000. for t in  [expectdata[2]] ]


   graphC.plot( norm_x, norm_y, label = d, c=Colors[n], ls="--"  )
   graphB.plot( energy_x, energy_y, label = d, c=Colors[n], ls=":" )

   n += 1

graphA.set_yscale( "log" )
graphB.ticklabel_format(style='sci', axis='y', scilimits=(-3,+3))
graphC.ticklabel_format(style='sci', axis='y', scilimits=(-3,+3))
graphB.set_xlabel( "t / fs" )

graphA.set_ylabel( "$\Delta$t / fs" )
graphB.set_ylabel( "$\Delta$ E / meV" )
graphC.set_ylabel( "$\Delta$ norm $\cdot$ 10$^3$" )

graphA.legend()

energy_line = matplotlib.lines.Line2D([], [], color='black', marker=None, ls=":", label='energy')
norm_line = matplotlib.lines.Line2D([], [], color='black', marker=None, ls="--", label='norm')
graphB.legend( handles=[energy_line, norm_line] )

plt.show()

plt.show()
