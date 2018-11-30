#! /usr/bin/env python

import numpy as np
import sys
import os
from scipy.interpolate import interp1d
from scipy.optimize import bisect

filtertime = 200.  # fs

# list of colors
Colors = [ "blue", "green", 'red', 'magenta', 'cyan', 'black' ]

# define the calculation directories
# when no command argument is given, use the working dir
if sys.argv[1:]:
   DirList = sys.argv[1:]
else:
   DirList = [ os.getcwd() ]
   
# define graphs
import matplotlib.pyplot as plt
graphA = plt.subplot("211")
graphB = plt.subplot("212")

# define function to compute the first zero of Re(Auto) which is used to shift the final spectrum
def FirstRootReAuto( TList, ReAutoList):
     # identify a suitable interval to look for the root
     firstNegative = 3+next(i for i in range(len(ReAutoList)) if ReAutoList[i] < 0.0)
     # interpolate with spline in the interval
     splineinterp = interp1d( TList[0:firstNegative], ReAutoList[0:firstNegative] )
     # find zero by bisection
     return bisect(splineinterp, 0.0, TList[firstNegative-1])

# cycle over directories
n = 0
for d in DirList:

   # check the presence of a file with the autocorrelation function:
   # "auto" MCTDH-style function with 1st:t,2nd:real(Auto),3rd:img(Auto)
   # "wf_prop" vMCG-style function with 2nd:t,5th:real(Auto),6th:img(Auto)
   # "expectations.dat" new vMCG-style function with 1st:t,4th:real(Auto),5th:img(Auto)

   #cat wf_prop | awk ' {print $2,$5,$6,sqrt($5**2+$6**2)}' > auto; PyPlot auto
   if os.path.isfile(os.path.join(d,"wf_prop")):
      # vMCG file is present
      f = os.path.join(d,"wf_prop")
      columns  = (1,4,5)
   elif os.path.isfile(os.path.join(d,"auto")):
      # MCTDH file is present
      f = os.path.join(d,"auto")
      columns  = (0,1,2)
   elif os.path.isfile(os.path.join(d,"expectations.dat")):
      # MCTDH file is present
      f = os.path.join(d,"expectations.dat")
      columns  = (0,3,4)
   else:
      print( "error reading directory "+str(d)+": autocorrelation file is missing!" )
      sys.exit()

   # read cross correlation 
   t,r,i = np.loadtxt(f, unpack=True, usecols=columns)
   if len(t)>10001:
      t = t[0:10001]
      r = r[0:10001]
      i = i[0:10001]
#   if len(t)<10001:
#      r = np.array( list(r) + [0.0]*(10001-len(r)) )
#      i = np.array( list(i) + [0.0]*(10001-len(i)) )
#      t = np.array(range(len(r)))*(t[1]-t[0])

   # Compute shift in energy
#   Eshift = np.pi/2.0/FirstRootReAuto( t/0.0241887, r ) 
 
   # scale auto correlation function by the norm
   Auto = list(r + 1j*i)

   # apply frequency filter
   for i in range(len(Auto)):
       Auto[i] = Auto[i]*np.exp(-t[i]/filtertime)

   tminus = [ -z for z in t[-1:0:-1]]
   Autominus = [ z.conjugate() for z in Auto[-1:0:-1]]
   t    = list(t[:])
   Auto = Auto[:]

   # compute fourier transform
   Spectra = np.fft.ifft(Auto[:]+Autominus[:])*t[1]/0.0241887
   Spectra *= len(Spectra)
   FreqSpacing = 2.0*np.pi/(t[1]/0.0241887*(len(t+tminus)))*219474.
   FreqGrid = np.array(range(len(Spectra)))*FreqSpacing 

   graphA.plot( FreqGrid, Spectra.real, label = d, ls="-", c=Colors[n] )
#   graphB.plot( tminus+t, np.imag(np.array(Autominus+Auto)), label = d, ls="-", c=Colors[n] )
   graphB.plot( tminus+t, abs(np.array(Autominus+Auto)), label = d, ls="-", c=Colors[n] )

   n += 1

#graphA.set_yscale( "log" )
# graphA.set_xlim( (3000.,17000.) )
graphB.set_xlim( (-10. ,1000.  ) )

graphA.legend()
plt.show()
