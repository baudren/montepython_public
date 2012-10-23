import pywlik
import numpy as nm

datadir = "FILLME"
#datadir = "/path/to/likelihood_v4p1"
ttmin = 2
ttmax = 32
temin = 2
temax = 32
use_gibbs = 0
use_lowlpol = 1

wmaplike = pywlik.wlik(datadir,ttmin,ttmax,temin,temax,use_gibbs,use_lowlpol)

cls = nm.loadtxt("wmap_7_2-32.cltest")

loglike = wmaplike(cls)

print "got %g expected %g"%(loglike,-845.483)