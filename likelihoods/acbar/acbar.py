# there is no specific likelihood code for this experiment, because it falls in the category of CMB experiments described in the "newdat" format. The class below inherits the properties of a general class "likelihood_newdat", which knows how to deal with all experiments in "newdat" format.
# however, the part related to nuisance parameters is specific to each experiment and is entirely defined below.

from likelihood_class import likelihood_newdat

#### remove these after defening function
import numpy as np
import math

class acbar(likelihood_newdat):

  def __init__(self,path,data,command_line=False):
    likelihood_newdat.__init__(self,path,data,command_line)
    self._need_nuisance_parameters(data,['A_SZ'])

    #####################
    # the following should become a function self._read_contamination_spectrum(self,['A_SZ'])
    # defined in generic class 'likelihood'

    # read spectrum contamination (so far, assumes only temperature contamination; will be trivial to generalize to polarization when such templates will become relevant)
    self.A_SZ_contamination=np.zeros(self.l_max+1,'float64')
    try:  
      for line in open(self.data_directory+self.A_SZ_file,'r'):
        l=int(float(line.split()[0]))
        if ((l >=2) and (l <= self.l_max)):
          self.A_SZ_contamination[l]=float(line.split()[1])/(l*(l+1.)/2./math.pi)
    except:
      print 'you must define a file name '+string(self)+'.A_SZ_file containing the contamination spectrum regulated by the nuisance parameter A_SZ'
      exit()

    # read renormalization factor
    # if it is not there, assume it is one, i.e. do not renormalize  
    try:
      self.A_SZ_contamination *= float(self.A_SZ_scale)
    except:
      pass

    # read central value of nuisance parameter
    # if it is not there, assume one by default
    try:
      exec self.A_SZ_bestfit
    except:
      self.A_SZ_bestfit=1.

    # read variance of nuisance parameter
    # if it are not there, assume flat prior (encoded through variance=0)
    try:
      exec self.A_SZ_variance
    except:
      self.A_SZ_variance=0.

    #############################

  def _loglkl(self,_cosmo,data):

    # get Cl's from CLASS
    cl = likelihood_newdat._get_cl(self,_cosmo)

    #####################
    # the following should become a function self._add_contamination_spectrum(self,['A_SZ'])
    # defined in generic class 'likelihood'
    # add contamination spectra multiplied by nuisance parameters
    # !!!!!!!!!!! This works only when nuisance parameter fixed
    for l in range(2,self.l_max):
      cl[0,l] += data.A_SZ*self.A_SZ_contamination[l]

    # get likelihood
    loglkl = likelihood_newdat._compute_loglkl(self,cl,_cosmo,data)

    #####################
    # the following should become a function self._add_nuisance_prior(self,['A_SZ'])
    # defined in generic class 'likelihood'
    # add prior on nuisance parameters
    if (self.A_SZ_variance>0):
      loglkl += -0.5*((data.A_SZ-self.A_SZ_bestfit)/self.A_SZ_variance)**2

    return loglkl

