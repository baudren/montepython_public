# there is no specific likelihood code for this experiment, because it falls in the category of CMB experiments described in the "newdat" format. The class below inherits the properties of a general class "likelihood_newdat", which knows how to deal with all experiments in "newdat" format.
# however, the part related to nuisance parameters is specific to each experiment and is entirely defined below.

from likelihood_class import likelihood_newdat

#### remove these after defening function
import numpy as np
import math

class spt(likelihood_newdat):

  def __init__(self,path,data,command_line=False):
    likelihood_newdat.__init__(self,path,data,command_line)
    self._need_nuisance_parameters(data,['SPT_SZ'])
    self._need_nuisance_parameters(data,['SPT_PS'])
    self._need_nuisance_parameters(data,['SPT_CL'])

    #####################
    # the following should become a function self._read_contamination_spectrum(self,['SPT_SZ'])

    # read spectrum contamination (so far, assumes only temperature contamination; will be trivial to generalize to polarization when such templates will become relevant)
    self.SPT_SZ_contamination=np.zeros(self.l_max+1,'float64')
    try:  
      for line in open(self.data_directory+self.SPT_SZ_file,'r'):
        l=int(float(line.split()[0]))
        if ((l >=2) and (l <= self.l_max)):
          self.SPT_SZ_contamination[l]=float(line.split()[1])/(l*(l+1.)/2./math.pi)
    except:
      print 'you must define a file name acbar.SPT_SZ_file containing the contamination spectrum regulated by the nuisance parameter SPT_SZ'
      exit()

    # read normalization factor  
    try:
      self.SPT_SZ_contamination *= float(self.SPT_SZ_scale)
    except:
      pass

    # read central value of nuisance parameter
    # if it is not there, assume one by default
    try:
      exec self.SPT_SZ_bestfit
    except:
      self.SPT_SZ_bestfit=1.

    # read variance of nuisance parameter
    # if it are not there, assume flat prior (encoded through variance=0)
    try:
      exec self.SPT_SZ_variance
    except:
      self.SPT_SZ_variance=0.

    #############################

    #####################
    # the following should become a function self._read_contamination_spectrum(self,['SPT_PS'])

    # read spectrum contamination (so far, assumes only temperature contamination; will be trivial to generalize to polarization when such templates will become relevant)
    self.SPT_PS_contamination=np.zeros(self.l_max+1,'float64')
    try:  
      for line in open(self.data_directory+self.SPT_PS_file,'r'):
        l=int(float(line.split()[0]))
        if ((l >=2) and (l <= self.l_max)):
          self.SPT_PS_contamination[l]=float(line.split()[1])/(l*(l+1.)/2./math.pi)
    except:
      print 'you must define a file name acbar.SPT_PS_file containing the contamination spectrum regulated by the nuisance parameter SPT_PS'

    # read normalization factor  
    try:
      self.SPT_PS_contamination *= float(self.SPT_PS_scale)
    except:
      pass

    # read central value of nuisance parameter
    # if it is not there, assume one by default
    try:
      exec self.SPT_PS_bestfit
    except:
      self.SPT_PS_bestfit=1.

    # read variance of nuisance parameter
    # if it are not there, assume flat prior (encoded through variance=0)
    try:
      exec self.SPT_PS_variance
    except:
      self.SPT_PS_variance=0.

    #############################

    #####################
    # the following should become a function self._read_contamination_spectrum(self,['SPT_CL'])

    # read spectrum contamination (so far, assumes only temperature contamination; will be trivial to generalize to polarization when such templates will become relevant)
    self.SPT_CL_contamination=np.zeros(self.l_max+1,'float64')
    try:  
      for line in open(self.data_directory+self.SPT_CL_file,'r'):
        l=int(float(line.split()[0]))
        if ((l >=2) and (l <= self.l_max)):
          self.SPT_CL_contamination[l]=float(line.split()[1])/(l*(l+1.)/2./math.pi)
    except:
      print 'you must define a file name acbar.SPT_CL_file containing the contamination spectrum regulated by the nuisance parameter SPT_CL'

    # read normalization factor  
    try:
      self.SPT_CL_contamination *= float(self.SPT_CL_scale)
    except:
      pass

    # read central value of nuisance parameter
    # if it is not there, assume one by default
    try:
      exec self.SPT_CL_bestfit
    except:
      self.SPT_CL_bestfit=1.

    # read variance of nuisance parameter
    # if it are not there, assume flat prior (encoded through variance=0)
    try:
      exec self.SPT_CL_variance
    except:
      self.SPT_CL_variance=0.

    #############################

  def _loglkl(self,_cosmo,data):

    # get Cl's from CLASS
    cl = likelihood_newdat._get_cl(self,_cosmo)

    # add contamination spectra multiplied by nuisance parameters
    for l in range(2,self.l_max):
      cl[0,l] += data.SPT_SZ*self.SPT_SZ_contamination[l]
      cl[0,l] += data.SPT_PS*self.SPT_PS_contamination[l]
      cl[0,l] += data.SPT_CL*self.SPT_CL_contamination[l]

    # get likelihood
    loglkl = likelihood_newdat._compute_loglkl(self,cl,_cosmo,data)

    # add prior on nuisance parameters
    if (self.SPT_SZ_variance>0):
      loglkl += -0.5*((data.SPT_SZ-self.SPT_SZ_bestfit)/self.SPT_SZ_variance)**2
    if (self.SPT_PS_variance>0):
      loglkl += -0.5*((data.SPT_PS-self.SPT_PS_bestfit)/self.SPT_PS_variance)**2
    if (self.SPT_CL_variance>0):
      loglkl += -0.5*((data.SPT_CL-self.SPT_CL_bestfit)/self.SPT_CL_variance)**2

    return loglkl
