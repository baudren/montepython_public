import os,sys
from collections import OrderedDict as od
import numpy as np
import math

# Definition of classes for the likelihoods to inherit, for every different
# types of expected likelihoods. Included so far, models for a Clik, simple
# prior, and newdat likelihoods.

# General class that will only define the store_lkl function, common for every
# single likelihood. It will copy the content of self.path, copied from
# initialization
class likelihood():

  def __init__(self,path,data,command_line=False):

    self.name = self.__class__.__name__
    self.folder = os.path.abspath(data.path['MontePython'])+'/../likelihoods/'+self.name+'/'
    self._read_from_file(path,data)

    if ((command_line is not False) and (os.path.exists(command_line.folder+'/log.dat') is False)):
      self._store_lkl_params(command_line)

    

  def _loglkl(self,_cosmo,data):
    raise NotImplementedError('Must implement method _loglkl() in your likelihood')

  def _store_lkl_params(self,command_line):
    log = open(command_line.folder+'log.param','a')
    tolog = open(self.path,'r')
    log.write("\n\n#-----Likelihood-{0}-----\n".format(self.name))
    for line in tolog:
      log.write(line)
    tolog.seek(0)
    tolog.close()
    log.close()

  def _read_from_file(self,path,data):
    self.path = path
    if os.path.isfile(path):
      data_file = open(path,'r')
      for line in data_file:
	if line.find('#')==-1:
	  if line.find(self.name+'.')!=-1:
	    exec(line.replace(self.name+'.','self.'))
      data_file.seek(0)
      data_file.close()

  def _get_cl(self,_cosmo):

    # get C_l^XX from CLASS
    cl = _cosmo.lensed_cl()

    # convert dimensionless C_l's to C_l in muK**2
    T = _cosmo._T_cmb()
    for key in cl.iterkeys():
      cl[key] *= (T*1.e6)**2

    return cl

  def _need_Class_args(self,data,dictionary):
    for key,value in dictionary.iteritems():
      try :
	data.Class_args[key]
	try:
	  float(data.Class_args[key])
	  num_flag = True
	except:
	  num_flag = False
      except KeyError:
	try:
	  float(value)
	  num_flag = True
	  data.Class_args[key] = 0
	except:
	  num_flag = False
	  data.Class_args[key] = ''
      if num_flag is False:
	if data.Class_args[key].find(value)==-1:
	  data.Class_args[key] += ' '+value+' '
      else:
	if float(data.Class_args[key])<value:
	  data.Class_args[key] = value


  def _read_contamination_spectra(self,data):
    failure = False
    for nuisance in self.use_nuisance:
      if nuisance not in data.nuisance_params.iterkeys():
        try:
          exec "data.%s" % nuisance
        except:
          print nuisance+' must be defined, either fixed or varying, for {0} likelihood'.format(self.name)
	  failure=True
    if failure:
      exit()

    for nuisance in self.use_nuisance:
      # read spectrum contamination (so far, assumes only temperature contamination; will be trivial to generalize to polarization when such templates will become relevant)
      exec "self.%s_contamination=np.zeros(self.l_max+1,'float64')" % nuisance
      try:  
        exec "File = open(self.data_directory+self.%s_file,'r')" % nuisance
      except:
        print 'you must define a file name '+self.name+'.'+nuisance+'_file containing the contamination spectrum regulated by the nuisance parameter '+nuisance
        exit()
   
      for line in File:
        l=int(float(line.split()[0]))
        if ((l >=2) and (l <= self.l_max)):
          exec "self.%s_contamination[l]=float(line.split()[1])/(l*(l+1.)/2./math.pi)" % nuisance

      # read renormalization factor
      # if it is not there, assume it is one, i.e. do not renormalize  
      try:
        exec "self.%s_contamination *= float(self.%s_scale)" % (nuisance,nuisance)
      except:
        pass

      # read central value of nuisance parameter
      # if it is not there, assume one by default
      try:
        exec "self.%s_prior_center" % nuisance
      except:
        exec "self.%s_prior_center=1." % nuisance

      # read variance of nuisance parameter
      # if it are not there, assume flat prior (encoded through variance=0)
      try:
        exec "self.%s_prior_variance" % nuisance
      except:
        exec "self.%s_prior_variance=0." % nuisance

  def _add_contamination_spectra(self,cl,data):

    # first, test whether the nuisance param is fixed or varying.
    for nuisance in self.use_nuisance:
      if nuisance in data.nuisance_param_names:
        nuisance_value = float(data.vector[np.where(data.param_names == nuisance)[0][0]])
      else:
        exec "nuisance_value = data.%s" % nuisance

      # add contamination spectra multiplied by nuisance parameters
      for l in range(2,self.l_max):
        exec "cl['tt'][l] += nuisance_value*self.%s_contamination[l]" % nuisance

    return cl

  def _add_nuisance_prior(self,loglkl,data):

    # first, test whether the nuisance param is fixed or varying.
    for nuisance in self.use_nuisance:
      if nuisance in data.nuisance_param_names:
        nuisance_value = float(data.vector[np.where(data.param_names == nuisance)[0][0]])
      else:
        exec "nuisance_value = data.%s" % nuisance

      # add prior on nuisance parameters
      exec "if (self.%s_prior_variance>0): loglkl += -0.5*((nuisance_value-self.%s_prior_center)/self.%s_prior_variance)**2" % (nuisance,nuisance,nuisance)

    return loglkl

# Likelihood type for prior
class likelihood_prior(likelihood):

  def _loglkl(self):
    raise NotImplementedError('Must implement method __loglkl() in your likelihood')




# Generic type for newdat likelihood (spt, boomerang, etc)
class likelihood_newdat(likelihood):

  def __init__(self,path,data,command_line=False):

    likelihood.__init__(self,path,data,command_line)

    self._need_Class_args(data,{'lensing':'yes', 'output':'tCl lCl pCl'})
      
    # open .newdat file
    newdatfile=open(self.data_directory+self.file,'r')

    # find beginning of window functions file names
    window_name=newdatfile.readline().strip('\n').replace(' ','')
    
    # initialize list of fist and last band for each type
    band_num = np.zeros(6,'int')
    band_min = np.zeros(6,'int')
    band_max = np.zeros(6,'int')

    # read number of bands for each of the six types TT, EE, BB, EB, TE, TB
    line = newdatfile.readline()
    for i in range(6):
      band_num[i] = int(line.split()[i])

    # read string equal to 'BAND_SELECTION' or not
    line = str(newdatfile.readline()).strip('\n').replace(' ','')

    # if yes, read 6 lines containing 'min, max'
    if (line=='BAND_SELECTION'):
      for i in range(6):
        line = newdatfile.readline()
        band_min[i] = int(line.split()[0])   
        band_max[i] = int(line.split()[1])

    # if no, set min to 1 and max to band_num (=use all bands)   
    else:
      band_min=[1 for i in range(6)]
      band_max=band_num

    # read line defining calibration uncertainty
    # contains: flag (=0 or 1), calib, calib_uncertainty
    line = newdatfile.readline()
    calib=float(line.split()[1])
    if (int(line.split()[0])==0):
      self.calib_uncertainty=0
    else:
      self.calib_uncertainty=float(line.split()[2])
 
    # read line defining beam uncertainty
    # contains: flag (=0, 1 or 2), beam_width, beam_sigma
    line = newdatfile.readline()
    beam_type = int(line.split()[0])
    if (beam_type > 0):
      self.has_beam_uncertainty = True
    else:
      self.has_beam_uncertainty = False
    beam_width = float(line.split()[1])
    beam_sigma = float(line.split()[2])
    
    # read flag (= 0, 1 or 2) for lognormal distributions and xfactors
    line = newdatfile.readline()
    likelihood_type = int(line.split()[0])
    if (likelihood_type > 0):
      self.has_xfactors = True
    else:
      self.has_xfactors = False

    # declare array of quantitites describing each point of measurement
    # size yet unknown, it will be found later and stored as self.num_points
    self.obs=np.array([],'float64')
    self.var=np.array([],'float64')
    self.beam_error=np.array([],'float64')
    self.has_xfactor=np.array([],'bool')
    self.xfactor=np.array([],'float64')

    # temporary array to know which bands are actually used
    used_index=np.array([],'int')

    index=-1

    # scan the lines describing each point of measurement
    for cltype in range(6):
      if (int(band_num[cltype]) != 0):
        # read name (but do not use it)
        newdatfile.readline()
        for band in range(int(band_num[cltype])):
          # read one line corresponding to one measurement
          line = newdatfile.readline()
          index += 1

          # if we wish to actually use this measurement
          if ((band >= band_min[cltype]-1) and (band <= band_max[cltype]-1)):

            used_index=np.append(used_index,index)

            self.obs=np.append(self.obs,float(line.split()[1])*calib**2)

            self.var=np.append(self.var,(0.5*(float(line.split()[2])+float(line.split()[3]))*calib**2)**2)

            self.xfactor=np.append(self.xfactor,float(line.split()[4])*calib**2)

            if ((likelihood_type == 0) or ((likelihood_type == 2) and (int(line.split()[7])==0))):
              self.has_xfactor=np.append(self.has_xfactor,[False])
            if ((likelihood_type == 1) or ((likelihood_type == 2) and (int(line.split()[7])==1))):
              self.has_xfactor=np.append(self.has_xfactor,[True])

            if (beam_type == 0):
              self.beam_error=np.append(self.beam_error,0.)
            if (beam_type == 1):
              l_mid=float(line.split()[5])+0.5*(float(line.split()[5])+float(line.split()[6]))
              self.beam_error=np.append(self.beam_error,abs(math.exp(-l_mid*(l_mid+1)*1.526e-8*2.*beam_sigma*beam_width)-1.))
            if (beam_type == 2):
              if (likelihood_type == 2):
                self.beam_error=np.append(self.beam_error,float(line.split()[8]))
              else:
                self.beam_error=np.append(self.beam_error,float(line.split()[7]))

        # now, skip and unused part of the file (with sub-correlation matrices)
        for band in range(int(band_num[cltype])):      
          newdatfile.readline()

    # number of points that we will actually use
    self.num_points=np.shape(self.obs)[0]

    # total number of points, including unused ones
    full_num_points=index+1

    # read full correlation matrix
    full_covmat=np.zeros((full_num_points,full_num_points),'float64')
    for point in range(full_num_points):
      full_covmat[point]=newdatfile.readline().split()

    # extract smaller correlation matrix for points actually used
    covmat=np.zeros((self.num_points,self.num_points),'float64')
    for point in range(self.num_points):
      covmat[point]=full_covmat[used_index[point],used_index]

    # recalibrate this correlation matrix
    covmat *= calib**4

    # redefine the correlation matrix, the observed points and their variance in case of lognormal likelihood
    if (self.has_xfactors):

      for i in range(self.num_points):

        for j in range(self.num_points):
          if (self.has_xfactor[i]):
            covmat[i,j] /= (self.obs[i]+self.xfactor[i])
          if (self.has_xfactor[j]):
            covmat[i,j] /= (self.obs[j]+self.xfactor[j])
            
      for i in range(self.num_points):
        if (self.has_xfactor[i]):
          self.var[i] /= (self.obs[i]+self.xfactor[i])**2
          self.obs[i] = math.log(self.obs[i]+self.xfactor[i])

    # invert correlation matrix
    self.inv_covmat=np.linalg.inv(covmat)

    # read window function files a first time, only for finding the smallest and largest l's for each point
    self.win_min=np.zeros(self.num_points,'int')
    self.win_max=np.zeros(self.num_points,'int')
    for point in range(self.num_points):
      for line in open(self.data_directory+'windows/'+window_name+str(used_index[point]+1),'r'):
        if any([float(line.split()[i]) != 0. for i in range(1,len(line.split()))]):
          if (self.win_min[point]==0):
            self.win_min[point]=int(line.split()[0])
          self.win_max[point]=int(line.split()[0])

    # infer from format of window function files whether we will use polarisation spectra or not 
    num_col=len(line.split())
    if (num_col == 2):
      self.has_pol=False
    else:
      if (num_col == 5):
        self.has_pol=True
      else:
        print 'window function files are understood of they contain 2 columns (l TT)'
        print 'or 5 columns (l TT TE EE BB)'
        print 'in this case the number of columns is',num_col
        exit()

    # define array of window functions
    self.window=np.zeros((self.num_points,max(self.win_max)+1,num_col-1),'float64')

    # go again through window function file, this time reading window functions; 
    # that are distributed as: l TT (TE EE BB)
    # where the last columns contaim W_l/l, not W_l
    # we mutiply by l in order to store the actual W_l
    for point in range(self.num_points):
      for line in open(self.data_directory+'windows/'+window_name+str(used_index[point]+1),'r'):
        l=int(line.split()[0])
        if (((self.has_pol==False) and (len(line.split()) !=2)) or ((self.has_pol==True) and (len(line.split()) !=5))):
           print 'for given experiment, all window functions should have the same number of columns, 2 or 5. This is not the case here.'
           exit()
        if ((l>=self.win_min[point]) and (l<=self.win_max[point])):
          self.window[point,l,:]=[float(line.split()[i]) for i in range(1,len(line.split()))]
          self.window[point,l,:]*=l

    # eventually, initialise quantitites used in the marginalization over nuisance parameters
    if ((self.has_xfactors == True) and ((self.calib_uncertainty > 1.e-4) or (self.has_beam_uncertainty == True))):
      self.halfsteps=5
      self.margeweights = np.zeros(2*self.halfsteps+1,'float64') 
      for i in range(-self.halfsteps,self.halfsteps+1):
        self.margeweights[i+self.halfsteps]=np.exp(-(float(i)*3./float(self.halfsteps))**2/2)
      self.margenorm=sum(self.margeweights)

    # store maximum value of l needed by window functions
    self.l_max=max(self.win_max)  

    # impose that CLASS computes Cl's up to maximum l needed by window function
    self._need_Class_args(data,{'l_max_scalars':self.l_max})

    # deal with nuisance parameters
    try:
      self.use_nuisance
    except:
      self.use_nuisance = []
    self._read_contamination_spectra(data)

    # end of initialisation

  def _loglkl(self,_cosmo,data):

    # get Cl's from CLASS
    cl = self._get_cl(_cosmo)

    # add contamination spectra multiplied by nuisance parameters
    cl = self._add_contamination_spectra(cl,data)

    # get likelihood
    loglkl = self._compute_loglkl(cl,_cosmo,data)

    # add prior on nuisance parameters
    loglkl = self._add_nuisance_prior(loglkl,data)

    return loglkl



  def _compute_loglkl(self,cl,_cosmo,data):

    # checks that Cl's have been computed up to high enough l given window function range. Normally this has been imposed before, so this test could even be supressed.
    if (np.shape(cl['tt'])[0]-1 < self.l_max):
      print 'CLASS computed Cls till l=',np.shape(cl['tt'])[0]-1,'while window functions need',self.l_max
      exit()

    # compute theoretical bandpowers, store them in theo[points]
    theo=np.zeros(self.num_points,'float64')

    for point in range(self.num_points):

      # find bandpowers B_l by convolving C_l's with [(l+1/2)/2pi W_l]
      for l in range(self.win_min[point],self.win_max[point]):

        theo[point] += cl['tt'][l]*self.window[point,l,0]*(l+0.5)/2./math.pi

        if (self.has_pol):
          theo[point] += (cl['te'][l]*self.window[point,l,1] + cl['ee'][l]*self.window[point,l,2] + cl['bb'][l]*self.window[point,l,3])*(l+0.5)/2./math.pi

    # allocate array for differencve between observed and theoretical bandpowers
    difference=np.zeros(self.num_points,'float64')

    # depending on the presence of lognormal likelihood, calibration uncertainty and beam uncertainity, use several methods for marginalising over nuisance parameters:

    # first method: numerical integration over calibration uncertainty:
    if (self.has_xfactors and ((self.calib_uncertainty > 1.e-4) or self.has_beam_uncertainty)):
    
      chisq_tmp=np.zeros(2*self.halfsteps+1,'float64')
      chisqcalib=np.zeros(2*self.halfsteps+1,'float64')
      beam_error=np.zeros(self.num_points,'float64')

      # loop over various beam errors
      for ibeam in range(2*self.halfsteps+1):

        # beam error
        for point in range(self.num_points):
          if (self.has_beam_uncertainty):
            beam_error[point]=1.+self.beam_error[point]*(ibeam-self.halfsteps)*3/float(self.halfsteps)
          else:
            beam_error[point]=1.

        # loop over various calibraion errors
        for icalib in range(2*self.halfsteps+1):

          # calibration error
          calib_error=1+self.calib_uncertainty*(icalib-self.halfsteps)*3/float(self.halfsteps)
         
          # compute difference between observed and theoretical points, after correcting the later for errors
          for point in range(self.num_points):

            # for lognormal likelihood, use log(B_l+X_l)
            if (self.has_xfactor[point]):
              difference[point]=self.obs[point]-math.log(theo[point]*beam_error[point]*calib_error+self.xfactor[point])
            # otherwise use B_l
            else:
              difference[point]=self.obs[point]-theo[point]*beam_error[point]*calib_error

              # find chisq with those corrections
         #chisq_tmp[icalib]=np.dot(np.transpose(difference),np.dot(self.inv_covmat,difference))
          chisq_tmp[icalib]=np.dot(difference,np.dot(self.inv_covmat,difference))

        minchisq=min(chisq_tmp)

       # find chisq marginalized over calibration uncertainty (if any)
        tot=0
        for icalib in range(2*self.halfsteps+1):
          tot += self.margeweights[icalib]*math.exp(max(-30.,-(chisq_tmp[icalib]-minchisq)/2.))

        chisqcalib[ibeam]=-2*math.log(tot/self.margenorm)+minchisq

     # find chisq marginalized over beam uncertainty (if any)
      if (self.has_beam_uncertainty):

        minchisq=min(chisqcalib)

        tot=0
        for ibeam in range(2*self.halfsteps+1):
          tot += self.margeweights[ibeam]*math.exp(max(-30.,-(chisqcalib[ibeam]-minchisq)/2.))

        chisq=-2*math.log(tot/self.margenorm)+minchisq

      else:
        chisq=chisqcalib[0]

   # second method: marginalize over nuisance parameters (if any) analytically
    else:

     # for lognormal likelihood, theo[point] should contain log(B_l+X_l)
      if (self.has_xfactors): 
        for point in range(self.num_points):
          if (self.has_xfactor[point]):
            theo[point]=math.log(theo[point]+self.xfactor[point])
 
     # find vector of difference between observed and theoretical bandpowers
      difference=self.obs-theo  

     # find chisq
      chisq=np.dot(np.transpose(difference),np.dot(self.inv_covmat,difference))

     # correct eventually for effect of analytic marginalization over nuisance parameters
      if ((self.calib_uncertainty > 1.e-4) or self.has_beam_uncertainty):
            
        denom=1.                  
        tmp= np.dot(self.inv_covmat,theo)                             
        chi2op=np.dot(np.transpose(difference),tmp)
        chi2pp=np.dot(np.transpose(theo),tmp)
         
        if (self.has_beam_uncertainty):
          for points in range(self.num_points):
            beam[point]=self.beam_error[point]*theo[point]
          tmp=np.dot(self.inv_covmat,beam)
          chi2dd=np.dot(np.transpose(beam),tmp)
          chi2pd=np.dot(np.transpose(theo),tmp)
          chi2od=np.dot(np.transpose(difference),tmp)
            
        if (self.calib_uncertainty > 1.e-4):
          wpp=1/(chi2pp+1/self.calib_uncertainty**2)
          chisq=chisq-wpp*chi2op**2
          denom = denom/wpp*self.calib_uncertainty**2
        else:
          wpp=0

        if (self.has_beam_uncertainty):
          wdd=1/(chi2dd-wpp*chi2pd**2+1)
          chisq=chisq-wdd*(chi2od-wpp*chi2op*chi2pd)**2
          denom=denom/wdd

        chisq+=log(denom)

   # finally, return ln(L)=-chi2/2

    self.loglkl = - 0.5 * chisq 
    return self.loglkl



class likelihood_clik(likelihood):

  

  def __init__(self,path,data,command_line=False):

    try:
      import clik
    except ImportError:
      print " /|\  You must first activate the binaries from the Clik distribution,"
      print "/_o_\ please run : source /path/to/clik/bin/clik_profile.sh"
      print "      and try again."
      exit()
    likelihood.__init__(self,path,data,command_line)
    self._need_Class_args(data,{'lensing':'yes', 'output':'tCl lCl pCl'})
    self.clik = clik.clik(self.path_clik)
    self.l_max = max(self.clik.get_lmax())
    self._need_Class_args(data,{'l_max_scalar':self.l_max})

    # deal with nuisance parameters
    try:
      self.use_nuisance
    except:
      self.use_nuisance = []
    self._read_contamination_spectra(data)

  def _loglkl(self,_cosmo,data):

    # get Cl's from CLASS
    cl = self._get_cl(_cosmo)
  
    # add contamination spectra multiplied by nuisance parameters
    cl = self._add_contamination_spectra(cl,data)

    # allocate array of Cl's and nuisance parameters
    tot=np.zeros(np.sum(self.clik.get_lmax())+6+len(self.clik.get_extra_parameter_names()))

    # fill with Cl's
    index=0
    for i in range(np.shape(self.clik.get_lmax())[0]):
      if (self.clik.get_lmax()[i] >-1):
        for j in range(self.clik.get_lmax()[i]+1):
          if (i==0):
            tot[index+j]=cl['tt'][j]        
          if (i==1):
            tot[index+j]=cl['ee'][j]
          if (i==2):
            tot[index+j]=cl['bb'][j]
          if (i==3):
            tot[index+j]=cl['te'][j]
          if (i==4):
            tot[index+j]=cl['tb'][j]
          if (i==5):
            tot[index+j]=cl['eb'][j]

        index += self.clik.get_lmax()[i]

    # fill with nuisance parameters
    for nuisance in self.clik.get_extra_parameter_names():
      if nuisance in data.nuisance_param_names:
        nuisance_value = float(data.vector[np.where(data.param_names == nuisance)[0][0]])
      else:
        try:
          exec "nuisance_value = data.%s" % nuisance
        except:
          print 'the likelihood needs a parameter '+nuisance
          print 'you must pass it through the input file (as a free nuisance parameter or a fixed parameter)'
          exit()
      index += 1
      tot[index]=nuisance_value

    # compute likelihood
    loglkl=self.clik(tot)[0]

    # add prior on nuisance parameters
    loglkl = self._add_nuisance_prior(loglkl,data)

    return loglkl
