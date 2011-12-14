import os,sys
from collections import OrderedDict as od

# Definition of classes for the likelihoods to inherit, for every different
# types of expected likelihoods. Included so far, models for a Clik, simple
# prior, and newdat likelihoods.

# General class that will only define the store_lkl function, common for every
# single likelihood. It will copy the content of self.path, copied from
# initialization
class likelihood():

  def __init__(self,path,data,command_line=False):

    self.name = self.__class__.__name__
    self.folder = sys.path[0]+'/'
    self._read_from_file(path,data)
    if ((command_line is not False) and (os.path.exists(command_line.folder+'/log.dat') is False)):
      self._store_lkl_params(command_line)
    

  def _loglkl(self,_cosmo,data):
    raise NotImplementedError('Must implement method _loglkl() in your likelihood')

  def _store_lkl_params(self,command_line):
    log = open(command_line.folder+'log.param','a')
    tolog = open(self.path,'r')
    log.write("\n#-----Likelihood-{0}-----\n\n".format(self.name))
    for line in tolog:
      log.write(line)
    tolog.seek(0)
    tolog.close()
    log.close()

  def _read_from_file(self,path,data):
    self.path = path
    self.nuisance = od()
    if os.path.isfile(path):
      data_file = open(path,'r')
      for line in data_file:
	if line.find('#')==-1:
	  if line.find(self.name+'.')!=-1:
	    exec(line.replace(self.name+'.','self.'))
      data_file.seek(0)
      data_file.close()
    for key,value in self.nuisance.iteritems():
      data.nuisance_params[self.name+'_'+key] = value


# Likelihood type for prior
class likelihood_prior(likelihood):

  def _loglkl(self):
    raise NotImplementedError('Must implement method __loglkl() in your likelihood')


class likelihood_newdat(likelihood):

  def __init__(self,path,data,command_line=False):

    likelihood.__init__(self,path,command_line)
    
    # open .newdat file
    newdatfile = open(self.data_directory+self.file,'r')

    # (can be moved later) define string for window function file name
    window_name = newdatfile.readline().strip('\n')
    
    # read number of bands for six modes
    # order fixed to be TT, EE, BB, EB, TE, TB
    band_num = newdatfile.readline().split()

    # probably useless
    #    self.has_pol=False
    #    if (any(self.band_num[2:6] != 0))
    #       self.has_pol=True

    # read string equal to 'BAND_SELECTION' or not
    # if yes, read 6 lines 'min, max'
    band_min = np.zeros(6,'int')
    band_max = np.zeros(6,'int')

    line = str(newdatfile.readline()).strip('\n')
    if (line=='BAND_SELECTION'):
      for i in range(6):
        # read line
        line = newdatfile.readline()
        band_min[i] = int(line.split()[0])   
        band_max[i] = int(line.split()[1])   
    else:
      band_min=[1 for i in range(6)]
      band_max=band_num
    #print band_min
    #print band_max

    # read line = flag, calib, calib_uncertainty
    line = newdatfile.readline()
    calib=float(line.split()[1])
    if (int(line.split()[0])==0):
      self.calib_uncertainty=0
    else:
      self.calib_uncertainty=float(line.split()[2])

    # read line = 0or1or2, beam_width, beam_sigma
    line = newdatfile.readline()
    beam_type = int(line.split()[0])
    if (beam_type > 0):
      self.has_beam_uncertainty = True
    else:
      self.has_beam_uncertainty = False
    beam_width = float(line.split()[1])
    beam_sigma = float(line.split()[2])
    
    # read integer = 0, 1 or 2 for xfactors
    likelihood_type = int(newdatfile.readline())
 
    if (likelihood_type > 0):
      self.has_xfactors = True
    else:
      self.has_xfactors = False

    # declare array of qunatitites describing each point of measurement
    # size yet unknown, will be decided later and kept as self.num_points
    #self.win_min=np.array([],'int')
    #self.win_max=np.array([],'int')
    #self.window=np.array(lmax,[],'float64')
    self.obs=np.array([],'float64')
    self.var=np.array([],'float64')
    self.beam_err=np.array([],'float64')
    #self.has_pol=np.array([],'bool')
    self.has_xfactor=np.array([],'bool')
    #self.xfactor=np.array([],'float64')
    used_index=np.array([],'int')

    index=-1

    for cltype in range(6):
      if (band_num[cltype] != 0):
        #read name
        newdatfile.readline()
        for band in range(int(band_num[cltype])):
          # read line
          line = newdatfile.readline()
          index += 1

          if ((band >= band_min[cltype]-1) and (band <= band_max[cltype]-1)):

            used_index=np.append(used_index,index)

            self.obs=np.append(self.obs,float(line.split()[1])*calib**2)

            self.var=np.append(self.var,(0.5*(float(line.split()[2])+float(line.split()[3]))*calib**2)**2)

            if (likelihood_type == 0):
              self.has_xfactor=np.append(self.has_xfactor,[False])
            if (likelihood_type == 1):
              self.has_xfactor=np.append(self.has_xfactor,[True])
            if (likelihood_type == 2):
              self.has_xfactor=np.append(self.has_xfactor,int(line.split()[7]))

            if (beam_type == 0):
              self.beam_error=0
            if (beam_type == 1):
              l_mid=int(line.split()[5])+0.5*(int(line.split()[5])+int(line.split()[6]))
              self.beam_error=abs(exp(-l_mid*(l_mid+1)*1.526e-8*2*beam_sigma*beam_width-1.))
            if (beam_type == 2):
              if (likelihood_type == 2):
                self.beam_error=np.append(self.beam_error,int(line.split()[8]))
              else:
                self.beam_error=np.append(self.beam_error,int(line.split()[7]))

        for band in range(int(band_num[cltype])):      
          newdatfile.readline()

    self.num_points=np.shape(self.obs)[0]
    full_num_points=index+1

    #read correlation matrix
    full_inv_covmat=np.zeros((full_num_points,full_num_points),'float64')
    self.inv_covmat=np.zeros((self.num_points,self.num_points),'float64')
    #print self.inv_covmat
    for point in range(full_num_points):
      full_inv_covmat[point]=newdatfile.readline().split()
    for point in range(self.num_points):
      self.inv_covmat[point]=full_inv_covmat[used_index[point],used_index]

    #read window functions
    self.win_min=np.zeros(self.num_points,'int')
    self.win_max=np.zeros(self.num_points,'int')
    
#    window=np.array([],'float64')

    for point in range(self.num_points):
      for line in open(abs_folder+self.data_directory+'windows/'+window_name+str(point+1),'r'):
        if any([float(line.split()[i]) != 0. for i in range(1,len(line.split()))]):
          if (self.win_min[point]==0):
            self.win_min[point]=line.split()[0]
          self.win_max[point]=line.split()[0]

    #print self.num_points,max(self.win_max)+1,len(line.split())-1
    self.window=np.zeros((self.num_points,max(self.win_max)+1,len(line.split())-1),'float64')

    for point in range(self.num_points):
      for line in open(abs_folder+self.data_directory+'windows/'+window_name+str(point+1),'r'):
        l=int(line.split()[0])
        if ((l>=self.win_min[point]) and (l<=self.win_max[point])):
          self.window[point,l,:]=[line.split()[i] for i in range(1,len(line.split()))]

    #for l in range (self.win_min[17],self.win_max[17]+1):
    #  print l,self.window[17,l,:]

    #eventually prepare marginalization oiver nuisance parameters
    if ((self.has_xfactors == True) and ((self.calib_uncertainty > 1.e-4) or (self.has_beam_uncertainty == True))):
      self.halfsteps=5
      self.margeweights = np.zeros(2*self.halfsteps+1,'float64') 
      for i in range(-self.halfsteps,self.halfsteps+1):
        self.margeweights[i+self.halfsteps]=np.exp(-(float(i)*3./float(self.halfsteps))**2/2)
      self.margenorm=sum(self.margeweights)

    #    covmat = np.zeros((np.shape(self.z)[0],np.shape(self.z)[0]),'float64')
    #    i=0
    #    for line in open(abs_folder+covmat_filename,'r'):
    #      covmat[i] = line.split()
    #      i+=1

    #    self.inv_covmat=np.linalg.inv(covmat)
    #    self.trace_inv_covmat=np.trace(self.inv_covmat)

  def _loglkl(self):
    raise NotImplementedError('Must implement method __loglkl() in your likelihood')
