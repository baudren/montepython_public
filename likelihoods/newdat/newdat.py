import os
import numpy as np

class newdat():
  
  def __init__(self,path):

    # path of .data file
    abs_folder = os.path.dirname(os.path.abspath(path))+'/'

    # parsing .data file
    for line in open(path,'r'):
      if line.find('#')==-1:
	if line.find('newdat.')!=-1:
	  exec(line.replace('newdat.','self.'))
      
    # open .newdat file
    newdatfile=open(abs_folder+self.directory+self.file,'r')

    # (can be moved later) define string for window function file name
    window_name=newdatfile.readline().strip('\n')
    #print window_name
    
    # read number of bands for six modes
    # order fixed to be TT, EE, BB, EB, TE, TB
    band_num = newdatfile.readline().split()
    #print band_num

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
    #print self.calib_uncertainty

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

    #print self.num_points
    #print full_num_points
    #print self.obs
    #print used_index

    #read correlation matrix
    full_inv_covmat=np.zeros((full_num_points,full_num_points),'float64')
    self.inv_covmat=np.zeros((self.num_points,self.num_points),'float64')
    #print self.inv_covmat
    for point in range(full_num_points):
      full_inv_covmat[point]=newdatfile.readline().split()
    for point in range(self.num_points):
      self.inv_covmat[point]=full_inv_covmat[used_index[point],used_index]

    #print full_inv_covmat
    #print 'inv_covmat:'
    #print self.inv_covmat

    #read window functions
    self.win_min=np.zeros(self.num_points,'int')
    self.win_max=np.zeros(self.num_points,'int')
    
#    window=np.array([],'float64')

    for point in range(self.num_points):
      for line in open(abs_folder+self.directory+'windows/'+window_name+str(point+1),'r'):
        if any([float(line.split()[i]) != 0. for i in range(1,len(line.split()))]):
          if (self.win_min[point]==0):
            self.win_min[point]=line.split()[0]
          self.win_max[point]=line.split()[0]

    #print self.win_min,self.win_max

    #print self.num_points,max(self.win_max)+1,len(line.split())-1
    self.window=np.zeros((self.num_points,max(self.win_max)+1,len(line.split())-1),'float64')

    for point in range(self.num_points):
      for line in open(abs_folder+self.directory+'windows/'+window_name+str(point+1),'r'):
        l=int(line.split()[0])
        if ((l>=self.win_min[point]) and (l<=self.win_max[point])):
          self.window[point,l,:]=[line.split()[i] for i in range(1,len(line.split()))]

    #print "point18"
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

  def _loglkl(self,_cosmo,data):

   #cl = nm.ndarray([4,lmax+1], dtype=nm.double)
   #cl = _cosmo.lensed_cl(-1,False)
   cl=np.array(_cosmo.lensed_cl(),'float64')
   
#   for point in range(self.num_points)
      # sum over l and types of from win_min[point] to w_max[point] of cl[l][type]*window[point][l][type], where type = TT or TT EE BB TE
#      theo[point]=np.sum(cl[l,:]*window[point,l])
                                       
      #if ((self.has_xfactors == True) && (self.has_xfactors[point] == True))
       #  theo[point]=log(theo[point]+xfactor[point])
 
  #    difference[point]=self.obs[point]-theo[point]                      
     
  # if (self.has_xfactors && (self.calib_uncertainty > 1.e-4 || self.has_beam_uncertainty == True))
      #GetCalibMargexChisq
   #   for j=-self.halfsteps,self.halfstep
    #     if (self.has_beam_uncertainty == True)
     #       theo[] *= (1+self.beam_error[]*j*3/float(self.halfsteps))
      #   for i=-self.halfsteps,self.halfstep
       #     calib=1+self.calib_uncertainty*i*3/float(halfsteps)
       #     chisq[i]=np.dot(np.transpose(difference),np.dot(self.inv_covmat,difference))
       #  minchisq=min(chisq[])
       #  chisqcalib[j]=-2*log(sum(self.mar...)

#         .....

#   else
#      chisq=np.dot(np.transpose(difference),np.dot(self.inv_covmat,difference))

#      if (self.calib_uncertainty > 1.e-4 || self.has_beam_uncertainty == True)
                              
#         tmp= np.dot(self.inv_covmat,theo)                             
#         chi2op=np.dot(np.transpose(difference),tmp)
 #        chi2pp=np.dot(np.transpose(theo),tmp)
         
#         if (self.has_beam_uncertainty == True)
#            beam[]=self.beam_error[]*theo[]
#            tmp=np.dot(self.inv_covmat,beam)
#            chi2dd=np.dot(np.transpose(beam),tmp)
#            chi2pd=np.dot(np.transpose(theo),tmp)
#            chi2od=np.dot(np.transpose(difference),tmp)
            
#         if (self.calib_uncertainty > 1.e-4)
#            wpp=1/(chi2pp+1/self.calib_uncertainty**2)
#            chisq=chisq-wpp*chi2op**2
#         else
#            wpp=0

#         if (self.has_beam_uncertainty == True)
#            wdd=1/(chi2dd-wpp*chi2pd**2+1)
#            chisq=chisq-wdd*(chi2od-wpp*chi2op*chi2pd**2)
#            denom=denom/wdd

#    self.loglkl = - 0.5 * chisq 
#    return self.loglkl


newdat = newdat('newdat.data')
