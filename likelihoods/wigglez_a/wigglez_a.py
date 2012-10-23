# simple likelihood for sdss_lrgDR4
from likelihood_class import likelihood
import numpy as np
import math

class wigglez_a(likelihood):

  def __init__(self,path,data,command_line,log_flag,default):

    likelihood.__init__(self,path,data,command_line,log_flag,default)

    # require P(k) from class
    self.need_cosmo_arguments(data,{'output':'mPk'})
    if self.use_halofit:
      self.need_cosmo_arguments(data,{'non linear':'halofit'}) 

    # read values of k (in h/Mpc)
    self.k_size=self.max_mpk_kbands_use-self.min_mpk_kbands_use+1
    self.mu_size=1
    self.k = np.zeros((self.k_size),'float64')
    self.kh = np.zeros((self.k_size),'float64')

    datafile = open(self.data_directory+self.kbands_file,'r')
    for i in range(self.num_mpk_kbands_full):
      line = datafile.readline()
      if i+2 > self.min_mpk_kbands_use and i < self.max_mpk_kbands_use:
        self.kh[i-self.min_mpk_kbands_use+1]=float(line.split()[0])
    datafile.close()      

    khmax = self.kh[-1]

    # check if need hight value of k for giggleZ
    if self.Use_giggleZ:
      datafile = open(self.data_directory+self.giggleZ_fidpk_file,'r')
      counter = 1
      line = datafile.readline()
      k=float(line.split()[0])
      while (k<khmax):
        counter += 1
        line = datafile.readline()
        k=float(line.split()[0])
      datafile.close()
      self.k_fid_size=counter  
      khmax=k

    # require k_max and z_max from the cosmological module
    self.need_cosmo_arguments(data,{'P_k_max_h/Mpc':khmax,'z_max_pk':self.redshift})

    # In case of a comparison, stop here
    if not default:
      return

    # read information on different regions in the sky
    if (self.has_regions):
      self.num_regions = len(self.used_region)
      self.num_regions_used = 0
      for i in range(self.num_regions):
        if (self.used_region[i]):
          self.num_regions_used += 1
      if (self.num_regions_used == 0):
        print 'mpk: no regions begin used in this data set'
        exit()
    else:
      self.num_regions = 1
      self.num_regions_used = 1
      self.used_region = [True]
 
    # read window functions
    self.n_size=self.max_mpk_points_use-self.min_mpk_points_use+1

    self.window=np.zeros((self.num_regions,self.n_size,self.k_size),'float64')

    datafile = open(self.data_directory+self.windows_file,'r')
    for i_region in range(self.num_regions):
      if self.num_regions>1:
        line = datafile.readline()
      for i in range(self.num_mpk_points_full):
        line = datafile.readline()
        if i+2 > self.min_mpk_points_use and i < self.max_mpk_points_use:
          for j in range(self.k_size):  
            self.window[i_region,i-self.min_mpk_points_use+1,j]=float(line.split()[j+self.min_mpk_kbands_use-1])
    datafile.close()    

    # read measurements
    self.P_obs=np.zeros((self.num_regions,self.n_size),'float64')
    self.P_err=np.zeros((self.num_regions,self.n_size),'float64')

    datafile = open(self.data_directory+self.measurements_file,'r')
    for i_region in range(self.num_regions):
      for i in range(2):
        line = datafile.readline()
      for i in range(self.num_mpk_points_full):
        line = datafile.readline()
        if i+2 > self.min_mpk_points_use and i < self.max_mpk_points_use:
          self.P_obs[i_region,i-self.min_mpk_points_use+1]=float(line.split()[3])  
          self.P_err[i_region,i-self.min_mpk_points_use+1]=float(line.split()[4])  
    datafile.close()

    # read covariance matrices
    self.invcov=np.zeros((self.num_regions,self.n_size,self.n_size),'float64')
    cov = np.zeros((self.n_size,self.n_size),'float64')
    invcov_tmp = np.zeros((self.n_size,self.n_size),'float64')

    datafile = open(self.data_directory+self.covmat_file,'r')
    for i_region in range(self.num_regions):
      for i in range(1):
        line = datafile.readline()
      for i in range(self.num_mpk_points_full):
        line = datafile.readline()
        if i+2 > self.min_mpk_points_use and i < self.max_mpk_points_use:
          for j in range(self.num_mpk_points_full):
            if j+2 > self.min_mpk_points_use and j < self.max_mpk_points_use:
              cov[i-self.min_mpk_points_use+1,j-self.min_mpk_points_use+1]=float(line.split()[j])
      invcov_tmp=np.linalg.inv(cov)
      for i in range(self.n_size):
        for j in range(self.n_size):
          self.invcov[i_region,i,j]=invcov_tmp[i,j]
    datafile.close()

    # read fiducial model
    if self.Use_giggleZ:
      self.P_fid=np.zeros((self.k_fid_size),'float64')
      self.k_fid=np.zeros((self.k_fid_size),'float64')
      datafile = open(self.data_directory+self.giggleZ_fidpk_file,'r')
      for i in range(self.k_fid_size):
        line = datafile.readline()
        self.k_fid[i]=float(line.split()[0])
        self.P_fid[i]=float(line.split()[1])
      datafile.close()

    return

  # compute likelihood

  def loglkl(self,_cosmo,data):

    print "get here"

    # reduced Hubble parameter
    h = _cosmo._h()
     
    if self.use_scaling:
      # angular diameter distance at this redshift, in Mpc/h
      d_angular = _cosmo._angular_distance(self.redshift)
      d_angular *= h  

      # radial distance at this redshift, in Mpc/h
      r,Hz = _cosmo.z_of_r([self.redhsift])
      d_radial = self.redshift*h/Hz[0]

      # scaling factor = (d_angular**2 * d_radial)^(1/3) relative
      # to a fiducial model
      scaling = pow((d_angular/self.d_angular_fid)**2*(d_radial/self.d_radial_fid),1./3.)
    else:
      scaling = 1

    # get rescaled values of k in 1/Mpc
    self.k = self.kh *h*scaling

    # get P(k) at right values of k, convert it to (Mpc/h)^3 and rescale it
    P_lin = np.zeros((self.k_size),'float64')

    if self.Use_giggleZ:

      P = np.zeros((self.k_fid_size),'float64')  
      
      for i in range(self.k_fid_size):

        P[i] = _cosmo._pk(self.k_fid[i]*h,self.redshift)

        power=0
        for j in range(6):
          power += self.giggleZ_fidpoly[j]*self.k_fid[i]**j

        # rescale P by fiducial model and get it in (Mpc/h)**3
        P[i] *= pow(10,power)/self.P_fid[i]*(h/scaling)**3

      # get rescaled values of k in 1/Mpc
      self.k=self.kh *h*scaling

      P_lin = np.interp(self.k,self.k_fid,P)

    else:
      # get rescaled values of k in 1/Mpc
      self.k=self.kh *h*scaling
      # get values of P(k) in Mpc**3
      for i in range(self.k_size):
        P_lin[i] = _cosmo._pk(self.k[i],self.redshift)
      # get rescaled values of P(k) in (Mpc/h)**3
      P_lin *= (h/scaling)**3

    do_marge = self.Q_marge

    W_P_th =  np.zeros((self.n_size),'float64')

    if do_marge and self.Q_flat:

      P_th =  np.zeros((self.k_size),'float64')
      for i in range(self.k_size):
        P_th[i] = P_lin[i]/(1.+self.Ag*self.k[i]) 

      k2 =  np.zeros((self.k_size),'float64')
      for i in range(self.k_size):
        k2[i] = P_th[i] * self.k[i]**2

      W_P_th = np.dot(self.window,P_th)

      W_P_th_k2 =  np.zeros((self.n_size),'float64')
      W_P_th_k2 = np.dot(self.window,k2)
                
      w = np.zeros((self.n_size),'float64') 
      w=1./(self.P_err*self.P_err)

      covdat = np.zeros((self.n_size),'float64') 
      covdat=self.P_obs*w  

      covth = np.zeros((self.n_size),'float64') 
      covth=W_P_th*w

      covth_k2 = np.zeros((self.n_size),'float64') 
      covth_k2=W_P_th_k2*w

      offdiag=sum(covth*W_P_th_k2)

      Mat=np.zeros((2,2),'float64')
      Mat=[[sum(covth*W_P_th),offdiag],[offdiag,sum(covth_k2*W_P_th_k2)]]

      Vec=np.zeros((2),'float64')
      Vec=[sum(covdat*W_P_th),sum(covdat*W_P_th_k2)]

      chi2=-sum(self.P_obs*covdat)+np.dot(Vec,np.dot(np.linalg.inv(Mat),Vec))-math.log(np.linalg.det(Mat))
      return -chi2/2

    else:  

      if (self.Q_sigma == 0):
        do_marge = False

      if (self.Use_jennings or self.Use_simpledamp):
        print "case with Use_jennings or Use_simpledamp not coded"
        exit()
      else :
        print "starting analytic marginalisation over bias"

        P_data_large =  np.zeros((self.n_size*self.num_regions_used),'float64')
        W_P_th_large =  np.zeros((self.n_size*self.num_regions_used),'float64')
        cov_dat_large =  np.zeros((self.n_size*self.num_regions_used),'float64')
        cov_th_large =  np.zeros((self.n_size*self.num_regions_used),'float64')

        normV = 0

        if do_marge:
          nQ = 6
          dQ = 0.4
        else:
          nQ=0
          dQ=0

        chisq = np.zeros((nQ*2+1),'float64')
        calweights = np.zeros((nQ*2+1),'float64')

        for iQ in range(-nQ,nQ+1):
          
          P_th =  np.zeros((self.k_size),'float64')
          for i in range(self.k_size):
            if self.Q_marge:
              Q = self.Q_mid +iQ*self.Q_sigma*dQ 
              P_th[i] = P_lin[i]*(1+Q*self.k[i]**2)/(1.+self.Ag*self.k[i]) 
            else:
              P_th[i] = P_lin[i]

          for i_region in range(self.num_regions):
            if self.used_region[i_region]:
              imin = i_region*self.n_size
              imax = (i_region+1)*self.n_size-1

              W_P_th = np.dot(self.window[i_region,:],P_th)
              for i in range(self.n_size):
                P_data_large[imin+i] = self.P_obs[i_region,i]
                W_P_th_large[imin+i] = W_P_th[i]
                cov_dat_large[imin+i] = np.sum(self.invcov[i_region,i,:]*self.P_obs[i_region,:])
                cov_th_large[imin+i] = np.sum(self.invcov[i_region,i,:]*W_P_th[:])
          normV += np.sum(W_P_th_large*cov_th_large)
          b_out = np.sum(W_P_th_large*cov_dat_large)/np.sum(W_P_th_large*cov_th_large)
          print "bias value",b_out
          chisq[iQ+nQ] = np.sum(P_data_large*cov_dat_large)  - np.sum(W_P_th_large*cov_dat_large)**2/normV
          
          if do_marge:
            calweights[iQ+nQ] = exp(-(iQ*dQ)**2/2)
          else: 
            return -chisq[iQ+nQ]/2

      if do_marge:
        if not self.Use_jennings:
          minchisq=np.min(chisq)
          lnlike = np.sum(exp(-(chisq[:]-minchisq)/2)*calweights[:])/np.sum(calweights[:])
          if (lnlike == 0):
            return data.boundary_loglike
          else:
            print "case with Use_jennings not coded"
            exit()

