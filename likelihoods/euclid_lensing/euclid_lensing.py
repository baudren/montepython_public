from likelihood_class import likelihood
import os
import numpy as np
import math
# Adapted from Julien Lesgourgues

class euclid_lensing(likelihood):

  def __init__(self, path, data, command_line):

    likelihood.__init__(self, path, data, command_line)

    self.need_cosmo_arguments(data,{'output':'mPk'})

    # Define array of l values, and initialize them
    self.l = np.zeros(self.nlmax,'float64')
    for nl in range(self.nlmax):
      self.l[nl] = 1.*math.exp(self.dlnl*nl)
    
    #print self.l[:]
    #exit()

    ########################################################
    # Find distribution of dn_dz (not normalized) in each bin
    ########################################################
    
    # Assuming each bin contains the same number of galaxies, we find the bin
    # limits in z space

    # Compute the total number of galaxies until zmax (no normalization yet)

    n_tot = 0.
    for z in np.arange(0,self.zmax+self.dz,self.dz):
      gd_1 = self.galaxy_distribution(z)
      gd_2 = self.galaxy_distribution(z+self.dz)
      n_tot += 0.5*(gd_1+gd_2)*self.dz

    # For each bin, compute the limit in z space
  
    # Create the array that will contain the z boundaries for each bin. The
    # first value is already correctly set to 0.
    self.z_bin_edge = np.zeros(self.nbin+1,'float64')

    for Bin in range(self.nbin-1):

      bin_count = 0.
      z =self.z_bin_edge[Bin]

      while (bin_count <= n_tot/self.nbin):
	gd_1 = self.galaxy_distribution(z)
	gd_2 = self.galaxy_distribution(z+self.dz)
	bin_count += 0.5*(gd_1+gd_2)*self.dz
	z += self.dz

      self.z_bin_edge[Bin+1] = z

    self.z_bin_edge[self.nbin] = self.zmax
    
    # Fill array of discrete z values
    self.z = np.zeros(self.nzmax,'float64')
    for nz in range(self.nzmax):
      self.z[nz] = (nz*1.0)/(self.nzmax-1.0)*self.zmax

    # Force the cosmological module to store Pk for redshifts up to max(self.z)
    self.need_cosmo_arguments(data,{'z_max_pk':self.z[-1]})
    # Force the cosmological module to store Pk for k up to an arbitrary number
    # (since self.r is not yet decided)... TODO
    self.need_cosmo_arguments(data,{'P_k_max_1/Mpc':self.k_max})

    # Fill distribution for each bin (convolving with photo_z distribution)
    self.eta_z = np.zeros((self.nzmax,self.nbin),'float64')
    for Bin in range(self.nbin):

      for nz in range(self.nzmax):
	z = self.z[nz]
	self.eta_z[nz,Bin] = 0.

	for nz2 in range(1,self.nzmax):

	  if ((self.z[nz2] >= self.z_bin_edge[Bin]) and (self.z[nz2] <= self.z_bin_edge[Bin+1])):
	    gd  = self.galaxy_distribution(self.z[nz2])
	    pzd = self.photo_z_distribution(z,self.z[nz2])
	    integrand_plus = gd*pzd
	  else:
	    integrand_plus = 0.

	  if ((self.z[nz2-1] >= self.z_bin_edge[Bin]) and (self.z[nz2-1] <= self.z_bin_edge[Bin+1])):
	    gd  = self.galaxy_distribution(self.z[nz2-1])
	    pzd = self.photo_z_distribution(z,self.z[nz2-1])
	    integrand_minus = gd*pzd
	  else:
	    integrand_minus = 0.
	  
	  self.eta_z[nz,Bin] += 0.5*(integrand_plus+integrand_minus)*(self.z[nz2]-self.z[nz2-1])

    #for nz in range(self.nzmax):
    #  print self.z[nz],self.eta_z[nz,0],self.eta_z[nz,1],self.eta_z[nz,2],self.eta_z[nz,3],self.eta_z[nz,4],self.galaxy_distribution(self.z[nz])
    #exit()  

    # integrate eta(z) over z (in view of normalizing it to one)
    self.eta_norm = np.zeros(self.nbin,'float64')
    for Bin in range(self.nbin):
      self.eta_norm[Bin] = np.sum(0.5*(self.eta_z[1:,Bin]+self.eta_z[:-1,Bin])*(self.z[1:]-self.z[:-1]))

    ################
    # Noise spectrum 
    ################
    
    # Number of galaxies per steradian
    self.noise = 3600.*self.gal_per_sqarcmn*(180./math.pi)**2

    # Number of galaxies per steradian per bin
    self.noise = self.noise/self.nbin
    
    # Noise spectrum (diagonal in bin*bin space, independent of l and Bin)
    self.noise = self.rms_shear**2/self.noise

    # TEST
    #self.noise = 0 
    
    
    ###########
    # Read data
    ###########

    # If the file exists, initialize the fiducial values
    self.Cl_fid = np.zeros((self.nlmax,self.nbin,self.nbin),'float64')
    self.fid_values_exist = False
    if os.path.exists(self.data_directory+'/'+self.fiducial_file):
      self.fid_values_exist = True
      fid_file = open(self.data_directory+'/'+self.fiducial_file,'r')
      line = fid_file.readline()
      while line.find('#')!=-1:
	line = fid_file.readline()
      while (line.find('\n')!=-1 and len(line)==1):
	line = fid_file.readline()
      for nl in range(self.nlmax):
	for Bin in range(self.nbin):
	  for Bin2 in range(self.nbin):
	    self.Cl_fid[nl,Bin,Bin2] = float(line)
	    line = fid_file.readline()
      
    # Else the file will be created in the loglkl() function. 
    return

  # Galaxy distribution, returns the function D(z) from the notes
  def galaxy_distribution(self,z):

    zmean = 0.9
    z0 = zmean/1.412

    galaxy_dist = z**2*math.exp(-(z/z0)**(1.5))
    
    return galaxy_dist

  # Photo z distribution
  def photo_z_distribution(self,z,zph):

    # Standard error on dz/(1+z)
    sigma_ph = 0.05

    # Note: you must normalize it yourself to one if you want to get nice plots
    # of the galaxy distribution function in each bin (otherwise, the spectra
    # will remain correct, but each D_i(x) will loot strangely normalized when
    # compared to the original D(z)

    photo_z_dist = math.exp(-0.5*((z-zph)/sigma_ph/(1.+z))**2)/sigma_ph/(1.+z)/math.sqrt(2.*math.pi)
    
    return photo_z_dist

  def loglkl(self, cosmo, data):

    # One wants to obtain here the relation between z and r, this is done by
    # asking the cosmological module with the function z_of_r
    self.r = np.zeros(self.nzmax,'float64')
    self.dzdr= np.zeros(self.nzmax,'float64')

    self.r,self.dzdr =  cosmo.z_of_r(self.z)

    # Compute now the selection function eta(r) = eta(z) dz/dr normalized to one.
    self.eta_r = np.zeros((self.nzmax,self.nbin),'float64')
    for Bin in range(self.nbin):
      self.eta_r[:,Bin] = self.eta_z[:,Bin]*self.dzdr[:]/self.eta_norm[Bin]

    # Compute function g_i(r), that depends on r and the bin
    # g_i(r) = 2r(1+z(r)) int_0^+\infty drs eta_r(rs) (rs-r)/rs
    g = np.zeros((self.nzmax,self.nbin),'float64')

    for Bin in range(self.nbin):
      for nr in range(1,self.nzmax-1):
        g[nr,Bin] = np.sum(0.5*(self.eta_r[nr+1:,Bin]*(self.r[nr+1:]-self.r[nr])/self.r[nr+1:] + self.eta_r[nr:-1,Bin]*(self.r[nr:-1]-self.r[nr])/self.r[nr:-1])*(self.r[nr+1:]-self.r[nr:-1]))

        g[nr,Bin] *= 2.*self.r[nr]*(1.+self.z[nr])

    # Get power spectrum P(k=l/r,z(r)) from cosmological module
    pk = np.zeros((self.nlmax,self.nzmax),'float64')
    for index_l in range(self.nlmax):
      for index_z in range(1,self.nzmax):
        if (self.l[index_l]/self.r[index_z] > self.k_max):
          print 'you should increase euclid_lensing.k_max up to at least',self.l[index_l]/self.r[index_z]
          exit()
	pk[index_l,index_z] = cosmo.pk(self.l[index_l]/self.r[index_z],self.z[index_z])

    # Recover the non_linear scale computed by halofit. If no scale was
    # affected, set the scale to one, and make sure that the nuisance parameter
    # epsilon is set to zero
    k_sigma = np.zeros(self.nzmax, 'float64')
    if (cosmo.nonlinear_method == 0):
      k_sigma[:]=1.e6
    else :
      k_sigma = cosmo.nonlinear_scale(self.z,self.nzmax)

    # Define the alpha function, that will characterize the theoretical
    # uncertainty. Chosen to be 0.001 at low k, raise between 0.1 and 0.2 to
    # self.theoretical_error
    alpha = np.zeros((self.nlmax,self.nzmax),'float64')
    #self.theoretical_error = 0.1
    for index_l in range(self.nlmax):
      k = self.l[index_l]/self.r[1:]
      #alpha[index_l,1:] = (np.tanh(60*(k[:] - 0.5*k_sigma[1:])) + 1.)/(2./(self.theoretical_error-0.001)) + 0.001
      alpha[index_l,1:] = np.log(1.+k[:]/k_sigma[1:])/(1.+np.log(1.+k[:]/k_sigma[1:]))*self.theoretical_error

      #if index_l == 20:
        #for index_k in range(1,self.nzmax):
          #print k[index_k-1],alpha[index_l,index_k],math.log(1.+k[index_k-1]/k_sigma[index_k])/(1.+math.log(1.+k[index_k-1]/k_sigma[index_k]))*0.5

    # recover the e_th_nu part of the error function
    e_th_nu = self.coefficient_f_nu*cosmo.Omega_nu/cosmo.Omega_m

    # Compute the Error E_th_nu function
    E_th_nu = np.zeros((self.nlmax,self.nzmax),'float64')
    for index_l in range(1,self.nlmax):
      E_th_nu[index_l,:] = np.log(1. + self.l[index_l]/k_sigma[:]*self.r[:]) / (1. + np.log(1. + self.l[index_l]/k_sigma[:]*self.r[:])) * e_th_nu

    # Add the error function, with the nuisance parameter, to P_nl_th
    for index_l in range(self.nlmax):
      pk[index_l,:] *= (1. + data.mcmc_parameters['epsilon']['current']*data.mcmc_parameters['epsilon']['scale']*E_th_nu[index_l,:])
      #pk[index_l,:] *= (1. + alpha[index_l,:])

    # Start loop over l for computation of C_l^shear
    Cl_integrand = np.zeros((self.nzmax,self.nbin,self.nbin),'float64')
    Cl = np.zeros((self.nlmax,self.nbin,self.nbin),'float64')
    # Start loop over l for computation of E_l
    El_integrand = np.zeros((self.nzmax,self.nbin,self.nbin),'float64')
    El = np.zeros((self.nlmax,self.nbin,self.nbin),'float64')

    #for index_z in range(self.nzmax):
    #  print self.z[index_z],self.r[index_z],g[index_z,0]
    #exit() 

    for nl in range(self.nlmax):
  
      # find Cl_integrand = (g(r) / r)**2 * P(l/r,z(r))
      for Bin1 in range(self.nbin):
        for Bin2 in range(self.nbin):
          Cl_integrand[1:,Bin1,Bin2] = g[1:,Bin1] * g[1:,Bin2]/(self.r[1:]**2) * pk[nl,1:] 
          El_integrand[1:,Bin1,Bin2] = g[1:,Bin1] * g[1:,Bin2]/(self.r[1:]**2) * pk[nl,1:] * alpha[nl,1:]
      
      # Integrate over r to get C_l^shear_ij = P_ij(l)
      # C_l^shear_ij = 9/16 Omega0_m^2 H_0^4 \sum_0^rmax dr (g_i(r) g_j(r) /r**2) P(k=l/r,z(r))
      for Bin1 in range(self.nbin):
	for Bin2 in range(self.nbin):
          Cl[nl,Bin1,Bin2] = np.sum(0.5*(Cl_integrand[1:,Bin1,Bin2]+Cl_integrand[:-1,Bin1,Bin2])*(self.r[1:]-self.r[:-1]))
	  Cl[nl,Bin1,Bin2] *= 9./16.*(cosmo.Omega0_m())**2 # in units of Mpc**4
	  Cl[nl,Bin1,Bin2] *= (cosmo.h()/2997.9)**4 # dimensionless

          El[nl,Bin1,Bin2] = np.sum(0.5*(El_integrand[1:,Bin1,Bin2]+El_integrand[:-1,Bin1,Bin2])*(self.r[1:]-self.r[:-1]))
	  El[nl,Bin1,Bin2] *= 9./16.*(cosmo.Omega0_m())**2 # in units of Mpc**4
	  El[nl,Bin1,Bin2] *= (cosmo.h()/2997.9)**4 # dimensionless
	  if Bin1 == Bin2:
	    Cl[nl,Bin1,Bin2] += self.noise

    # Write fiducial model spectra if needed (exit in that case)
    if self.fid_values_exist is False:
      # Store the values now, and exit.
      fid_file = open(self.data_directory+'/'+self.fiducial_file,'w')
      fid_file.write('# Fiducial parameters')
      for key,value in data.mcmc_parameters.iteritems():
	fid_file.write(', %s = %.5g' % (key,value['current']*value['scale']))
      fid_file.write('\n')
      for nl in range(self.nlmax):
	for Bin1 in range(self.nbin):
	  for Bin2 in range(self.nbin):
	    fid_file.write("%.8g\n" % Cl[nl,Bin1,Bin2])
      print '\n\n /|\  Writting fiducial model in {0}'.format(self.data_directory+self.fiducial_file)
      print '/_o_\ for {0} likelihood'.format(self.name)
      return +1

    # Now that the fiducial model is stored, we add the El to both Cl and
    # Cl_fid (we create a new array, otherwise we would modify the self.Cl_fid
    # from one step to the other)
    #Cl          += El

    # Spline Cl[nl,Bin1,Bin2] along l
    ddCl     = np.zeros((self.nlmax,self.nbin,self.nbin),'float64')
    u_spline = np.zeros(self.nlmax,'float64')
    for Bin1 in range(self.nbin):
      for Bin2 in range(self.nbin):
        for nl in range(1,self.nlmax-1):
          sig_spline = (self.l[nl]-self.l[nl-1]) / (self.l[nl+1] - self.l[nl])
          p_spline   = sig_spline*ddCl[nl-1,Bin1,Bin2]+2.
          ddCl[nl,Bin1,Bin2] = (sig_spline-1.)/p_spline
          u_spline[nl] = (6.*((Cl[nl+1,Bin1,Bin2] - Cl[nl,Bin1,Bin2])/(self.l[nl+1]-self.l[nl]) - (Cl[nl,Bin1,Bin2]-Cl[nl-1,Bin1,Bin2])/(self.l[nl]-self.l[nl-1]))/(self.l[nl+1]-self.l[nl-1]) - sig_spline*u_spline[nl-1])/p_spline
        for nl in range(self.nlmax-2,-1,-1):
          ddCl[nl,Bin1,Bin2] = ddCl[nl,Bin1,Bin2]*ddCl[nl+1,Bin1,Bin2] + u_spline[nl]

    # Spline El[nl,Bin1,Bin2] along l
    ddEl     = np.zeros((self.nlmax,self.nbin,self.nbin),'float64')
    u_spline = np.zeros(self.nlmax,'float64')
    for Bin1 in range(self.nbin):
      for Bin2 in range(self.nbin):
        for nl in range(1,self.nlmax-1):
          sig_spline = (self.l[nl]-self.l[nl-1]) / (self.l[nl+1] - self.l[nl])
          p_spline   = sig_spline*ddEl[nl-1,Bin1,Bin2]+2.
          ddEl[nl,Bin1,Bin2] = (sig_spline-1.)/p_spline
          u_spline[nl] = (6.*((El[nl+1,Bin1,Bin2] - El[nl,Bin1,Bin2])/(self.l[nl+1]-self.l[nl]) - (El[nl,Bin1,Bin2]-El[nl-1,Bin1,Bin2])/(self.l[nl]-self.l[nl-1]))/(self.l[nl+1]-self.l[nl-1]) - sig_spline*u_spline[nl-1])/p_spline
        for nl in range(self.nlmax-2,-1,-1):
          ddEl[nl,Bin1,Bin2] = ddEl[nl,Bin1,Bin2]*ddEl[nl+1,Bin1,Bin2] + u_spline[nl]

    # Spline Cl_fid[nl,Bin1,Bin2]  along l
    ddCl_fid = np.zeros((self.nlmax,self.nbin,self.nbin),'float64')
    for Bin1 in range(self.nbin):
      for Bin2 in range(self.nbin):
	for nl in range(1,self.nlmax-1):
	  sig_spline = (self.l[nl]-self.l[nl-1]) / (self.l[nl+1] - self.l[nl])
	  p_spline   = sig_spline*ddCl_fid[nl-1,Bin1,Bin2]+2.
	  ddCl_fid[nl,Bin1,Bin2] = (sig_spline-1.)/p_spline
	  u_spline[nl] = (6.*((self.Cl_fid[nl+1,Bin1,Bin2] - self.Cl_fid[nl,Bin1,Bin2])/(self.l[nl+1]-self.l[nl]) - (self.Cl_fid[nl,Bin1,Bin2]-self.Cl_fid[nl-1,Bin1,Bin2])/(self.l[nl]-self.l[nl-1]))/(self.l[nl+1]-self.l[nl-1]) - sig_spline*u_spline[nl-1])/p_spline
	for nl in range(self.nlmax-2,-1,-1):
	  ddCl_fid[nl,Bin1,Bin2] = ddCl_fid[nl,Bin1,Bin2]*ddCl_fid[nl+1,Bin1,Bin2] + u_spline[nl]

    # Compute likelihood
    chi2 = 0.
    Cov_theory = np.zeros((self.nbin,self.nbin),'float64')
    Cov_observ = np.zeros((self.nbin,self.nbin),'float64')
    Cov_error  = np.zeros((self.nbin,self.nbin),'float64')

    # Prepare interpolation for every integer value of l, from the array
    # self.l, to finally compute the likelihood (sum over all l's)
    #dof = self.nbin*(self.nbin+1)/2.
    #dof = 1
    #dof = self.nbin**2
    #dof = self.nbin
    dof = 1./(int(self.l[-1])-int(self.l[0])+1)
    #print int(self.l[0]),int(self.l[-1]),int(self.l[-1])-int(self.l[0])+1

    for lll in range(int(self.l[0]),int(self.l[-1])+1):

      # Determine the closest non integer value.
      klo = 1
      khi = self.nlmax
      while (khi-klo > 1):
	k = (khi+klo)/2
	if (self.l[k] > lll):
	  khi = k
	else:
	  klo = k
      h = self.l[khi]-self.l[klo]
      if (h == 0.):
	print 'Problem: h=0 in splint of C[l,Bin]'
	exit()
      a = ( self.l[khi] - lll)/h
      b = (lll - self.l[klo])/h

      for Bin1 in range(self.nbin):
        Cov_theory[Bin1,:] = a*Cl[klo,Bin1,:] + b*Cl[khi,Bin1,:] + ((a**3-a)*ddCl[klo,Bin1,:] + (b**3-b)*ddCl[khi,Bin1,:])*(h**2)/6.
        Cov_observ[Bin1,:] = a*self.Cl_fid[klo,Bin1,:] + b*self.Cl_fid[khi,Bin1,:] + ((a**3-a)*ddCl_fid[klo,Bin1,:] + (b**3-b)*ddCl_fid[khi,Bin1,:])*(h**2)/6.
        Cov_error[Bin1,:]  = a*El[klo,Bin1,:] + b*El[khi,Bin1,:] + ((a**3-a)*ddEl[klo,Bin1,:] + (b**3-b)*ddEl[khi,Bin1,:])*(h**2)/6.

      #print lll,Cov_theory[0,0]-self.noise,Cov_theory[4,4]-self.noise,self.noise,Cov_error[0,0],Cov_error[4,4],math.sqrt(2./(2.*lll+1)*dof)*Cov_theory[0,0],math.sqrt(2./(2.*lll+1)*dof)*Cov_theory[4,4]

      #Debugging
      #Cov_theory *= (1e8)
      #Cov_observ *= (1e8)
      #Cov_error  *= (1e8)
      #Cov_error = np.zeros((self.nbin,self.nbin),'float64')
      #Cov_error = Cov_theory / 5
      #Cov_observ =  Cov_theory
      #if (lll == 300):
      #Cov_observ = Cov_theory+4.*Cov_error
      #else:
      #  Cov_observ = np.copy(Cov_theory)
      #if (lll==300):
      #Cov_observ += Cov_error
      #Cov_observ = Cov_theory
      #Cov_observ = Cov_theory*(1.+1e-2)

      det_theory = np.linalg.det(Cov_theory)
      det_observ = np.linalg.det(Cov_observ)

      if (self.theoretical_error > 0):

      #det_cross_obs = 0
      #for i in range(self.nbin):
      #	newCov = np.copy(Cov_theory)
      #	newCov[:,i] = Cov_observ[:,i]
      #	det_cross_obs += np.linalg.det(newCov)

        det_cross_err = 0
        for i in range(self.nbin):
          newCov = np.copy(Cov_theory)
          newCov[:,i] = Cov_error[:,i]
          det_cross_err += np.linalg.det(newCov)

      #det_cross_obs_err = 0
      #if (self.nbin>1):
      #  for i in range(self.nbin):
      #    for j in range(self.nbin):
      #      if (i!=j):
      #        newCov = np.copy(Cov_theory)
      #        newCov[:,i] = Cov_observ[:,i]
      #        newCov[:,j] = Cov_error[:,j]
      #        det_cross_obs_err += np.linalg.det(newCov)


      #epsilon_l = (2.*lll+1.)*self.fsky/2.*(det_cross_err*det_cross_obs/(det_theory)**2-(det_cross_err+det_cross_obs_err)/det_theory)/dof

        # Newton method
        # Find starting point for the method:
        #start = (Cov_observ[0,0] - Cov_theory[0,0])/Cov_error[0,0]*(1.+0.3)
        start = 0
        step = 0.001*det_theory/det_cross_err
        error = 1
        old_chi2 = -1.*data.boundary_loglike
        error_tol = 0.01
        epsilon_l = start
        #print 'starting'
        #print 'start is',start
        while error > error_tol:
          vector = np.array([epsilon_l-step,epsilon_l,epsilon_l+step])
        # Computing the function on three neighbouring points
          function_vector = np.zeros(3,'float64') 
          for k in range(3):
            Cov_theory_plus_error = Cov_theory + vector[k] * Cov_error

            det_theory_plus_error = np.linalg.det(Cov_theory_plus_error)

            det_theory_plus_error_cross_obs = 0
            for i in range(self.nbin):
              newCov = np.copy(Cov_theory_plus_error)
              newCov[:,i] = Cov_observ[:,i]
              det_theory_plus_error_cross_obs += np.linalg.det(newCov)
            function_vector[k] = (2.*lll+1.)*self.fsky*(det_theory_plus_error_cross_obs/det_theory_plus_error + math.log(det_theory_plus_error/det_observ) - self.nbin ) + dof*vector[k]**2 

          # Computing first
          first_d  = (function_vector[2]-function_vector[0]) / (vector[2]-vector[0])
          second_d = (function_vector[2]+function_vector[0]-2*function_vector[1]) / (vector[2]-vector[1])**2

          # Updating point and error
          epsilon_l = vector[1] - first_d/second_d
          error = abs(function_vector[1] - old_chi2)
          old_chi2 = function_vector[1]
        #print old_chi2,epsilon_l
      #print
        #print 'error and epsilon',error,epsilon_l
        #print epsilon_l,function_vector,error
      

      # End Newton
      
      #epsilon_l = 0
      #print det_cross_obs,self.nbin*det_theory
      #print det_cross_obs_err,det_cross_err*(self.nbin-1)

        Cov_theory_plus_error = Cov_theory + epsilon_l * Cov_error

        det_theory_plus_error = np.linalg.det(Cov_theory_plus_error)

        det_theory_plus_error_cross_obs = 0
        for i in range(self.nbin):
          newCov = np.copy(Cov_theory_plus_error)
          newCov[:,i] = Cov_observ[:,i]
          det_theory_plus_error_cross_obs += np.linalg.det(newCov)

        chi2 += (2.*lll+1.)*self.fsky*(det_theory_plus_error_cross_obs/det_theory_plus_error + math.log(det_theory_plus_error/det_observ) - self.nbin ) + dof*epsilon_l**2 

      else:

        #det_theory = np.linalg.det(Cov_theory)
        #det_observ = np.linalg.det(Cov_observ)
        det_cross  = 0.
        for i in range(self.nbin):
          newCov = np.copy(Cov_theory)
          newCov[:,i] = Cov_observ[:,i]
          det_cross += np.linalg.det(newCov)

        chi2 += (2.*lll+1.)*self.fsky*(det_cross/det_theory + math.log(det_theory/det_observ) - self.nbin)
      #print epsilon_l,chi2

    
    # Finally adding a gaussian prior on the epsilon nuisance parameter
    #print
    #print chi2
    #exit()
    chi2+=(data.mcmc_parameters['epsilon']['current']*data.mcmc_parameters['epsilon']['scale'])**2

    return -chi2/2.
