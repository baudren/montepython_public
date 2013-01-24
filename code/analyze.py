# Written by Benjamin Audren
# Plotting routine adapted from Karim Benabed, and the pmc code

import os,sys
import io
import math
import numpy as np

# Module for handling display
import matplotlib.pyplot as plt
# The root plotting module, to change options like font sizes, etc...
import matplotlib

# Module to handle warnings from matplotlib
import warnings

class info:

  def __init__(self,command_line):
  
    warnings.filterwarnings("error")
    # Check if the scipy module has the interpolate method correctly installed
    # (should be the case on every linux distribution with standard numpy)
    try:
      from scipy.interpolate import interp1d
      self.has_interpolate_module = True
    except ImportError:
      self.has_interpolate_module = False
      print('No cubic interpolation done (no interpolate method found in scipy), only linear')

    # At this points, Files could contain either a list of files (that could be
    # only one) or a folder.
    Files     = command_line.files
    binnumber = command_line.bins

    # Save the extension to output files
    self.extension = command_line.extension

    # Read a potential file describing changes to be done for the parameter
    # names, and number of paramaters plotted (can be let empty, all will then
    # be plotted). Initialized at empty structures.
    self.to_change  = {}
    self.to_plot    = []
    self.new_scales = {}

    if command_line.optional_plot_file is not None:
      for line in open(command_line.optional_plot_file[0],'r'):
	exec(line.replace('info.','self.'))

    # Prepare the files, according to the case, load the log.param, and
    # prepare the output (plots folder, .covmat, .info and .log files). After
    # this step, self.files will contain all chains.
    self.prepare(Files)

    # Compute the mean, maximum of likelihood, 1-sigma variance for this
    # main folder. This will create the self.spam chain
    self.convergence()

    # Create the main chain, which consists in all elements of self.spam put
    # together. This will serve for the  plotting.
    chain = np.copy(self.spam[0])
    for i in range(len(self.spam)-1):
      chain = np.append(chain,self.spam[i+1],axis=0)

    # In case of comparison, launch the prepare and convergence methods, with
    # an additional flag: is_main_chain=False. This will ensure that all the
    # information will not be stored to self.files, self.covmat... but to the
    # specified output, respectively comp_Files, comp_spam... 
    if command_line.comp is not None:
      comp_Files,comp_folder,comp_param = self.prepare(command_line.comp,is_main_chain = False)
      comp_spam,comp_ref_names,comp_tex_names,comp_backup_names,comp_plotted_parameters,comp_boundaries,comp_mean = self.convergence(is_main_chain = False,Files = comp_Files,param = comp_param)
      comp_mean = comp_mean[0]

      # Create comp_chain
      comp_chain = np.copy(comp_spam[0])
      for i in range(len(comp_spam)-1):
	comp_chain = np.append(comp_chain,comp_spam[i+1],axis=0)

    # Total number of steps.
    weight = sum(chain)[0] 

    # Covariance matrix computation (for the whole chain)
    self.mean   = self.mean[0]
    self.var    = self.var[0]
    self.covar  = np.zeros((len(self.ref_names),len(self.ref_names)))
    
    print('--> Computing covariance matrix')
    for i in range(len(self.ref_names)):
      for j in range(i,len(self.ref_names)):
        self.covar[i,j] = np.sum( chain[:,0]*((chain[:,i+2]-self.mean[i])*(chain[:,j+2]-self.mean[j])))/weight
	if i!=j:
	  self.covar[j,i]=self.covar[i,j]

    # Writing it out in name_of_folder.covmat
    self.cov.write('# ')
    for i in range(len(self.ref_names)):
      string = self.backup_names[i]
      if i != len(self.ref_names)-1:
        string+=','
      self.cov.write('%-16s' % string )
    self.cov.write('\n')
    # Removing scale factors in order to store true parameter covariance
    self.covar = np.dot(self.scales.T,np.dot(self.covar,self.scales))
    for i in range(len(self.ref_names)):
      for j in range(len(self.ref_names)):
        if self.covar[i][j]>0:
          self.cov.write(' %.5e\t' % self.covar[i][j])
        else:
          self.cov.write('%.5e\t' % self.covar[i][j])
      self.cov.write('\n')

    # Sorting by likelihood: a will hold the list of indices where the points
    # are sorted with increasing likelihood.
    a=chain[:,1].argsort(0)
    total=chain[:,0].sum()

    # Writing the best-fit model in name_of_folder.bestfit
    self.bf.write('# ')
    for i in range(len(self.ref_names)):
      string = self.backup_names[i]
      if i != len(self.ref_names)-1:
        string+=','
      self.bf.write('%-16s' % string )
    self.bf.write('\n')
    # Removing scale factors in order to store true parameter values
    for i in range(len(self.ref_names)):
      bfvalue = chain[a[0],2+i]*self.scales[i,i]
      if bfvalue>0:
          self.bf.write(' %.5e\t' % bfvalue)
      else:
          self.bf.write('%.5e\t' % bfvalue)
    self.bf.write('\n')                                       
                                          
    # Defining the sigma contours (1, 2 and 3-sigma)
    self.lvls = (68.26,95.4,99.7)

    # Computing 1,2 and 3-sigma errors, and plot. This will create the triangle
    # and 1d plot by default. 
    # If you also specified a comparison folder, it will create a versus plot
    # with the 1d comparison of all the common parameters, plus the 1d
    # distibutions for the others.
    self.bounds = np.zeros((len(self.ref_names),len(self.lvls),2))
    if command_line.plot == True:
      if command_line.comp is None:
	self.plot_triangle(chain,command_line,bin_number=binnumber,levels=self.lvls)
      else:
	self.plot_triangle(chain,command_line,bin_number=binnumber,levels=self.lvls,comp_chain=comp_chain,comp_ref_names = comp_ref_names,comp_tex_names = comp_tex_names,comp_backup_names = comp_backup_names,comp_plotted_parameters=comp_plotted_parameters,comp_folder = comp_folder,comp_boundaries = comp_boundaries,comp_mean = comp_mean)


    # Creating the array indices to hold the proper ordered list of plotted
    # parameters
    indices = []
    for i in range(len(self.plotted_parameters)):
      indices.append(self.ref_names.index(self.plotted_parameters[i]))

    print('--> Writing .info and .tex files')
    # Write down to the .h_info file all necessary information
    self.h_info.write(' param names\t:\t')
    self.v_info_names = []
    for i in indices:
      if self.scales[i,i] != 1: 
	if (float(self.scales[i,i]) > 100. or (self.scales[i,i]) < 0.01):
	  string = ' %0.e%s' % (1./self.scales[i,i],self.ref_names[i])
	elif (float(self.scales[i,i] < 1)):
	  string = ' %2d%s' % (1./self.scales[i,i],self.ref_names[i])
	else :
	  string = ' %2g%s' % (1./self.scales[i,i],self.ref_names[i])
      else:
	string = ' %s' % self.ref_names[i]
      self.v_info_names.append(string)
      self.h_info.write("%-16s" % string)

    self.write_h(self.h_info,indices,'R-1 values','%.6f',self.R)
    self.write_h(self.h_info,indices,'Best Fit  ','%.6e',chain[a[0],2:])
    self.write_h(self.h_info,indices,'mean      ','%.6e',self.mean)
    self.write_h(self.h_info,indices,'sigma     ','%.6e',(self.bounds[:,0,1]-self.bounds[:,0,0])/2.)
    self.h_info.write('\n')
    self.write_h(self.h_info,indices,'1-sigma - ','%.6e',self.bounds[:,0,0])
    self.write_h(self.h_info,indices,'1-sigma + ','%.6e',self.bounds[:,0,1])
    self.write_h(self.h_info,indices,'2-sigma - ','%.6e',self.bounds[:,1,0])
    self.write_h(self.h_info,indices,'2-sigma + ','%.6e',self.bounds[:,1,1])
    self.write_h(self.h_info,indices,'3-sigma - ','%.6e',self.bounds[:,2,0])
    self.write_h(self.h_info,indices,'3-sigma + ','%.6e',self.bounds[:,2,1])
    # bounds 
    self.h_info.write('\n')
    self.write_h(self.h_info,indices,'1-sigma > ','%.6e',self.mean+self.bounds[:,0,0])
    self.write_h(self.h_info,indices,'1-sigma < ','%.6e',self.mean+self.bounds[:,0,1])
    self.write_h(self.h_info,indices,'2-sigma > ','%.6e',self.mean+self.bounds[:,1,0])
    self.write_h(self.h_info,indices,'2-sigma < ','%.6e',self.mean+self.bounds[:,1,1])
    self.write_h(self.h_info,indices,'3-sigma > ','%.6e',self.mean+self.bounds[:,2,0])
    self.write_h(self.h_info,indices,'3-sigma < ','%.6e',self.mean+self.bounds[:,2,1])

    # Write vertical info file
    self.v_info.write('%-15s\t: %-6s %-10s %-10s %-10s %-11s %-10s %-11s %-10s %-10s %-10s %-10s %-10s' % ('param names','R-1','Best fit','mean','sigma','1-sigma -','1-sigma +','2-sigma -','2-sigma +','1-sigma >','1-sigma <','2-sigma >','2-sigma <'))
    for i in range(len(self.v_info_names)):
      name = self.v_info_names[i]
      #index = self.v_info_names.index(name)
      index = indices[i]
      self.v_info.write('\n%-15s\t: %.4f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e' % (name,
        self.R[index],chain[a[0],2:][index],
        self.mean[index],
        (self.bounds[:,0,1][index]-self.bounds[:,0,0][index])/2.,
        self.bounds[:,0,0][index],self.bounds[:,0,1][index],
        self.bounds[:,1,0][index],self.bounds[:,1,1][index],
        self.mean[index]+self.bounds[:,0,0][index],self.mean[index]+self.bounds[:,0,1][index],
        self.mean[index]+self.bounds[:,1,0][index],self.mean[index]+self.bounds[:,1,1][index]))

    # Writing the .tex file that will have this table prepared to be imported
    # in a tex document.
    self.write_tex(indices)

  # Scan the whole input folder, and include all chains in it (ending with
  # .txt to keep compatibility with CosmoMC)
  def prepare(self,Files,is_main_chain=True):

    # If the input command was an entire folder, then grab everything in it
    if os.path.isdir(Files[0]): 
      if Files[0][-1]!='/':
	Files[0]+='/'
      folder = Files[0]
      Files = [folder+elem for elem in os.listdir(folder) if (elem.find('.txt')!=-1 and elem.find('__')!=-1)]

    # Else, one needs to recover the folder, depending on the case
    else:
      if (len(Files[0].split('/'))==0 or (Files[0].split('/')[0]=='.')):
	folder = './'
      else:
	folder = ''
	for i in range(len(Files[0].split('/')[:-1])):
	  folder += Files[0].split('/')[i]+'/'

    # Remove too small files to potentially eliminate any problems of chains
    # being too short, and sub-folders (as the ./plots/ one that will be
    # created after the first run anyway
    for elem in np.copy(Files):
      if os.path.isdir('{0}'.format(elem)) is True:
        Files.remove(elem)
      # Note, this limit with the size is taylored for not too huge number of
      # parameters. Be aware that it might not work when having more than,
      # say, 20 free parameters.
      elif os.path.getsize(elem) < 600:   	  
        Files.remove(elem)

    # Check if the log.param file exists
    if os.path.isfile(folder+'log.param') is True:
      if os.path.getsize(folder+'log.param')>0:
	param = open(folder+'log.param','r')
      else:
	print('\n\n  The log param file {0} seems empty'.format(folder+'log.param'))
	exit()
    else:
      print('\n\n  The log param file {0} is absent ?'.format(folder+'log.param'))
      exit()

    # If the folder has no subdirectory, then go for a simple infoname,
    # otherwise, call it with the last name
    if (len(folder.split('/')) <= 2 and folder.split('/')[-1] == ''):
      v_infoname = folder+folder.rstrip('/')+'.v_info'
      h_infoname = folder+folder.rstrip('/')+'.h_info'
      texname  = folder+folder.rstrip('/')+'.tex'
      covname  = folder+folder.rstrip('/')+'.covmat'
      logname  = folder+folder.rstrip('/')+'.log'
      bfname  = folder+folder.rstrip('/')+'.bestfit'
    else:
      v_infoname = folder+folder.split('/')[-2]+'.v_info'
      h_infoname = folder+folder.split('/')[-2]+'.h_info'
      texname  = folder+folder.split('/')[-2]+'.tex'
      covname  = folder+folder.split('/')[-2]+'.covmat'
      logname  = folder+folder.split('/')[-2]+'.log'
      bfname  = folder+folder.split('/')[-2]+'.bestfit'

    # Distinction between the main chain and the comparative one, instead of
    # storing everything into the class, return it
    if is_main_chain:
      self.v_info  = open(v_infoname,'w')
      self.h_info  = open(h_infoname,'w')
      self.tex   = open(texname,'w')
      self.cov   = open(covname,'w')
      self.log   = open(logname,'w')
      self.bf   = open(bfname,'w')
      self.param = param

      self.Files = Files
      self.folder= folder
      return True
    else:
      return Files,folder,param


  def convergence(self,is_main_chain=True,Files=None,param=None):

    # We here have a list of files, that may be of length 1. If this
    # happens, then we split the only file in 3 subchains, otherwise we compute
    # normally the R coefficient. 
    spam = list()

    # Recovering default ordering of parameters 
    ref_names = []
    tex_names = []
    boundaries= []
    scales    = []

    # Backup names
    backup_names = []

    # Derived parameters
    derived_names = []
    derived_tex_names = []

    plotted_parameters = []

    if is_main_chain:
      Files = self.Files
      param = self.param

    # Recovering parameter names and scales, creating tex names,
    for line in param:
      if line.find('#')==-1:
	if line.find('data.experiments')!=-1:
	  self.experiments = line.split('=')[-1].replace('[','').replace(']','').replace('\n','').replace("'","").split(',')
	if line.find('data.parameters')!=-1:
	  name = line.split("'")[1]
          backup = name
          # Rename the names according the .extra file (opt)
	  if name in self.to_change.iterkeys():
	    name = self.to_change[name]
	  if (float(line.split('=')[-1].split(',')[-3].replace(' ','')) != 0 or str(line.split('=')[-1].split(',')[-1].replace(' ','').replace(']','').replace('\n','').replace("'","").replace("\t",'')) == 'derived' ):
            # The real name is always kept, to have still the class names in
            # the covmat
            backup_names.append(backup)
            if self.to_plot==[]:
              plotted_parameters.append(name)
	    else:
              if name in self.to_plot:
                plotted_parameters.append(name)
	    temp = [float(elem) for elem in line.split(",")[1:3]]
	    boundaries.append(temp)
	    ref_names.append(name)
	    scales.append(float(line.split('=')[-1].split(",")[4].replace(' ','')))
            if name in self.new_scales.iterkeys():
              scales[-1] = self.new_scales[name]
	    number = 1./scales[-1]
	    tex_names.append(io.get_tex_name(name,number=number))
    scales = np.diag(scales)
    param.seek(0)

    # Log param names for the main chain
    if is_main_chain:
      for elem in ref_names:
	self.log.write("%s   " % elem)
      self.log.write("\n")

    # Total number of steps done:
    total_number_of_steps = 0
    total_number_of_accepted_steps = 0

    # max_lkl will be appended with all the maximum likelihoods of files, then
    # will be replaced by its own maximum. This way, the global maximum
    # likelihood will be used as a reference, and not each chain's maximum.
    max_lkl = []

    # Circle through all files to find the maximum (and largest filename)
    length_of_largest_filename = 0
    print('--> Finding global maximum of likelihood')
    for File in Files:
      i=Files.index(File)
      if len(File.split('/')[-1])>length_of_largest_filename:
        length_of_largest_filename = len(File.split('/')[-1])
      # cheese will brutally contain everything in the chain File being scanned
      # Small trick, to analyze CosmoMC files directly, since the convention of
      # spacing is different, we have to test for the configuration of the
      # line. If it starts with three blanck spaces, it will be a CosmoMC file,
      # so every element will be separated with three spaces
      if line.startswith("   "):
        cheese = (np.array([[float(elem) for elem in line[4:].split()] for line in open(File,'r')]))
      # else it is the normal Monte Python convention
      else:
        cheese = (np.array([[float(elem) for elem in line.split()] for line in open(File,'r')]))

      max_lkl.append(min(cheese[:,1])) # beware, it is the min because we are talking about '- log likelihood'

    # Selecting only the true maximum.
    try:
      max_lkl = min(max_lkl)
    except ValueError:
      print('No decently sized chain was found. Please wait a bit to analyze this folder')
      exit()

    # Restarting the circling through files
    for File in Files:
      i=Files.index(File)
      # To improve presentation, and print only once the full path of the
      # analyzed folder, we recover the length of the path name, and create an
      # empty complementary string of this length
      index_slash = File.rfind('/')
      complementary_string = ''
      for j in range(index_slash+2):
        complementary_string+=' '
      if i ==0:
        exec "print '--> Scanning file %-{0}s' % File,".format(length_of_largest_filename)
      else:
        exec "print '                 %s%-{0}s' % (complementary_string,File.split('/')[-1]),".format(length_of_largest_filename)
      # cheese will brutally contain everything in the chain File being scanned
      cheese = (np.array([[float(elem) for elem in line.split()] for line in open(File,'r')]))
      local_max_lkl = min(cheese[:,1]) # beware, it is the min because we are talking about '- log likelihood'
      line_count = 0
      for line in open(File,'r'):
	line_count+=1
      if is_main_chain:
        self.log.write("%s\t Number of steps:%d\tSteps accepted:%d\tacc = %.2g\tmin(-loglike) = %.2f " % (File,sum(cheese[:,0]),line_count,line_count*1.0/sum(cheese[:,0]),local_max_lkl))
        self.log.write("\n")
        total_number_of_steps += sum(cheese[:,0])
        total_number_of_accepted_steps += line_count

      # Removing burn-in
      start = 0
      try:
	while cheese[start,1]>max_lkl+3:
	  start+=1
        print('  \t: Removed {0}\t points of burn-in'.format(start))
      except IndexError:
        print('  \t: Removed everything: chain not converged')


      # ham contains cheese without the burn-in, if there are any points left (more than 5)
      if np.shape(cheese)[0] > start+5:
	ham = np.copy(cheese[start::])

	# Deal with single file case
	if len(Files) == 1:
	  print('  Beware, convergence computed for a single file')
	  bacon   = np.copy(cheese[::3,:])
	  egg     = np.copy(cheese[1::3,:])
	  sausage = np.copy(cheese[2::3,:])

	  spam.append(bacon)
	  spam.append(egg)
	  spam.append(sausage)
	  continue
	
	# Adding resulting table to spam
	spam.append(ham)

    # Applying now new rules for scales
    for name in self.new_scales.iterkeys():
      index = ref_names.index(name)
      for i in range(len(spam)):
        spam[i][:,index+2] *= 1./scales[index,index]

    # Now that the list spam contains all the different chains removed of their
    # respective burn-in, proceed to the convergence computation

    # Test the length of the list
    if len(spam) == 0:
      print('No decently sized chain was found. Please wait a bit to analyze this folder')
      exit()

    # 2D arrays for mean and var, one column will contain the total (over all
    # chains) mean (resp.  variance), and each other column the respective
    # chain mean (resp. chain variance).  R only contains the values for each
    # parameter
    mean     = np.zeros((len(spam)+1,np.shape(spam[0])[1]-2))
    var      = np.zeros((len(spam)+1,np.shape(spam[0])[1]-2))
    R 	     = np.zeros(np.shape(spam[0])[1]-2)

    # Store the total number of points, and the total in each chain
    total    = np.zeros(len(spam)+1)
    for j in range(len(spam)):
      total[j+1] = np.sum(spam[j][:,0])
    total[0] = np.sum(total[1:])

    # Compute mean and variance for each chain
    print('--> Computing mean values')
    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
        submean     = np.sum(spam[j][:,0]*spam[j][:,i+2])
        mean[j+1,i] = submean/total[j+1]
        mean[0,i]  += submean
      mean[0,i] /= total[0]
    
    print('--> Computing variance')
    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
	var[0,i]     += np.sum(spam[j][:,0]*(spam[j][:,i+2]-mean[0,i]  )**2)
	var[j+1,i]    = np.sum(spam[j][:,0]*(spam[j][:,i+2]-mean[j+1,i])**2)/(total[j+1]-1)
      var[0,i]  /= (total[0]-1)
    
    # Gelman Rubin Diagnostic:
    # Computes a quantity linked to the ratio of the mean of the variances of
    # the different chains (within), and the variance of the means (between)
    # Note: This is not strictly speaking the Gelman Rubin test, defined for
    # same-length MC chains. Our quantity is defined without the square root,
    # which should not change much the result: a small sqrt(R) will still be a
    # small R. The same convention is used in CosmoMC, except for the weighted
    # average: we decided to do the average taking into account that longer
    # chains should count more
    within  = 0
    between = 0

    print('--> Computing convergence')
    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
        within  += total[j+1]*var[j+1,i]
        between += total[j+1]*(mean[j+1,i]-mean[0,i])**2
      within  /= total[0]
      between /= (total[0]-1)

      R[i] = between/within
      #if i == 0:
        #print ' -> R is ',R[i],'\tfor ',ref_names[i]
      #else:
        #print '         ',R[i],'\tfor ',ref_names[i]
      if i == 0:
        print ' -> R is %.6f' % R[i],'\tfor ',ref_names[i]
      else:
        print '         %.6f' % R[i],'\tfor ',ref_names[i]
    
    # Log finally the total number of steps, and absolute loglikelihood
    self.log.write("--> Total    number    of    steps: %d\n" % total_number_of_steps)
    self.log.write("--> Total number of accepted steps: %d\n" % total_number_of_accepted_steps)
    self.log.write("--> Minimum of -logLike           : %.2f" % max_lkl)

    # If the analysis is done with the main folder (and not the comparison
    # one), store all relevant quantities in the class.
    if is_main_chain:
      self.spam = spam

      self.ref_names    = ref_names
      self.tex_names    = tex_names
      self.boundaries   = boundaries

      self.backup_names = backup_names

      self.mean = mean
      self.var  = var
      self.R    = R

      self.scales = scales

      self.plotted_parameters = plotted_parameters

      return True
    else:
      return spam,ref_names,tex_names,backup_names,plotted_parameters,boundaries,mean

  # Plotting routine, also computes the sigma errors. Partly imported from
  # Karim Benabed in pmc.
  def plot_triangle(self,chain,command_line,bin_number=20,scales=(),levels=(68.26,95.4,99.7),aspect=(16,16),fig=None,tick_at_peak=False,comp_chain = None,comp_ref_names = None,comp_tex_names = None, comp_backup_names = None, comp_plotted_parameters = None, comp_folder = None,comp_boundaries = None,comp_mean = None):

    # If comparison is asked, don't plot 2d levels
    if command_line.comp is not None:
      plot_2d   = False
      comp      = True
      comp_done = False
    else:
      plot_2d = True
      comp    = False
      comp_done = False

    # Pre configuration of the output, note that changes to the font size will
    # occur later on as well, to obtain a nice scaling.
    matplotlib.rc('text',usetex = True)
    matplotlib.rc('font',size=11)
    matplotlib.rc('xtick',labelsize='8')
    matplotlib.rc('ytick',labelsize='8')
    lvls = np.array(levels)/100.

    # Create the figures
    if plot_2d:
      if fig:
	fig2d = plt.figure(fig,aspect)
      else:
	fig2d = plt.figure(1,figsize=aspect)

    exps = ''
    for exp in self.experiments:
      exps+=exp
      exps+=', '
    exps = exps[:-2].replace('_',' ')

    #plt.figtext(0.4,0.95,'Experiments: '+exps,fontsize=40,alpha=0.6)
    #plt.figtext(0.9, 0.7,'Monte Python',fontsize=70, rotation=90,alpha=0.15)
    fig1d = plt.figure(2,figsize=aspect)

    # clear figure
    plt.clf()

    # Recover the total number of parameters to potentially plot
    n = np.shape(chain)[1]-2
    if not scales:
     scales = np.ones(n)
    scales=np.array(scales)
    
    mean = self.mean*scales
    var  = self.var*scales**2

    # 1D plot
    max_values = np.max(chain[:,2:],axis=0)*scales
    min_values = np.min(chain[:,2:],axis=0)*scales
    span = (max_values-min_values)

    best_minus_lkl = np.min(chain[:,1],axis=0)

    if comp:
      comp_max_values = np.max(comp_chain[:,2:],axis=0)
      comp_min_values = np.min(comp_chain[:,2:],axis=0)
      comp_span = (comp_max_values-comp_min_values)

    # Define the place of ticks
    if tick_at_peak:
      pass
    else:
      ticks = np.array((min_values+span*0.1,(max_values+min_values)/2.,max_values-span*0.1)).T
      x_range = np.array((min_values,max_values)).T
      if comp:
	comp_ticks = np.array((comp_min_values+comp_span*0.1,(comp_max_values+comp_min_values)/2.,comp_max_values-comp_span*0.1)).T
	comp_x_range = np.array((comp_min_values,comp_max_values)).T
	for i in range(np.shape(comp_ticks)[0]):
	  if abs(comp_x_range[i][0]-comp_boundaries[i][0]) < comp_span[i]/bin_number :
	    comp_ticks[i][0] = comp_boundaries[i][0]
	    comp_x_range[i][0] = comp_boundaries[i][0]
	  if abs(comp_x_range[i][1]-comp_boundaries[i][1]) < comp_span[i]/bin_number :
	    comp_ticks[i][2] = comp_boundaries[i][1]
	    comp_x_range[i][1] = comp_boundaries[i][1]

    for i in range(np.shape(ticks)[0]):
      if abs(x_range[i][0]-self.boundaries[i][0]) < span[i]/bin_number :
	ticks[i][0] = self.boundaries[i][0]
	x_range[i][0] = self.boundaries[i][0]
      if abs(x_range[i][1]-self.boundaries[i][1]) < span[i]/bin_number :
	ticks[i][2] = self.boundaries[i][1]
	x_range[i][1] = self.boundaries[i][1]
      

    # Borders stuff, might need adjustement for printing on paper.
    fig1d.subplots_adjust(bottom=0.03, left=.07, right=0.98, top=0.93, hspace=.35)
    if plot_2d:
      fig2d.subplots_adjust(bottom=0.03, left=.07, right=0.98, top=0.93, hspace=.35)

    # In case of a comparison, figure out which names are shared, which are
    # unique and thus require a simple treatment.
    if comp:
      backup_comp_names = np.copy(comp_plotted_parameters)
      #print 'backup_comp_names is',backup_comp_names
      #print 'comp_backup_names is',comp_backup_names
      #print 'comp_ref_names is',comp_ref_names
      #print 'comp_tex_names is',comp_tex_names
      #print 'comp_plotted_parameters is',comp_plotted_parameters

      for i in range(len(self.plotted_parameters)):
	if self.plotted_parameters[i] in comp_plotted_parameters:
	  comp_plotted_parameters.remove(self.plotted_parameters[i])

      num_columns = int(round(math.sqrt(len(self.plotted_parameters) + len(comp_plotted_parameters)))) 
      num_lines   = int(math.ceil((len(self.plotted_parameters)+len(comp_plotted_parameters))*1.0/num_columns))
    else:
      num_columns = int(round(math.sqrt(len(self.plotted_parameters))))
      num_lines   = int(math.ceil(len(self.plotted_parameters)*1.0/num_columns))

    # Actual plotting
    print('-----------------------------------------------')
    for i in range(len(self.plotted_parameters)):

        print ' -> Computing histograms for ',self.plotted_parameters[i]

	index = self.ref_names.index(self.plotted_parameters[i])
	# Adding the subplots to the respective figures, this will be the diagonal for the triangle plot.
	if plot_2d:
	  ax2d=fig2d.add_subplot(len(self.plotted_parameters),len(self.plotted_parameters),i*(len(self.plotted_parameters)+1)+1,yticks=[])
        ax1d = fig1d.add_subplot(num_lines,num_columns,i+1,yticks=[])

	# normalized histogram
	hist,bin_edges=np.histogram(chain[:,index+2],bins=bin_number,weights=chain[:,0],normed=False)
	bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])

	# interpolated histogram (if available)
	interp_hist,interp_grid = self.cubic_interpolation(hist,bincenters)
	interp_hist /= np.max(interp_hist)

	if comp:
	  try:
	    # For the names in common, the following line will not output an
	    # error. Then compute the comparative histogram
	    #ii = np.where( comp_ref_names == self.plotted_parameters[i] )[0][0]
	    ii = comp_ref_names.index(self.plotted_parameters[i] )
	    comp_hist,comp_bin_edges = np.histogram(comp_chain[:,ii+2],bins=bin_number,weights=comp_chain[:,0],normed=False)
	    comp_bincenters = 0.5*(comp_bin_edges[1:]+comp_bin_edges[:-1])
	    interp_comp_hist,interp_comp_grid = self.cubic_interpolation(comp_hist,comp_bincenters)
	    interp_comp_hist /= np.max(interp_comp_hist)
	    comp_done = True
	  except ValueError : # If the name was not found, return the error. This will be then plotted at the end
	    comp_done = False
	if comp:
	  if not comp_done:
	    print('{0} was not found in the second folder'.format(self.plotted_parameters[i]))

	# minimum credible interval (method by Jan Haman). Fails for multimodal histograms
	bounds = self.minimum_credible_intervals(hist,bincenters,lvls)
	if bounds is False: # print out the faulty histogram (try reducing the binnumber to avoir this)
	  print(hist)
	else:
	  for elem in bounds:
	    for j in (0,1):
	      elem[j] -= self.mean[index]
	  self.bounds[index] = bounds

	if comp_done:
	  comp_bounds = self.minimum_credible_intervals(comp_hist,comp_bincenters,lvls)
	  if comp_bounds is False:
	    print(comp_hist)
	  else:
	    for elem in comp_bounds:
	      for j in (0,1):
		elem[j] -= comp_mean[ii]

	# plotting
	if plot_2d:
	  ax2d.set_xticks(ticks[index])
	  fontsize2d,ticksize2d = self.get_fontsize( len(self.tex_names) )
	  ax2d.set_xticklabels(['%.3g' % s for s in ticks[index]],fontsize=ticksize2d)
	  ax2d.set_title('%s= $%.3g^{+%.3g}_{%.3g}$' % (self.tex_names[index],self.mean[index],bounds[0][1],bounds[0][0]),fontsize=fontsize2d)
	  ax2d.plot(interp_grid,interp_hist,color='red',linewidth=2,ls='-')
	  ax2d.axis([x_range[index][0], x_range[index][1],0,1.05])

	fontsize1d,ticksize1d = self.get_fontsize( max(num_columns,num_lines))
	ax1d.set_title('%s= $%.3g^{+%.3g}_{%.3g}$' % (self.tex_names[index],self.mean[index],bounds[0][1],bounds[0][0]),fontsize=fontsize1d)
	ax1d.set_xticks(ticks[index])
	ax1d.set_xticklabels(['%.3g' % s for s in ticks[index]],fontsize=ticksize1d)
	ax1d.axis([x_range[index][0], x_range[index][1],0,1.05])
	
	if comp_done:
	  # complex variation of intervals
          ii = comp_ref_names.index(self.plotted_parameters[i])
	  if comp_x_range[ii][0] > x_range[index][0]:
	    comp_ticks[ii][0] = ticks[index][0]
	    comp_x_range[ii][0] = x_range[index][0]
	  if comp_x_range[ii][1] < x_range[index][1]:
	    comp_ticks[ii][2] = ticks[index][2]
	    comp_x_range[ii][1] = x_range[index][1]
	  comp_ticks[ii][1] = (comp_x_range[ii][1]+comp_x_range[ii][0])/2.
	  ax1d.set_xticks(comp_ticks[ii])
	  ax1d.set_xticklabels(['%.3g' % s for s in comp_ticks[ii]],fontsize=ticksize1d)
	  ax1d.axis([comp_x_range[ii][0], comp_x_range[ii][1],0,1.05])


	ax1d.plot(interp_grid,interp_hist,color='black',linewidth=2,ls='-')
	if comp_done:
	  ax1d.plot(interp_comp_grid,interp_comp_hist,color='red',linewidth=2,ls='-')

	# mean likelihood (optional, if comparison, it will not be printed)
        if (plot_2d and command_line.mean_likelihood):
          try:
            lkl_mean=np.zeros(len(bincenters),'float64')
            norm=np.zeros(len(bincenters),'float64')
            for j in range(len(bin_edges)-1):
              tmp = np.array([elem for elem in chain[:,:] if (elem[index+2]>=bin_edges[j] and elem[index+2]<=bin_edges[j+1])],'float')
              lkl_mean[j] += np.sum(np.exp( best_minus_lkl - tmp[:,1])*tmp[:,0])
              norm[j] += np.sum(tmp,axis=0)[0]
            lkl_mean /= norm
            #lkl_mean *= max(hist)/max(lkl_mean)
            lkl_mean /= max(lkl_mean)
            interp_lkl_mean,interp_grid = self.cubic_interpolation(lkl_mean,bincenters)
            ax2d.plot(interp_grid,interp_lkl_mean,color='red',ls='--',lw=2)
            ax1d.plot(interp_grid,interp_lkl_mean,color='red',ls='--',lw=4)
          except:
            print 'could not find likelihood contour for ',self.plotted_parameters[i]

	if command_line.subplot is True:
	  if not comp:
	    extent2d = ax2d.get_window_extent().transformed(fig2d.dpi_scale_trans.inverted())
	    fig2d.savefig(self.folder+'plots/{0}_{1}.{2}'.format(self.folder.split('/')[-2],self.plotted_parameters[index],self.extension), bbox_inches=extent2d.expanded(1.1, 1.4))
	  else:
	    extent1d = ax1d.get_window_extent().transformed(fig1d.dpi_scale_trans.inverted())
	    fig1d.savefig(self.folder+'plots/{0}_{1}.{2}'.format(self.folder.split('/')[-2],self.plotted_parameters[index],self.extension), bbox_inches=extent1d.expanded(1.1, 1.4))

	# Now do the rest of the triangle plot
	if plot_2d:
	  for j in range(i):
	    second_index = self.ref_names.index(self.plotted_parameters[j])
	    ax2dsub=fig2d.add_subplot(len(self.plotted_parameters),len(self.plotted_parameters),(i)*len(self.plotted_parameters)+j+1)
	    n,xedges,yedges=np.histogram2d(chain[:,index+2],chain[:,second_index+2],weights=chain[:,0],bins=(bin_number,bin_number),normed=False)
	    extent = [x_range[second_index][0], x_range[second_index][1], x_range[index][0], x_range[index][1]]
	    x_centers = 0.5*(xedges[1:]+xedges[:-1])
	    y_centers = 0.5*(yedges[1:]+yedges[:-1])

	    ax2dsub.set_xticks(ticks[second_index])
	    if i == len(self.plotted_parameters)-1:
	      ax2dsub.set_xticklabels(['%.3g' % s for s in ticks[second_index]],fontsize=ticksize2d)
	    else:
	      ax2dsub.set_xticklabels([''])

	    ax2dsub.set_yticks(ticks[index])
	    if j == 0:
	      ax2dsub.set_yticklabels(['%.3g' % s for s in ticks[index]],fontsize=ticksize2d)
	    else:
	      ax2dsub.set_yticklabels([''])
	    ax2dsub.imshow(n, extent=extent, aspect='auto',interpolation='gaussian',origin='lower',cmap=matplotlib.cm.Reds)

	    # plotting contours, using the ctr_level method (from Karim Benabed)
            try:
              cs = ax2dsub.contour(y_centers,x_centers,n,extent=extent,levels=self.ctr_level(n,lvls),colors="k",zorder=5)
              #ax2dsub.clabel(cs, cs.levels[:-1], inline=True,inline_spacing=0, fmt=dict(zip(cs.levels[:-1],[r"%d \%%"%int(l*100) for l in lvls[::-1]])), fontsize=19)
            except Warning:
              print('     /!\  The routine could not find the contour of the "%s-%s" 2d-plot'% (self.plotted_parameters[i],self.plotted_parameters[j]))
              pass

	    if command_line.subplot is True:
	      # Store the individual 2d plots
	      fig_temp=plt.figure(3,figsize=(6,6))
	      fig_temp.clf()
	      ax_temp=fig_temp.add_subplot(111)
	      ax_temp.set_xticks(ticks[second_index])
	      ax_temp.set_yticks(ticks[index])
	      ax_temp.imshow(n, extent=extent, aspect='auto',interpolation='gaussian',origin='lower',cmap=matplotlib.cm.Reds)
	      ax_temp.set_xticklabels(['%.3g' % s for s in ticks[second_index]],fontsize=ticksize2d)
	      ax_temp.set_yticklabels(['%.3g' % s for s in ticks[index]],fontsize=ticksize2d)
	      ax_temp.set_title('%s vs %s' % (self.tex_names[index],self.tex_names[second_index]),fontsize=fontsize1d)
              try:
                cs = ax_temp.contour(y_centers,x_centers,n,extent=extent,levels=self.ctr_level(n,lvls),colors="k",zorder=5)
                #ax_temp.clabel(cs, cs.levels[:-1], inline=True,inline_spacing=0, fmt=dict(zip(cs.levels[:-1],[r"%d \%%"%int(l*100) for l in lvls[::-1]])), fontsize=19)
              except Warning:
                print('     /!\  The routine could not find the contour of the "%s-%s" 2d-plot'% (self.plotted_parameters[i],self.plotted_parameters[j]))
                pass

	      fig_temp.savefig(self.folder+'plots/{0}_2d_{1}-{2}.{3}'.format(self.folder.split('/')[-2],self.ref_names[index],self.ref_names[second_index],self.extension))

	      # store the coordinates of the points for further plotting.
	      plot_file = open(self.folder+'plots/{0}_2d_{1}-{2}.dat'.format(self.folder.split('/')[-2],self.ref_names[index],self.ref_names[second_index]),'w')
	      plot_file.write('# contour for confidence level {0}\n'.format(levels[2]))
	      for elem in cs.collections[0].get_paths():
		points = elem.vertices
		for k in range(np.shape(points)[0]):
		  plot_file.write("%.8g\t %.8g\n" % (points[k,0],points[k,1]))
	      plot_file.write("\n\n")

	      plot_file.write('# contour for confidence level {0}\n'.format(levels[1]))
	      for elem in cs.collections[1].get_paths():
		points = elem.vertices
		for k in range(np.shape(points)[0]):
		  plot_file.write("%.8g\t %.8g\n" % (points[k,0],points[k,1]))
	      plot_file.write("\n\n")

	      plot_file.write('# contour for confidence level {0}\n'.format(levels[0]))
	      for elem in cs.collections[2].get_paths():
		points = elem.vertices
		for k in range(np.shape(points)[0]):
		  plot_file.write("%.8g\t %.8g\n" % (points[k,0],points[k,1]))
	      plot_file.close()

    # Plot the remaining 1d diagram for the parameters only in the comp folder
    if comp:
      #if len(self.plotted_parameters) == len(self.ref_names):
      for i in range(len(self.plotted_parameters),len(self.plotted_parameters)+len(comp_plotted_parameters)):

        ax1d = fig1d.add_subplot(num_lines,num_columns,i+1,yticks=[])
        ii   = comp_ref_names.index(comp_plotted_parameters[i-len(self.plotted_parameters)])

        comp_hist,comp_bin_edges = np.histogram(comp_chain[:,ii+2],bins=bin_number,weights=comp_chain[:,0],normed=False)
        comp_bincenters = 0.5*(comp_bin_edges[1:]+comp_bin_edges[:-1])
        interp_comp_hist,interp_comp_grid = self.cubic_interpolation(comp_hist,comp_bincenters)
        interp_comp_hist /= np.max(interp_comp_hist)

        comp_bounds = self.minimum_credible_intervals(comp_hist,comp_bincenters,lvls)
        if comp_bounds is False:
          print(comp_hist)
        else:
          for elem in comp_bounds:
            for j in (0,1):
              elem[j] -= comp_mean[ii]
        ax1d.set_xticks(comp_ticks[ii])
        ax1d.set_xticklabels(['%.3g' % s for s in comp_ticks[ii]],fontsize=ticksize1d)
        ax1d.axis([comp_x_range[ii][0], comp_x_range[ii][1],0,1.05])
        ax1d.set_title('%s= $%.3g^{+%.3g}_{%.3g}$' % (comp_tex_names[ii],comp_mean[ii],comp_bounds[0][1],comp_bounds[0][0]),fontsize=fontsize1d)
        ax1d.plot(interp_comp_grid,interp_comp_hist,color='red',linewidth=2,ls='-')
        if command_line.subplot is True:
          extent1d = ax1d.get_window_extent().transformed(fig1d.dpi_scale_trans.inverted())
          fig1d.savefig(self.folder+'plots/{0}_{1}.{2}'.format(self.folder.split('/')[-2],self.ref_names[i],self.extension), bbox_inches=extent1d.expanded(1.1, 1.4))
	
    # If plots/ folder in output folder does not exist, create it
    if os.path.isdir(self.folder+'plots') is False:
      os.mkdir(self.folder+'plots')
    print('-----------------------------------------------')
    print('--> Saving figures to .{0} files'.format(self.extension))
    if plot_2d:
      fig2d.savefig(self.folder+'plots/{0}_triangle.{1}'.format(self.folder.split('/')[-2],self.extension),  bbox_inches=0)
    if comp:
      fig1d.savefig(self.folder+'plots/{0}-vs-{1}.{2}'.format(self.folder.split('/')[-2],comp_folder.split('/')[-2],self.extension), bbox_inches=0)
    else:
      fig1d.savefig(self.folder+'plots/{0}_1d.{1}'.format(self.folder.split('/')[-2],self.extension), bbox_inches=0)


  # Extract the contours for the 2d plots (KB)
  def ctr_level(self,histogram2d,lvl,infinite = False):

    hist=histogram2d.flatten()*1.
    hist.sort()
    cum_hist=np.cumsum(hist[::-1])
    cum_hist/=cum_hist[-1]

    alvl=np.searchsorted(cum_hist,lvl)[::-1]
    clist=[0]+[hist[-ii] for ii in alvl]+[np.max(hist)]
    if not infinite:
      return clist[1:]
    return clist

  # Extract minimum credible intervals (method from Jan Haman)
  def minimum_credible_intervals(self,histogram,bincenters,levels):
    bounds = np.zeros((len(levels),2))
    j = 0
    delta = bincenters[1]-bincenters[0]
    left_edge = max(histogram[0] - 0.5*(histogram[1]-histogram[0]),0.)
    right_edge = max(histogram[-1] + 0.5*(histogram[-1]-histogram[-2]),0.)
    failed = False
    for level in levels:
      norm = float((sum(histogram)-0.5*(histogram[0]+histogram[-1]))*delta)
      norm += 0.25*(left_edge+histogram[0])*delta
      norm += 0.25*(right_edge+histogram[-1])*delta
      water_level_up   = max(histogram)*1.0
      water_level_down = min(histogram)*1.0
      top = 0.
      
      ii=0
      while ((abs((top/norm)-level) > 0.0001) and not failed):
	top=0.
	water_level = (water_level_up + water_level_down)/2.
	ontop = [elem for elem in histogram if elem>water_level]
	indices = [i for i in range(len(histogram)) if histogram[i]>water_level]
	# check for multimodal posteriors
	if ((indices[-1]-indices[0]+1)!=len(indices)):
	  print('    /!\ Warning: could not derive minimum credible intervals (multimodal posterior)')
	  failed = True
	  break
	top = (sum(histogram[indices])-0.5*(histogram[indices[0]]+histogram[indices[-1]]))*(delta)

	# left
	if indices[0]>0:
	  top += 0.5 * (water_level + histogram[indices[0]]) * delta *(histogram[indices[0]]-water_level)/(histogram[indices[0]]-histogram[indices[0]-1])
	else:
	  if (left_edge > water_level):
	    top += 0.25*(left_edge+histogram[indices[0]])*delta
	  else:
	    top += 0.25 * (water_level + histogram[indices[0]]) * delta * (histogram[indices[0]]-water_level)/(histogram[indices[0]]-left_edge)

	# right
	if indices[-1]<(len(histogram)-1) :
	  top += 0.5 * (water_level + histogram[indices[-1]]) * (delta)*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-histogram[indices[-1]+1])
	else:
	  if (right_edge > water_level):
	    top += 0.25*(right_edge+histogram[indices[-1]])*delta
	  else:
	    top += 0.25 * (water_level + histogram[indices[-1]]) * delta * (histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-right_edge)

	if top/norm >= level:
	  water_level_down = water_level
	else:
	  water_level_up = water_level
	# safeguard, just in case
	ii+=1
	if (ii>1000):
          print('    /!\ Warning: the loop to check for sigma deviations was too long to converge')
	  break

      #print top,norm,abs(top/norm)
      #print bincenters[indices]
      #print histogram[indices],water_level
      #print

      # min
      if indices[0]>0:
	bounds[j][0] = bincenters[indices[0]] - delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-histogram[indices[0]-1])
      else:
	if (left_edge > water_level):          
	  bounds[j][0] = bincenters[0]-0.5*delta
	else:
	  bounds[j][0] = bincenters[indices[0]] - 0.5*delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-left_edge)

      # max
      if indices[-1]<(len(histogram)-1):
	bounds[j][1] = bincenters[indices[-1]]+ delta*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-histogram[indices[-1]+1])
      else:
	if (right_edge > water_level):
	  bounds[j][1] = bincenters[-1]+0.5*delta
	else:
	  bounds[j][1] = bincenters[indices[-1]]+ 0.5*delta*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-right_edge)

      j+=1
	
    #print
    return bounds

  def write_h(self,file,indices,name,string,quantity,modifiers=None):
    file.write('\n '+name+'\t:\t')
    for i in indices:
      if quantity[i] >= 0:
        space_string = ' '
      else:
        space_string = ''
      file.write(space_string+string % quantity[i]+'\t')

  def write_tex(self,indices):
    
    self.tex.write("\documentclass{article}\n")
    self.tex.write("\\begin{document}\n")
    self.tex.write("\\begin{tabular}{|")
    for i in indices:
      self.tex.write("c")
    self.tex.write("|}\n \hline")

    for i in indices:
      self.tex.write("%s " % self.tex_names[i])
      if i!=indices[-1]:
        self.tex.write(" & ")
      else:
        self.tex.write(" \\\\ \hline\n")

    for i in indices:
      self.tex.write("$%.4g_{%.2g}^{+%.2g}$" % (self.mean[i],self.bounds[i,0,0],self.bounds[i,0,1]))
      if i!=indices[-1]:
        self.tex.write(" & ")
      else:
        self.tex.write(" \\\\ \hline\n")

    self.tex.write("\\end{tabular}\n")
    self.tex.write("\\end{document}")

  def cubic_interpolation(self,hist,bincenters):
    if self.has_interpolate_module:
      interp_grid = np.linspace(bincenters[0],bincenters[-1],len(bincenters)*10)
      from scipy.interpolate import interp1d
      f = interp1d(bincenters,hist,kind='cubic')
      interp_hist = f(interp_grid)
      return interp_hist,interp_grid
    else:
      return hist,bincenters

  # Empirical method to adjust font size on the plots to fit the number of
  # parameters. Feel free to modify to your needs.
  def get_fontsize(self,diag_length):
    # Approximate values to have roughly nice displays font size
    #fontsize = round( 19 - (diag_length-5)*1.38)
    #ticksize = round( 14 - (diag_length-5)*1)
    # If the above does not work, please fix the values with the following two
    # lines (and commenting the above)
    fontsize = 15
    ticksize = 14
    return fontsize,ticksize
