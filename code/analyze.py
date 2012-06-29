import os,sys
import re
import io
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


class info:

  def __init__(self,command_line):
  
    # Check if the scipy module has the interpolate method correctly installed
    # (should be the case on every computer...)
    try:
      from scipy.interpolate import interp1d
      self.has_interpolate_module = True
    except ImportError:
      self.has_interpolate_module = False
      print 'No cubic interpolation done (no interpolate method found in scipy), only linear'

    # At this points, Files could contain either a list of files (that could be
    # only one) or a folder.
    Files     = command_line.files
    binnumber = command_line.bins

    # Prepare the files, according to the case, load the log.param, and
    # prepare the output (plots folder, .covmat, .info and .log files). After
    # this step, self.files will contain all chains.
    self.prepare(Files)
    # And compute the mean, maximum of likelihood, 1-sigma variance for this
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
      comp_spam,comp_ref_names,comp_tex_names,comp_boundaries,comp_mean = self.convergence(is_main_chain = False,Files = comp_Files,param = comp_param)
      comp_mean = comp_mean[0]
      # Create comp_chain
      comp_chain = np.copy(comp_spam[0])
      for i in range(len(comp_spam)-1):
	comp_chain = np.append(comp_chain,comp_spam[i+1],axis=0)


    weight = sum(chain)[0] # Total number of steps.

    # Covariance matrix computation (for the whole chain)
    self.mean   = self.mean[0]
    self.var    = self.var[0]
    self.covar  = np.zeros((len(self.ref_names),len(self.ref_names)))
    
    for i in range(len(self.ref_names)):
      for j in range(i,len(self.ref_names)):
	for k in range(np.shape(chain)[0]):
	  self.covar[i][j]+=chain[k][0]*((chain[k][i+2]-self.mean[i])*(chain[k][j+2]-self.mean[j]))
	self.covar[i][j]  /= weight
	if i!=j:
	  self.covar[j][i]=self.covar[i][j]

    self.cov.write('# ')
    for i in range(len(self.ref_names)):
      string = self.ref_names[i]
      if i != len(self.ref_names)-1:
        string+=','
      self.cov.write('%-16s' % string )
    self.cov.write('\n')
    self.covar = np.dot(self.scales.T,np.dot(self.covar,self.scales))
    for i in range(len(self.ref_names)):
      for j in range(len(self.ref_names)):
        if self.covar[i][j]>0:
          self.cov.write(' %.5e\t' % self.covar[i][j])
        else:
          self.cov.write('%.5e\t' % self.covar[i][j])
      self.cov.write('\n')

    # Sorting by likelihood
    a=chain[:,1].argsort(0)
    total=chain[:,0].sum()
    self.lvls = [total*0.6826,total*0.954,total*0.997]

    # Computing 1,2 and 3-sigma errors, and plot. This will create the triangle
    # plot and the 1d by default. If you also specified a comparison folder, it
    # will create a versus plot with the 1d comparison of all the common
    # parameters, plus the 1d distibutions for the others.
    self.bounds = np.zeros((len(self.ref_names),len(self.lvls),2))
    if command_line.plot == True:
      if command_line.comp is None:
	self.plot_triangle(chain,command_line,bin_number=binnumber)
      else:
	self.plot_triangle(chain,command_line,bin_number=binnumber,comp_chain=comp_chain,comp_ref_names = comp_ref_names,comp_tex_names = comp_tex_names,comp_folder = comp_folder,comp_boundaries = comp_boundaries,comp_mean = comp_mean)

    # Write down to the .info file all necessary information
    self.info.write('\n param names:\t')
    for i in range(len(self.ref_names)):
      if self.scales[i,i] != 1: 
	if (float(self.scales[i,i]) > 100. or (self.scales[i,i]) < 0.01):
	  string = '%0.e%s' % (1./self.scales[i,i],self.ref_names[i])
	elif (float(self.scales[i,i] < 1)):
	  string = '%2d%s' % (1./self.scales[i,i],self.ref_names[i])
	else :
	  string = '%2g%s' % (1./self.scales[i,i],self.ref_names[i])
      else:
	string = '%s' % self.ref_names[i]
      self.info.write("%-16s" % string)

    #for elem in self.ref_names:
      #self.info.write('%s\t' % (elem))
    self.info.write('\n R-1 values:\t')
    for elem in self.R:
      self.info.write('%.6f\t' % (elem-1.))
    self.info.write('\n Best Fit:\t')
    for elem in chain[a[0],2:]:
      self.info.write('%.6e\t' % (elem))
    self.info.write('\n mean    :\t')
    for elem in self.mean:
      self.info.write('%.6e\t' % (elem))
    self.info.write('\n\n 1-sigma - :\t')
    for elem in self.bounds:
      self.info.write('%.6e\t' % (elem[0][0]))
    self.info.write('\n 1-sigma + :\t')
    for elem in self.bounds:
      self.info.write(' %.6e\t' % (elem[0][1]))
    self.info.write('\n 2-sigma - :\t')
    for elem in self.bounds:
      self.info.write('%.6e\t' % (elem[1][0]))
    self.info.write('\n 2-sigma + :\t')
    for elem in self.bounds:
      self.info.write(' %.6e\t' % (elem[1][1]))
    self.info.write('\n 3-sigma - :\t')
    for elem in self.bounds:
      self.info.write('%.6e\t' % (elem[2][0]))
    self.info.write('\n 3-sigma + :\t')
    for elem in self.bounds:
      self.info.write(' %.6e\t' % (elem[2][1]))

    # bounds 
    self.info.write('\n\n 1-sigma > :\t')
    for i in range(np.shape(self.bounds)[0]):
      self.info.write('%.6e\t' % (self.mean[i]+self.bounds[i,0,0]))
    self.info.write('\n 1-sigma < :\t')
    for i in range(np.shape(self.bounds)[0]):
      self.info.write('%.6e\t' % (self.mean[i]+self.bounds[i,0,1]))
    self.info.write('\n 2-sigma > :\t')
    for i in range(np.shape(self.bounds)[0]):
      self.info.write('%.6e\t' % (self.mean[i]+self.bounds[i,1,0]))
    self.info.write('\n 2-sigma < :\t')
    for i in range(np.shape(self.bounds)[0]):
      self.info.write('%.6e\t' % (self.mean[i]+self.bounds[i,1,1]))
    self.info.write('\n 3-sigma > :\t')
    for i in range(np.shape(self.bounds)[0]):
      self.info.write('%.6e\t' % (self.mean[i]+self.bounds[i,2,0]))
    self.info.write('\n 3-sigma < :\t')
    for i in range(np.shape(self.bounds)[0]):
      self.info.write('%.6e\t' % (self.mean[i]+self.bounds[i,2,1]))

  def prepare(self,Files,is_main_chain=True):

    # Scan the whole input folder, and include all chains in it (all files with
    # a '.' in their name will be discarded as not mcmc chains (If you want to
    # have additionnal files than the .param in this folder, please keep this
    # in  mind

    if os.path.isdir(Files[0]): 
      if Files[0][-1]!='/':
	Files[0]+='/'
      folder = Files[0]
      Files = [folder+elem for elem in os.listdir(folder) if (elem.find('.txt')!=-1 and elem.find('__')!=-1)]
      for elem in np.copy(Files):
	if os.path.isdir('{0}'.format(elem)) is True:
	  Files.remove(elem)
	elif os.path.getsize(elem) < 600: # If the file is too short in lines, remove it from the list
	  Files.remove(elem)

    # Recover the folder, depending on the case
    else:
      if (len(Files[0].split('/'))==0 or (Files[0].split('/')[0]=='.')):
	folder = './'
      else:
	folder = ''
	for i in range(len(Files[0].split('/')[:-1])):
	  folder += Files[0].split('/')[i]+'/'
	for elem in np.copy(Files):
	  if os.path.isdir('{0}'.format(elem)) is True:
	    Files.remove(elem)
	  elif os.path.getsize(elem) < 600:
	    Files.remove(elem)

    # Check if the log.param file exists
    if os.path.isfile(folder+'log.param') is True:
      if os.path.getsize(folder+'log.param')>0:
	param = open(folder+'log.param','r')
      else:
	print '\n\n  The log param file {0} seems empty'.format(folder+'log.param')
	exit()
    else:
      print '\n\n  The log param file {0} is absent ?'.format(folder+'log.param')
      exit()

    # If the folder has no subdirectory, then go for a simple infoname,
    # otherwise, call it with the last name
    if (len(folder.split('/')) <= 2 and folder.split('/')[-1] == ''):
      infoname = folder+folder.rstrip('/')+'.info'
      covname  = folder+folder.rstrip('/')+'.covmat'
      logname  = folder+folder.rstrip('/')+'.log'
    else:
      infoname = folder+folder.split('/')[-2]+'.info'
      covname  = folder+folder.split('/')[-2]+'.covmat'
      logname  = folder+folder.split('/')[-2]+'.log'

    # Distinction between the main chain and the comparative one, instead of
    # storing everything into the class, return it
    if is_main_chain:
      self.info  = open(infoname,'w')
      self.cov   = open(covname,'w')
      self.log   = open(logname,'w')
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

    derived_names = []
    derived_tex_names = []

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
	  if (float(line.split('=')[-1].split(',')[-3].replace(' ','')) != 0 or str(line.split('=')[-1].split(',')[-1].replace(' ','').replace(']','').replace('\n','').replace("'","").replace("\t",'')) == 'derived' ):
	    temp = [float(elem) for elem in line.split(",")[1:3]]
	    boundaries.append(temp)
	    ref_names.append(name)
	    scales.append(float(line.split('=')[-1].split(",")[4].replace(' ','')))
	    number = 1./scales[-1]
	    tex_names.append(io.get_tex_name(name,number=number))
	#if line.find('data.derived_parameters_list')!=-1:
	  #for elem in line.split('=')[-1].strip(' ').replace('[','').replace(']','').replace('\n','').replace("'","").split(','):
	    #ref_names.append(elem)
	    #tex_names.append(io.get_tex_name(elem,1))
	    #boundaries.append([-1.,-1.])
	    #scales.append(1.)
	  #print derived_names,derived_tex_names
	  #exit()
    scales = np.diag(scales)
    param.seek(0)

    # Log param names for the main chain
    if is_main_chain:
      for elem in ref_names:
	self.log.write("%s   " % elem)
      self.log.write("\n")

    # Total number of steps done:
    total = 0

    # max_lkl will be appended with all the maximum likelihoods of files, then
    # will be replaced by its own maximum. This way, the global maximum
    # likelihood will be used as a reference, and not each chain's maximum.
    max_lkl = []

    # Circle through all files to find the maximum
    print 'Finding global maximum of likelihood'
    for File in Files:
      i=Files.index(File)
      # cheese will brutally contain everything in the chain File being scanned
      cheese = (np.array([[float(elem) for elem in line.split()] for line in open(File,'r')]))
      max_lkl.append(min(cheese[:,1])) # beware, it is the min because we are talking about '- log likelihood'

    # Selecting only the true maximum.
    max_lkl = min(max_lkl)

    # Restarting the circling through files
    for File in Files:
      i=Files.index(File)
      print 'scanning file {0}'.format(File)
      # cheese will brutally contain everything in the chain File being scanned
      cheese = (np.array([[float(elem) for elem in line.split()] for line in open(File,'r')]))
      line_count = 0
      for line in open(File,'r'):
	line_count+=1
      self.log.write("%s\t Number of steps:%d\tSteps accepted:%d\tacc = %.2g\tmin(-loglike) = %.5g " % (File,sum(cheese[:,0]),line_count,line_count*1.0/sum(cheese[:,0]),max_lkl))
      self.log.write("\n")
      total += sum(cheese[:,0])

      # Removing burn-in
      start = 0
      while cheese[start,1]>max_lkl+2:
	start+=1
      print '  Removed {0} points of burn-in'.format(start)

      # ham contains cheese without the burn-in
      ham = np.copy(cheese[start::])

      # Deal with single file case
      if len(Files) == 1:
	print '  Beware, convergence computed for a single file'
	bacon   = np.copy(cheese[::3,:])
	egg     = np.copy(cheese[1::3,:])
	sausage = np.copy(cheese[2::3,:])

	spam.append(bacon)
	spam.append(egg)
	spam.append(sausage)
	continue
      
      # Adding resulting table to spam
      spam.append(ham)


    # Now that the list spam contains all the different chains removed of their
    # respective burn-in, proceed to the convergence computation
    length  = np.mean([sum(elem[:])[0] for elem in spam])

    # 2D arrays for mean and var, one column will contain the total (over all
    # chains) mean (resp.  variance), and each other column the respective
    # chain mean (resp. chain variance).  R only contains the values for each
    # parameter
    mean     = np.zeros((len(spam)+1,np.shape(spam[0])[1]-2))
    var      = np.zeros((len(spam)+1,np.shape(spam[0])[1]-2))
    R 	     = np.zeros(np.shape(spam[0])[1]-2)
    
    # Quantities for the Gelman Rubin diagnostic
    within  = 0
    between = 0

    # Compute mean and variance for each chain
    print 'Computing mean values'
    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
	for k in range(np.shape(spam[j])[0]):
	  mean[j+1,i] += spam[j][k,0]*spam[j][k,i+2]
	mean[j+1,i] = mean[j+1,i]/sum(spam[j])[0]
	mean[0,i]  += mean[j+1,i]/len(spam)
    
    print 'Computing variance'
    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
	subvar = 0
	for k in range(np.shape(spam[j])[0]):
	  subvar 	   += spam[j][k,0]*(spam[j][k,i+2]-mean[0,i]  )**2
	  var[j+1,i] += spam[j][k,0]*(spam[j][k,i+2]-mean[j+1,i])**2

	var[0,i]   += subvar/(sum(spam[j])[0]-1)/len(spam)
	var[j+1,i]  = var[j+1,i]/(sum(spam[j])[0]-1)

    
    # Gelman Rubin Diagnostic:
    # Computes a quantity linked to the ratio of the mean of the variances of
    # the different chains (within), and the variance of the means (between)
    # TODO: verify the validity of this computation for chains of different
    # length !
    print 'Computing convergence'
    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
	within  += var[j+1,i]/len(spam)
	between += sum(spam[j])[0]*(mean[j+1,i]-mean[0,i])**2 / (len(spam)-1)

      R[i] = math.sqrt(((1-1/length)*within+(len(spam)+1)/(len(spam)*length)*between)/within)
    
    # Log finally the total number of steps
    self.log.write("--> Total number of steps:%d" % total)
    if is_main_chain:
      self.spam = spam

      self.ref_names = ref_names
      self.tex_names = tex_names
      self.boundaries = boundaries

      self.mean = mean
      self.var  = var
      self.R = R

      self.scales = scales

      return True
    else:
      return spam,ref_names,tex_names,boundaries,mean

  # Plotting routine, also computes the sigma errors. Partly imported from
  # Karim Benabed in pmc. However, many options from this method (if not all of
  # them) are unused, meaning : select, legend, show_prop, show_peak,
  # show_extra, add_legend, tick_at_peak, convolve
  def plot_triangle(self,chain,command_line,select=None,bin_number=20,scales=(),legend=(),levels=(68.26,95.4,99.7),show_prop=True,fill=68.26,show_mean=True,show_peak=True,show_extra=None,add_legend=r"$=%(peak).4g^{+%(up).3g}_{-%(down).3g}$",aspect=(16,16),fig=None,tick_at_peak=False,convolve=True,comp_chain = None,comp_ref_names = None,comp_tex_names = None, comp_folder = None,comp_boundaries = None,comp_mean = None):

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

    plt.figtext(0.4,0.95,'Experiments: '+exps,fontsize=40,alpha=0.6)
    plt.figtext(0.9, 0.7,'Monte Python',fontsize=70, rotation=90,alpha=0.15)
    fig1d = plt.figure(2,figsize=aspect)

    # clear figure
    plt.clf()

    #fig2d.savefig(self.folder+'plots/{0}_triangle.pdf'.format(self.folder.split('/')[-2]))
    #exit()
    n = np.shape(chain)[1]-2
    if not scales:
     scales = np.ones(n)
    scales=np.array(scales)

    #################
    # Beginning of unused stuff
    if select==None:
      select = range(n)

    if show_mean:
      if not isinstance(show_mean,(list,tuple)): 
	show_mean=(True,)*(n)
    else:
      if not isinstance(show_mean,(list,tuple)): 
	show_mean=(False,)*(n)

    if show_peak:
      if not isinstance(show_peak,(list,tuple)): 
	show_peak=(True,)*(n)
    else:
      if not isinstance(show_peak,(list,tuple)): 
	show_peak=(False,)*(n)

    n = len(select)
    # End of unused stuff
    ################
    
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
    fig1d.subplots_adjust(bottom=0.03, left=.04, right=0.98, top=0.95, hspace=.35)
    if plot_2d:
      fig2d.subplots_adjust(bottom=0.03, left=.04, right=0.98, top=0.98, hspace=.35)
    

    # In case of a comparison, figure out which names are shared, which are
    # unique and thus require a simple treatment.
    if comp:
      index = 1
      backup_comp_names = np.copy(comp_ref_names)
      for i in range(len(self.ref_names)):
	index +=1
	if self.ref_names[i] in comp_ref_names:
	  comp_ref_names.remove(self.ref_names[i])
	  comp_tex_names.remove(self.tex_names[i])
      for name in comp_ref_names:
	index +=1
      num_columns = round(math.sqrt(index)) 
      num_lines   = math.ceil((len(self.ref_names)+len(comp_ref_names))*1.0/num_columns)
    else:
      num_columns = round(math.sqrt(len(self.ref_names)))
      num_lines   = math.ceil(len(self.ref_names)*1.0/num_columns)

    # Actual plotting
    for i in range(len(self.ref_names)):

      # Adding the subplots to the respective figures, this will be the diagonal for the triangle plot.
      if plot_2d:
	ax2d=fig2d.add_subplot(len(self.ref_names),len(self.ref_names),i*(len(self.ref_names)+1)+1,yticks=[])
      ax1d = fig1d.add_subplot(num_lines,num_columns,i+1,yticks=[])

      # normalized histogram
      hist,bin_edges=np.histogram(chain[:,i+2],bins=bin_number,weights=chain[:,0],normed=False)
      bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])

      # interpolated histogram (if available)
      interp_hist,interp_grid = self.cubic_interpolation(hist,bincenters)
      interp_hist /= np.max(interp_hist)

      if comp:
	try:
	  # For the names in common, the following line will not output an
	  # error. Then compute the comparative histogram
	  ii = np.where( backup_comp_names == self.ref_names[i] )[0][0]
	  comp_hist,comp_bin_edges = np.histogram(comp_chain[:,ii+2],bins=bin_number,weights=comp_chain[:,0],normed=False)
	  comp_bincenters = 0.5*(comp_bin_edges[1:]+comp_bin_edges[:-1])
	  interp_comp_hist,interp_comp_grid = self.cubic_interpolation(comp_hist,comp_bincenters)
	  interp_comp_hist /= np.max(interp_comp_hist)
	  comp_done = True
	except IndexError : # If the name was not found, return the error. This will be then plotted at the end
	  comp_done = False
      if comp:
	if not comp_done:
	  print '{0} was not found in the second folder'.format(self.ref_names[i])

      # minimum credible interval (method by Jan Haman). Fails for multimodal histograms
      bounds = self.minimum_credible_intervals(hist,bincenters,lvls)
      if bounds is False: # print out the faulty histogram (try reducing the binnumber to avoir this)
	print hist
      else:
	for elem in bounds:
	  for j in (0,1):
	    elem[j] -= self.mean[i]
	self.bounds[i] = bounds

      if comp_done:
	comp_bounds = self.minimum_credible_intervals(comp_hist,comp_bincenters,lvls)
	if comp_bounds is False:
	  print comp_hist
	else:
	  for elem in comp_bounds:
	    for j in (0,1):
	      elem[j] -= comp_mean[ii]

      # plotting
      if plot_2d:
	ax2d.set_xticks(ticks[i])
	fontsize2d,ticksize2d = self.get_fontsize( len(self.tex_names) )
	ax2d.set_xticklabels(['%.4g' % s for s in ticks[i]],fontsize=ticksize2d)
	ax2d.set_title('%s= $%.4g^{+%.4g}_{%.4g}$' % (self.tex_names[i],self.mean[i],bounds[0][1],bounds[0][0]),fontsize=fontsize2d)
	ax2d.plot(interp_grid,interp_hist,color='red',linewidth=2,ls='-')
	ax2d.axis([x_range[i][0], x_range[i][1],0,1.05])

      fontsize1d,ticksize1d = self.get_fontsize( max(num_columns,num_lines))
      ax1d.set_title('%s= $%.4g^{+%.4g}_{%.4g}$' % (self.tex_names[i],self.mean[i],bounds[0][1],bounds[0][0]),fontsize=fontsize1d)
      ax1d.set_xticks(ticks[i])
      ax1d.set_xticklabels(['%.4g' % s for s in ticks[i]],fontsize=ticksize1d)
      ax1d.axis([x_range[i][0], x_range[i][1],0,1.05])
      
      if comp_done:
	# complex variation of intervals
	if comp_x_range[ii][0] > x_range[i][0]:
	  comp_ticks[ii][0] = ticks[i][0]
	  comp_x_range[ii][0] = x_range[i][0]
	if comp_x_range[ii][1] < x_range[i][1]:
	  comp_ticks[ii][2] = ticks[i][2]
	  comp_x_range[ii][1] = x_range[i][1]
	comp_ticks[ii][1] = (comp_x_range[ii][1]+comp_x_range[ii][0])/2.
	ax1d.set_xticks(comp_ticks[ii])
	ax1d.set_xticklabels(['%.4g' % s for s in comp_ticks[ii]],fontsize=ticksize1d)
	ax1d.axis([comp_x_range[ii][0], comp_x_range[ii][1],0,1.05])


      ax1d.plot(interp_grid,interp_hist,color='black',linewidth=2,ls='-')
      if comp_done:
	ax1d.plot(interp_comp_grid,interp_comp_hist,color='red',linewidth=2,ls='-')

      # mean likelihood (optional, if comparison, it will not be printed)
      if plot_2d:
	lkl_mean=np.zeros(len(bincenters),'float64')
	norm=np.zeros(len(bincenters),'float64')
	for j in range(len(bin_edges)-1):
	  for k in range(np.shape(chain)[0]):
	    if (chain[k,i+2]>=bin_edges[j] and chain[k,i+2]<=bin_edges[j+1]):
	      lkl_mean[j] += math.exp( best_minus_lkl - chain[k,1])*chain[k,0]
	      norm[j] += chain[k,0]
	  lkl_mean[j] /= norm[j]
	#lkl_mean *= max(hist)/max(lkl_mean)
	lkl_mean /= max(lkl_mean)
	interp_lkl_mean,interp_grid = self.cubic_interpolation(lkl_mean,bincenters)
	ax2d.plot(interp_grid,interp_lkl_mean,color='red',ls='--',lw=2)
	ax1d.plot(interp_grid,interp_lkl_mean,color='red',ls='--',lw=4)

      if command_line.subplot is True:
	if not comp_done:
	  extent2d = ax2d.get_window_extent().transformed(fig2d.dpi_scale_trans.inverted())
	  fig2d.savefig(self.folder+'plots/{0}_{1}.pdf'.format(self.folder.split('/')[-2],self.ref_names[i]), bbox_inches=extent2d.expanded(1.1, 1.4))
	else:
	  extent1d = ax1d.get_window_extent().transformed(fig1d.dpi_scale_trans.inverted())
	  fig1d.savefig(self.folder+'plots/{0}_{1}.pdf'.format(self.folder.split('/')[-2],self.ref_names[i]), bbox_inches=extent1d.expanded(1.1, 1.4))

      # Now do the rest of the triangle plot
      if plot_2d:
	for j in range(i):
	  ax2dsub=fig2d.add_subplot(len(self.ref_names),len(self.ref_names),(i)*len(self.ref_names)+j+1)
	  n,xedges,yedges=np.histogram2d(chain[:,i+2],chain[:,j+2],weights=chain[:,0],bins=(bin_number,bin_number),normed=False)
	  extent = [x_range[j][0], x_range[j][1], x_range[i][0], x_range[i][1]]
	  x_centers = 0.5*(xedges[1:]+xedges[:-1])
	  y_centers = 0.5*(yedges[1:]+yedges[:-1])

	  ax2dsub.set_xticks(ticks[j])
	  if i == len(self.ref_names)-1:
	    ax2dsub.set_xticklabels(['%.4g' % s for s in ticks[j]],fontsize=ticksize2d)
	  else:
	    ax2dsub.set_xticklabels([''])

	  ax2dsub.set_yticks(ticks[i])
	  if j == 0:
	    ax2dsub.set_yticklabels(['%.4g' % s for s in ticks[i]],fontsize=ticksize2d)
	  else:
	    ax2dsub.set_yticklabels([''])
	  ax2dsub.imshow(n, extent=extent, aspect='auto',interpolation='gaussian',origin='lower',cmap=matplotlib.cm.Reds)

	  # plotting contours, using the ctr_level method (from Karim Benabed)
	  cs = ax2dsub.contour(y_centers,x_centers,n,extent=extent,levels=self.ctr_level(n,lvls),colors="k",zorder=5)
	  ax2dsub.clabel(cs, cs.levels[:-1], inline=True,inline_spacing=0, fmt=dict(zip(cs.levels[:-1],[r"%d \%%"%int(l*100) for l in lvls[::-1]])), fontsize=19)

	  if command_line.subplot is True:
	    # Store the individual 2d plots
	    fig_temp=plt.figure(3,figsize=(6,6))
	    fig_temp.clf()
	    ax_temp=fig_temp.add_subplot(111)
	    ax_temp.set_xticks(ticks[j])
	    ax_temp.set_yticks(ticks[i])
	    ax_temp.imshow(n, extent=extent, aspect='auto',interpolation='gaussian',origin='lower',cmap=matplotlib.cm.Reds)
	    ax_temp.set_xticklabels(['%.4g' % s for s in ticks[j]],fontsize=ticksize2d)
	    ax_temp.set_yticklabels(['%.4g' % s for s in ticks[i]],fontsize=ticksize2d)
	    ax_temp.set_title('%s vs %s' % (self.tex_names[i],self.tex_names[j]),fontsize=fontsize1d)
	    cs = ax_temp.contour(y_centers,x_centers,n,extent=extent,levels=self.ctr_level(n,lvls),colors="k",zorder=5)
	    ax_temp.clabel(cs, cs.levels[:-1], inline=True,inline_spacing=0, fmt=dict(zip(cs.levels[:-1],[r"%d \%%"%int(l*100) for l in lvls[::-1]])), fontsize=19)
	    fig_temp.savefig(self.folder+'plots/{0}_2d_{1}-{2}.pdf'.format(self.folder.split('/')[-2],self.ref_names[i],self.ref_names[j]))

	    # store the coordinates of the points for further plotting.
	    plot_file = open(self.folder+'plots/{0}_2d_{1}-{2}.dat'.format(self.folder.split('/')[-2],self.ref_names[i],self.ref_names[j]),'w')
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
      for i in range(len(self.ref_names),len(self.ref_names)+len(comp_ref_names)):

	ax1d = fig1d.add_subplot(num_lines,num_columns,i+1,yticks=[])
	ii = np.where(backup_comp_names == comp_ref_names[i-len(self.ref_names)])[0][0]

	comp_hist,comp_bin_edges = np.histogram(comp_chain[:,ii+2],bins=bin_number,weights=comp_chain[:,0],normed=False)
	comp_bincenters = 0.5*(comp_bin_edges[1:]+comp_bin_edges[:-1])
	interp_comp_hist,interp_comp_grid = self.cubic_interpolation(comp_hist,comp_bincenters)
	interp_comp_hist /= np.max(interp_comp_hist)

	comp_bounds = self.minimum_credible_intervals(comp_hist,comp_bincenters,lvls)
	if comp_bounds is False:
	  print comp_hist
	else:
	  for elem in comp_bounds:
	    for j in (0,1):
	      elem[j] -= comp_mean[ii]
	ax1d.set_xticks(comp_ticks[ii])
	ax1d.set_xticklabels(['%.4g' % s for s in comp_ticks[ii]],fontsize=ticksize1d)
	ax1d.axis([comp_x_range[ii][0], comp_x_range[ii][1],0,1.05])
	ax1d.set_title('%s= $%.4g^{+%.4g}_{%.4g}$' % (comp_tex_names[i-len(self.ref_names)],comp_mean[ii],comp_bounds[0][1],comp_bounds[0][0]),fontsize=fontsize1d)
	ax1d.plot(interp_comp_grid,interp_comp_hist,color='red',linewidth=2,ls='-')
	if command_line.subplot is True:
	  extent1d = ax1d.get_window_extent().transformed(fig1d.dpi_scale_trans.inverted())
	  fig1d.savefig(self.folder+'plots/{0}_{1}.pdf'.format(self.folder.split('/')[-2],self.ref_names[i]), bbox_inches=extent1d.expanded(1.1, 1.4))
	
    # If plots/ folder in output folder does not exist, create it
    if os.path.isdir(self.folder+'plots') is False:
      os.mkdir(self.folder+'plots')
    if plot_2d:
      fig2d.savefig(self.folder+'plots/{0}_triangle.pdf'.format(self.folder.split('/')[-2]))
    if comp:
      fig1d.savefig(self.folder+'plots/{0}-vs-{1}.pdf'.format(self.folder.split('/')[-2],comp_folder.split('/')[-2]))
    else:
      fig1d.savefig(self.folder+'plots/{0}_1d.pdf'.format(self.folder.split('/')[-2]))


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
	  print '  Warning : Can not derive minimum credible intervals for this multimodal posterior'
	  print histogram
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
	  print '\n\n  the loop to check for sigma deviations was too long to converge'
	  break

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
	
    return bounds

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
    # for a diagonal of 5, fontsize of 19, for a diagonal of 13, fontsize of 8
    fontsize = 19
#round( 19 - (diag_length-5)*1.38)
    ticksize = 19
#round( 14 - (diag_length-5)*1)
    return fontsize,ticksize
