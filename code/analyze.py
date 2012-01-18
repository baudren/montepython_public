import os,sys
import re
import io
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib
import math
import numpy as np
from scipy.interpolate import interp1d


class info:

  def __init__(self,command_line):
  
    Files     = command_line.files
    binnumber = command_line.bins
    self.prepare(Files)
    self.convergence()

    chain = np.copy(self.spam[0])
    for i in range(len(self.spam)-1):
     chain = np.append(chain,self.spam[i+1],axis=0)

    if command_line.comp is not None:
      comp_Files,comp_folder,comp_log = self.prepare(command_line.comp,is_main_chain = False)
      comp_spam,comp_ref_names,comp_tex_names,comp_boundaries,comp_mean = self.convergence(is_main_chain = False,Files = comp_Files,log = comp_log)
      comp_mean = comp_mean[0]
      comp_chain = np.copy(comp_spam[0])
      for i in range(len(comp_spam)-1):
	comp_chain = np.append(comp_chain,comp_spam[i+1],axis=0)


    weight = sum(chain)[0]

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

    self.cov.write('#{0}\n'.format(self.ref_names))
    for i in range(len(self.ref_names)):
      for j in range(len(self.ref_names)):
	self.cov.write('{0} '.format(self.covar[i][j]))
      self.cov.write('\n')

    a=chain[:,1].argsort(0)
    total=chain[:,0].sum()
    self.lvls = [total*0.68,total*0.95,total*0.99]

    self.bounds = np.zeros((len(self.ref_names),len(self.lvls),2))
    if command_line.comp is None:
      self.plot_triangle(chain,command_line,bin_number=binnumber)
    else:
      self.plot_triangle(chain,command_line,bin_number=binnumber,comp_chain=comp_chain,comp_ref_names = comp_ref_names,comp_tex_names = comp_tex_names,comp_folder = comp_folder,comp_boundaries = comp_boundaries)

    self.info.write('\n param names:\t')
    for elem in self.ref_names:
      self.info.write('%s\t' % (elem))
    self.info.write('\n R-1 values:\t')
    for elem in self.R:
      self.info.write('%.7f\t' % (elem-1.))
    self.info.write('\n Best Fit:\t')
    for elem in chain[a[0],2:]:
      self.info.write('%.7f\t' % (elem))
    self.info.write('\n mean    :\t')
    for elem in self.mean:
      self.info.write('%.7f\t' % (elem))
    self.info.write('\n 1-sigma - :\t')
    for elem in self.bounds:
      self.info.write('%.7f\t' % (elem[0][0]))
    self.info.write('\n 1-sigma + :\t')
    for elem in self.bounds:
      self.info.write('%.7f\t' % (elem[0][1]))
    #self.info.write('\n 2-sigma - :\t')
    #for elem in self.bounds:
      #self.info.write('%.7f\t' % (elem[1][0]))
    #self.info.write('\n 2-sigma + :\t')
    #for elem in self.bounds:
      #self.info.write('%.7f\t' % (elem[1][1]))

  def prepare(self,Files,is_main_chain=True):

    # Scan the whole input folder, and include all chains in it (all files with a
    # '.' in their name will be discarded as not mcmc chains (If you want to have
    # additionnal files than the log.dat and .param in this folder, please keep
    # this in  mind

    if os.path.isdir(Files[0]): 
      if Files[0][-1]!='/':
	Files[0]+='/'
      folder = Files[0]
      Files = [folder+elem for elem in os.listdir(folder) if elem.find('.')==-1]
      for elem in Files:
	if os.path.isdir('{0}'.format(elem)) is True:
	  Files.remove(elem)

    # Recover the folder, depending on the case
    else:
      if (len(Files[0].split('/'))==0 or (Files[0].split('/')[0]=='.')):
	folder = './'
      else:
	folder = ''
	for i in range(len(Files[0].split('/')[:-1])):
	  folder += Files[0].split('/')[i]+'/'

    # Check if the log.dat file exists
    if os.path.isfile(folder+'log.param') is True:
      if os.path.getsize(folder+'log.param')>0:
	log = open(folder+'log.param','r')
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
    else:
      infoname = folder+folder.split('/')[-2]+'.info'
      covname  = folder+folder.split('/')[-2]+'.covmat'

    if is_main_chain:
      self.info  = open(infoname,'w')
      self.cov   = open(covname,'w')
      self.log   = log

      self.Files = Files
      self.folder= folder
      return True
    else:
      return Files,folder,log

  def convergence(self,is_main_chain=True,Files=None,log=None):

    # We here have a list of files, that may be of length 1. If this
    # happens, then we split the only file in 3 subchains, otherwise we compute
    # normally the R coefficient. 
    spam = list()

    # Recovering default ordering of parameters 
    ref_names = []
    tex_names = []
    boundaries= []

    if is_main_chain:
      Files = self.Files
      log   = self.log

    # Defining a list of names that should have a \ in their tex names
    tex_greek = ['omega','tau','alpha','beta','delta','nu','Omega']
    for line in log:
      if line.find('#')==-1:
	if (line.find('data.Class_params')!=-1 or line.find('data.nuisance_params')!=-1):
	  name = line.split("'")[1]
	  if len(line.split('=')[-1].split(',')) == 4:
	    if line.split('=')[-1].split(',')[-1].replace(']\n','').replace(' ','') != '0':
	      temp = [float(elem) for elem in line.split(",")[1:3]]
	      boundaries.append(temp)
	      ref_names.append(name)
	      for elem in tex_greek:
		if elem in name:
		  name="""\\"""+name
	      if name.find('_')!=-1:
		temp_name = name.split('_')[0]+'_{'
		for i in range(1,len(name.split('_'))):
		  temp_name += name.split('_')[i]
		temp_name += '}'
		name = temp_name
	      tex_names.append('${0}$'.format(name))
	  elif len(line.split('=')[-1].split(',')) == 5:
	    if line.split('=')[-1].split(',')[-2].replace(' ','') != 0:
	      temp = [float(elem) for elem in line.split(",")[1:3]]
	      boundaries.append(temp)
	      ref_names.append(name)
	      number = 1./float(line.split('=')[-1].split(',')[-1].replace(']\n','').replace(' ',''))
	      if number < 1000:
		tex_names.append("$%0.d~%s$" % (number,name))
	      else:
		tex_names.append("$%0.e%s$" % (number,name))
		m = re.search(r'(?:\$[0-9]*e\+[0]*)([0-9]*)(.*)',tex_names[-1])
		tex_names[-1] = '$10^{'+m.groups()[0]+'}'+m.groups()[1]
    log.seek(0)

    for File in Files:
      i=Files.index(File)
      print 'scanning file {0}'.format(File)
      cheese = (np.array([[float(elem) for elem in line.split()] for line in open(File,'r')]))
      max_lkl = min(cheese[:,1]) # beware, it is the min because we are talking about '- log likelihood'

      start = 0
      while cheese[start,1]>max_lkl+2:
	start+=1
      print '  Removed {0} points of burn-in'.format(start)
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


    # Convergence computation
    length  = np.mean([sum(elem[:])[0] for elem in spam])

    # 2 D arrays for mean and var, one column will contain the total (over all
    # chains) mean (resp.  variance), and each other column the respective chain
    # mean (resp. chain variance).  R only contains the values for each parameter
    mean     = np.zeros((len(spam)+1,np.shape(spam[0])[1]-2))
    var      = np.zeros((len(spam)+1,np.shape(spam[0])[1]-2))
    R 	     = np.zeros(np.shape(spam[0])[1]-2)
    within  = 0
    between = 0

    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
	for k in range(np.shape(spam[j])[0]):
	  mean[j+1,i] += spam[j][k,0]*spam[j][k,i+2]
	mean[j+1,i] = mean[j+1,i]/sum(spam[j])[0]
	mean[0,i]  += mean[j+1,i]/len(spam)
    
    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
	subvar = 0
	for k in range(np.shape(spam[j])[0]):
	  subvar 	   += spam[j][k,0]*(spam[j][k,i+2]-mean[0,i]  )**2
	  var[j+1,i] += spam[j][k,0]*(spam[j][k,i+2]-mean[j+1,i])**2

	var[0,i]   += subvar/(sum(spam[j])[0]-1)/len(spam)
	var[j+1,i]  = var[j+1,i]/(sum(spam[j])[0]-1)

    
    # Gelman Rubin Diagnostic
    # verify the validity of this computation for different length !
    for i in range(np.shape(mean)[1]):
      for j in range(len(spam)):
	within  += var[j+1,i]/len(spam)
	between += sum(spam[j])[0]*(mean[j+1,i]-mean[0,i])**2 / (len(spam)-1)

      R[i] = math.sqrt(((1-1/length)*within+(len(spam)+1)/(len(spam)*length)*between)/within)
    
    if is_main_chain:
      self.spam = spam

      self.ref_names = ref_names
      self.tex_names = tex_names
      self.boundaries = boundaries

      self.mean = mean
      self.var  = var
      self.R = R

      return True
    else:
      return spam,ref_names,tex_names,boundaries,mean

  def plot_triangle(self,chain,command_line,select=None,bin_number=20,scales=(),legend=(),levels=(68.26,95.4,99.7),show_prop=True,fill=68.26,show_mean=True,show_peak=True,show_extra=None,add_legend=r"$=%(peak).4g^{+%(up).3g}_{-%(down).3g}$",aspect=(16,16),fig=None,tick_at_peak=False,convolve=True,comp_chain = None,comp_ref_names = None,comp_tex_names = None, comp_folder = None,comp_boundaries = None):

    # If comparison is asked, don't plot 2d levels
    if command_line.comp is not None:
      plot_2d = False
      comp    = True
    else:
      plot_2d = True
      comp    = False

    matplotlib.rc('text',usetex = True)
    matplotlib.rc('font',size=11)
    matplotlib.rc('xtick',labelsize='8')
    matplotlib.rc('ytick',labelsize='8')
    lvls = np.array(levels)/100.

    if plot_2d:
      if fig:
	fig2d = plt.figure(fig,aspect)
      else:
	fig2d = plt.figure(1,figsize=aspect)

    # TEST
    fig1d = plt.figure(2,figsize=aspect)
    # clear figure
    plt.clf()

    n = np.shape(chain)[1]-2
    # to be modified by data.param.itervalues()[4]
    if select==None:
      select = range(n)
    if not scales:
     scales = np.ones(n)
    scales=np.array(scales)


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
    mean = self.mean*scales
    var  = self.var*scales**2
    #pmax = 
    # 1D plot
    max_values = np.max(chain[:,2:],axis=0)*scales
    min_values = np.min(chain[:,2:],axis=0)*scales
    span = (max_values-min_values)

    best_minus_lkl = np.min(chain[:,1],axis=0)

    if comp:
      comp_max_values = np.max(comp_chain[:,2:],axis=0)
      comp_min_values = np.min(comp_chain[:,2:],axis=0)
      comp_span = (comp_max_values-comp_min_values)

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
      
    if plot_2d:
      fig2d.subplots_adjust(bottom=0.03, left=.04, right=0.98, top=0.98, hspace=.35)

    if comp:
      index = 0
      backup_comp_names = np.copy(comp_ref_names)
      for name in self.ref_names:
	index +=1
	if name in comp_ref_names:
	  comp_ref_names.remove(name)
      for name in comp_ref_names:
	index +=1
      num_columns = round(math.sqrt(index)) 
      num_lines   = round((len(self.ref_names)+len(comp_ref_names))*1.0/num_columns)
    else:
      num_columns = round(math.sqrt(len(self.ref_names)))
      num_lines   = round(len(self.ref_names)*1.0/num_columns)

    for i in range(len(self.ref_names)):

      if plot_2d:
	ax2d=fig2d.add_subplot(len(self.ref_names),len(self.ref_names),i*(len(self.ref_names)+1)+1,yticks=[])

      ax1d = fig1d.add_subplot(num_columns,num_lines,i,yticks=[])

      # histogram
      hist,bin_edges=np.histogram(chain[:,i+2],bins=bin_number,weights=chain[:,0],normed=False)
      bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])

      interp_hist,interp_grid = self.cubic_interpolation(hist,bincenters)

      if comp:
	try:
	  #ii = backup_comp_names.index(self.ref_names[i])
	  ii = np.where( backup_comp_names == self.ref_names[i] )[0][0]
	  comp_hist,comp_bin_edges = np.histogram(comp_chain[:,ii+2],bins=bin_number,weights=comp_chain[:,0],normed=False)
	  comp_hist *= max(hist)/max(comp_hist)
	  comp_bincenters = 0.5*(comp_bin_edges[1:]+comp_bin_edges[:-1])
	  comp_done = True
	except ValueError :
	  print 'ooh'
	  comp_done = False
      if comp:
	if not comp_done:
	  print '{0} was not found in the comparison folder'.format(self.ref_names[i])

      # minimum credible interval
      bounds = self.minimum_credible_intervals(hist,bincenters,lvls)
      if bounds is False:
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
	      elem[j] -= self.mean[i]

      # plotting
      if plot_2d:
	ax2d.set_xticks(ticks[i])
	ax2d.set_xticklabels(['%.4g' % s for s in ticks[i]])
	ax2d.set_title('%s= %.4g' % (self.tex_names[i],mean[i]))
	ax2d.plot(bincenters,hist,color='red',linewidth=2,ls='-')
	ax2d.axis([x_range[i][0], x_range[i][1],0,np.max(hist)])

      ax1d.set_title('%s= %.4g' % (self.tex_names[i],mean[i]))
      ax1d.set_xticks(ticks[i])
      ax1d.set_xticklabels(['%.4g' % s for s in ticks[i]])
      ax1d.axis([x_range[i][0], x_range[i][1],0,np.max(hist)])
      
      if comp_done:
	if comp_span[ii] >= span[i]:
	  ax1d.set_xticks(comp_ticks[ii])
	  ax1d.set_xticklabels(['%.4g' % s for s in comp_ticks[ii]])
	  ax1d.axis([comp_x_range[i][0], comp_x_range[i][1],0,np.max(comp_hist)])

      ax1d.plot(interp_grid,interp_hist,color='black',linewidth=2,ls='-')
      if comp_done:
	ax1d.plot(comp_bincenters,comp_hist,color='red',linewidth=2,ls='-')


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
	lkl_mean *= max(hist)/max(lkl_mean)
	interp_lkl_mean,interp_grid = self.cubic_interpolation(lkl_mean,bincenter)
	ax2d.plot(interp_grid,interp_lkl_mean,color='red',ls='--',lw=2)
	ax1d.plot(interp_grid,interp_lkl_mean,color='red',ls='--',lw=4)

      if plot_2d:
	for j in range(i):
	  ax2dsub=fig2d.add_subplot(len(self.ref_names),len(self.ref_names),(i)*len(self.ref_names)+j+1)
	  n,xedges,yedges=np.histogram2d(chain[:,i+2],chain[:,j+2],weights=chain[:,0],bins=(bin_number,bin_number),normed=False)
	  extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
	  x_centers = 0.5*(xedges[1:]+xedges[:-1])
	  y_centers = 0.5*(yedges[1:]+yedges[:-1])

	  ax2dsub.set_xticks(ticks[j])
	  if i == len(self.ref_names)-1:
	    ax2dsub.set_xticklabels(['%.4g' % s for s in ticks[j]])
	  else:
	    ax2dsub.set_xticklabels([''])

	  ax2dsub.set_yticks(ticks[i])
	  if j == 0:
	    ax2dsub.set_yticklabels(['%.4g' % s for s in ticks[i]])
	  else:
	    ax2dsub.set_yticklabels([''])
	  ax2dsub.imshow(n.T, extent=extent, aspect='auto',interpolation='gaussian',origin='lower',cmap=matplotlib.cm.Reds)

	  # plotting contours
	  cs = ax2dsub.contour(y_centers,x_centers,n.T,extent=extent,levels=self.ctr_level(n.T,lvls),colors="k",zorder=5)
	  ax2dsub.clabel(cs, cs.levels[:-1], inline=True,inline_spacing=0, fmt=dict(zip(cs.levels[:-1],[r"%d \%%"%int(l*100) for l in lvls[::-1]])), fontsize=6)

    if comp:
      for i in range(len(self.ref_names),len(self.ref_names)+len(comp_ref_names)):

	ax1d = fig1d.add_subplot(num_columns,num_lines,i+len(self.ref_names),yticks=[])

	ii = backup_comp_names.index(comp_ref_names[i-len(self.ref_names)])
	comp_hist,comp_bin_edges = np.histogram(comp_chain[:,ii+2],bins=bin_number,weights=comp_chain[:,0],normed=False)
	comp_bincenters = 0.5*(comp_bin_edges[1:]+comp_bin_edges[:-1])

	comp_bounds = self.minimum_credible_intervals(comp_hist,comp_bincenters,lvls)
	if comp_bounds is False:
	  print comp_hist
	else:
	  for elem in comp_bounds:
	    for j in (0,1):
	      elem[j] -= self.mean[i-len(self.ref_names)]
	ax1d.plot(comp_bincenters,comp_hist,color='green',linewidth=2,ls='-')
	
    # If plots/ folder in output folder does not exist, create it
    if os.path.isdir(self.folder+'plots') is False:
      os.mkdir(self.folder+'plots')
    if plot_2d:
      fig2d.savefig(self.folder+'plots/{0}_triangle.pdf'.format(self.folder.split('/')[-2]))
    if comp:
      fig1d.savefig(self.folder+'plots/{0}-vs-{1}.pdf'.format(self.folder.split('/')[-2],comp_folder.split('/')[-2]))
    else:
      fig1d.savefig(self.folder+'plots/{0}_1d.pdf'.format(self.folder.split('/')[-2]))






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

  def minimum_credible_intervals(self,histogram,bincenter,levels):
    bounds = np.zeros((len(levels),2))
    j = 0
    delta = bincenter[1]-bincenter[0]
    for level in levels:
      norm = float((sum(histogram)-0.5*(histogram[0]+histogram[-1]))*delta)
      water_level_up   = max(histogram)*1.0
      water_level_down = min(histogram)*1.0
      top = 0.
      
      ii=0
      while (abs((top/norm)-level) > 0.0001):
	top=0.
	water_level = (water_level_up + water_level_down)/2.
	ontop = [elem for elem in histogram if elem>water_level]
	indices = [i for i in range(len(histogram)) if histogram[i]>water_level]
	# check for multimodal posteriors
	if ((indices[-1]-indices[0]+1)!=len(indices)):
	  print '\n\n  Can not derive minimum credible intervals for this multimodal posterior'
	  print histogram
	  break
	top = (sum(histogram[indices])-0.5*(histogram[indices[0]]+histogram[indices[-1]]))*(delta)

	# left
	if indices[0]>0:
	  top += 0.5 * (water_level + histogram[indices[0]]) * (delta)*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-histogram[indices[0]-1])

	# right
	if indices[-1]<(len(histogram)-1) :
	  top += 0.5 * (water_level + histogram[indices[-1]]) * (delta)*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-histogram[indices[-1]+1])

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
	bounds[j][0] = bincenter[indices[0]] - delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-histogram[indices[0]-1])
      else:
	bounds[j][0] = bincenter[0]
      # max
      if indices[-1]<(len(histogram)-1):
	bounds[j][1] = bincenter[indices[-1]]+ delta*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-histogram[indices[-1]+1])
      else:
	bounds[j][1] = bincenter[-1]
      j+=1
	
    return bounds

  def cubic_interpolation(self,hist,bincenter):
    interp_grid = np.linspace(bincenter[0],bincenter[-1],len(bincenter)*10)
    interp_hist = interp1d(bincenter,hist,interp_grid)
    return interp_hist,interp_grid

