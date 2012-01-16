import os,sys
import re
import io
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib
import math
import numpy as np


class info:

  def __init__(self,Files,binnumber):
  
    self.prepare(Files)
    self.convergence()

    chain = np.copy(self.spam[0])
    for i in range(len(self.spam)-1):
     chain = np.append(chain,self.spam[i+1],axis=0)

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

    # sorting
    a=chain[:,1].argsort(0)
    total=chain[:,0].sum()
    Sum=0
    for i in range(len(a)):
      Sum+=chain[a[i],0]
      if Sum<total*0.99:
	if Sum<total*0.95:
	  if Sum<total*0.68:
	    a68=i
	  else:
	    a95= i
	else:
	  a99=i

    self.lvls = [total*0.68,total*0.95,total*0.99]
    self.plot_triangle(chain,bin_number=binnumber)

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
    self.info.write('\n 1-sigma :\t')
    for elem in self.var:
      self.info.write('%.7f\t' % (math.sqrt(elem) ))

  def prepare(self,Files):

    # Scan the whole input folder, and include all chains in it (all files with a
    # '.' in their name will be discarded as not mcmc chains (If you want to have
    # additionnal files than the log.dat and .param in this folder, please keep
    # this in  mind

    if os.path.isdir(Files[0]): 
      if Files[0][-1]!='/':
	Files[0]+='/'
      self.folder = Files[0]
      Files = [Files[0]+elem for elem in os.listdir(Files[0]) if elem.find('.')==-1]

    # Recover the folder, depending on the case
    else:
      if (len(Files[0].split('/'))==0 or (Files[0].split('/')[0]=='.')):
	self.folder = './'
      else:
	self.folder = ''
	for i in range(len(Files[0].split('/')[:-1])):
	  self.folder += Files[0].split('/')[i]+'/'

    # Check if the log.dat file exists
    if os.path.isfile(self.folder+'log.param') is True:
      if os.path.getsize(self.folder+'log.param')>0:
	self.log = open(self.folder+'log.param','r')
      else:
	print '\n\n  The log param file {0} seems empty'.format(self.folder+'log.dat')
	exit()
    else:
      print '\n\n  The log param file {0} is absent ?'.format(self.folder+'log.dat')
      exit()

    # If the folder has no subdirectory, then go for a simple infoname,
    # otherwise, call it with the last name
    if (len(self.folder.split('/')) <= 2 and self.folder.split('/')[-1] == ''):
      infoname = self.folder+self.folder.rstrip('/')+'.info'
      covname  = self.folder+self.folder.rstrip('/')+'.covmat'
    else:
      infoname = self.folder+self.folder.split('/')[-2]+'.info'
      covname  = self.folder+self.folder.split('/')[-2]+'.covmat'

    self.info  = open(infoname,'w')
    self.cov   = open(covname,'w')

    # Updating the files names, taking into account the removed ones from
    # cleaning
    if len(Files)>=2:
      Files = [ self.folder+elem for elem in os.listdir(self.folder) if elem.find('.')==-1 ]

    self.Files = Files
    return True

  def convergence(self):

    # We here have a list of files, that may be of length 1. If this
    # happens, then we split the only file in 3 subchains, otherwise we compute
    # normally the R coefficient. 
    self.spam = list()

    # Recovering default ordering of parameters 
    self.ref_names = []
    self.tex_names = []
    self.boundaries= []

    # Defining a list of names that should have a \ in their tex names
    tex_greek = ['omega','tau','alpha','beta','delta','nu','Omega']
    for line in self.log:
      if line.find('#')==-1:
	if (line.find('data.Class_params')!=-1 or line.find('data.nuisance_params')!=-1):
	  name = line.split("'")[1]
	  if len(line.split('=')[-1].split(',')) == 4:
	    if line.split('=')[-1].split(',')[-1].replace(']\n','').replace(' ','') != '0':
	      temp = [float(elem) for elem in line.split(",")[1:3]]
	      self.boundaries.append(temp)
	      self.ref_names.append(name)
	      for elem in tex_greek:
		if elem in name:
		  name="""\\"""+name
	      if name.find('_')!=-1:
		temp_name = name.split('_')[0]+'_{'
		for i in range(1,len(name.split('_'))):
		  temp_name += name.split('_')[i]
		temp_name += '}'
		name = temp_name
	      self.tex_names.append('${0}$'.format(name))
	  elif len(line.split('=')[-1].split(',')) == 5:
	    if line.split('=')[-1].split(',')[-2].replace(' ','') != 0:
	      temp = [float(elem) for elem in line.split(",")[1:3]]
	      self.boundaries.append(temp)
	      self.ref_names.append(name)
	      number = 1./float(line.split('=')[-1].split(',')[-1].replace(']\n','').replace(' ',''))
	      if number < 1000:
		self.tex_names.append("$%0.d~%s$" % (number,name))
	      else:
		self.tex_names.append("$%0.e%s$" % (number,name))
		m = re.search(r'(?:\$[0-9]*e\+[0]*)([0-9]*)(.*)',self.tex_names[-1])
		self.tex_names[-1] = '$10^{'+m.groups()[0]+'}'+m.groups()[1]
    self.log.seek(0)

    for File in self.Files:
      i=self.Files.index(File)
      print 'scanning file {0}'.format(File)
      cheese = (np.array([[float(elem) for elem in line.split()] for line in open(File,'r')]))
      max_lkl = min(cheese[:,1]) # beware, it is the min because we are talking about '- log likelihood'

      start = 0
      while cheese[start,1]>max_lkl+2:
	start+=1
      print '  Removed {0} points of burn-in'.format(start)
      ham = np.copy(cheese[start::])

      # Deal with single file case
      if len(self.Files) == 1:
	print '  Beware, convergence computed for a single file'
	bacon   = np.copy(cheese[::3,:])
	egg     = np.copy(cheese[1::3,:])
	sausage = np.copy(cheese[2::3,:])

	self.spam.append(bacon)
	self.spam.append(egg)
	self.spam.append(sausage)
	continue
      
      # Adding resulting table to spam
      self.spam.append(ham)


    # Convergence computation
    length  = np.mean([sum(elem[:])[0] for elem in self.spam])

    # 2 D arrays for mean and var, one column will contain the total (over all
    # chains) mean (resp.  variance), and each other column the respective chain
    # mean (resp. chain variance).  R only contains the values for each parameter
    self.mean     = np.zeros((len(self.spam)+1,np.shape(self.spam[0])[1]-2))
    self.var      = np.zeros((len(self.spam)+1,np.shape(self.spam[0])[1]-2))
    self.R 	  = np.zeros(np.shape(self.spam[0])[1]-2)
    within  = 0
    between = 0

    for i in range(np.shape(self.mean)[1]):
      for j in range(len(self.spam)):
	for k in range(np.shape(self.spam[j])[0]):
	  self.mean[j+1,i] += self.spam[j][k,0]*self.spam[j][k,i+2]
	self.mean[j+1,i] = self.mean[j+1,i]/sum(self.spam[j])[0]
	self.mean[0,i]  += self.mean[j+1,i]/len(self.spam)
    
    for i in range(np.shape(self.mean)[1]):
      for j in range(len(self.spam)):
	subvar = 0
	for k in range(np.shape(self.spam[j])[0]):
	  subvar 	   += self.spam[j][k,0]*(self.spam[j][k,i+2]-self.mean[0,i]  )**2
	  self.var[j+1,i] += self.spam[j][k,0]*(self.spam[j][k,i+2]-self.mean[j+1,i])**2

	self.var[0,i]   += subvar/(sum(self.spam[j])[0]-1)/len(self.spam)
	self.var[j+1,i]  = self.var[j+1,i]/(sum(self.spam[j])[0]-1)

    
    # Gelman Rubin Diagnostic
    # verify the validity of this computation for different length !
    for i in range(np.shape(self.mean)[1]):
      for j in range(len(self.spam)):
	within  += self.var[j+1,i]/len(self.spam)
	between += sum(self.spam[j])[0]*(self.mean[j+1,i]-self.mean[0,i])**2 / (len(self.spam)-1)

      self.R[i] = math.sqrt(((1-1/length)*within+(len(self.spam)+1)/(len(self.spam)*length)*between)/within)
    return True

  def plot_triangle(self,chain,select=None,bin_number=20,scales=(),legend=(),levels=(68.26,95.4,99.7),show_prop=True,fill=68.26,show_mean=True,show_peak=True,show_extra=None,add_legend=r"$=%(peak).4g^{+%(up).3g}_{-%(down).3g}$",aspect=(16,16),fig=None,tick_at_peak=False,convolve=True):

    matplotlib.rc('text',usetex = True)
    matplotlib.rc('font',size=11)
    matplotlib.rc('xtick',labelsize='8')
    matplotlib.rc('ytick',labelsize='8')
    lvls = np.array(levels)/100.

    if fig:
      fig = plt.figure(fig,aspect)
    else:
      fig = plt.figure(1,figsize=aspect)

    # TEST
    fig2 = plt.figure(2,figsize=aspect)
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

    if tick_at_peak:
      pass
    else:
      ticks = np.array((min_values+span*0.1,(max_values+min_values)/2.,max_values-span*0.1)).T
      x_range = np.array((min_values,max_values)).T
      #ticks = np.array((min_values,(max_values+min_values)/2.,max_values)).T

    for i in range(np.shape(ticks)[0]):
      if abs(x_range[i][0]-self.boundaries[i][0]) < span[i]/bin_number :
	ticks[i][0] = self.boundaries[i][0]
	x_range[i][0] = self.boundaries[i][0]
      if abs(x_range[i][1]-self.boundaries[i][1]) < span[i]/bin_number :
	ticks[i][2] = self.boundaries[i][1]
	x_range[i][1] = self.boundaries[i][1]
      
    #fig.subplots_adjust(bottom=0.03, left=.02, right=0.98, top=0.98, hspace=.35)
    fig.subplots_adjust(bottom=0.03, left=.04, right=0.98, top=0.98, hspace=.35)

    for i in range(len(self.ref_names)):

      ax=fig.add_subplot(len(self.ref_names),len(self.ref_names),i*(len(self.ref_names)+1)+1,yticks=[])

      num_column = round(math.sqrt(len(self.ref_names)))
      ax1d = fig2.add_subplot(num_column,round(len(self.ref_names)*1.0/num_column),i,yticks=[])

      # histogram
      n,bin_edges=np.histogram(chain[:,i+2],bins=bin_number,weights=chain[:,0],normed=False)
      bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])
      #self.minimum_credible_intervals(n,bincenters,lvls)
      #exit()
      ax.set_xticks(ticks[i])
      ax.set_xticklabels(['%.4g' % s for s in ticks[i]])
      ax.set_title('%s= %.4g' % (self.tex_names[i],mean[i]))
      #ax.plot(bincenters,n,color='red',linewidth=2,ls='steps')
      ax.plot(bincenters,n,color='red',linewidth=2,ls='-')
      ax.axis([x_range[i][0], x_range[i][1],0,np.max(n)])

      ax1d.set_xticks(ticks[i])
      ax1d.set_xticklabels(['%.4g' % s for s in ticks[i]])
      ax1d.set_title('%s= %.4g' % (self.tex_names[i],mean[i]))
      ax1d.plot(bincenters,n,color='red',linewidth=2,ls='-')
      ax1d.axis([x_range[i][0], x_range[i][1],0,np.max(n)])
      # mean likelihood (optional)
      
      mean=np.zeros(len(bincenters),'float64')
      norm=np.zeros(len(bincenters),'float64')
      for j in range(len(bin_edges)-1):
	for k in range(np.shape(chain)[0]):
	  if (chain[k,i+2]>bin_edges[j] and chain[k,i+2]<bin_edges[j+1]):
	    mean[j] += chain[k,1]*chain[k,0]
	    norm[j] += chain[k,0]
	mean[j] /= norm[j]
      mean *= max(n)/max(mean)
      ax.plot(bincenters,mean,color='red',ls='--')
      ax1d.plot(bincenters,mean,color='red',ls='--')

      for j in range(i):
	ax1=fig.add_subplot(len(self.ref_names),len(self.ref_names),(i)*len(self.ref_names)+j+1)
	n,xedges,yedges=np.histogram2d(chain[:,i+2],chain[:,j+2],weights=chain[:,0],bins=(bin_number,bin_number),normed=False)
	extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
	x_centers = 0.5*(xedges[1:]+xedges[:-1])
	y_centers = 0.5*(yedges[1:]+yedges[:-1])

	ax1.set_xticks(ticks[j])
	if i == len(self.ref_names)-1:
	  ax1.set_xticklabels(['%.4g' % s for s in ticks[j]])
	else:
	  ax1.set_xticklabels([''])

	ax1.set_yticks(ticks[i])
	if j == 0:
	  ax1.set_yticklabels(['%.4g' % s for s in ticks[i]])
	else:
	  ax1.set_yticklabels([''])
	ax1.imshow(n.T, extent=extent, aspect='auto',interpolation='gaussian',origin='lower',cmap=matplotlib.cm.Reds)

	# smoothing the histogram, to have nicer contours
	#n = self.smoothing_hist(n,200)
	# plotting contours
	cs = ax1.contour(y_centers,x_centers,n.T,extent=extent,levels=self.ctr_level(n.T,lvls),colors="k",zorder=5)
	ax1.clabel(cs, cs.levels[:-1], inline=True,inline_spacing=0, fmt=dict(zip(cs.levels[:-1],[r"%d \%%"%int(l*100) for l in lvls[::-1]])), fontsize=6)


    fig.savefig(self.folder+'{0}_triangle.pdf'.format(self.folder.split('/')[-2]))
    fig2.savefig(self.folder+'{0}_1d.pdf'.format(self.folder.split('/')[-2]))

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
    norm = float((sum(histogram)-0.5*(histogram[0]+histogram[-1]))*(bincenter[-1]-bincenter[0]))
    print histogram
    print histogram[-1]
    print bincenter
    print norm
    for level in levels:
      water_level_up   = max(histogram)
      water_level_down = 0
      
      while ((water_level_up-water_level_down > max(histogram)/100000.)):
	top=0
	water_level = (water_level_up + water_level_down)/2.
	indices = [i for i in range(len(histogram)) if histogram[i]>water_level]
	# check for multimodal posteriors
	if ((indices[-1]-indices[0]+1)!=len(indices)):
	  print '\n\n  Can not derive minimum credible intervals for this multimodal posterior'
	  return False
	top = (sum(histogram[indices])-0.5*(histogram[indices[0]]+histogram[indices[-1]]))*(bincenter[indices[-1]]-bincenter[indices[0]])

	# left
	if indices[0]>0:
	  top += 0.5 * (water_level + histogram[indices[0]]) * (bincenter[indices[0]] - bincenter[indices[0]-1] )*(histogram[indices[0]] - water_level)/(histogram[indices[0]]-histogram[indices[0]-1]) 

	# right
	if indices[-1]<(len(histogram)-1) :
	  top += 0.5 * (water_level + histogram[indices[-1]]) * (bincenter[indices[-1]+1] - bincenter[indices[-1]]) * (histogram[indices[-1]] - water_level) / (histogram[indices[-1]]-histogram[indices[-1]+1])
	
	if top/norm > level:
	  water_level_down = water_level
	else:
	  water_level_up = water_level
	print top/norm,level,water_level_down,water_level_up
	

    exit()
    return sigmas
