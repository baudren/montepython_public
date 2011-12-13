import os,sys
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
	  self.covar[i][j]+=chain[k,0]*((chain[k,i+2]-self.mean[i])*(chain[k,j+2]-self.mean[j]))
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
	    a95=i
	else:
	  a99=i

    # Plotting parameters
    font = {'size'   : 10}
    matplotlib.rc('font', **font)

    # 1D and 2D plots
    fig = plt.figure(1,(24,16))
    fig.subplots_adjust(bottom=0.03, left=.02, right=0.98, top=0.98, hspace=.35)

    for i in range(len(self.ref_names)):
      ax=fig.add_subplot(len(self.ref_names),len(self.ref_names),i*(len(self.ref_names)+1)+1,yticks=[])

      n,bins,patches=plt.hist(chain[:,i+2],bins=binnumber,weights=chain[:,0],normed=True,color='red')
      ax.set_xticks(np.linspace(round(min(chain[:,i+2]),3),round(max(chain[:,i+2]),3),5))
      ax.set_xticklabels(['%1.3f' % s for s in np.linspace(round(min(chain[:,i+2]),3),round(max(chain[:,i+2]),3),5)])
      ax.set_xlabel('{0}'.format(self.ref_names[i]))
      #y = mlab.normpdf( bins, mean[i], var[i])
      #y = Max*exp(-(mean[i], var[i])
      #print mean[i],var[i]
      #plt.plot(bins,1./np.sqrt(2*np.pi*np.sqrt(var[i]))*np.exp(-(bins-mean[i])**2/(2.*var[i])),color='blue')
      bins = np.linspace(min(bins),max(bins),1000)
      ax.plot(bins,1./np.sqrt(2*np.pi*self.var[i])*np.exp(-(bins-self.mean[i])**2/(2.*self.var[i])),color='blue',linewidth=2)

      #l = plt.plot(bins,y, 'r--', linewidth=1)
      #print y

      for j in range(i):
	ax1=fig.add_subplot(len(self.ref_names),len(self.ref_names),(i)*len(self.ref_names)+j+1)
	n,xedges,yedges=np.histogram2d(chain[:,i+2],chain[:,j+2],bins=(binnumber,binnumber))
	extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
	plt.imshow(n, extent=extent, aspect='auto',interpolation='nearest')

    fig.savefig(self.folder+'hist.pdf')
    #plt.show()
    self.info.write('\n parameter names:\t')
    for elem in self.ref_names:
      self.info.write('{0}\t'.format(elem))
    self.info.write('\n R values:\t')
    for elem in self.R:
      self.info.write('{0}\t'.format(elem))
    self.info.write('\n mean    :\t')
    for elem in self.mean:
      self.info.write('{0}\t'.format(elem))

  def prepare(self,Files):

    # Scan the whole input folder, and include all chains in it (all files with a
    # '.' in their name will be discarded as not mcmc chains (If you want to have
    # additionnal files than the log.dat and .param in this folder, please keep
    # this in  mind
    if os.path.isdir(Files[0]): 
      if Files[0].find('/')==-1:
	Files[0]+='/'
      Files = [Files[0]+elem for elem in os.listdir(Files[0]) if elem.find('.')==-1]

    # Recover the folder, depending on the case
    if (len(Files[0].split('/'))==0 or (Files[0].split('/')[0]=='.')):
      self.folder = './'
    else:
      self.folder = Files[0].split('/')[0]+'/'

    # Check if the log.dat file exists
    if os.path.isfile(self.folder+'log.dat') is True:
      if os.path.getsize(self.folder+'log.dat')>0:
	self.log = open(self.folder+'log.dat','r')
      else:
	print '\n\n  The companion log file {0} seems empty'.format(self.folder+'log.dat')
	exit()
    else:
      print '\n\n  The companion log file {0} is absent ?'.format(self.folder+'log.dat')
      exit()

    # Cleaning the folder if necessary
    numlines=sum(1 for line in self.log)
    self.log.seek(0)
    if numlines < len(Files):
      io.clean(self.folder)

    infoname = self.folder+self.folder.rstrip('/')+'.info'
    covname  = self.folder+self.folder.rstrip('/')+'.covmat'

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

    # Recovering default ordering of parameters (first line in log file)
    self.ref_names = self.log.readline().split('\t')[1].strip('[').strip(']').split()

    self.log.seek(0)

    for File in self.Files:
      i=self.Files.index(File)
      print 'scanning file {0}'.format(File)
      cheese = (np.array([[float(elem) for elem in line.split()] for line in open(File,'r')]))
      max_lkl = min(cheese[:,1]) # beware, it is the min because we are talking about '- log likelihood'

      start = 0
      while cheese[start,1]>max_lkl+4:
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
      
      # Recover names for this chain
      for line in self.log:
	if line.find(File.split('/')[-1])!=-1:
	  names = line.split('\t')[1].strip('[').strip(']').split()
      self.log.seek(0)

      # Verify the order is in agreement, otherwise, proceed to a swapping
      for name in names:
	if name != self.ref_names[names.index(name)]:
	  ham,names = swap(ham,names,names.index(name),self.ref_names.index(name))

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

  def swap(ham,names,i,j):
    # first swap the things in names
    temp 	   = names[i]
    names[i] = names[j]
    names[j] = temp

    # then swap the columns
    temp 	   = np.copy(ham[:,i])
    ham[:,i] = ham[:,j]
    ham[:,j] = temp
    return ham,names
