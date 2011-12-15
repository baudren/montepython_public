import os,sys
import random as rd
import numpy  as np
from collections import OrderedDict as od
from datetime    import date

import mcmc
import io

class data:

  def __init__(self,command_line,path,default=True):
    # First of all, distinguish between the two cases, genuine initialization
    # or comparison with existing data
    if default:
      self.param = command_line.par
    else:
      self.param = command_line.folder+'/log.param'

    data.path = path

    # Initialisation of all the data needed to run the Boltzmann code
    # First apply default parameters,
    # Then try and modify according to command_line.param custom parameter file
    # Finally overwrites with special command_line arguments.
    rd.seed()
    self.Class_params=od()
    self.Class_param_names=[]
    self.Class=[]    # Create the Class param vector

    self.nuisance_params=od()
    self.nuisance_param_names=[]
    self.nuisance=[] # Create the nuisance vector

    self.params=od()
    self.param_names=[]

    self.vector=[]   # Will contain the merging of the two above

    self.Class_args={} # Contains the arguments of the Class instance

    try:
      param_file = open(self.param,'r')
    except IOError:
      print "\n /|\  Error in initializing the data class,\n/_o_\ parameter file {0} does not point to a file".format(self.param)
      exit()
    self._read_file(param_file)

    # Recover Class version and subversion
    if default:
      svn_file = open(path['class']+'/include/svnversion.h','r')
      self.subversion = svn_file.readline().split()[-1].replace('"','')
      svn_file.close()
      for line in open(path['class']+'/include/common.h','r'):
	if line.find('_VERSION_')!=-1:
	  self.version = line.split()[-1].replace('"','')
	  break
    else:
      self._read_version(param_file)
    param_file.close()

    self._parameters()

    if default:
      sys.stdout.write('testing likelihoods for:\n')
      for i in range(len(self.exp)):
	sys.stdout.write(self.exp[i]+',\t')
      sys.stdout.write('\n')
      
      # logging the parameter file (only if folder does not exist !)
      if command_line.folder[-1]!='/':
	command_line.folder+='/'
      if not os.path.exists(command_line.folder):
	os.mkdir(command_line.folder)
        io.log_parameters(self,command_line)
      
      self.lkl=dict()
	
	# adding the likelihood directory to the path, to import the module
	# then, for each library, calling an instance of the likelihood.
	# Beware, though, if you add new likelihoods, they should go to the
	# folder likelihoods/yourlike/yourlike.py, and contain a yourlike.data,
	# otherwise the following set of commands will not work anymore.

      for elem in self.exp:

	folder = os.path.abspath(path['MontePython'])+ "/../likelihoods/%s" % elem
	if folder not in sys.path:
	  sys.path.insert(0, folder)
	exec "import %s" % elem
	if self.param.find('log.param')==-1:
	  exec "self.lkl['%s'] = %s.%s('%s/%s.data',self,command_line)"% (elem,elem,elem,folder,elem)
	else:
	  exec "self.lkl['%s'] = %s.%s(self.param,self)"% (elem,elem,elem)
    
      # log Class_args, and eventually all fixed parameters
      if self.param.find('log.param')==-1:
	io.log_Class_args(self,command_line)
  
    for i in range(len(self.Class)): # Initialize the arguments
      mcmc.jump(self,self.Class_param_names[i],self.Class[i])

  def __cmp__(self,other):
    self.uo_params  = {}
    # comparing Class versions
    if self.version != other.version:
      print '/!\ Warning, you are running with a different version of Class'
    if len(list(set(other.exp).symmetric_difference(set(self.exp))))==0: # if all the experiments are tested again
      for key,elem in self.params.iteritems():
	self.uo_params[key]=elem
      other.uo_params = {}
      for key,elem in other.params.iteritems():
	other.uo_params[key]=elem
      return cmp(self.uo_params,other.uo_params) # and if all the unordered parameters have the same value
    else:
      return -1
    
  def _read_file(self,_file):
    for line in _file:
      if line.find('#')==-1:
	if line.find('data.')!=-1:
	  exec(line.replace('data.','self.'))
    _file.seek(0)

  def _read_version(self,_file):
    # Read the first line (Class version)
    first_line = _file.readline()
    self.version = first_line.split()[1]
    self.subversion = first_line.split()[-1].replace(')','').replace('-','')
    _file.seek(0)

  def _parameters(self):
    failure=False
    for key,value in self.Class_params.iteritems():
      if value[3] == 0:
	print 'nope'
      self.Class_param_names.append(key)
      self.Class_params[key] = value
      temp=rd.gauss(value[0],value[3])
      while (value[1]!=-1 and temp<value[1] and failure==False):
	if value[1]>value[0]:
	  print('  Warning: you might have inconsistently set the min boundary for {0} parameter'.format(key))
	  failure=True
	temp=rd.gauss(value[0],value[3])
      failure=False
      while (value[2]!=-1 and temp>value[2] and failure==False):
	if value[2]<value[0]:
	  print '  Warning: you might have inconsistently set the max boundary for {0} parameter'.format(key)
	  failure=True
	temp=rd.gauss(value[0],value[3])
      self.Class.append(temp)
    
    for key,value in self.nuisance_params.iteritems():
      self.nuisance_param_names.append(key)
      self.nuisance_params[key] = value
      temp=rd.gauss(value[0],value[3])
      while (value[1]!=-1 and temp<value[1] and failure==False):
	if value[1]>value[0]:
	  print '  Warning: you might have inconsistently set the min boundary for {0} parameter'.format(key)
	  failure=True
	temp=rd.gauss(value[0],value[3])
      failure=False
      while (value[2]!=-1 and temp>value[2] and failure==False):
	if value[2]<value[0]:
	  print '  Warning: you might have inconsistently set the max boundary for {0} parameter'.format(key)
	  failure=True
	temp=rd.gauss(value[0],value[3])
      self.nuisance.append(temp)
      


    self.params.update(self.Class_params)
    self.params.update(self.nuisance_params)
    self.param_names = np.append(self.Class_param_names,self.nuisance_param_names)
    self._update_vector()

  def _nuisance_parameters(self):
    for key,value in self.nuisance_params.iteritems():
      self.nuisance_param_names.append(key)
      temp=rd.gauss(value[0],value[3])
      while (value[1]!=-1 and temp<value[1] and failure==False):
	if value[1]>value[0]:
	  print '  Warning: you might have inconsistently set the min boundary for {0} parameter'.format(key)
	  failure=True
	temp=rd.gauss(value[0],value[3])
      failure=False
      while (value[2]!=-1 and temp>value[2] and failure==False):
	if value[2]<value[0]:
	  print '  Warning: you might have inconsistently set the max boundary for {0} parameter'.format(key)
	  failure=True
	temp=rd.gauss(value[0],value[3])
      self.nuisance.append(temp)

  def _update_vector(self):
    self.vector = np.append(self.Class,self.nuisance)

  def _transmit_vector(self,vector):
    for i in range(len(self.Class)):
      self.Class[i]    = vector[i]
    for i in range(len(self.Class),len(self.Class)+len(self.nuisance)):
      self.nuisance[i-len(self.Class)] = vector[i]
