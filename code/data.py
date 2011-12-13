import os,sys
import random as rd
import numpy  as np
from collections import OrderedDict as od
from datetime    import date

import mcmc

class data:

  def __init__(self,param,path,default=True):
    # Initialisation of all the data needed to run the Boltzmann code
    # First apply default parameters,
    # Then try and modify according to command_line.param custom parameter file
    # Finally overwrites with special command_line arguments.

    rd.seed()
    self.params=od()
    self.param_names=[]
    self.theta=[] # Create the theta vector
    self.args={}

    try:
      param_file = open(param,'r')
    except IOError:
      print "\n /|\  Error in initializing the data class,\n/_o_\ parameter file {0} does not point to a file".format(param)
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
	exec "self.lkl['%s'] = %s.%s('%s/%s.data')"% (elem,elem,elem,folder,elem)

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
      return False
    
  def _read_file(self,_file):
    for line in _file:
      if line.find('#')==-1:
	if line.find('data')!=-1:
	  exec(line.replace('data','self'))
    _file.seek(0)

  def _read_version(self,_file):
    # Read the first line (Class version)
    first_line = _file.readline()
    self.version = first_line.split()[1]
    self.subversion = first_line.split()[-1].replace(')','').replace('-','')
    _file.seek(0)

  def _parameters(self):
    failure=False
    for key,value in self.params.iteritems():
      self.param_names.append(key)
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
      self.theta.append(temp)

    for i in range(len(self.theta)): # Initialize the arguments
      mcmc.jump(self,self.param_names[i],self.theta[i])
