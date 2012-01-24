import os,sys
import random as rd
import numpy  as np
from collections import OrderedDict as od
from datetime    import date

import mcmc
import io

class data:

  def __init__(self,command_line,path,default=True):

    # Initialisation of the random seed
    rd.seed()

    # Distinguish between the two cases, genuine initialization or comparison
    # with existing data
    if default:
      self.param = command_line.par
    else:
      self.param = command_line.folder+'/log.param'

    # Recover jumping method from command_line
    self.jumping = command_line.jumping
    self.path = path

    # Creation of the two main dictionnaries:

    # -- Class_arguments, that will contain everything for the cosmological code
    # Class, and will be updated from:

    # -- mcmc_parameters, an ordered dictionary of dictionaries that will contain
    # everything needed by the Monte-Carlo procedure. Every parameter name will
    # be the key of a dictionary, containing: -initial configuration, role,
    # status last_accepted point, and current point

    self.Class_arguments = {}
    self.mcmc_parameters = od()

    # Read from the parameter file to fill properly the mcmc_parameters dictionary.
    self.fill_mcmc_parameters()

    self.update_Class_arguments()

    # Recover Class version and subversion
    if default:
      svn_file = open(path['class']+'/include/svnversion.h','r')
      self.subversion = svn_file.readline().split()[-1].replace('"','')
      svn_file.close()
      for line in open(path['class']+'/include/common.h','r'):
	if line.find('_VERSION_')!=-1:
	  self.version = line.split()[-1].replace('"','')
	  break
    else: # read in the parameter file
      self.read_version(self.param_file)

    # End of initialisation with the parameter file
    self.param_file.close()

    log_flag = False
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
	log_flag = True

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
	  exec "self.lkl['%s'] = %s.%s('%s/%s.data',self,command_line,log_flag)"% (elem,elem,elem,folder,elem)
	else:
	  exec "self.lkl['%s'] = %s.%s(self.param,self)"% (elem,elem,elem)
    
      if log_flag:
	io.log_Class_arguments(self,command_line)

  def __cmp__(self,other):
    self.uo_params  = {}
    other.uo_params = {}

    # Comparing Class versions
    if self.version != other.version:
      print '/!\ Warning, you are running with a different version of Class'

    # Check if all the experiments are tested again,
    if len(list(set(other.exp).symmetric_difference(set(self.exp))))==0: 
      # and if all the unordered parameters,fixed or varying, have the same value
      for key,elem in self.mcmc_parameters.iteritems():
	self.uo_params[key]=elem['initial']
      for key,elem in other.mcmc_parameters.iteritems():
	other.uo_params[key]=elem['initial']
      return cmp(self.uo_params,other.uo_params) 
    else:
      return -1
    
  def read_file(self,_file):
    for line in _file:
      if line.find('#')==-1:
	if line.split('=')[0].find('data.')!=-1:
	  exec(line.replace('data.','self.'))
    _file.seek(0)

  def read_version(self,_file):
    # Read the first line (Class version)
    first_line = _file.readline()
    self.version = first_line.split()[1]
    self.subversion = first_line.split()[-1].replace(')','').replace('-','')
    _file.seek(0)

  def fill_mcmc_parameters(self):

    # Define temporary quantities, only to simplify the input in the parameter
    # file
    self.params=od()

    # Read from the parameter file everything
    try:
      self.param_file = open(self.param,'r')
    except IOError:
      print "\n /|\  Error in initializing the data class,\n/_o_\ parameter file {0} does not point to a file".format(self.param)
      exit()
    self.read_file(self.param_file)

    self.from_input_to_mcmc_parameters(self.params)


  def from_input_to_mcmc_parameters(self,dictionary):
    for key,value in dictionary.iteritems():
      self.mcmc_parameters[key] = od()
      self.mcmc_parameters[key]['initial'] = value[0:5]
      self.mcmc_parameters[key]['role']    = value[-1]
      self.mcmc_parameters[key]['tex_name']= io.get_tex_name(key)
      if value[3] == 0:
	self.mcmc_parameters[key]['status']    = 'fixed'
	self.mcmc_parameters[key]['current']   = value[0]
      else:
	self.mcmc_parameters[key]['status']    = 'varying'

  def get_mcmc_parameters(self,table_of_strings):
    table = []
    for key,value in self.mcmc_parameters.iteritems():
      number = 0
      for subkey,subvalue in value.iteritems():
	for string in table_of_strings:
	  if subvalue == string:
	    number += 1
      if number == len(table_of_strings):
	table.append(key)
    return table

  def update_Class_arguments(self):
    for elem in self.get_mcmc_parameters(['cosmo']):
      try:
	self.Class_arguments[elem]   = self.mcmc_parameters[elem]['current']*self.mcmc_parameters[elem]['initial'][4]
      except KeyError:
	pass
