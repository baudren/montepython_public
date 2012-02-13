import os,sys                             # Basic python module to handle file writing, and obtaining system information
import math
import random as rd 
import numpy  as np                       # Numerical Python module
from collections import OrderedDict as od # A modified version of Python dictionary, 
					  # in order to keep track of the order in it (much the same as in an array)
from datetime    import date              

import io # Needs to talk to io.py file for the logging of parameters

class data:

  # - path contains the configuration you inputed in your .conf file,
  # - the default flag is here to distinguish whether it is a genuine
  # initialization of data (default=True), or whether it is a comparison with
  # the log.param of an existing folder (default=False)
  def __init__(self,command_line,path,default=True):

    # Initialisation of the random seed
    rd.seed()

    # Distinguish between the two cases, genuine initialization or comparison
    # with existing data (in this case, grab the log.param)
    if default:
      self.param = command_line.param
    else:
      self.param = command_line.folder+'/log.param'

    # Recover jumping method from command_line
    self.jumping = command_line.jumping
    self.path = path

    # Creation of the two main dictionnaries:

    # -- Class_arguments, that will contain everything for the cosmological code
    # Class (no need to keep track of the order), and will be updated from:

    # -- mcmc_parameters, an ordered dictionary of dictionaries that will contain
    # everything needed by the Monte-Carlo procedure. Every parameter name will
    # be the key of a dictionary, containing: -initial configuration, role,
    # status last_accepted point, and current point

    self.Class_arguments = {}
    self.mcmc_parameters = od()

    # Read from the parameter file to fill properly the mcmc_parameters dictionary.
    self.fill_mcmc_parameters()

    # From mcmc_parameters, update Class_arguments (Class_arguments will always
    # contain the latest position, and forget anything about previous points,
    # whereas in mcmc_parameters, we will keep track of the last accepted point
    self.update_Class_arguments()

    # Recover Class version and subversion,
    if default:
      svn_file = open(path['class']+'/include/svnversion.h','r')
      self.subversion = svn_file.readline().split()[-1].replace('"','')
      svn_file.close()
      for line in open(path['class']+'/include/common.h','r'):
	if line.find('_VERSION_')!=-1:
	  self.version = line.split()[-1].replace('"','')
	  break
    else: # read in the existing parameter file
      self.read_version(self.param_file)

    # End of initialisation with the parameter file
    self.param_file.close()

    # log_flag, initially at False, will help determine if the code should
    # log the parameter file in the folder
    log_flag = False

    # For a true initialization, one should then initialize the likelihoods.
    # This step is obviously skipped for a comparison
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
	# Logging of parameters
        io.log_parameters(self,command_line)
	log_flag = True

    self.lkl=dict()
      
    # adding the likelihood directory to the path, to import the module
    # then, for each library, calling an instance of the likelihood.
    # Beware, though, if you add new likelihoods, they should go to the
    # folder likelihoods/yourlike/yourlike.py, and contain a yourlike.data,
    # otherwise the following set of commands will not work anymore.

    # For the logging if log_flag is True, each likelihood will log its
    # parameters

    for elem in self.exp:

      folder = os.path.abspath(path['MontePython'])+ "/../likelihoods/%s" % elem
      # add the folder of the likelihood to the path of libraries to...
      if folder not in sys.path:
	sys.path.insert(0, folder)
      # ... import easily the likelihood.py program
      exec "import %s" % elem
      # Initialize the likelihoods. Depending on the values of command_line,
      # log_flag and default, the routine will call slightly different things.
      # If log_flag, log.param will be appended. If default, some
      # precomputation will be made. Finally if not default, only a dictionary
      # containing the .data file will be created, for comparison purpose.
      exec "self.lkl['%s'] = %s.%s('%s/%s.data',self,command_line,log_flag,default)"% (elem,elem,elem,folder,elem)

    # Finally, log the Class_arguments used. This comes in the end, because
    # it can be modified inside the likelihoods init functions
    if log_flag:
      io.log_Class_arguments(self,command_line)
      io.log_default_configuration(self,command_line)

  # Redefinition of the 'compare' method for two instances of this data class.
  # It will decide which basic operations to perform when the code asked if two
  # instances are the same (in case you want to launch a new chain in an
  # existing folder, with your own parameter file)
  def __cmp__(self,other):

    # Comparing Class versions (warning only, will not fail the comparison)
    if self.version != other.version:
      print '/!\ Warning, you are running with a different version of Class'

    # Defines unordered version of the dictionaries of parameters
    self.uo_params  = {}
    other.uo_params = {}

    # Check if all the experiments are tested again,
    if len(list(set(other.exp).symmetric_difference(set(self.exp))))==0: 
      # Check that they have been called with the same .data file, stored in
      # dictionary when initializing.
      for exp in self.exp:
	for elem in self.lkl[exp].dictionary:
	  if self.lkl[exp].dictionary[elem]!=other.lkl[exp].dictionary[elem]:
	    print 'in your parameter file: ',self.lkl[exp].dictionary
	    print 'in log.param:           ',other.lkl[exp].dictionary
	    return -1

      # Fill in the unordered version of dictionaries
      for key,elem in self.mcmc_parameters.iteritems():
	self.uo_params[key]=elem['initial']
      for key,elem in other.mcmc_parameters.iteritems():
	other.uo_params[key]=elem['initial']


      # And finally compare them (standard comparison between dictionnaries,
      # will return True if both have the same keys and values associated to
      # them.
      return cmp(self.uo_params,other.uo_params) 
    else:
      return -1
    
  # Method defined to read the parameter file
  def read_file(self,_file):
    for line in _file:
      if line.find('#')==-1:
	if line.split('=')[0].find('data.')!=-1:
	  exec(line.replace('data.','self.'))
    _file.seek(0)

  # Extract version and subversion from an existing log.param
  def read_version(self,_file):
    # Read the first line (Class version)
    first_line = _file.readline()
    self.version = first_line.split()[1]
    self.subversion = first_line.split()[-1].replace(')','').replace('-','')
    _file.seek(0)

  # Initializes the ordered dictionary mcmc_parameters
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

    # Transform from params dictionnary to mcmc_parameters dictionary of
    # dictionaries, method defined just below
    self.from_input_to_mcmc_parameters(self.params)


  def from_input_to_mcmc_parameters(self,dictionary):
    # At the end of this initialization, every field but one is filled for
    # every parameter, be it fixed or varying. The missing field is the
    # 'last_accepted' one, that will be filled in in the mcmc part.
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

  # Method that returns a convenient, ordered array of parameter names that
  # correspond of the table_of_strings argument. 
  
  # For instance, if table_of_strings=['varying'], this routine will return all
  # the varying parameters in mcmc_parameters, cosmological or nuisance
  # parameters indifferently. If asked with ['varying','nuisance'], only the
  # nuisance parameters of the above will be returned.
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

  # Put in Class_arguments the current values of mcmc_parameters
  def update_Class_arguments(self):
    # For all elements in any cosmological parameters
    for elem in self.get_mcmc_parameters(['cosmo']):
      try:
	# If they have a current value already, use it
	self.Class_arguments[elem]   = self.mcmc_parameters[elem]['current']*self.mcmc_parameters[elem]['initial'][4]
      except KeyError: # it will go there if there is not yet a 'current' field,
	pass           # In this case, nothing to do.

    for elem in self.get_mcmc_parameters(['cosmo']):
      if elem == 'Omega_Lambda':
        try:
          omega_b      = self.Class_arguments['omega_b']
          omega_cdm    = self.Class_arguments['omega_cdm']
          Omega_Lambda = self.Class_arguments['Omega_Lambda']
          self.Class_arguments['h']   = math.sqrt( (omega_b+omega_cdm) / (1.-Omega_Lambda) )
          del self.Class_arguments[elem]
        except (KeyError):
          pass
