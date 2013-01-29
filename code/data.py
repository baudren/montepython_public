import os,sys                             # Basic python module to handle file writing, and obtaining system information
import math
import random as rd 
import numpy  as np                       # Numerical Python module

# A modified version of Python dictionary
# in order to keep track of the order in it (much the same as in an array)
try:
  from collections import OrderedDict as od 
except: 
  # in case an older version of Python is used, this module does not belong to
  # collections. Please remind to put that into your PYTHONPATH variable.
  from ordereddict import OrderedDict as od
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
      self.param = command_line.folder+'log.param'

    # Recover jumping method from command_line
    self.jumping 	= command_line.jumping
    self.jumping_factor = command_line.jumping_factor
    self.path           = path

    # Define the boundary loglike, the value used to defined a loglike that is
    # out of bounds. If a loglike is affected to this value, it will
    # automatically rejected
    self.boundary_loglike = -1e30

    # Creation of the two main dictionnaries:

    # -- cosmo_arguments, that will contain everything for the cosmological
    # code (so far, only Class), and will be updated from:

    # -- mcmc_parameters, an ordered dictionary of dictionaries that will contain
    # everything needed by the Monte-Carlo procedure. Every parameter name will
    # be the key of a dictionary, containing: -initial configuration, role,
    # status last_accepted point, and current point

    self.cosmo_arguments = {}
    self.mcmc_parameters = od()

    # Read from the parameter file to fill properly the mcmc_parameters dictionary.
    self.fill_mcmc_parameters()

    # Determine which cosmological code is in use
    if path['cosmo'].find('class') != -1:
      self.cosmological_module_name = 'Class'
    else:
      self.cosmological_module_name = None

    # Recover the cosmological code version (and subversion if relevant). To
    # implement a new cosmological code, please add another case to the test
    # below.
    if default:
      if self.cosmological_module_name == 'Class':
        svn_file = open(path['cosmo']+'/include/svnversion.h','r')
        self.subversion = svn_file.readline().split()[-1].replace('"','')
        svn_file.close()
        for line in open(path['cosmo']+'/include/common.h','r'):
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

    # Record if the cosmological parameters were changed (slow step)
    self.need_cosmo_update = True

    # For a true initialization, one should then initialize the likelihoods.
    # This step is obviously skipped for a comparison
    if default:

      sys.stdout.write('Testing likelihoods for:\n -> ')
      for i in range(len(self.experiments)):
	sys.stdout.write(self.experiments[i]+', ')
      sys.stdout.write('\n')
      
      # logging the parameter file (only if folder does not exist !)
      if command_line.folder[-1]!='/':
	command_line.folder+='/'
      if (os.path.exists(command_line.folder) and not os.path.exists(command_line.folder+'log.param')):
        if command_line.param != None:
          print '/!\   Detecting empty folder, logging the parameter file'
          io.log_parameters(self,command_line)
          log_flag = True
      if not os.path.exists(command_line.folder) :
	os.mkdir(command_line.folder)
	# Logging of parameters
        io.log_parameters(self,command_line)
	log_flag = True

    self.lkl=od()
      
    # adding the likelihood directory to the path, to import the module
    # then, for each library, calling an instance of the likelihood.
    # Beware, though, if you add new likelihoods, they should go to the
    # folder likelihoods/yourlike/yourlike.py, and contain a yourlike.data,
    # otherwise the following set of commands will not work anymore.

    # For the logging if log_flag is True, each likelihood will log its
    # parameters

    for elem in self.experiments:

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
      
    # Storing parameters by blocks of speed
    self.group_parameters_in_blocks()

    # Finally, log the cosmo_arguments used. This comes in the end, because
    # it can be modified inside the likelihoods init functions
    if log_flag:
      io.log_cosmo_arguments(self,command_line)
      io.log_default_configuration(self,command_line)


  # Block mcmc parameters in terms of speed (this will make as many blocks as
  # they are likelihoods with nuisance parameters plus one (the slow block of
  # cosmology))
  # IMPORTANT NOTE: this routine assumes that the nuisance parameters are
  # already written sequentially, and grouped together (not necessarilly in the
  # order described in data.experiments). If you mix up the different nuisance
  # parameters in the .param file, this routine will not function as intended.
  def group_parameters_in_blocks(self):
    array = []
    # First obvious block is all cosmological parameters
    array.append(len(self.get_mcmc_parameters(['varying','cosmo'])))

    # Then, store all nuisance parameters
    nuisance = self.get_mcmc_parameters(['varying','nuisance'])

    # Then circle through them
    index = 0
    while index < len(nuisance):
      elem = nuisance[index]
      flag = False
      # For each one, check if they belong to a likelihood
      for likelihood in self.lkl.itervalues():
        if elem in likelihood.use_nuisance:
          # If yes, store the number of nuisance parameters needed for this likelihood.
          flag = True
          array.append(len(likelihood.use_nuisance)+array[-1])
          index += len(likelihood.use_nuisance)
          continue
      if not flag:
        # If the loop reaches this part, it means this nuisance parameter was
        # associated with no likelihood: this should not happen
        print '/!\ nuisance parameter {0} is associated to no likelihood\n'.format(elem)
        exit()


    # Store the result
    self.blocks_parameters = array


  # Redefinition of the 'compare' method for two instances of this data class.
  # It will decide which basic operations to perform when the code asked if two
  # instances are the same (in case you want to launch a new chain in an
  # existing folder, with your own parameter file)
  def __cmp__(self,other):

    # Comparing cosmological code versions (warning only, will not fail the comparison)
    if self.version != other.version:
      print '/!\ Warning, you are running with a different version of your cosmological code'

    # Defines unordered version of the dictionaries of parameters
    self.uo_parameters  = {}
    other.uo_parameters = {}

    # Check if all the experiments are tested again,
    if len(list(set(other.experiments).symmetric_difference(set(self.experiments))))==0: 
      # Check that they have been called with the same .data file, stored in
      # dictionary when initializing.
      for experiment in self.experiments:
	for elem in self.lkl[experiment].dictionary:
	  if self.lkl[experiment].dictionary[elem]!=other.lkl[experiment].dictionary[elem]:
	    print 'in your parameter file: ',self.lkl[experiment].dictionary
	    print 'in log.param:           ',other.lkl[experiment].dictionary
	    return -1

      # Fill in the unordered version of dictionaries
      for key,elem in self.mcmc_parameters.iteritems():
	self.uo_parameters[key]=elem['initial']
      for key,elem in other.mcmc_parameters.iteritems():
	other.uo_parameters[key]=elem['initial']


      # And finally compare them (standard comparison between dictionnaries,
      # will return True if both have the same keys and values associated to
      # them.
      return cmp(self.uo_parameters,other.uo_parameters) 
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
    # Read the first line (cosmological code version)
    first_line = _file.readline()
    self.version = first_line.split()[1]
    self.subversion = first_line.split()[-1].replace(')','').replace('-','')
    _file.seek(0)

  # Initializes the ordered dictionary mcmc_parameters
  def fill_mcmc_parameters(self):

    # Define temporary quantities, only to simplify the input in the parameter
    # file
    self.parameters=od()

    # Read from the parameter file everything
    try:
      self.param_file = open(self.param,'r')
    except IOError:
      print "\n /|\  Error in initializing the data class,\n/_o_\ parameter file {0} does not point to a file".format(self.param)
      exit()
    self.read_file(self.param_file)

    # Transform from parameters dictionnary to mcmc_parameters dictionary of
    # dictionaries, method defined just below
    self.from_input_to_mcmc_parameters(self.parameters)


  def from_input_to_mcmc_parameters(self,dictionary):
    # At the end of this initialization, every field but one is filled for
    # every parameter, be it fixed or varying. The missing field is the
    # 'last_accepted' one, that will be filled in in the mcmc part.
    for key,value in dictionary.iteritems():
      self.mcmc_parameters[key] = od()
      self.mcmc_parameters[key]['initial'] = value[0:4]
      self.mcmc_parameters[key]['scale']   = value[4]
      self.mcmc_parameters[key]['role']    = value[-1]
      self.mcmc_parameters[key]['tex_name']= io.get_tex_name(key)
      if value[3] == 0:
	self.mcmc_parameters[key]['status']    = 'fixed'
	self.mcmc_parameters[key]['current']   = value[0]
      else :
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

  # Routine to determine whether the value of cosmological parameters were
  # changed, and if no, to skip recomputation.
  def check_for_slow_step(self,new_step):

    parameter_names = self.get_mcmc_parameters(['varying'])
    cosmo_names     = self.get_mcmc_parameters(['cosmo'])

    need_change = 0

    # For all elements in the varying parameters:
    for elem in parameter_names:
      i = parameter_names.index(elem)
      # If it is a cosmological parameter
      if elem in cosmo_names:
        if self.mcmc_parameters[elem]['current'] != new_step[i]:
          need_change += 1

    # If any cosmological value was changed, 
    if need_change > 0:
      self.need_cosmo_update = True
    else:
      self.need_cosmo_update = False

    
    for likelihood in self.lkl.itervalues():
      # If the cosmology changed, you need to recompute the likelihood anyway
      if self.need_cosmo_update:
        likelihood.need_update = True
        continue
      # Otherwise, check if the nuisance parameters of this likelihood were
      # changed
      need_change = 0
      for elem in parameter_names:
        i = parameter_names.index(elem)
        if elem in likelihood.use_nuisance:
          if self.mcmc_parameters[elem]['current'] != new_step[i]:
            need_change += 1
      if need_change>0:
        likelihood.need_update = True
      else:
        likelihood.need_update = False


  # Put in cosmo_arguments the current values of mcmc_parameters
  def update_cosmo_arguments(self):
    # For all elements in any cosmological parameters
    for elem in self.get_mcmc_parameters(['cosmo']):
      # Fill in the dictionnary with the current value of parameters
      self.cosmo_arguments[elem]   = self.mcmc_parameters[elem]['current']*self.mcmc_parameters[elem]['scale']

    # For all elements in the cosmological parameters from the mcmc list,
    # translate any-one that is not directly a Class parameter into one.
    # The try: except: syntax ensures that the first call 
    for elem in self.get_mcmc_parameters(['cosmo']):
      # infer h from Omega_Lambda and delete Omega_Lambda
      if elem == 'Omega_Lambda':
        omega_b      = self.cosmo_arguments['omega_b']
        omega_cdm    = self.cosmo_arguments['omega_cdm']
        Omega_Lambda = self.cosmo_arguments['Omega_Lambda']
        self.cosmo_arguments['h']   = math.sqrt( (omega_b+omega_cdm) / (1.-Omega_Lambda) )
        del self.cosmo_arguments[elem]
      # infer omega_cdm from Omega_L and delete Omega_L
      if elem == 'Omega_L':
        omega_b      = self.cosmo_arguments['omega_b']
        h = self.cosmo_arguments['h']
        Omega_L = self.cosmo_arguments['Omega_L']
        self.cosmo_arguments['omega_cdm']  = (1.-Omega_L)*h*h-omega_b
        del self.cosmo_arguments[elem]
      if elem == 'ln10^{10}A_s':
        self.cosmo_arguments['A_s']   = math.exp(self.cosmo_arguments[elem])/1.e10
        del self.cosmo_arguments[elem]
      if elem == 'exp_m_2_tau_As':
        tau_reio = self.cosmo_arguments['tau_reio']
        self.cosmo_arguments['A_s']   = self.cosmo_arguments[elem]*math.exp(2.*tau_reio)
        del self.cosmo_arguments[elem]
      if elem == 'f_cdi':
        self.cosmo_arguments['n_cdi']   =self.cosmo_arguments['n_s']
      if elem == 'beta':
        self.cosmo_arguments['alpha'] = 2.*self.cosmo_arguments['beta']
      # We only do that on xe_1, for there is at least one of them.
      if elem.find('xe_1') != -1:
        # To pass this option, you must have set a number of cosmological settings
        # reio_parametrization to reio_bins_tanh, binned_reio_z set, and binned_reio_num
        # First, you need to set reio_parametrization to reio_bins_tanh
        if (self.cosmo_arguments['reio_parametrization'] != 'reio_bins_tanh'):
          print ' /|\  Warning, you set binned_reio_xe to some values'
          print '/_o_\ without setting reio_parametrization to reio_bins_tanh'
          exit()
        else:
          try:
            size = self.cosmo_arguments['binned_reio_num']
          except (KeyError):
            print ' /|\  You need to set reio_binnumber to the value corresponding to'
            print '/_o_\ the one in binned_reio_xe'
            exit()
        string = ''
        for i in range(1,size+1):
          string += '%.4g' % self.cosmo_arguments['xe_%d' % i]
          del self.cosmo_arguments['xe_%d' % i]
          if i != size:
            string += ','
        self.cosmo_arguments['binned_reio_xe'] = string
