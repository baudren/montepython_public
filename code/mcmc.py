from classy import Class
import os,sys
import math
import random as rd
import numpy as np
import io
import data

# COMPUTE LKL
def compute_lkl(_cosmo,data):
  # Prepare the cosmological module with the new set of parameters
  _cosmo.set(data.Class_arguments)

  # Compute the model
  failure=False
  try:
    _cosmo._compute(["lensing"])
  except NameError :
    return True,0 
  except KeyboardInterrupt:
    exit()

  # For each desired likelihood, compute its value against the theoretical model
  loglike=0
  for likelihood in data.lkl.itervalues():
    loglike+=likelihood.loglkl(_cosmo,data)

  # Clean the cosmological strucutre
  _cosmo._struct_cleanup(set(["lensing","nonlinear","spectra","primordial","transfer","perturb","thermodynamics","background","bessel"]))

  return failure,loglike

def read_args_from_chain(data,chain):
  # Chain is defined as a special File (that inherits from File, and has two
  # new method, head and tail). The class is defined in code/io.py
  Chain = io.File(chain,'r')
  parameter_names = data.get_mcmc_parameters(['varying'])

  # BE CAREFUL: Here it works because of the particular presentation of the
  # chain, and the use of tabbings. Please keep this in mind if having
  # difficulties
  i = 1
  for elem in parameter_names:
    data.mcmc_parameters[elem]['last_accepted'] = float(Chain.tail(1)[0].split('\t')[i])
    i+=1

def get_cov(data,command_line):
  np.set_printoptions(precision=2,linewidth=150)
  parameter_names = data.get_mcmc_parameters(['varying'])
  i=0

  # if the user wants to use a covmat file
  if command_line.cov is not None:
    cov=open('{0}'.format(command_line.cov),'r')
    for line in cov:
      if line.find('#')!=-1:
	covnames = line.strip('#[').strip(']\n').replace("'","").split(', ')
	M=np.zeros((len(covnames),len(covnames)),'float64')
	rot=np.zeros((len(covnames),len(covnames)))
      else:
	line=line.split()
	for j in range(len(line)):
	  M[i][j]=np.array(line[j],'float64')
	i+=1

    # Deal with the all problematic cases. 

    # First, rotate M for the parameters to be well ordered, even if some names
    # are missing or some are in extra.
    temp_names = []
    for elem in parameter_names:
      if elem in covnames:
	temp_names.append(elem)
    for k in range(len(covnames)):
      for h in range(len(covnames)):
	try:
	  if covnames[k]==temp_names[h]:
	    rot[h][k] = 1.
	  else:
	    rot[h][k] = 0.
	except IndexError:
	  rot[h][k] = 0.
    print rot
    print covnames,parameter_names
    print '\nInput covariance matrix:'
    print covnames
    print M
    M=np.dot(rot,np.dot(M,rot))

    print '\nFirst treatment'
    print M
    
    M_temp    = np.ones((len(parameter_names),len(parameter_names)),'float64')
    indices_1 = np.zeros(len(parameter_names))
    indices_2 = np.zeros(len(covnames))
    #Remove names that are in param_names but not in covnames
    for k in range(len(parameter_names)):
      if parameter_names[k] in covnames:
	indices_1[k]=1
    for zeros in np.where(indices_1 == 0)[0]:
      M_temp[zeros,:] = 0
      M_temp[:,zeros] = 0
    #Remove names that are in covnames but not in param_names
    for h in range(len(covnames)):
      if covnames[h] in parameter_names:
	indices_2[h]=1
    print indices_2
    for zeros in np.where(indices_2 == 0)[0]:
      M[zeros,:] = 0
      M[:,zeros] = 0
    # super trick
    print M
    print M_temp
    M_temp[M_temp == 1]=M[M!=0]
    M = np.copy(M_temp)
    # on all other lines, just use sigma^2
    for zeros in np.where(indices_1 == 0)[0]:
      M[zeros,zeros] = np.array(data.mcmc_parameters[parameter_names[zeros]]['initial'][3],'float64')**2

  # else, take sigmas^2.
  else:
    M = np.identity(len(parameter_names),'float64')
    for elem in parameter_names:
      M[i][i]=np.array(data.mcmc_parameters[elem]['initial'][3],'float64')**2
      i+=1

  print '\nDeduced starting covariance matrix:'
  print parameter_names
  print M

  #inverse, and diagonalization
  eigv,eigV=np.linalg.eig(np.linalg.inv(M))
  return eigv,eigV

def get_new_pos(data,eigv,U,k):
  
  parameter_names = data.get_mcmc_parameters(['varying'])
  vector_new=np.zeros(len(parameter_names),'float64')
  sigmas=np.zeros(len(parameter_names),'float64')

  # Write the vector of last accepted points, 
  vector = np.zeros(len(parameter_names),'float64')
  try:
    for elem in parameter_names:
      vector[parameter_names.index(elem)] = data.mcmc_parameters[elem]['last_accepted']
  except KeyError: # Else take the mean value
    for elem in parameter_names:
      vector[parameter_names.index(elem)] = data.mcmc_parameters[elem]['initial'][0]

  flag=1
  while flag!=0:
    flag=0 # initialize: there are no problems at first
    rd.seed()

    # Choice here between sequential and global change of direction
    if data.jumping == 'global':
      for i in range(len(vector)):
	sigmas[i]=(math.sqrt(1/eigv[i]/len(vector)))*rd.gauss(0,1)*2.4
    elif data.jumping == 'sequential':
      i = k%len(vector)
      sigmas[i] = (math.sqrt(1/eigv[i]))*rd.gauss(0,1)*2.4
    else:
      print '\n\n  Jumping method unknown (accepted : global, sequential)'

    vector_new = vector + np.dot(U,sigmas)
    i=0
    for elem in parameter_names:
      value = data.mcmc_parameters[elem]['initial']
      if(value[1]!=-1 and vector_new[i]<value[1]):
	flag+=1 # if a boundary value is reached, increment
      elif(value[2]!=-1 and vector_new[i]>value[2]):
	flag+=1 # same
      i+=1

  i=0
  for elem in parameter_names:
    data.mcmc_parameters[elem]['current'] = vector_new[i]
    i+=1
    
  # Propagate the information towards the Class arguments
  data.update_Class_arguments()


def accept_step(data):
  for elem in data.get_mcmc_parameters(['varying']):
    data.mcmc_parameters[elem]['last_accepted'] = data.mcmc_parameters[elem]['current']

#---------------MCMC-CHAIN----------------------------------------------
def chain(_cosmo,data,command_line):
  # initialisation
  num_failure=10
  loglike=0
  failure=False
  sigma_eig,U=get_cov(data,command_line)
  failed=0

  # if restart wanted, pick initial value for arguments
  if command_line.restart is not None:
    read_args_from_chain(data,command_line.restart)

  # Pick a position (from last accepted point if restart, from the mean value
  # else)
  get_new_pos(data,sigma_eig,U,failed)
  # Compute the starting Likelihood
  failure,loglike=compute_lkl(_cosmo,data)

  while ((failure is True) and (failed<=num_failure)):
    failed +=1
    get_new_pos(data,sigma_eig,U,failed)
    failure,loglike=compute_lkl(_cosmo,data)

  if failure is True:
    print ' /|\   Class tried {0} times to initialize with given parameters, and failed...'.format(num_failure+2)
    print '/_o_\  You might want to change your starting values, or pick default ones!'
    exit()

  max_loglike = loglike

  acc,rej=0.0,0.0	# acceptance and rejection number count
  N=1			# number of time the system stayed in the current position
  io.print_parameters(sys.stdout,data)

  k = 1
  while (k <= command_line.N and failed <= num_failure):

    get_new_pos(data,sigma_eig,U,failed)
    failure,newloglike=compute_lkl(_cosmo,data)
    
    if(failure==True):
      failed += 1
      print 'Warning: Class failed due to choice of parameters, picking up new values'
      print data.Class_arguments
      continue

    # Harmless trick to avoid exponentiating large numbers
    if newloglike >= loglike:
      alpha = 1.
    else:
      alpha=np.exp(newloglike-loglike)

    if ((alpha == 1.) or (rd.uniform(0,1) < alpha)): #accept step

      accept_step(data)
      io.print_vector([data.out,sys.stdout],N,loglike,data)
      loglike=newloglike
      if loglike > max_loglike:
	max_loglike = loglike
      acc+=1.0
      N=1
      
    else:
      rej+=1.0
      N+=1
      
    if acc % data.write_step ==0:
      io.refresh_file(data)
    k += 1

  if (failed == num_failure):
    print ' /|\   The computation failed {0} times, \n'.format(failed)
    print '/_o_\  Please check the values of your parameters'

  rate=acc/(acc+rej)
  print '#  {0} steps done, acceptance rate:'.format(command_line.N),rate
  if command_line.restart is not None:
    os.remove(command_line.restart)
    print '  deleting starting point of the chain {0}'.format(command_line.restart)
  return rate,max_loglike
