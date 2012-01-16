from classy import Class
import os,sys
import math
import random as rd
import numpy as np
import io
import data

# COMPUTE LKL
def compute_lkl(_cosmo,Data,args):
  # Prepare the cosmological module with the new set of parameters
  _cosmo.set(args)

  # Compute the model
  failure=False
  try:
    _cosmo._compute(["lensing"])
  except:
    return True,0 

  # For each desired likelihood, compute its value against the theoretical model
  loglike=0
  for likelihood in Data.lkl.itervalues():
    loglike+=likelihood._loglkl(_cosmo,Data)

  # Clean the cosmological strucutre
  _cosmo._struct_cleanup(set(["lensing","nonlinear","spectra","primordial","transfer","perturb","thermodynamics","background","bessel"]))

  return failure,loglike

def read_args_from_chain(data,chain):
  Chain = io.File(chain,'r')
  i = 1
  for elem in data.Class_param_names:
    data.Class_args[elem] = float(Chain.tail(1)[0].split('\t')[i])
    data.Class[i-1] = data.Class_args[elem]
    data.vector[i-1] = data.Class_args[elem]
    i+=1
  for elem in data.nuisance_param_names:
    data.vector[i-1] = float(Chain.tail(1)[0].split('\t')[i])
    data.nuisance[i-1-len(data.Class)] = data.vector[i-1]
    i+=1

def get_cov(data,command_line):
  np.set_printoptions(precision=2,linewidth=150)
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
    for elem in data.param_names:
      if elem in covnames:
	temp_names.append(elem)
    for k in range(len(covnames)):
      for h in range(len(covnames)):
	if covnames[k]==temp_names[h]:
	  rot[k][h] = 1.
	else:
	  rot[k][h] = 0.
    print M
    M=np.dot(rot,np.dot(M,rot))
    print M
    
    M_temp    = np.ones((len(data.param_names),len(data.param_names)),'float64')
    indices_1 = np.zeros(len(data.param_names))
    indices_2 = np.zeros(len(covnames))
    #Remove names that are in param_names but not in covnames
    for k in range(len(data.param_names)):
      if data.param_names[k] in covnames:
	indices_1[k]=1
    for zeros in np.where(indices_1 == 0)[0]:
      M_temp[zeros,:] = 0
      M_temp[:,zeros] = 0
    #Remove names that are in covnames but not in param_names
    for h in range(len(covnames)):
      if covnames[h] in data.param_names:
	indices_2[h]=1
    for zeros in np.where(indices_2 == 0)[0]:
      M[zeros,:] = 0
      M[:,zeros] = 0
    # super trick
    M_temp[M_temp == 1]=M[M!=0]
    M = np.copy(M_temp)
    # on all other lines, just use sigma^2
    for zeros in np.where(indices_1 == 0)[0]:
      M[zeros,zeros] = np.array(data.params[data.param_names[zeros]][3],'float64')**2
    print M



  # else, take sigmas^2.
  else:
    M = np.identity(len(data.vector),'float64')
    for elem in data.params:
      M[i][i]=np.array(data.params[elem][3],'float64')**2
      i+=1

  #inverse, and diagonalization
  eigv,eigV=np.linalg.eig(np.linalg.inv(M))
  return eigv,eigV

def get_new_pos(data,eigv,U,k):
  vector_new=np.zeros(len(data.vector),'float64')
  sigmas=np.zeros(len(data.vector),'float64')
  flag=1
  while flag!=0:
    flag=0 # initialize: there are no problems at first
    rd.seed()

    # Choice here between sequential and global change of direction
    if data.jumping == 'global':
      for i in range(len(data.vector)):
	sigmas[i]=(math.sqrt(1/eigv[i]/len(data.vector)))*rd.gauss(0,1)*2.4
    elif data.jumping == 'sequential':
      i = k%len(data.vector)
      sigmas[i] = (math.sqrt(1/eigv[i]))*rd.gauss(0,1)*2.4
    else:
      print '\n\n  Jumping method unknown (accepted : global, sequential)'

    vector_new=data.vector+np.dot(U,sigmas)
    i=0
    for value in data.params.itervalues():
      if(value[1]!=-1 and vector_new[i]<value[1]):
	flag+=1 # if a boundary value is reached, increment
      elif(value[2]!=-1 and vector_new[i]>value[2]):
	flag+=1 # same
      else:
	flag+=0 # else keep it at zero
      i+=1
  data._transmit_vector(vector_new)
  return vector_new

def jump(data,Direction,Range):
  # check if there is a multiplier for the number
  multiplier = 1.
  if Direction in data.Class_params.iterkeys():
    if len(data.Class_params[Direction])==5:
      multiplier = data.Class_params[Direction][4]

  data.Class_args[Direction]=Range*multiplier
  return data.Class_args


#---------------MCMC-CHAIN----------------------------------------------
def chain(_cosmo,data,command_line):
  # initialisation
  num_failure=2
  loglike=0
  failure=False
  sigma_eig,U=get_cov(data,command_line)
  failed=0

  # if restart wanted, change initial value for arguments
  if command_line.restart is not None:
    read_args_from_chain(data,command_line.restart)
    for i in range(len(data.Class)):
      data.Class_args = jump(data,data.Class_param_names[i],data.Class[i])
  failure,loglike=compute_lkl(_cosmo,data,data.Class_args)
  while ((failure is True) and (failed<=num_failure)):
    failed +=1
    data.vector=get_new_pos(data,sigma_eig,U,failed)
    for i in range(len(data.Class)):
      data.Class_args = jump(data,data.Class_param_names[i],data.Class[i])
    for i in range(len(data.Class),len(data.Class)+len(data.nuisance)):
      data.nuisance[i-len(data.Class)] = data.vector[i]
    failure,loglike=compute_lkl(_cosmo,data,data.Class_args)
  if failure is True:
    print ' /|\   Class tried {0} times to initialize with given parameters, and failed...'.format(num_failure+2)
    print '/_o_\  You might want to change your starting values, or pick default ones!'
    exit()

  max_loglike = loglike

  acc,rej=0.0,0.0	#acceptance and rejection number count
  io.print_parameters(sys.stdout,data)
  N=1			#number of steps
  for k in range(command_line.N):
    vector_new=get_new_pos(data,sigma_eig,U,failed+k)
    for i in range(len(data.Class)):
      newargs = jump(data,data.Class_param_names[i],vector_new[i])
    for i in range(len(data.Class),len(data.Class)+len(data.nuisance)):
      data.nuisance[i-len(data.Class)] = data.vector[i]
    failure,newloglike=compute_lkl(_cosmo,data,newargs)
    if(failure==True):
      failed += 1
      k-=2
      print 'Warning: Class failed due to choice of parameters, picking up new values'
      print newargs
      continue
    # Harmless trick to avoid exponentiating large numbers
    if newloglike >= loglike:
      alpha = 1.
    else:
      alpha=np.exp(newloglike-loglike)
    if ((alpha == 1.) or (rd.uniform(0,1) < alpha)): #accept step
      io.print_vector([data.out,sys.stdout],N,loglike,data)
      data.Class_args=newargs
      loglike=newloglike
      if loglike > max_loglike:
	max_loglike = loglike
      data.vector = vector_new
      acc+=1.0
      N=1
    else:
      rej+=1.0
      N+=1
    if acc % data.write_step ==0:
      io.refresh_file(data)
    if (failed >= num_failure):
      print ' /|\   The computation failed too many times, \n'
      print '/_o_\  Please check the values of your parameters'
  rate=acc/(acc+rej)
  print '#  {0} steps done, acceptance rate:'.format(command_line.N),rate
  if command_line.restart is not None:
    os.remove(command_line.restart)
    print '  deleting starting point of the chain {0}'.format(command_line.restart)
  return rate,max_loglike
