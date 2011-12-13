from classy import Class
import sys
import math
import random as rd
import numpy as np
import io
import data

# COMPUTE LKL
def compute_lkl(_cosmo,Data):
  # Prepare the cosmological module with the new set of parameters
  _cosmo.set(Data.args)

  # Compute the model
  failure=False
  try:
    _cosmo._compute(["lensing"])
  except:
    return True,0 

  # For each desired likelihood, compute its value against the theoretical model
  loglike=0
  for likelihood in Data.lkl.itervalues():
    try:
      loglike+=likelihood._loglkl(_cosmo,Data)
    except:
      print '\n  Failure in likelihood {0}'.format(likelihood)
      return True,0
  
  # Clean the cosmological strucutre
  _cosmo._struct_cleanup(set(["lensing","nonlinear","spectra","primordial","transfer","perturb","thermodynamics","background","bessel"]))

  return failure,loglike

def read_args_from_chain(data,chain):
  Chain = io.File(chain,'r')
  i = 1
  for elem in data.param_names:
    data.args[elem] = float(Chain.tail(1)[0].split('\t')[i])
    data.theta[i-1] = data.args[elem]
    i+=1

def get_cov(data,command_line):
  M=np.identity(len(data.theta),'float64')
  np.set_printoptions(precision=2,linewidth=150)
  rot=np.zeros((len(data.theta),len(data.theta)))
  i=0
  # if the user wants to use a covmat file
  if command_line.cov is not None:
    cov=open('{0}'.format(command_line.cov),'r')
    for line in cov:
      if line.find('#')!=-1:
	covnames = line.strip('#[').strip(']\n').replace("'","").split(', ')
	for k in range(len(data.theta)):
	  if covnames[k] not in data.param_names:
	    print ' /|\  Error, you do not have the same parameters\n/_o_\ in your cov matrix and chain!' 
	    print data.param_names,covnames
	    exit()
	  for h in range(len(data.theta)):
	    if covnames[k]==data.param_names[h]:
	      rot[k][h]=1
	    else:
	      rot[k][h]=0
      else:
	line=line.split()
	for j in range(len(line)):
	  M[i][j]=np.array(line[j],'float64')
	i+=1
    #Apply rotation
    M=np.dot(rot,np.dot(M,rot))

  # else, take sigmas^2.
  else:
    for elem in data.params:
      M[i][i]=np.array(data.params[elem][3],'float64')**2
      i+=1

  #inverse, and diagonalization
  eigv,eigV=np.linalg.eig(np.linalg.inv(M))
  return eigv,eigV

def get_new_pos(data,eigv,U):
  theta_new=np.zeros(len(data.theta),'float64')
  sigmas=np.zeros(len(data.theta),'float64')
  flag=1
  while flag!=0:
    flag=0 # initialize: there are no problems at first
    rd.seed()
    for i in range(len(data.theta)):
      sigmas[i]=(math.sqrt(1/eigv[i]))*rd.gauss(0,1)/2.4
    theta_new=data.theta+np.dot(U,sigmas)
    i=0
    for value in data.params.itervalues():
      if(value[1]!=-1 and theta_new[i]<value[1]):
	flag+=1 # if a boundary value is reached, increment
      elif(value[2]!=-1 and theta_new[i]>value[2]):
	flag+=1 # same
      else:
	flag+=0 # else keep it at zero
    if flag!=0:
      print 'one turn for nothing'
  return theta_new

def jump(data,Direction,Range):
  data.args[Direction]=Range
  if Direction.find('A_s')!=-1:
    data.args[Direction]=Range*1e-9


#---------------MCMC-CHAIN----------------------------------------------
def chain(_cosmo,data,command_line,out):
  # initialisation
  loglike=0
  failure=False
  sigma_eig,U=get_cov(data,command_line)
  failed=0

  # if restart wanted, change initial value for arguments
  if command_line.restart is not None:
    read_args_from_chain(data,command_line.restart)
    for i in range(len(data.theta)):
      jump(data,data.param_names[i],data.theta[i])

  failure,loglike=compute_lkl(_cosmo,data)
  while ((failure is True) and (failed<=98)):
    failed +=1
    data.theta=get_new_pos(data,sigma_eig,U)
    for i in range(len(data.theta)):
      jump(data,data.param_names[i],data.theta[i])
    failure,loglike=compute_lkl(_cosmo,data)
  if failure is True:
    print ' /|\   Class tried 100 times to initialize with given parameters, and failed...'
    print '/_o_\  You might want to change your starting values, or pick default ones!'
    exit()

  min_loglike = loglike

  acc,rej=0.0,0.0	#acceptance and rejection number count
  io.print_parameters(sys.stdout,data.param_names)
  N=1			#number of steps
  for k in range(command_line.N):
    theta_new=get_new_pos(data,sigma_eig,U)
    for i in range(len(theta_new)):
      newargs=jump(data,data.param_names[i],theta_new[i])
    failure,newloglike=compute_lkl(_cosmo,data)
    if(failure==True):
      failed += 1
      k-=2
      print 'Warning: Class failed due to choice of parameters, picking up new values'
      print newargs
      continue
    alpha=np.exp(newloglike-loglike)
    if ((alpha>1) or (rd.uniform(0,1) < alpha)): #accept step
      io.print_theta([out,sys.stdout],N,loglike,data)
      args=newargs
      loglike=newloglike
      if loglike > min_loglike:
	min_loglike = loglike
      data.theta=theta_new
      acc+=1.0
      N=1
    else:
      rej+=1.0
      N+=1
    if acc % data.write_step ==0:
      out = io.refresh_file(out,data.out_name)
    if (failed >= 98):
      print ' /|\   The computation failed too many times, \n'
      print '/_o_\  Please check the values of your parameters'
  if N>1:
    io.print_theta([out,sys.stdout],N-1,loglike,data)
  rate=acc/(acc+rej)
  print '#{0} steps done, acceptance rate:'.format(command_line.N),rate
  return rate,min_loglike

def clik_loglikelihood(_clik,cl):
  loglkl=_clik(cl)
  return loglkl
