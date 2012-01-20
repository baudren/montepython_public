import os
import re
import sys
import mcmc
import random as rd
import numpy as np
from collections import OrderedDict as od
from datetime import date


def log_parameters(data,command_line):
  log     = open(command_line.folder+'/log.param','w')
  param_file = open(command_line.par,'r')
  log.write("#-----Class {0} (subversion {1})-----\n\n".format(data.version,data.subversion))
  for line in param_file:
    log.write(line)
  param_file.close()
  log.close()

def log_Class_arguments(data,command_line):
  if len(data.Class_arguments) >= 1:
    log     = open(command_line.folder+'/log.param','a')
    log.write('\n\n#-----------Class-arguments---------\n')
    log.write('data.Class_arguments.update({0})\n'.format(data.Class_arguments))
    log.close()

def print_parameters(out,data):
  param = data.get_mcmc_parameters(['varying'])
  out.write('#  -LogLkl\t')
  for i in range(len(param)):
    if len(data.mcmc_parameters[param[i]]['initial'])==5:
      number = 1./(data.mcmc_parameters[param[i]]['initial' ][4])
      if number < 1000:
	out.write('%0.d%s\t' % (number,param[i]))
      else:
	out.write('%0.e%s\t' % (number,param[i]))
    else:
      out.write('{0}\t'.format(param[i]))
  out.write('\n')

def print_vector(out,N,loglkl,data):
  for j in range(len(out)):
    out[j].write('%d  %.3f\t' % (N,-loglkl))
    for elem in data.get_mcmc_parameters(['varying']):
      out[j].write('%.6f\t' % data.mcmc_parameters[elem]['last_accepted'])
    out[j].write('\n')

def refresh_file(data):
  data.out.close()
  data.out=open(data.out_name,'a')

def create_output_files(command_line,data):
  if command_line.restart is None:
    number = command_line.N
  else:
    number = int(command_line.restart.split('/')[-1].split('__')[0].split('_')[1]) + command_line.N

  # output file
  outname_base='{0}_{1}__'.format(date.today(),number)
  suffix=0
  for files in os.listdir(command_line.folder):
    if files.find(outname_base)!=-1:
      if int(files.split('__')[-1])>suffix:
        suffix=int(files.split('__')[-1])
  suffix+=1
  data.out=open(command_line.folder+outname_base+str(suffix),'w')
  print '  Creating {0}{1}{2}'.format(command_line.folder,outname_base,suffix)
  data.out_name='{0}{1}{2}'.format(command_line.folder,outname_base,suffix)
  # in case of a restart, copying the whole thing in the new file
  if command_line.restart is not None:
    for line in open(command_line.restart,'r'):
      data.out.write(line)

def get_tex_name(name,number=0):
  tex_greek = ['omega','tau','alpha','beta','delta','nu','Omega']
  for elem in tex_greek:
    if elem in name:
      name="""\\"""+name
  if number==0: 
    if name.find('_')!=-1:
      temp_name = name.split('_')[0]+'_{'
      for i in range(1,len(name.split('_'))):
	temp_name += name.split('_')[i]
      temp_name += '}'
      name = temp_name
    name = "${0}$".format(name)
    return name
  elif number < 1000:
    name = "$%0.d~%s$" % (number,name)
  else:
    temp_name = "$%0.e%s$" % (number,name)
    m = re.search(r'(?:\$[0-9]*e\+[0]*)([0-9]*)(.*)',temp_name)
    name = '$10^{'+m.groups()[0]+'}'+m.groups()[1]
  return name

# New class
class File(file):
  def head(self, lines_2find=1):
    self.seek(0)                            #Rewind file
    return [self.next() for x in xrange(lines_2find)]

  def tail(self, lines_2find=1):  
    self.seek(0, 2)                         #go to end of file
    bytes_in_file = self.tell()             
    lines_found, total_bytes_scanned = 0, 0
    while (lines_2find+1 > lines_found and
       bytes_in_file > total_bytes_scanned): 
      byte_block = min(1024, bytes_in_file-total_bytes_scanned)
      self.seek(-(byte_block+total_bytes_scanned), 2)
      total_bytes_scanned += byte_block
      lines_found += self.read(1024).count('\n')
    self.seek(-total_bytes_scanned, 2)
    line_list = list(self.readlines())
    return line_list[-lines_2find:]

