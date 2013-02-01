import numpy as np
import sys

mat1_name = sys.argv[-2]
mat2_name = sys.argv[-1]

file1 = open(mat1_name,'r')
file2 = open(mat2_name,'r')

i=0
for line in file1:
  if line.find('#')!=-1:
    list1 = line.strip('#').replace(' ','').replace('\n','').split(',')
    mat1  = np.zeros((len(list1),len(list1)),'float64')
  else:
    line=line.split()
    for j in range(len(line)):
      mat1[i][j]=np.array(line[j],'float64')
    i+=1

i=0
for line in file2:
  if line.find('#')!=-1:
    list2 = line.strip('#').replace(' ','').replace('\n','').split(',')
    mat2  = np.zeros((len(list2),len(list2)),'float64')
  else:
    line=line.split()
    for j in range(len(line)):
      mat2[i][j]=np.array(line[j],'float64')
    i+=1

print
print 'Here are the two matrices to merge: ',mat1_name,mat2_name
print 'I will keep cosmology from ',mat1_name, 'and nuisance from ',mat2_name,'. If you want the contrary, change the order of calling'
raw_input('press enter to go on')
print 
print 'List of parameters from ',mat1_name
for elem in list1:
  i = list1.index(elem)
  print i+1,elem
print 
print 'I will cut the N first parameters of ',mat1_name
print '\n'
N = int(raw_input('Specify N: '))
print 'I will then include the matrix elements of the first matrix up to ',list1[N-1]
print 'List of parameters from ',mat2_name
for elem in list2:
  i = list2.index(elem)
  print i+1,elem
print 
print 'I will append the ending parameters of ',mat2_name,' starting from the Mth one'
print '\n'
M = int(raw_input('Specify M: '))
print 'I will then include the matrix elements of the second matrix from ',list2[M-1]

submat1 = np.zeros((N,N),'float64')
submat2 = np.zeros((len(list2)-M+1,len(list2)-M+1),'float64')

submat1 = mat1[:N,:N]
submat2 = mat2[M-1:,M-1:]

final_mat = np.zeros((N+len(list2)-M+1,N+len(list2)-M+1),'float64')
final_mat[:N,:N] = submat1
final_mat[N:,N:] = submat2

final_list = list1[:N]+list2[M-1:]

print final_list
result = raw_input('writing the resulting covariance matrix to\n')

output = open(result,'w')
output.write('# ')
for i in range(len(final_list)):
  string = final_list[i]
  if i != len(final_list)-1:
    string+=','
  output.write('%-16s' % string )
output.write('\n')
for i in range(len(final_list)):
  for j in range(len(final_list)):
    if final_mat[i][j]>=0:
      output.write(' %.5e\t' % final_mat[i][j])
    else:
      output.write('%.5e\t' % final_mat[i][j])
  output.write('\n')

print 'done !'
