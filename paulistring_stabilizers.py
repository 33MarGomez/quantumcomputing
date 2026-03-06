#basic automator seeking stabilizers in similarity transform of Pauli matrices
#via D. Gottesman,
#https://arxiv.org/abs/quant-ph/9807006

import numpy as np

strings = ['PauliX','PauliY','PauliZ']
paulix = np.array([[0,1],[1,0]],dtype=complex)
pauliy = np.array([[0,-1j],[1j,0]],dtype=complex)
pauliz = np.array([[1,0],[0,-1]],dtype=complex)
identity = np.array([[1,0],[0,1]])

task_list = np.array([paulix,pauliy,pauliz])
indexes = [0,1,2]
print('Format = Transformation, Target, eigenvalue')

for i in indexes:
  i_2 = i - 1
  i_3 = i - 2
  dump = np.einsum('ij,jk->ik',task_list[i_2,:],task_list[i,:])
  final = np.einsum('ij,jk->ik',dump,task_list[i_2,:])
  if np.all(final == task_list[i,:]):
    print(strings[i_2] + ', ' + strings[i] + ', +1')
  if np.all(final == -1.0 * task_list[i,:]):
    print(strings[i_2] + ', ' + strings[i] + ', -1')
  dump = np.einsum('ij,jk->ik',task_list[i_3,:],task_list[i,:])
  final = np.einsum('ij,jk->ik',dump,task_list[i_3,:])
  if np.all(final == task_list[i,:]):
    print(strings[i_3] + ', ' + strings[i] + ', +1')
  if np.all(final == -1.0 * task_list[i,:]):
    print(strings[i_3] + ', ' + strings[i] + ', -1')

dump = np.einsum('ij,kl->ijkl',identity,pauliy)
row0 = np.asarray([dump[0,0,0],dump[0,1,0]]).flatten()
row1 = np.asarray([dump[0,0,1],dump[0,1,1]]).flatten()
row2 = np.asarray([dump[1,0,0],dump[1,1,0]]).flatten()
row3 = np.asarray([dump[1,0,1],dump[1,1,1]]).flatten()
final = np.asarray([row0,row1,row2,row3])
print(final)

strings = ['Identity', 'PauliX','PauliY','PauliZ']
paulix = np.array([[0,1],[1,0]],dtype=complex)
pauliy = np.array([[0,-1j],[1j,0]],dtype=complex)
pauliz = np.array([[1,0],[0,-1]],dtype=complex)
identity = np.array([[1,0],[0,1]])

task_list = np.array([identity,paulix,pauliy,pauliz])
indexes = [0,1,2,3]
for i in indexes:
  otimes = [i]
  otimes.extend([i-1,i-2,i-3])
  j = [0,1,2,3]
  for r in otimes:
    target = np.einsum('ij,kl->ijkl',task_list[i,:],task_list[r,:])
    row0 = np.asarray([target[0,0,0],target[0,1,0]]).flatten()
    row1 = np.asarray([target[0,0,1],target[0,1,1]]).flatten()
    row2 = np.asarray([target[1,0,0],target[1,1,0]]).flatten()
    row3 = np.asarray([target[1,0,1],target[1,1,1]]).flatten()
    target = np.asarray([row0,row1,row2,row3])
    for e in j:
      q = [e, e-1, e-2, e-3]
      for t in q:
        transform = np.einsum('ij,kl->ijkl',task_list[e,:],task_list[t,:])
        row0 = np.asarray([transform[0,0,0],transform[0,1,0]]).flatten()
        row1 = np.asarray([transform[0,0,1],transform[0,1,1]]).flatten()
        row2 = np.asarray([transform[1,0,0],transform[1,1,0]]).flatten()
        row3 = np.asarray([transform[1,0,1],transform[1,1,1]]).flatten()
        transform = np.asarray([row0,row1,row2,row3])
        dump = np.einsum('ij,jk->ik',transform,target)
        final = np.einsum('ij,jk->ik',dump,transform)
        if np.all(final == target):
          print(strings[e] + ' otimes ' + strings[t] + ', ' + strings[i] + ' otimes ' + strings[r] + ', +1')
        if np.all(final == -1.0 * target):
          print(strings[e] + ' otimes ' + strings[t] + ', ' + strings[i] + ' otimes ' + strings[r] + ', -1')
