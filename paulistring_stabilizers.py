#basic automator seeking stabilizers in similarity transform of Pauli matrices
#via D. Gottesman,
#https://arxiv.org/abs/quant-ph/9807006

import numpy as np
import matplotlib.pyplot as plt

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

paulix = np.array([[0,1],[1,0]],dtype=complex)
pauliy = np.array([[0,-1j],[1j,0]],dtype=complex)
pauliz = np.array([[1,0],[0,-1]],dtype=complex)
identity = np.array([[1,0],[0,1]])

task_list = np.array([pauliz, paulix, pauliy, identity])

master_graph = []
master_iterator = [0,1,2,3]

for target_x in master_iterator: #this is the outer per-target basis
  outer_write = [0,1,2,3] #this cycles outer per-transform basis. Row-wise write goes here
  for transform_x in outer_write:
    row = []
    outer_interest = [0,1,2,3] #this cycles over the inner matrix of the target product
    for target_y in outer_interest:
      inner_transform = [0,1,2,3] #this cycles over the inner matrix of the transform
      for transform_y in inner_transform:
        target = np.einsum('ij,kl->ijkl',task_list[target_x,:],task_list[target_y,:])
        row0 = np.asarray([target[0,0,0],target[0,1,0]]).flatten()
        row1 = np.asarray([target[0,0,1],target[0,1,1]]).flatten()
        row2 = np.asarray([target[1,0,0],target[1,1,0]]).flatten()
        row3 = np.asarray([target[1,0,1],target[1,1,1]]).flatten()
        target = np.asarray([row0,row1,row2,row3])
        transform = np.einsum('ij,kl->ijkl',task_list[transform_x,:],task_list[transform_y,:])
        row0 = np.asarray([transform[0,0,0],transform[0,1,0]]).flatten()
        row1 = np.asarray([transform[0,0,1],transform[0,1,1]]).flatten()
        row2 = np.asarray([transform[1,0,0],transform[1,1,0]]).flatten()
        row3 = np.asarray([transform[1,0,1],transform[1,1,1]]).flatten()
        transform = np.asarray([row0,row1,row2,row3])
        dump = np.einsum('ij,jk->ik',transform,target)
        final = np.einsum('ij,jk->ik',dump,transform)
        if np.all(final == target):
          row.append(1)
        if np.all(final == -1.0 * target):
          row.append(-1)
    master_graph.append(row)

master_graph = np.asarray(master_graph)
data = master_graph
labels = ['Z', 'X', 'Y', 'I'] * 4

fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(data, cmap='RdBu', vmin=-1, vmax=1)


ax.set_xticks(range(16), labels=labels)
ax.set_yticks(range(16), labels=labels)
outer_pos = [1.5, 5.5, 9.5, 13.5]
outer_labels = ['Z', 'X', 'Y', 'I']
ax.xaxis.tick_top()

# Secondary X (Bottom)
secax_x = ax.secondary_xaxis('top', functions=(lambda x: x, lambda x: x))
secax_x.set_xticks(outer_pos, labels=outer_labels)
secax_x.tick_params(axis='x', which='both', length=0, pad=20) # Push labels down
# Secondary Y (Left)
secax_y = ax.secondary_yaxis('left', functions=(lambda y: y, lambda y: y))
secax_y.set_yticks(outer_pos, labels=outer_labels)
secax_y.tick_params(axis='y', which='both', length=0, pad=25) # Push labels left

ax.set_xlabel("Lower Qubit (Target Operator by/ Similarity Transform)", fontsize=12, labelpad=35)
ax.xaxis.set_label_position('top')
ax.set_ylabel("Upper Qubit (Target Operator by/ Similarity Transform)", fontsize=12, labelpad=40)

plt.tight_layout()
plt.show()

task_list = np.array([pauliz, paulix, pauliy, identity])
strings = ['PauliZ', 'PauliX', 'PauliY', 'Identity']
hh = np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])
cnot12 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
master_iterator = [0,1,2,3]

#HH and CNOT if you're curious
for n in master_iterator:
  outer_write = [0,1,2,3]
  for y in outer_write:
    target = np.einsum('ij,kl->ijkl',task_list[n,:],task_list[y,:])
    row0 = np.asarray([target[0,0,0],target[0,1,0]]).flatten()
    row1 = np.asarray([target[0,0,1],target[0,1,1]]).flatten()
    row2 = np.asarray([target[1,0,0],target[1,1,0]]).flatten()
    row3 = np.asarray([target[1,0,1],target[1,1,1]]).flatten()
    target = np.asarray([row0,row1,row2,row3])
    dump = np.einsum('ij,jk->ik',hh,target)
    final = np.einsum('ij,jk->ik',dump,hh)
    final = final/4
    if np.all(final == target):
      print('HH is +1 stabilizer to ' + strings[n] + ' otimes ' + strings[y])
    if np.all(final == -1.0 * target):
      print('HH is -1 stabilizer to ' + strings[n] + ' otimes ' + strings[y])
    dump = np.einsum('ij,jk->ik',cnot12,target)
    final = np.einsum('ij,jk->ik',dump,cnot12)
    if np.all(final == target):
      print('CNOT12 is +1 stabilizer to ' + strings[n] + ' otimes ' + strings[y])
    if np.all(final == -1.0 * target):
      print('CNOT12 is -1 stabilizer to ' + strings[n] + ' otimes ' + strings[y])
    dump = np.einsum('ij,jk->ik',target,hh)
    final = np.einsum('ij,jk->ik',dump,target)
    if np.all(final == hh): #the sqrt2 **-1 scalar seperates out
      print(strings[n] + ' otimes ' + strings[y] + ' is a +1 stabilizer to HH')
    if np.all(final == -1.0 * hh): #the sqrt2 **-1 scalar seperates out
      print(strings[n] + ' otimes ' + strings[y] + ' is a -1 stabilizer to HH')
    dump = np.einsum('ij,jk->ik',target,cnot12)
    final = np.einsum('ij,jk->ik',dump,target)
    if np.all(final == cnot12):
      print(strings[n] + ' otimes ' + strings[y] + ' is a +1 stabilizer to CNOT12')
    if np.all(final == -1.0 * cnot12):
      print(strings[n] + ' otimes ' + strings[y] + ' is a +1 stabilizer to CNOT12')
