pip install pennylane --upgrade # run this in its own cell, colab

import pennylane as qml

dev0 = qml.device('default.qubit', wires=5)

def brachialplexusfunction():
  '''
  The measurements are not underbaked, there's a spectroscopic reason no co-dependent statistics are used;
  it is population intensities, why take a weaker signal if the symmetries are conserved in the body?
  '''
  qml.Hadamard(wires=0)
  qml.Hadamard(wires=4)
  qml.CNOT(wires=[0,1])
  qml.CNOT(wires=[4,3])
  qml.SWAP(wires=[1,2])
  qml.Toffoli(wires=[1,2,3])
  qml.Toffoli(wires=[0,4,1])
  return qml.expval(qml.PauliZ(wires=1)), qml.expval(qml.PauliZ(wires=2)),qml.expval(qml.PauliZ(wires=3))

qnode_0 = qml.QNode(brachialplexusfunction, dev0)
print('The Brachial Nexus as it is:')
print(qml.draw(qnode_0)())

truth_guide = ['wire 0 = 0','w0=1','w1=0','w1=1','w2=0','w2=1','w3=0','w3=1','w4=0','w4=1'] #[None]
inputs_guide = ['C5','C6','C7','C8','T1']
outputs_guide = ['C5-C6_superposition','Median','Lateral_anterior','Post. Cord','C8-T1_superposition']
median_is_wire1= [None, None, 0,'w0=w4=1',
                  'see w1','see w1, but mentally invert', None, None, None, None]
lateral_ant_is_wire2 = ['see w2 = 0', 'see w2 = 1', 'check on w0', 'check on w0 !but mentally invert!',
                        0, 1, None, None, None, None]
Posterior_cord_is_wire3 = [None, None,'nothing','see if w2 is also 1','nothing','only flip if w1 is also 1',
                         'check on w1,w2','check on w1,w2','do nothing','is now flipped']

scratchworks_guide = ['starting state', 'chronological action 1 explanation', 'state update', 'chronological explanation 2', 'accompanying state']
scratch_work_w3 = [0,'T1 is 1', 1, 'w1 and w2 both 1', 0, 'read']
scratch_work_w3 = [0,'T1 is 1', 1, 'Either of w1, w2 has a zero or both', 1, 'read']
scratch_work_master_canonical_state = [0, 'C5 is 1', 1, 'if Median w1 is also 1, flip Posterior Cord w3', 1, 'read w2']
scratch_work_master_canonical_state = [1, 'C5 is 1', 0, 'Posterior Cord w3 guaranteed not to flip', 0, 'read w2']
scratch_work_master_canonical_state = [0, 'C5 flipped the upper trunk and T1 flipped the lower trunk', 1, 'read w1']
scratch_work_master_canonical_state = [0, 'C5 or T1 flipped one or no measurements between them but not both' , 0, 'read w1']
scratch_work_master_canonical_state = [0, 'irradiate', 1, 'C5 is 0', 1,
                                       'read w2. If w3 flipped twice, w1 =1, otherwise w1 =0']
scratch_work_master_canonical_state = [0, 'do not irradiate', 0, 'C5 is 1', 1, 'read w2']
so_what_is_c7 =  ['interfere w4, interfere w0, later in the computer and run the experiment multiple times to converge the answer, lots of tracer',
                  'Based on C6 irradiation state, and w3 the exact w0 and w4 states and carry on computing with them.']

dev1 = qml.device('default.qubit', wires=5)

def brachialplexus_convoludenexus():
  '''
  If the joining is better represented by state read, a CNOT exists. This means the brachial plexus crosses over itself to mess with med students.
  '''
  qml.Hadamard(wires=0)
  qml.Hadamard(wires=4)
  qml.CNOT(wires=[0,1])
  qml.CNOT(wires=[4,3])
  qml.Toffoli(wires=[1,2,3])
  qml.CNOT(wires=[1,2])
  qml.Toffoli(wires=[0,4,2])
  return qml.expval(qml.PauliZ(wires=1)), qml.expval(qml.PauliZ(wires=2)),qml.expval(qml.PauliZ(wires=3))

qnode_1 = qml.QNode(brachialplexus_convoludenexus, dev1)
print('An arbitrary Brachial Nexus cross-over is just med-school torture:')
print(qml.draw(qnode_1)())
