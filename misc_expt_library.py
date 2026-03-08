import pennylane as qml
import numpy as np

def arbitrary1qubit_Rot_from_U(unitary):
  """
  Applies Rot gate after finding its arguments. Aids in development by giving a maximum-sized decomposition RzRYRz to a choice unitary, use equivalent QubitUnitary function otherwise
  Make sure you rewrite the wires argument to work on the desired level as it has a placeholder value.
  u11 * u22 = cos**2(theta/2) meaning cos(theta) = cos**2(theta/2) - sin**2(theta/2)
  u21 * u22 = e**-2omega(i)/2 cos(theta/2)sin(theta/2)

  Arguments:
  unitary = single qubit unitary in array of rows [[a,b],[c,d]], dtype = complex

  Returns
  qml.state() = state vector
  phi,theta,omega = in order fed to Rot, so they can be reused.
  """
  unitary = np.array([[0,-1j],[1j,0]])
  theta = np.arccos(2.0 * unitary[0,0].item() * unitary[1,1].item() - 1.0)
  omega = np.arcsin(np.imag(unitary[1,0].item() * unitary[1,1].item() / (np.sin(0.5*theta)*np.cos(0.5*theta))))
  phi = np.arcsin(np.imag(-1.0 * unitary[0,1].item() * unitary[1,1].item() / (np.sin(0.5*theta)*np.cos(0.5*theta))))
  qml.Rot(phi,theta,omega, wires=0)

  return qml.state(), phi, theta, omega

def ibm2001_nmr_expt():
  '''
  The IBM 2001 NMR quantum computation paper uploaded to Pennylane. Initialize with wires = 7
  Inputs:
  None
  
  Outputs:
  qml.Probs() = probabilities of each qubit state in this seven-wire system.
  '''

  U1 = np.array([[1,0,0,0],[0,1,0,0],[0,0,1/np.sqrt(2.0)-(1/np.sqrt(2.0))*1j,0],[0,0,0,1/np.sqrt(2.0)+(1/np.sqrt(2.0))*1j]],dtype=complex)
  U2 = np.array([[1,0,0,0],[0,1,0,0],[0,0,np.cos(np.pi/8.0)-np.sin(np.pi/8.0)*1j,0],
                [0,0,0,np.cos(np.pi/8.0)+np.sin(np.pi/8.0)*1j]],dtype=complex)


  qml.Hadamard(wires=0)
  qml.Hadamard(wires=1)
  qml.Hadamard(wires=2)
  qml.CNOT(wires=[2,4])
  qml.CNOT(wires=[2,5])
  qml.CNOT(wires=[3,5])
  qml.Toffoli(wires=[1,5,3])
  qml.CNOT(wires=[3,5])
  qml.CNOT(wires=[6,4])
  qml.Toffoli(wires=[1,4,6])
  qml.CNOT(wires=[6,4])
  qml.Hadamard(wires=0)
  qml.QubitUnitary(U1,wires=[0,1])
  qml.Hadamard(wires=1)
  qml.SWAP(wires=[1,2])
  qml.QubitUnitary(U1,wires=[0,1])
  qml.QubitUnitary(U2,wires=[0,1])
  qml.SWAP(wires=[1,2])
  qml.Hadamard(wires=2)
  return qml.probs()
