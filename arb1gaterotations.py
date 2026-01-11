import pennylane as qml
import numpy as np

unitary = np.array([[0,-1j],[1j,0]])

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
  theta = np.arccos(2.0 * unitary[0,0].item() * unitary[1,1].item() - 1.0)
  omega = np.arcsin(np.imag(unitary[1,0].item() * unitary[1,1].item() / (np.sin(0.5*theta)*np.cos(0.5*theta))))
  phi = np.arcsin(np.imag(-1.0 * unitary[0,1].item() * unitary[1,1].item() / (np.sin(0.5*theta)*np.cos(0.5*theta))))
  qml.Rot(phi,theta,omega, wires=0)

  return qml.state(), phi, theta, omega
