import numpy as np
from scipy import stats as stat
import pywt

def db4errorvector(wavepacket,probability):
  """ 
  This divides by the number of std expected to have made up the geometric factor scaling coefficients to fourier mapping.
  """
  z_score = stat.norm.ppf(probability,loc = 0, scale = 1)
  scalar = (np.float64(1.0) + np.sqrt(np.float64(3.0)))/(np.float64(4.0)*(np.sqrt(np.float64(2.0)))) / z_score
  unscaledvector = pywt.dwt(data = wavepacket,wavelet = 'db4')
  scaledvector = np.array(unscaledvector)
  scaledvector = scalar * scaledvector

  return scalar, scaledvector

def db4fourierskimmer(signal,probability):
  """
  This takes the confidence of qft adherence off the top and leaves behind deviant behaviour.
  Add them back together for the original, fourier-amplified signal.
  If d scales the fourth entry to the expected coefficient, and z is number of standard deviations causing the signal, and q is fourier conforming
  dz*(signal)-q / (dz) = signal - (1/dz)*q

  variables:
  signal = signal vector
  probability = 0 < p < 1 of observing the signal. Used for n_std.

  returns:
  flowpast = filter coefficients for the modulus
  delta = surviving vector after losing fourier components

  """
  z_score = stat.norm.ppf(probability,loc = 0, scale = 1)
  scalar =  (1.0 + np.sqrt(3.0))/(4.0*(np.sqrt(2.0))) / signal[3]
  element = np.array([[(1.0 - np.sqrt(3.0))/(4.0 * np.sqrt(2.0))],
                          [(3.0 - np.sqrt(3.0))/(4.0* np.sqrt(2.0))],
                          [(3.0 + np.sqrt(3.0))/(4.0* np.sqrt(2.0))],
                          [(1.0 + np.sqrt(3.0))/(4.0*(np.sqrt(2.0)))]
                          ]
                         )
  dim = np.size(signal)
  subtraction = np.zeros(shape=(dim,dim))
  z = 0
  while z < dim-2:
    subtraction[z,z] = element[0].item()
    subtraction[z+1,z+3] = -1.0 * element[0].item()
    subtraction[z,z+2] = element[2].item()
    subtraction[z+1,z+1] = -1.0 * element[2].item()
    subtraction[z,z+3] = element[3].item()
    subtraction[z+1,z] = element[3].item()
    subtraction[z,z+1] = element[1].item()
    subtraction[z+1,z+2] = element[1].item()
    z += 2
  subtraction[-2,-2] = element[0].item()
  subtraction[-1,1] = -1.0 * element[0].item()
  subtraction[-2,0] = element[2].item()
  subtraction[-1,-1] = -1.0 * element[2].item()
  subtraction[-2,1] = element[3].item()
  subtraction[-1,-2] = element[3].item()
  subtraction[-1,0] = element[1].item()
  subtraction[-2,-1] = element[1].item()

  delta = signal -  np.matmul(subtraction,signal) * ((z_score*scalar)**(-1.0))
  flowpast = pywt.dwt(data = delta,wavelet = 'db4')

  return flowpast, delta

def coefficientconstructor(c_3):
  """
  Here is the derivation for the fourier wavelets in db4. On Wikipedia, and across 2 different books, the stated solution does not respect the symmetry
  of the fourier transform but does satisfy all the requirements because it is only incorrect when substituted at an intermediate step:

  (2*sqrt(3) -3)(1-sqrt(3) = 5 sqrt(3) - 15 NOT 3 + sqrt(3). It does work if you flip all the signs on the irrational, c_3 = 1 + sqrt(3) & c_1 = 3 - sqrt(3)

  to prove one of the moments is wrong. Here is the corrected coefficient order, despite being cyclical, order does matter here.
  This is just an efficient note-taking technique, this function essentially does nothing. The derivation works based on =0 not changing lin. comb.
  """
  #c_2 = c_0 * c_2 + c_1 * c_3 - c_2 * (c_0 -1.0 * c_1 + c_2 - c_3) + c_2 * (-1.0 * c_1 + 2* c_2 - 3 * c_3) + c_3 * (-1.0 * c_1 + 2 * c_2 - 3 * c_3)
  c_2 = np.sqrt(3.0)j * c_3
  c_1 = (2.0 * np.sqrt(3.0) - 3.0) * c_3
  c_0 = np.sqrt(1.0 - c_1**2 - c_2**2 - c_3**2)
  if c_0 - c_1 + c_2 - c_3 == 0:
    print('Zeroeth moment is true')
    return c_2
  if 2.0 * c_2 - c_1 - 3.0 * c_3 == 0:
    print('First moment is true')
    return c_1
  if c_0 == -1.0 * c_1 * c_3 / c_2:
    print('Inverse matrix is transposition')
    return c_0

signal = [0.0,0.5,0.9,0.11,0,0,0.03,0.2] #sample signal vector
p = 0.41 #sample probability
scalar,vector = db4fouriererrorvector(signal,p)
coef,propagator = db4fourierskimmer(signal,p)
print(str(scalar))
print(str(vector))
print('Modulus Filter coefficients = '+ str(coef))
print('Modulus = ' + str(propagator))
