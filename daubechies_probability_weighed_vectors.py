import numpy as np
from scipy import stats as stat
import pywt

def db4errorvector(wavepacket,probability):
  """ This divides by the number of std expected to have made up the geometric factor scaling coefficients to fourier mapping.
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
  flowpast = filter coefficients for the remainer
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

signal = [0.0,0.5,0.9,0.11,0,0,0.03,0.2] #sample signal vector
p = 0.41 #sample probability
scalar,vector = db4fouriererrorvector(signal,p)
coef,propagator = db4fourierskimmer(signal,p)
print(str(scalar))
print(str(vector))
print('Modulus Filter coefficients = '+ str(coef))
print('Modulus = ' + str(propagator))

#generates wavelet filter coefficients for db4 in fourier mapping
c_0 = 0
c_1 = 0
c_2 = 0.0
c_3 = (np.float64(1.0) + np.sqrt(np.float64(3.0)))/(np.float64(4.0)*(np.sqrt(np.float64(2.0))))
#occupy in reverse order

known_amplitudes = 1 #how many coefficients do you already know the answer to.

c_0 = np.float64(c_0)
c_1 = np.float64(c_1)
c_2 = np.float64(c_2)
c_3 = np.float64(c_3)

if known_amplitudes == 0:
  print("Standard solution. This one's free chump!")
  c_0 = (np.float64(1.0) - np.sqrt(np.float64(3.0)))/(np.float64(4.0)* np.sqrt(np.float64(2.0)))
  c_1 = (np.float64(3.0) - np.sqrt(np.float64(3.0)))/(np.float(4.0)* np.sqrt(np.float64(2.0)))
  c_2 = (np.float64(3.0) + np.sqrt(np.float64(3.0)))/(np.float64(4.0)* np.sqrt(np.float64(2.0)))
  c_3 = (np.float64(1.0) + np.sqrt(np.float64(3.0)))/(np.float64(4.0)*(np.sqrt(np.float64(2.0))))
  finalkey = np.array([c_0, c_1, c_2, c_3])
  print(str(finalkey))

if known_amplitudes == 3:
  if c_0 == np.float64(0.0):
    answer = np.float64(-1.0)* c_1 * c_3 /(c_2)
    finalkey = np.array([answer, c_1, c_2, c_3])
    print(str(finalkey))

if known_amplitudes == 2:
  c_1 = np.float64(2.0) * c_2 - np.float64(3.0) * c_3
  c_0 = c_1 - c_2 + c_3
  finalkey = np.array([c_0, c_1, c_2, c_3])

if known_amplitudes == 1:
  c_2 = np.sqrt(np.float64(3.0)) * c_3
  c_1 = (np.float64(2.0) * np.sqrt(np.float64(3.0)) - np.float64(3.0)) * c_3
  c_0 = c_1 - c_2 + c_3
  N = np.sqrt(np.float64(1.0)/((c_0**2) + (c_1**2) + (c_2**2) + (c_3**2)))
  c_0 = c_0 / N
  c_1 = c_1 / N
  c_2 = c_2 / N
  c_3 = c_3 / N
  finalkey = np.array([c_0, c_1, c_2, c_3])
  print(finalkey)
  check = finalkey[0] **2 + finalkey[3]**2 + finalkey[1]**2 + finalkey[2]**2
  print(str(check) + ' should be =1')
