pip install pennylane --upgrade #this must run in its own Jupyter (colab) cell

import pennylane as qml
import matplotlib.pyplot as plt
import math

counts = []
trials = 
swap_fail_decimal =
irradiation_trials_decimal =

dev0 = qml.device('default.qubit', wires=5)


def brachialplexusfunction(counts,trials,irradiation_trials_decimal,swap_fail_decimal):
  '''
  This will not calculate the right number counts * trials and swap_fail_decimal to run, that must happen outside the circuit.
  irradiation_trials_decimal will affect the ancilla (wire 2 excited state) by automatically detecting if it is minority or majority case,
  |1> is the assumed minority case.
  This function implements a successful swap case of irradiation trials while documenting the swap_fail_decimal outside itself.
  Its purpose is to eventually run deeper into the computer for extra calculations that can rebuild the nerve junctions.

  Inputs:
  counts = Comma-seperated List. Number of times the circuit is initialized to collect an experimental sample of outcomes.
  trials = Integer. Collects a new set of initializations for a fresh experimental sample
  irradiation_trials_decimal = float, decimal. Fraction of trials using an excited wire 1 ancilla excited scheme (irradiation, wire 2).
                               Will not admit equal numbers of irradiated trials (because then it is balanced!)
                               Will not admit more than half excited trials (because then it is just a flipped measurement!)
  swap_fail_decimal = float, decimal. Number of times SWAP is not a feature of the circuit. Note that this part is converted
                      to shots outside the circuit, as this function handles the SWAP-success cases.

  Output: Qnode callable. Will be in the qnode list form, call individual elements using qnode[0], qnode[1], qnode[2].
  qml.expval(qml.PauliZ(wires=1)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 1 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  qml.expval(qml.PauliZ(wires=2)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 2 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  qml.expval(qml.PauliZ(wires=3)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 2 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  OPTIONAL: Write a print(number_of_times_run_message) to track your rate in a print-out.
  These have a corresponding next pair entry in axis_1 corresponding to total shots to be plotted. Each is an expectation value in a wire
  individually plotted against those trials. This data set is red.
  '''
  if irradiation_trials_decimal < 0.5:
    qml.PauliX(wires=1)

  qml.Hadamard(wires=0)
  qml.Hadamard(wires=4)
  qml.CNOT(wires=[0,1])
  qml.CNOT(wires=[4,3])
  qml.SWAP(wires=[1,2])
  qml.Toffoli(wires=[1,2,3])
  qml.Toffoli(wires=[0,4,1])

  number_of_times_run_message = str(1.0 - swap_fail_decimal) + ' (/1.0) of all SWAP-affected cases have a case of SWAP.'
  return qml.expval(qml.PauliZ(wires=1)), qml.expval(qml.PauliZ(wires=2)),qml.expval(qml.PauliZ(wires=3))


def brachialplexusfunction_fail(counts,trials,irradiation_trials_decimal,swap_fail_decimal):
  '''
  This will not calculate the right number counts * trials and swap_fail_decimal to run, that must happen outside the circuit.
  irradiation_trials_decimal will affect the ancilla (wire 2 excited state) by automatically detecting if it is a minority or majority case,
  |1> is the assumed minority case.
  |0> and |1> are perfect in these simulations, you will get ==0 expectation values.
  This function implements a failed swap case of irradiation trials while documenting the swap_fail_decimal and other variables outside itself.
  Its purpose is to eventually run deeper into the computer for extra calculations that can rebuild the nerve junctions.
  Some extra theory is found in plexustheory_fail. Please do it, it teaches you how to read CNOT as a differential diagnosis and peeks under the
  hood to explain why you see final expectation values ==0.

  Inputs:
  counts = Comma-seperated List. Number of times the circuit is initialized to collect an experimental sample of outcomes.
  trials = Integer. Collects a new set of initializations for a fresh experimental sample
  irradiation_trials_decimal = float, decimal. Fraction of trials using an excited wire 1 ancilla excited scheme (irradiation, wire 2).
                               Will not admit equal numbers of irradiated trials (because then it is balanced!)
                               Will not admit more than half excited trials (because then it is just a flipped measurement!)
  swap_fail_decimal = float, decimal. Number of times SWAP is not a feature of the circuit. Note that this part is converted
                      to shots outside the circuit, as this function handles the SWAP-success cases.

  Output: Qnode callable. Will be in the qnode list form, call individual elements using qnode[0], qnode[1], qnode[2].
  qml.expval(qml.PauliZ(wires=1)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 1 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  qml.expval(qml.PauliZ(wires=2)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 2 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  qml.expval(qml.PauliZ(wires=3)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 2 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  OPTIONAL: Write a print(number_of_times_run_message) to track your rate in a print-out.
  These have a corresponding next pair entry in axis_1 corresponding to total shots to be plotted. Each is an expectation value in a wire
  individually plotted against those trials. This data set is red.
  '''
  if irradiation_trials_decimal < 0.5:
    qml.PauliX(wires=1)

  qml.Hadamard(wires=0)
  qml.Hadamard(wires=4)
  qml.CNOT(wires=[0,1])
  qml.CNOT(wires=[4,3])
  qml.Toffoli(wires=[1,2,3])
  qml.Toffoli(wires=[0,4,2])

  number_of_times_run_message = str(swap_fail_decimal) + ' (/1.0) of all SWAP-affected cases fail a case of SWAP.'
  return qml.expval(qml.PauliZ(wires=1)), qml.expval(qml.PauliZ(wires=2)),qml.expval(qml.PauliZ(wires=3))

def brachialplexus_binaryexchange(counts,trials,irradiation_trials_decimal):
  '''
  This will not calculate the right number counts * trials and swap_fail_decimal to run, that must happen outside the circuit.
  irradiation_trials_decimal will affect the ancilla (wire 2 excited state) by automatically detecting if it is a minority or majority case,
  |1> is the assumed minority case.
  This is going to give the outcomes of the wires in the absence of any SWAP structure, representing pure if-this-then-that logic in nerve joining.

  Inputs:
  counts = Comma-seperated List. Number of times the circuit is initialized to collect an experimental sample of outcomes.
  trials = Integer. Collects a new set of initializations for a fresh experimental sample
  irradiation_trials_decimal = float, decimal. Fraction of trials using an excited wire 1 ancilla excited scheme (irradiation, wire 2).
                               Will not admit equal numbers of irradiated trials (because then it is balanced!)
                               Will not admit more than half excited trials (because then it is just a flipped measurement!)

  Output: Qnode callable. Will be in the qnode list form, call individual elements using qnode[0], qnode[1], qnode[2].
  qml.expval(qml.PauliZ(wires=1)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 1 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  qml.expval(qml.PauliZ(wires=2)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 2 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  qml.expval(qml.PauliZ(wires=3)) = float, -1.0 =< number =< 1.0. Inner (dot) product of the wire 2 vector with z-axis.
                                    Components in 0- are multiplied by +1, components in 1- are multiplied by -1
                                    This is the final value it takes on after a number of shots contribute their outcomes.
  These have a corresponding previous pair in axis_1 corresponding to total shots to be plotted. Each is an expectation value in a wire
  individually plotted against those trials. This data set is black.
  '''
  if irradiation_trials_decimal < 0.5:
    qml.PauliX(wires=1)

  qml.Hadamard(wires=0)
  qml.Hadamard(wires=4)
  qml.CNOT(wires=[0,1])
  qml.CNOT(wires=[4,3])
  qml.Toffoli(wires=[1,2,3])
  qml.CNOT(wires=[1,2])
  qml.Toffoli(wires=[0,4,2])
  return qml.expval(qml.PauliZ(wires=1)), qml.expval(qml.PauliZ(wires=2)),qml.expval(qml.PauliZ(wires=3))

axis_1 = []
axis_2_wire1 = []
axis_2_wire2 = []
axis_2_wire3 = []

for i in counts:
  iterations_passed = 0
  while iterations_passed < trials:
    shot_counter = 0
    shots = i * irradiation_trials_decimal * (1.0-swap_fail_decimal)
    shots = int(shots)
    qnode0 = qml.QNode(brachialplexusfunction, dev0, shots = shots)(counts,trials,irradiation_trials_decimal,swap_fail_decimal)
    shot_counter = shot_counter + shots
    shots = i * irradiation_trials_decimal * swap_fail_decimal
    shots = int(shots)
    qnode1 = qml.QNode(brachialplexusfunction_fail, dev0, shots = shots)(counts,trials,irradiation_trials_decimal,swap_fail_decimal)
    shot_counter = shot_counter + shots
    shots = i * (1.0-irradiation_trials_decimal) * (1.0-swap_fail_decimal)
    shots = int(shots)
    qnode2 = qml.QNode(brachialplexusfunction, dev0, shots = shots)(counts,trials,(1.0-irradiation_trials_decimal),swap_fail_decimal)
    shot_counter = shot_counter + shots
    shots = i * (1.0-irradiation_trials_decimal) * swap_fail_decimal
    shots = int(shots)
    qnode3 = qml.QNode(brachialplexusfunction_fail, dev0, shots = shots)(counts,trials,(1.0-irradiation_trials_decimal),swap_fail_decimal)
    shot_counter = shot_counter + shots
    axis_1.append(shot_counter)
    axis_2_wire1.append(qnode0[0] * irradiation_trials_decimal * (1.0-swap_fail_decimal) +
                        qnode1[0] * irradiation_trials_decimal * swap_fail_decimal +
                        qnode2[0] * (1.0-irradiation_trials_decimal) * (1.0-swap_fail_decimal) +
                        qnode3[0] * (1.0-irradiation_trials_decimal) * swap_fail_decimal)
    axis_2_wire2.append(qnode0[1] * irradiation_trials_decimal * (1.0-swap_fail_decimal) +
                        qnode1[1] * irradiation_trials_decimal * swap_fail_decimal +
                        qnode2[1] * (1.0-irradiation_trials_decimal) * (1.0-swap_fail_decimal) +
                        qnode3[1] * (1.0-irradiation_trials_decimal) * swap_fail_decimal)
    axis_2_wire3.append(qnode0[2] * irradiation_trials_decimal * (1.0-swap_fail_decimal) +
                        qnode1[2] * irradiation_trials_decimal * swap_fail_decimal +
                        qnode2[2] * (1.0-irradiation_trials_decimal) * (1.0-swap_fail_decimal) +
                        qnode3[2] * (1.0-irradiation_trials_decimal) * swap_fail_decimal)
    shot_counter = 0
    shots = i * irradiation_trials_decimal
    shots = int(shots)
    qnode4 = qml.QNode(brachialplexus_binaryexchange, dev0, shots = shots)(counts,trials,irradiation_trials_decimal)
    shot_counter = shot_counter + shots
    shots = i * (1.0 - irradiation_trials_decimal)
    shots = int(shots)
    qnode5 = qml.QNode(brachialplexus_binaryexchange, dev0, shots = shots)(counts,trials,(1.0-irradiation_trials_decimal))
    shot_counter = shot_counter + shots
    axis_1.append(shot_counter)
    axis_2_wire1.append(qnode4[0] * irradiation_trials_decimal + qnode5[0] * (1.0-irradiation_trials_decimal))
    axis_2_wire2.append(qnode4[1] * irradiation_trials_decimal + qnode5[1] * (1.0-irradiation_trials_decimal))
    axis_2_wire3.append(qnode4[2] * irradiation_trials_decimal + qnode5[2] * (1.0-irradiation_trials_decimal))
    iterations_passed = iterations_passed + 1

plotter_iterator = 0
data_collected = len(axis_1)

graph1_obj, graph1_axisobj  = plt.subplots()
graph2_obj, graph2_axisobj  = plt.subplots()
graph3_obj, graph3_axisobj = plt.subplots()

while plotter_iterator <= (data_collected-1):
  x = math.log10(axis_1[plotter_iterator])
  y = axis_2_wire1[plotter_iterator]
  graph1_axisobj.scatter(x,y,marker='o', color='r')
  plotter_iterator = plotter_iterator + 2

plotter_iterator = 0

while plotter_iterator <= (data_collected-1):
  x = math.log10(axis_1[plotter_iterator])
  y = axis_2_wire2[plotter_iterator]
  graph2_axisobj.scatter(x,y,marker='o', color='r')
  plotter_iterator = plotter_iterator + 2

plotter_iterator = 0

while plotter_iterator <= (data_collected-1):
  x = math.log10(axis_1[plotter_iterator])
  y = axis_2_wire3[plotter_iterator]
  graph3_axisobj.scatter(x,y,marker='o', color='r')
  plotter_iterator = plotter_iterator + 2

plotter_iterator = 1

while plotter_iterator <= (data_collected-1):
  x = math.log10(axis_1[plotter_iterator])
  y = axis_2_wire1[plotter_iterator]
  graph1_axisobj.scatter(x,y,marker='o', color='k')
  plotter_iterator = plotter_iterator + 2

plotter_iterator = 1

while plotter_iterator <= (data_collected-1):
  x = math.log10(axis_1[plotter_iterator])
  y = axis_2_wire2[plotter_iterator]
  graph2_axisobj.scatter(x,y,marker='o', color='k')
  plotter_iterator = plotter_iterator + 2

plotter_iterator = 1

while plotter_iterator <= (data_collected-1):
  x = math.log10(axis_1[plotter_iterator])
  y = axis_2_wire3[plotter_iterator]
  graph3_axisobj.scatter(x,y,marker='o', color='k')
  plotter_iterator = plotter_iterator + 2

graph1_axisobj.set_xlim(0.0,math.log10(max(axis_1)))
graph2_axisobj.set_xlim(0.0,math.log10(max(axis_1)))
graph3_axisobj.set_xlim(0.0,math.log10(max(axis_1)))
graph1_axisobj.set_ylim(-1.0,1.0)
graph2_axisobj.set_ylim(-1.0,1.0)
graph3_axisobj.set_ylim(-1.0,1.0)
graph1_axisobj.set_xlabel("Log10 Number of Shots (individual)")
graph2_axisobj.set_xlabel("Log10 Number of Shots (individual)")
graph3_axisobj.set_xlabel("Log10 Number of Shots (individual)")
graph1_axisobj.set_ylabel("Wire 1 Signal (/arbitrary return)")
graph2_axisobj.set_ylabel("Wire 2 Signal (/arbitrary return)")
graph3_axisobj.set_ylabel("Wire 3 Signal (/arbitrary return)")
graph1_axisobj.set_title("Expectation number of " + str(trials) + " shot-weighed instances of computer", loc='left', pad=20)
graph2_axisobj.set_title("Expectation number of " + str(trials) + " shot-weighed instances of computer", loc='left', pad=20)
graph3_axisobj.set_title("Expectation number of " + str(trials) + " shot-weighed instances of computer", loc='left', pad=20)

print('Red sets have swap-fail conditions. Black sets are binary exchange, the nerve joining simply flips another signal as false')
plt.show()

print('Your circuits: \n SWAP-success spectroscopy')
qnode_0 = qml.QNode(brachialplexusfunction, dev0)
print(qml.draw(qnode_0)(1.0,1.0,1.0,1.0))

qnode_2 = qml.QNode(brachialplexusfunction_fail, dev0)
print(' SWAP-fail Spectroscopy')
print(qml.draw(qnode_2)(1.0,1.0,1.0,1.0))

print(' SWAP-indifferent spectroscopy')
qnode_1 = qml.QNode(brachialplexus_binaryexchange, dev0)
print(qml.draw(qnode_1)(1.0,1.0,1.0))

def plexusfail_theory(irradiation_trials_decimal):
  '''
  This is the probabilities of every bit combination in the absence of SWAP.
  Previously, I said to think about 1s and 0s in every part of the circuit, reflecting experimental signals, appropriately called the
  vector bundle model. This function peeks under the hood, Pennylane has perfect 0 probability for |1> in any |0> starting state, which is
  all of them except H and wire 1 |1> excitations. In the H case, a perfect 50/50 forms, the twist is that it is a product of psi^2, because
  the un-squared probabilities the computer sees are an addition, |0> + |1> mapping |0>, and subtraction |0> - |1> mapping |1>. There is a
  reciprocal root-2 I hid for brevity, responsible for the 1-total probability. The point is that wire 1 and wire 2 see 50/50 at all times
  of a descriptive signal subtraction or addition later on, which can be interfered conventionally like a polynomial:

  difference of squares = (a+b)(a-b) = a**2 - b**2
  parabola = (a+b)(a+b) = a**2 + 2ab + b**2
  parabola = (a-b)(a-b) = a**2 - 2ab + b**2

  The purpose of these functions is to practice reading CNOT as differential diagnoses, while returning the probabilities responsible for
  the graphs above. Try getting comfortable with this before incorporating the wire 1 excitation scheme. Read it as if this, then that. For example
  if C5 were 1 then C6 was not 1, it was 0 giving the top wire =1 (C5/C6 superposition) between them and the Posterior Trunk (wire 2) a depletion
  of a 1 (=0). Tiffoli gates are read identically, the median branch has not more than two zeroes between them, or else 1 and 0 reflect a C7 start,
  but not either of T1 or C5, so on. I find using refute-to in your notes will help you organize your diagnosis: C5 = 1 refute to C6 =0 from =1.

  Inputs:
  irradiation_trials_decimal = Determines whether minority case wire 1 |1> is active if <0.5 or majority case wire 1 |0>

  Outputs:
  qml.probs = List of probabilities adding up to one. Read in the conditional probability convention, wire 5 = 0 given 0000 = |00000> = entry[0].
              |00001> = entry[1], |00010> = entry[2], |00011> = entry[3]. There is a physical reason it is laid out like binary, irrelevant to
              the experiment.
  '''
  if irradiation_trials_decimal < 0.5:
    qml.PauliX(wires=1)

  qml.Hadamard(wires=0)
  qml.Hadamard(wires=4)
  qml.CNOT(wires=[0,1])
  qml.CNOT(wires=[4,3])
  qml.Toffoli(wires=[1,2,3])
  qml.Toffoli(wires=[0,4,2])

  return qml.probs()

qnode1 = qml.QNode(plexusfail_theory, dev0)

psi0 = qnode1(0.1)
psi1 = qnode1(0.6)

print('===\nProbabilities Breakdown\n===')
print('Bit Combination Probabilities of the circuits above are (normalized to 1):')
print('SWAP fails in presence of irradiation:')
print(psi0)
print('SWAP fails in absence of irradiation:')
print(psi1)

def plexus_theory(irradiation_trials_decimal):
  '''
  This is the probabilities of every bit combination in the presence of SWAP.
  Previously, I said to think about 1s and 0s in every part of the circuit, reflecting experimental signals, appropriately called the
  vector bundle model. This function peeks under the hood, Pennylane has perfect 0 probability for |1> in any |0> starting state, which is
  all of them except H and wire 1 |1> excitations. In the H case, a perfect 50/50 forms. For more details, see plexusfail_theory

  The purpose of these functions is to practice reading CNOT as differential diagnoses, while returning the probabilities responsible for
  the graphs above. Try getting comfortable with this before incorporating wire 1 (wire 2 dump) excitation schema. Read it as if this, then that.
  For example, if C5 were 1 then C6 was not 1, it was 0 giving the top wire =1 (C5/C6 superposition) between them and the Posterior Trunk (wire 2)
  a depletion of a 1 (=0). Tiffoli gates are read identically, the median branch has not more than two zeroes, or else 1 and 0 reflect a C7 start,
  but not either of T1 or C5, so on. I find using refute-to in your notes will help you organize your diagnosis: C5 = 1 refute to C6 =0 from =1.

  Inputs:
  irradiation_trials_decimal = Determines whether minority case wire 1 |1> is active if <0.5 or majority case wire 1 |0>

  Outputs:
  qml.probs = List of probabilities adding up to one. Read in the conditional probability convention, wire 5 = 0 given 0000 = |00000> = entry[0].
              |00001> = entry[1], |00010> = entry[2], |00011> = entry[3]. There is a physical reason it is laid out like binary, irrelevant to
              the experiment.
  '''
  if irradiation_trials_decimal < 0.5:
    qml.PauliX(wires=1)

  qml.Hadamard(wires=0)
  qml.Hadamard(wires=4)
  qml.CNOT(wires=[0,1])
  qml.CNOT(wires=[4,3])
  qml.SWAP(wires=[1,2])
  qml.Toffoli(wires=[1,2,3])
  qml.Toffoli(wires=[0,4,1])

  return qml.probs()

qnode2 = qml.QNode(plexus_theory, dev0)
psi2 = qnode2(0.1)
psi3 = qnode2(0.6)

print('SWAP succeeds in presence of irradiation:')
print(psi2)
print('SWAP succeeds in absence of irradiation:')
print(psi3)

def binaryexchange_theory(irradiation_trials_decimal):
  '''
  This is the probabilities of every bit combination in the case of CNOT targeting lateral cord integration.
  Previously, I said to think about 1s and 0s in every part of the circuit, reflecting experimental signals, appropriately called the
  vector bundle model. This function peeks under the hood, Pennylane has perfect 0 probability for |1> in any |0> starting state, which is
  all of them except H and wire 2 |1> excitations. In the H case, a perfect 50/50 forms. For more details, see plexusfail_theory

  The purpose of these functions is to practice reading CNOT as differential diagnoses, while returning the probabilities responsible for
  the graphs above. Try getting comfortable with this before incorporating wire 1 (wire 2 dump) excitation schema. Read it as if this, then that.
  For example, if C5 were 1 then C6 was not 1, it was 0 giving the top wire =1 (C5/C6 superposition) between them and the Posterior Trunk (wire 2)
  a depletion of a 1 (=0). Tiffoli gates are read identically, the median branch has not more than two zeroes, or else 1 and 0 reflect a C7 start,
  but not either of T1 or C5, so on. I find using refute-to in your notes will help you organize your diagnosis: C5 = 1 refute to C6 =0 from =1.

  Inputs:
  irradiation_trials_decimal = Determines whether minority case wire 1 |1> is active if <0.5 or majority case wire 1 |0>

  Outputs:
  qml.probs = List of probabilities adding up to one. Read in the conditional probability convention, wire 5 = 0 given 0000 = |00000> = entry[0].
              |00001> = entry[1], |00010> = entry[2], |00011> = entry[3]. There is a physical reason it is laid out like binary, irrelevant to
              the experiment.
  '''
  if irradiation_trials_decimal < 0.5:
    qml.PauliX(wires=1)

  qml.Hadamard(wires=0)
  qml.Hadamard(wires=4)
  qml.CNOT(wires=[0,1])
  qml.CNOT(wires=[4,3])
  qml.Toffoli(wires=[1,2,3])
  qml.CNOT(wires=[1,2])
  qml.Toffoli(wires=[0,4,2])
  return qml.probs()

qnode3 = qml.QNode(plexus_theory, dev0)

psi4 = qnode3(0.1)
psi5 = qnode3(0.6)

print('Binary nerve integration in presence of irradiation:')
print(psi4)
print('Binary nerve integration in absence of irradiation:')
print(psi5)

print('Just the SWAP-fail experiments:')
uppercase_psi = (swap_fail_decimal * irradiation_trials_decimal * psi0 +
                 swap_fail_decimal * (1.0-irradiation_trials_decimal)* psi1 +
                 (1.0-swap_fail_decimal) * irradiation_trials_decimal * psi2 +
                 (1.0-swap_fail_decimal) * (1.0 - irradiation_trials_decimal)* psi3)
##via m+n = 1, a+b = 1 of swap-fail and irradiation schemes on normalization q probability.
##On maq + mbq + naq + nbq = (m+n)aq + (m+n)bq = q
print(uppercase_psi)

print('Just the binary exchange experiment:')
uppercase_psi = ((1.0-irradiation_trials_decimal) * psi4 +
                 irradiation_trials_decimal * psi5)
##aq + bq = (a+b)q on a+b = 1 = q
print(uppercase_psi)

print('The equally-weighed bit contributions in a joining set:')
uppercase_psi = (psi0 + psi1 + psi2 + psi3 + psi4 + psi5) / 6.0
print(uppercase_psi)

print('Incorporating a joining set of experimental probability:')
uppercase_psi = (swap_fail_decimal * irradiation_trials_decimal * psi0 +
                 swap_fail_decimal * (1.0-irradiation_trials_decimal) * psi1 +
                 (1.0-swap_fail_decimal) * irradiation_trials_decimal * psi2 +
                 (1.0-swap_fail_decimal) * (1.0 - irradiation_trials_decimal) * psi3 +
                 (1.0-irradiation_trials_decimal) * psi4 +
                 irradiation_trials_decimal * psi5) / 2.0
##via m+n = 1, a+b = 1 for swap-fail and irradiation respectively and q normal state.
##on maq + mbq + naq + nbq + aq + bq and apply (m+n)*(aq + bq) on those last terms.
##gather m+n then a+b then just 2q.
print(uppercase_psi)
