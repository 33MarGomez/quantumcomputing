"""
This script runs preliminary tests and collection on basic building blocks to make quantum
random data go far. These preliminary tests use aer, classical simulation, to check for
desireable properties before permitting an API pass to IBM quantum hardware. These are all
the experiments I have gathered thus far and have uses for.

If you'd like to get involved, run this code and send me the .csv. Your generated values
will be used internally, but due to data integrity standards, will not be published.
This work is strictly for research purposes on the applicability of qubit-aligned
electronics and represents no dual-use cases.

-Marco
"""
!pip install 'qiskit[visualization]' qiskit-ibm-runtime qiskit-aer
!pip install pennylane --upgrade

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# IBM Runtime specific imports
from qiskit_ibm_runtime import SamplerV2 as Sampler, QiskitRuntimeService

import sys

import qiskit
import qiskit_aer
import qiskit_ibm_runtime

print("Python:", sys.version.split()[0])
print("qiskit:", qiskit.__version__)
print("qiskit-aer:", qiskit_aer.__version__)
print("qiskit-ibm-runtime:", qiskit_ibm_runtime.__version__)

import pandas as pd
import numpy as np
import pennylane as qpl

#run in a cell if in colab and delete to prevent accidental key disclosure
QiskitRuntimeService.save_account(
    channel="ibm_quantum_platform",
    token="YOUR_TOKEN_HERE",
    instance="YOUR_CDR",
    overwrite=True,
    set_as_default=True,
)

#script
three = qpl.device('default.qubit', wires=3)

@qpl.qnode(three)
def fabric_spin_echo():
    qpl.RX(phi=np.pi*0.5,wires=0)
    qpl.IsingZZ(phi=np.pi*0.5,wires=[0,1])
    qpl.IsingYY(phi=np.pi*-0.5,wires=[0,1])
    qpl.RY(phi=np.pi*2.0,wires=0)
    qpl.IsingZZ(phi=np.pi*0.5,wires=[1,2])
    qpl.IsingYY(phi=np.pi*-0.5,wires=[1,2])
    qpl.IsingXX(phi=np.pi,wires=[0,1])
    qpl.IsingXX(phi=np.pi,wires=[1,2])
    return qpl.state(),qpl.probs(wires=[0,1,2])

def check_sechos_expt():
    """
    Verifies aspects of the inner product with (s)pin (echos) applied to
    prevent qubit evolution during collection.
    
    returns
    Bool - Whether sechos has the correct properties of:
         |101> & |111> are conj. amplitudes of same magnitude
         |000> & |010> are real, same-magnitude amplitudes
         |001> and |011> are the same imaginary amplitude
         |100> and |110> are the same real amplitude
         all probabilities are equal
    
    Theory- Exchanging upper 1 and 0 amplitudes exchanges complex and real
    axis anti-symmetry on upper and lower registers. Because register
    position cannot be exchanged, and complex phasing is a laboratory-frame
    construct, no outstanding degree of freedom exists.
    In the first (np.pi/4.0) ZZ and RX application, cos |0> takes on the same
    value, the only such solution to not changing qubit amplitudes. As RX
    (np.pi/2.0) is -i|1> and this projection is shared with RZ|1>, applying
    YY towards -1|1> establishes chirality. The shared -1|1> Z cosine, -i|1>
    RZ sine, RY cos |0> establish this. A spin echo is then applied as it
    commutes with YY. Updates to the next qubit are expected to not impact
    evolution as any interference in the qubit should evolve it to its starting
    state in x-freedom. Therefore, XX should appear as the final ingredient in
    the inner product.
    """
    vector,squared = fabric_spin_echo()
    if (round(np.imag(vector[5]),4)
        == -1.0*round(np.imag(vector[7]),4)):
        pass
    else:
        print("|101> and |111> are not complex conjugate amplitudes\n"+
              "as expected. Pass to API aborted.\nPlease check your gates"
              +"\nas an undue property was found in qpl check")
        return False
    if (round(np.real(vector[0]),4)
        == -1.0*round(np.real(vector[2]),4)):
        pass
    else:
        print("|000> and |010> are not oppositely-signed amplitudes\n"+
              "as expected. Pass to API aborted.\nPlease check your gates"
              +"\nas an undue property was found in qpl check")
        return False
    if (round(np.real(vector[4]),4) == round(np.real(vector[6]),4) and
        round(np.imag(vector[1]),4) == round(np.imag(vector[3]),4) and
        abs(round(np.real(vector[4]),2)) > 0.1 and
        abs(round(np.imag(vector[1]),2)) > 0.1 and
        abs(round(np.real(vector[4]),2)) < 0.5 and
        abs(round(np.real(vector[4]),2)) < 0.5
    ):
        pass
    else:
        print("Failed to meet any of the following conditions:-\n" +
              "|001> and |011> are imaginary, equal and less than the "+
              "magnitude of the RX(ZZ) evolution. \n|100> and |110> are "+
              "real, equal and less than the magnitude of the RX(ZZ)"+
              "evolution\n-as expected. "+
              "Pass to API aborted.\nPlease check your gates"
              +"\nas an undue property was found in qpl check")
        return False
    squared = np.round(squared,3)
    if np.all(squared==squared[0]):
        return True
    else:
        print("Probabilities not uniform. Please review gate set\n"+
              "API call aborted for gate revision")
        return False

def check_ghz(counts):
    """
    Test function, all GHZ states satisfy expected properties. In this case,
    check only 000 and 111 occured from the simulator

    Args:
    counts:list[str] = data.meas.get_bitstrings() return to be processed.
    """
    num = range(len(counts))
    series = pd.Series(counts,index=num)
    series = series.astype(str)
    unexpecteds = series[(series!='000') & (series!='111')]
    if unexpecteds.size==0:
        return True
    else:
        print("Pass to API aborted.\nPlease check your gates"
              +"\nas an undue property was found in AER check"
        )
        return False

def AER_count_check(circuit, backend, shots=1024):
    """
    Strictly to visualize the outcome before proceding.

    Args:
        circuit (QuantumCircuit): The quantum circuit to run.
        backend: The Qiskit backend (real device or simulator).
        shots (int): The number of shots to run the circuit.

    Returns:
        dict = A count of each qubit combination across shots.
    """
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(circuit)
    sampler = Sampler(mode=backend)
    job = sampler.run([isa_circuit], shots=shots)
    result = job.result()
    return result[0].data.meas.get_counts()

def preset_run_collect(circuit, backend, shots=1024):
    """
    Runs a quantum circuit on a specified backend and returns the measurement counts.

    Args:
        circuit (QuantumCircuit): The quantum circuit to run.
        backend: The Qiskit backend (real device or simulator).
        shots (int): The number of shots to run the circuit.

    Returns:
        list[str] = A list of measurement bitstrings in 0 or 1.
    """
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(circuit)

    sampler = Sampler(mode=backend)

    job = sampler.run([isa_circuit], shots=shots)
    result = job.result()

    return result[0].data.meas.get_bitstrings()

def make_csv(counts,experiment_string,trial_string):
    """
    Creates a .csv using Pandas to save results. This also prints the full
    file as a back-up if the runtime is prematurely aborted by mistake, so
    it make be retrieved at a later time.

    Args
    counts:list[str] = data.meas.get_bitstrings() return to be processed.
    experiment_string = name of the phenomenon
    trial_string = current trial, not to overwrite on your disk

    Returns:
    None

    Writes:
    csv_{experiment_string}{trial_string_expt}.csv
    prints index,bitstring pairs of the experiment
    """
    num = range(len(counts))
    series = pd.Series(result,index=num)
    series = series.astype(str)
    series.index.name = "index"
    series = series.astype(str)
    series.name = "bitstring"
    file_name = f"csv_{experiment_string}{trial_string}.csv"
    series.to_csv(file_name)
    print(series.to_string())

ghz = QuantumCircuit(3)
ghz.h(0)
ghz.cx(0, 1)
ghz.cx(0, 2)
ghz.measure_all()

sechos = QuantumCircuit(3)
sechos.rx(np.pi*0.5,0)
sechos.rzz(np.pi*0.5,0,1)
sechos.ryy(np.pi*-0.5,0,1)
sechos.ry(np.pi*2.0,0)
sechos.rzz(np.pi*0.5,1,2)
sechos.ryy(np.pi*-0.5,1,2)
sechos.rxx(np.pi,0,1)
sechos.rxx(np.pi,1,2)
sechos.measure_all()

service = QiskitRuntimeService(channel="ibm_quantum_platform")
service = QiskitRuntimeService()
backend = AerSimulator()
counts = preset_run_collect(ghz, backend, shots=1024)
good_go = check_ghz(counts)
if good_go==True:
    expt_str = 'GHZ'
    trial_str = '1'
    backend = service.least_busy(operational=True, simulator=False, min_num_qubits=3)
    backend = AerSimulator() #trip fuse against run-all
    result = preset_run_collect(ghz, backend, shots=1024)
    make_csv(result,expt_str,trial_str)
good_go = check_sechos_expt()
if good_go==True:
    expt_str = 'fabric_spin_echo'
    trial_str = '1'
    backend = service.least_busy(operational=True, simulator=False, min_num_qubits=3)
    backend = AerSimulator() #trip fuse against run-all
    result = preset_run_collect(sechos, backend, shots=1024)
    make_csv(result,expt_str,trial_str)
