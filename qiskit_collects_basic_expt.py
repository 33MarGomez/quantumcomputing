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

#run in a cell if in colab and delete to prevent accidental key disclosure
QiskitRuntimeService.save_account(
    channel="ibm_quantum_platform",
    token="YOUR_TOKEN_HERE",
    instance="YOUR_CDR"
    overwrite=True,
    set_as_default=True,
)

#script
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

service = QiskitRuntimeService(channel="ibm_quantum_platform")
service = QiskitRuntimeService()
backend = AerSimulator()
counts = preset_run_collect(ghz, backend, shots=1024)
good_go = check_ghz(counts)
if good_go==True:
    expt_str = 'GHZ'
    trial_str = '1'
    backend = service.least_busy(operational=True, simulator=False, min_num_qubits=3)
    backend = AerSimulator()
    result = preset_run_collect(ghz, backend, shots=1024)
    make_csv(result,expt_str,trial_str)
