"""
Stores values in variable storage to randomly work through
bitstrings generated in a quantum computer.

Variables:
seed:int = choose a number, ideally itself picked from default_rng or
       random number strategy of your choice.
filename:str = for python to pick out the file to turn from .csv to
           a workable format.
values_grabbed:int = the total number of numbers you want 

returns:
storage:list[str] = the bitstrings in order of their choosing
"""

import pandas as pd
import numpy as np

SEED =
FILENAME =
VALUES_GRABBED =

rng = np.random.default_rng(SEED)

series = pd.read_csv(FILENAME,index_col="index",header=0,dtype=str)
rng_iter = series.shape[0]
values_iter = 0
storage = []

while values_iter < VALUES_GRABBED:
    current = rng.integers(0,rng_iter,dtype=int)
    value = series["bitstring"].pop(current)
    series.reset_index(drop=True)
    storage.append(value)
    rng_iter-=1
    values_iter+=1
