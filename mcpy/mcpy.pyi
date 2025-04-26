from __future__ import annotations
from typing import Iterable, Callable, overload
import numpy as np

'''
Stub file to assist type hinting of python front end.
All operations are compiled in c++ code and have been exposed to python.

The source code is located in the header files accompanied,
and the implementations in the cpp files.


The compile command is in mcpyext_main.cpp,
however the pybind developer header files need to be installed,
and the command needs to be altered for different platforms (e.g. Mac OS)
or different python version.

The compile command (and the corresponding .so binary file)
supports linux with python version 3.12
'''

OBSERVABLE = Callable[[STATE], float]


class Sample:

    '''
    Sample class is a collection of data of an observable (array of numbers)
    Its methods can be used to perform statistics
    '''

    def __init__(self, data: Iterable[float]):...

    @property
    def data(self)->np.ndarray:...

    @property
    def N(self)->int:...

    @property
    def stat(self)->str:...

    def mean(self)->float:...

    def std(self)->float:...

    def popul_std(self)->float:...

    def error(self)->float:...

    def bin_it(self)->BinnedSample:...


class BinnedSample:

    '''
    When initialized, binning analysis is performed,
    and each binning level is inside .binned_samples.
    '''

    def __init__(self, data: Iterable[float]):...

    @property
    def max_level(self)->int:...

    @property
    def converged(self)->bool:...

    @property
    def tau_estimate(self)->float:...

    @property
    def binned_samples(self)->list[Sample]:...



class STATE:
    '''
    Base class for any Markov Chain State.
    
    Although empty, it exists to assist abstract behaviors like isinstance(State_subclass, STATE) == True
    '''
    pass


class SpinState(STATE):

    @overload
    def __call__(self, i: int)->int:...

    @overload
    def __call__(self, i: int, j: int)->int:...

    @property
    def spins(self)->np.ndarray[int]:... #2D lattice of spins +1 / -1

    @property
    def sites(self)->int:... #number of sites

    @property
    def M(self)->float:... #magnetization

    @property
    def energy(self)->float:... #total energy



class MarkovChain:

    @property
    def state(self)->STATE:...

    def update(self, method: str, steps=1)->None:...


class IsingModel2DMarkovChain(MarkovChain):

    def __init__(self, T: float, Lx: int, Ly: int):...

    @property
    def state(self)->SpinState:...

    def ssf_update(self)->None:...

    def wolff_update(self)->None:...



class MonteCarlo:

    def __init__(self, markov_chain: MarkovChain): ...

    @property
    def data(self)->list[STATE]:...

    @property
    def N(self)->int:...

    def sample(self, A: OBSERVABLE)->Sample:...

    def update(self, method: str, steps: int, sweeps=0):...

    def thermalize(self, method: str, sweeps: int):...


class IsingModel2D(MonteCarlo):

    def __init__(self, T: float, Lx: int, Ly: int):...

    @property
    def data(self)->list[SpinState]:...

    @property
    def Temp(self)->float:...

    def sample(self, A: Callable[[SpinState], float])->Sample:...

    def energy_sample(self)->Sample:...

    def ssf_update(self, steps: int, sweeps=0)->None:...

    def wolff_update(self, steps: int, sweeps=0)->None:...

    def ssf_thermalize(self, sweeps: int)->None:...

    def wolff_thermalize(self, sweeps: int)->None:...

#perform many Monte Carlo simulations in parallel
def update_all(sims: Iterable[MonteCarlo], method: str, steps: int, sweeps=0, threads=-1)->None:...
