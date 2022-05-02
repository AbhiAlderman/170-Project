import argparse
from pathlib import Path
from typing import Callable, Dict

from instance import Instance
from solution import Solution
from file_wrappers import StdinFileWrapper, StdoutFileWrapper

#my imports
import cvxpy as cp
import numpy as np
import pickle

#load in helpers
#region
with open('smallpenalty.pkl', 'rb') as f:
    sp = pickle.load(f)
with open('medpenalty.pkl', 'rb') as f:
    mp = pickle.load(f)
with open('largepenalty.pkl', 'rb') as f:
    lp = pickle.load(f)
with open('smallsignal.pkl', 'rb') as f:
    ss = pickle.load(f)
with open('medsignal.pkl', 'rb') as f:
    ms = pickle.load(f)
with open('largesignal.pkl', 'rb') as f:
    ls = pickle.load(f)
#endregion

def solve_my(instance: Instance) -> Solution:
    n = instance.grid_side_length
    X = cp.Variable(shape=(n, n), boolean=True) #tower pleacemenat variable

    obj = 0
    for i in range(n):
        for j in range(n):
            numnearby = 0
            ##get all the points in a circle of radius r_p
            ##assume all small for now, change later!
            w = 0
            for (a,b) in sp[i,j]:
                w += X[a][b]
            obj += X[i][j] * 170 * np.exp(0.17 * w)
    constraints = []
    for city in instance.cities:
        i, j = city.x, city.y 
        
    prob = cp.Problem()
    return Solution(
        instance=instance,
        towers = []
    )