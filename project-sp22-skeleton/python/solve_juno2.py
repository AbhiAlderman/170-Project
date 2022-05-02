import argparse
from pathlib import Path
from typing import Callable, Dict

from instance import Instance
from solution import Solution
from file_wrappers import StdinFileWrapper, StdoutFileWrapper

#my imports
import cvxpy as cp
import cvxopt
import numpy as np
import pickle
from point import Point

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
##
with open('testpenalty.pkl', 'rb') as f:
    tp = pickle.load(f)
with open('testsignal.pkl', 'rb') as f:
    ts = pickle.load(f)
##
#endregion

def dstsq(x1, y1, x2, y2):
    return (x1 - x2)**2 + (y1 - y2)**2

global testvar

#region debug
def debm(matvar, d):
    for i in range(d):
        str = ''
        for j in range(d):
            str += '{} '.format(matvar[i][j].value)
        print(str)

#endregion

def solve_j2(instance: Instance) -> Solution:
    answer = []
    d = instance.grid_side_length
    dictp = {10: tp, 30: sp, 50: mp, 100: lp}[d] #penaltydict
    dicts = {10: ts, 30: ss, 50: ms, 100: ls}[d] #penaltydict
    rsig = instance.coverage_radius
    rpen = instance.penalty_radius
    Xs = []
    for i in range(d):
        Xs.append([] * d)
    print(Xs)
    for i in range(d):
        for j in range(d):
            Xs[i].append( cp.Variable(boolean=True) ) #tower pleacemenat variable

    obj = 0
    for i in range(d):
        for j in range(d):
            numnearby = 0
            ##get all the points in a circle of radius r_p
            ##assume all small for now, change later!
            w = 0
            for (a,b) in dictp[i,j]:
                w += Xs[a][b]
            #obj += Xs[i][j] * 170 * cp.exp(0.17 * w)
            obj += cp.exp( Xs[i][j] + cp.exp(0.17 * w) ) #black magic time
    obj = cp.Minimize(obj)
    constraints = []
    #selectors = []
    radcov = instance.coverage_radius
    #for city in instance.cities:
    for count in range(len(instance.cities)):
        city = instance.cities[count]
        i, j = city.x, city.y
        sumtows = 0
        for (a, b) in dicts[i, j]: #for each possible tower covering this city
            sumtows += Xs[a][b]
        constraints.append(sumtows >= 1)


    prob = cp.Problem(obj, constraints)
    prob.solve()
    #print('---------------------')
    #debm(Xs, d)
    towercoords = []
    for i in range(d):
        for j in range(d):
            print(Xs[i][j].value)
            if Xs[i][j].value > 0.1:
                towercoords.append( (i,j) )
                answer.append(Point(i, j))
    print('---')
    print(towercoords)
    return Solution(
        instance=instance,
        towers = answer
    )

def main(): #debug driver
    '''solve_j(Instance(
            grid_side_length=10,
            coverage_radius=1,
            penalty_radius=2,
            cities=[
                Point(x=5, y=9),
                Point(x=2, y=2),
                Point(x=1, y=1),
                Point(x=2, y=3),
                Point(x=0, y=1),
            ]))'''
    with open('inputs/small/001.in') as f:
        instance = Instance.parse(f.readlines())
        solver = solve_j2
        solution = solver(instance)
        assert solution.valid()
        with open('outputj/small/001.out') as g:
            print("# Penalty: ", solution.penalty(), file=g)
            solution.serialize(g)
    
if __name__ == "__main__":
    main()
