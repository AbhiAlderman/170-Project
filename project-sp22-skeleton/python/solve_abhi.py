"""Solves an instance.

Modify this file to implement your own solvers.

For usage, run `python3 solve.py --help`.
"""

import argparse
import numpy
import math

from pathlib import Path
from typing import Callable, Dict

from instance import Instance
from solution import Solution
from point import Point
from file_wrappers import StdinFileWrapper, StdoutFileWrapper


def solve_naive(instance: Instance) -> Solution:
    return Solution(
        instance=instance,
        towers=instance.cities,
    )

def solve_dp(instance: Instance) -> Solution:
    """
    DP approach:
    get list of cities
    get bounds for towers (0 -> instance.grid_size_length)
    pass these into a recursive helper func
    recursive func:
    choose tower that covers most (greedy)
    if no penalty, just choose it
    if penalty, choose min between choosing it or
        choosing the next best without penalty
    return recursive func, removing that tower from list
    
    PROBLEM: want to return the actual list of towers we want,
        but we need to return the penalty for min to work
    ^ it should be fine to pass in a list and edit it
        while still returning penalty, but check this 
        if it doesnt work
    """
    cities = instance.cities
    solution = recursiveHelper(cities) #solution = [penalty, tower list]
    return Solution(
        instance=instance,
        towers=solution[1],
    )


def recursiveHelper(cities, towers, penalty_radius):
    """
    if there are no cities left then just return the current penalty
    next, find greedy tower aka tower that covers the most
    """
    
    #case where there are no cities left
    if not cities:
        return [findPenalty(towers, penalty_radius), towers]
    

    #return list(penalty, tower list)
    

def findGreedyTower(cities, towers, penalty_radius, length):
    """
    find the next tower that covers the most cities
    returns a [towerP, towerMP]
    towerP = greediest tower w/ or w/o penalty
    towerMP = greediest tower w/ minimum penalty

    Psuedocode:
        double for loop to get each possible city coordinate
        for each possible tower:
            if tower exists in towers, skip
            find # of cities covered
            if not greater than towerP or towerNP move on
            if greater than, then:
                if has penalty, replace towerP
                if min penalty, replace towerMP and maybe towerP
            
    """
    #brute force method :(
        for i in range(length):
            for j in range(length):

def findPenalty(towers, penalty_radius):
    """
    This function finds the penalty with current towers
    does NOT check if towers are valid
    """
    penalty_radius = penalty_radius * penalty_radius
    total = 0
    for i in range(len(towers)):
        currentTower = towers[i]
        wj = 0
        for j in range(len(towers)):
            if i == j:
                continue
            if towers[i].distance_sq(towers[j]) <= penalty_radius:
                wj += 1
        total += 170 * math.exp(0.17 * wj)
    return total




SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive
}


# You shouldn't need to modify anything below this line.
def infile(args):
    if args.input == "-":
        return StdinFileWrapper()

    return Path(args.input).open("r")


def outfile(args):
    if args.output == "-":
        return StdoutFileWrapper()

    return Path(args.output).open("w")


def main(args):
    with infile(args) as f:
        instance = Instance.parse(f.readlines())
        solver = SOLVERS[args.solver]
        solution = solver(instance)
        assert solution.valid()
        with outfile(args) as g:
            print("# Penalty: ", solution.penalty(), file=g)
            solution.serialize(g)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve a problem instance.")
    parser.add_argument("input", type=str, help="The input instance file to "
                        "read an instance from. Use - for stdin.")
    parser.add_argument("--solver", required=True, type=str,
                        help="The solver type.", choices=SOLVERS.keys())
    parser.add_argument("output", type=str,
                        help="The output file. Use - for stdout.",
                        default="-")
    main(parser.parse_args())
