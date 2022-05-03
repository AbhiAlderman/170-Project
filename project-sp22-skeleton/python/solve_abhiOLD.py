"""Solves an instance.

Modify this file to implement your own solvers.

For usage, run `python3 solve.py --help`.
"""

import argparse
from ast import List
import math

from pathlib import Path
from tkinter import N
from typing import Callable, Dict

from instance import Instance
from solution import Solution
from point import Point
from file_wrappers import StdinFileWrapper, StdoutFileWrapper

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

def solve_naive(instance: Instance) -> Solution:
    return Solution(
        instance=instance,
        towers=instance.cities,
    )
"""
HELPER FUNCTIONS:
    solve_dp = dp solution to solve problem
    recursive helper = recursive step for solve dp
    findGreedyTower = returns [towerP, towerNP]
        towerP = best tower regardless of penalty
        towerNP = best tower w/ min penalty
        penalty defined as how many other cities will get penalized,
            not defined as the actual penalty value
    getPointsInRadius = get list of possible points from a radius
        3 helpers
        sp = for small penalty radius
        mp = for medium penalty radius
        lp = for large penalty radius
    findPenalty = finds the current penalty given current towers and penalty radius

"""
def solve_dp_small(instance: Instance) -> Solution:
    return solve_dp(instance, 8)
def solve_dp_medium(instance: Instance) -> Solution:
    return solve_dp(instance, 10)
def solve_dp_big(instance: Instance) -> Solution:
    return solve_dp(instance, 14)
    

def solve_dp(instance: Instance, penalty_radius) -> Solution:
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
    #solution = recursiveHelper(cities, [], penalty_radius) #solution = [penalty, tower list]
    solution = justGreedyRecursive(cities, [], penalty_radius)
    return Solution(solution, instance)

def justGreedyRecursive(cities: List, towers: List, penalty_radius: int):
    if not cities:
        return towers
    greedySolution = findGreedyTowers(cities, towers, penalty_radius)
    penaltyTower = greedySolution[0][0]
    penaltyTowerCities = greedySolution[0][1]
    newCities = [city for city in cities if city not in penaltyTowerCities]
    towers.append(penaltyTower)
    return justGreedyRecursive(newCities, towers, penalty_radius)


def recursiveHelper(cities: List, towers: List, penalty_radius: int):
    """
    if there are no cities left then just return the current penalty
    next, find greedy tower aka tower that covers the most
    """
    
    #case where there are no cities left
    if not cities:
        return [findPenalty(towers, penalty_radius), towers]
    towersList = findGreedyTowers(cities, towers, penalty_radius)
    penaltyTower = towersList[0][0]
    penaltyTowerCities = towersList[0][1]
    lowPenaltyTower = towersList[1][0]
    lowPenaltyTowerCities = towersList[1][1]
    if penaltyTower == lowPenaltyTower:
        newCities = [city for city in cities if city not in penaltyTowerCities]
        towers.append(penaltyTower)
        return recursiveHelper(newCities, towers, penalty_radius)
    else:
        highPenaltyCities = [city for city in cities if city not in penaltyTowerCities]
        lowPenaltyCities = [city for city in cities if city not in lowPenaltyTowerCities]
        highPenaltyTowers = towers
        highPenaltyTowers.append(penaltyTower)
        lowPenaltyTowers = towers
        lowPenaltyTowers.append(lowPenaltyTower)
        greedyPath = recursiveHelper(highPenaltyCities, highPenaltyTowers, penalty_radius)
        lowPenaltyPath = recursiveHelper(lowPenaltyCities, lowPenaltyTowers, penalty_radius)
        greedyPathCost = greedyPath[0]
        lowPenaltyPathCost = lowPenaltyPath[0]
        if greedyPathCost <= lowPenaltyPathCost:
            print("Next choice is greedy because " + str(greedyPathCost) + " <= " + str(lowPenaltyPathCost))
            return greedyPath
        else:
            print("Next choice is kind")
            return lowPenaltyPath 

    #return list(penalty, tower list)
    

#returns [[greedy tower w/ penalty, [cities covered by this tower]], [greedy tower w/ min penalty, [cities covered by this tower]]]
def findGreedyTowers(cities: List, towers: List, penalty_radius: int):
    """
    find the next tower that covers the most cities
    returns a [towerP, towerMP]
    towerP = greediest tower w/ or w/o penalty
    towerMP = greediest tower w/ minimum penalty

    Psuedocode:
        get each possible tower coordinate
        for each possible tower:
            if tower exists in towers, skip
            find # of cities covered
            if not greater than towerP or towerNP move on
            if greater than, then:
                if has penalty, replace towerP
                if min penalty, replace towerMP and maybe towerP
            
    """
    #best tower regardless of penalty
    towerP = Point(None, None)

    #best tower of lowest penalty
    towerMP = Point(None, None)

    #how much each tower covers
    towerPCoverage = 0
    towerMPCoverage = 0

    #the current smallest penalty
    currentMinPenalty = 205

    #list of towers with lowest penalty
    #[[cities covered, tower point, list of cities covered]]
    minPenaltyTowers = [[None, None, None]]

    #list of cities covered by towerP and towerMP
    towerPCities = []
    towerMPCities = []

    #get list of possible towers
    possibleTowers = getPossibleTowers(cities, penalty_radius)
    #each element is (x, y) not a point
    for newTower in possibleTowers:
        #make the element a point
        #check if point is already a tower
        if newTower in towers:
            continue
        
        #get all points within radius of this tower
        towerRadius = getPenaltyRadius(newTower, penalty_radius)
        
        #keep track of the current penalty and coverage of tower
        currentTowerCoverage = 0
        currentTowerPenalty = 0

        #keep track of cities within tower radius
        currentCitiesCovered = []
        #for each point, check if it is a tower or city
        for possibleCoordinate in towerRadius:
            currentPoint = Point(possibleCoordinate[0], possibleCoordinate[1])
            #if a city, add 1 to current coverage
            if currentPoint in cities and currentPoint.distance_sq(newTower) <= 9:
                currentTowerCoverage += 1
                currentCitiesCovered.append(currentPoint)
            #if a tower, add 1 to current penalty
            if currentPoint in towers:
                currentTowerPenalty += 1
        #now that we know this points coverage and penalty, check if its the best so far

        #if this tower doesnt cover any cities, we dont care about it
        if currentTowerCoverage == 0:
            continue

        #if covers more cities than our max, make our max
        if currentTowerCoverage > towerPCoverage:
            towerPCoverage = currentTowerCoverage
            towerP = newTower
            towerPCities = currentCitiesCovered
        
        if currentTowerPenalty > currentMinPenalty:
            continue
        #if covers same cities for same min penalty, add to list
        if currentTowerPenalty == currentMinPenalty and currentTowerCoverage == towerMPCoverage:
            minPenaltyTowers.append([currentTowerCoverage, newTower, currentCitiesCovered])
            

        #if covers more cities for same min penalty, make new list
        if currentTowerPenalty == currentMinPenalty and currentTowerCoverage > towerMPCoverage:
            minPenaltyTowers.clear()
            minPenaltyTowers.append([currentTowerCoverage, newTower, currentCitiesCovered])
            towerMPCoverage = currentTowerCoverage
            
        #if has lower min penalty, make new list
        if currentTowerPenalty < currentMinPenalty:
            minPenaltyTowers.clear()
            minPenaltyTowers.append([currentTowerCoverage, newTower, currentCitiesCovered])
            towerMPCoverage = currentTowerCoverage
            currentMinPenalty = currentTowerPenalty
            
    
    #get the best tower from the min penalty towers
    towerMPCoverage = 0
    minTowerCoverage = 0
    minTowerPoint = Point
    for i in range(len(minPenaltyTowers)):
        minTower = minPenaltyTowers[i]
        minTowerCoverage = minTower[0]
        minTowerPoint = minTower[1]
        minTowerList = minTower[2]
        if minTowerCoverage is None:
            continue
        if minTowerCoverage > towerMPCoverage:
            towerMPCoverage = minTowerCoverage
            towerMP = minTowerPoint
            towerMPCities = minTowerList
    if towerMP == Point(None, None):
        
        towerMP = towerP
        towerMPCities = towerPCities
    towersList = [[towerP, towerPCities], [towerMP, towerMPCities]]
    return towersList

    """
    print everything!
    return a list of cities that are now gonna be covered as well
    print("TowerP Coverage: " + str(towerPCoverage))
    print("TowerMP Coverage: " + str(towerMPCoverage))
    print("Current Min Penalty: " + str(currentMinPenalty))
    print("TowerP: " + str(towerP))
    print("TowerMP: " + str(towerMP))
    cities = [Point(0, 2), Point(1, 2), Point(2, 2), Point(3, 2), Point(4, 2), Point(2, 0), Point(2, 1), Point(2, 3), Point(2, 4), 
    Point(29, 27), Point(28, 27), Point(26, 27), Point(25, 27), Point(27, 25), Point(27, 26), Point(27, 28), Point(27, 29)]
    """

def findPenalty(towers: List, penalty_radius: int):
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

#get list of possible coordinates within penalty radius of tower
def getPenaltyRadius(tower: Point, radius: int):
    if radius == 8:
        return sp[tower.x, tower.y]
    elif radius == 10:
        return mp[tower.x, tower.y]
    else:
        return lp[tower.x, tower.y]

#get list of possible coordinates within service of cities
def getPossibleTowers(cities: List, radius: int):
    if radius == 8:
        return getPossibleTowersHelperSS(cities)
    elif radius == 10:
        return getPossibleTowersHelperMS(cities)
    else:
        return getPossibleTowersHelperLS(cities)

def getPossibleTowersHelperSS(cities: List):
    possibleTowers = []
    for city in cities:
        radiusCoordinates = ss[city.x, city.y]
        for coordinate in radiusCoordinates:
            newPoint = Point(coordinate[0], coordinate[1])
            if newPoint not in possibleTowers:
                possibleTowers.append(newPoint)
    return possibleTowers

def getPossibleTowersHelperMS(cities: List):
    possibleTowers = []
    for city in cities:
        radiusCoordinates = ms[city.x, city.y]
        for coordinate in radiusCoordinates:
            newPoint = Point(coordinate[0], coordinate[1])
            if newPoint not in possibleTowers:
                possibleTowers.append(newPoint)
    return possibleTowers

def getPossibleTowersHelperLS(cities: List):
    possibleTowers = []
    for city in cities:
        radiusCoordinates = ls[city.x, city.y]
        for coordinate in radiusCoordinates:
            newPoint = Point(coordinate[0], coordinate[1])
            if newPoint not in possibleTowers:
                possibleTowers.append(newPoint)
    return possibleTowers

SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive,
    "solve_dp_small": solve_dp_small,
    "solve_dp_medium": solve_dp_medium,
    "solve_dp_big": solve_dp_big 
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
