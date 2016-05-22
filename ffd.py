#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import operator
from functools import reduce
import random
from shapes import shape_candidates, char_range

class ConvexSpace:
    def __init__(self, coordinates, boundaries, tag=0):
        self.coordinates = coordinates
        self.boundaries = boundaries
        self.tag = tag

    def canFit(self, boundariesToFit):
        for thisS, s in zip(self.boundaries, boundariesToFit):
            if (thisS < s):
                return False
        return True

    def contains(self, other):
        for thisCoord, thisD, otherCoord, otherD in zip(self.coordinates, self.boundaries,
                                                other.coordinates, other.boundaries):
            if (otherCoord < thisCoord or otherCoord > thisCoord + thisD
                or otherCoord + otherD > thisCoord + thisD):
                return False
        return True

    # in principle, each dimension can be cut in 2 new ConvexSpaces
    # this is unless there is no space left at any or both sides of this ConvexSpace
    def minus(self, other):
        newSpaces = []
        for (i, thisCoord), thisD, otherCoord, otherD in zip(enumerate(self.coordinates),
                                                              self.boundaries, other.coordinates, other.boundaries):
            if (otherCoord > thisCoord and otherCoord < thisCoord + thisD):
                newBoundaries = list(self.boundaries)
                newBoundaries[i] = otherCoord - thisCoord
                newSpaces.append(ConvexSpace(list(self.coordinates), newBoundaries))
            if (otherCoord + otherD > thisCoord and otherCoord + otherD < thisCoord + thisD):
                newCoordinates = list(self.coordinates)
                newCoordinates[i] = otherCoord + otherD
                newBoundaries = list(self.boundaries)
                newBoundaries[i] = (thisCoord + thisD) - (otherCoord + otherD)
                newSpaces.append(ConvexSpace(newCoordinates, newBoundaries))
        return newSpaces

    # return intersection between this space and another as
    # a list of intervals in all dimensions or an empty list (if no intersection)
    def intersection(self, other):
        intersection = []
        for thisCoord, thisD, otherCoord, otherD in zip(self.coordinates, self.boundaries,
                                                        other.coordinates, other.boundaries):
            if (otherCoord + otherD >= thisCoord and otherCoord + otherD < thisCoord + thisD):
                newCoord = max(thisCoord, otherCoord)
                intersection.append((newCoord, otherCoord + otherD - newCoord))
            elif (otherCoord >= thisCoord and otherCoord <= thisCoord + thisD):
                intersection.append((otherCoord, min(otherCoord + otherD, thisCoord + thisD)))
            else:
                intersection = []
                break
        return intersection

    # try to join two spaces over the i-th dimension
    def join(self, other, i):
        # if one space is contained in the other, return the greater one
        if (self.contains(other)):
            return [self]
        if (other.contains(self)):
            return [other]
        # check if all dimensions overlap (we can only join in that case)
        intersection = self.intersection(other)
        if (not intersection):
            return []
        # if any dimension of the intersection is 0 (except for the i-th),
        # return (can not join in that case)
        if (any(map(lambda x: False if x[0] == i else x[1] == 0, enumerate(intersection)))):
            return []

        # take full interval of both spaces in i-th dimension
        # and all other intervals as is from intersection to create joined space
        joinedCoordinates = [x[0] for x in intersection]
        joinedCoordinates[i] = min(self.coordinates[i], other.coordinates[i])
        joinedBoundaries = [x[1]-x[0] for x in intersection]
        joinedBoundaries[i] = max(self.coordinates[i] + self.boundaries[i],
                                  other.coordinates[i] + other.boundaries[i]) - joinedCoordinates[i]
        joinedSpace = ConvexSpace(joinedCoordinates, joinedBoundaries)

        # subtract the original spaces from the joined space and return all resulting spaces
        # which are not contained in the joined space already
        return [s for s in [joinedSpace] + joinedSpace.minus(self) +
        joinedSpace.minus(other) if s == joinedSpace or not joinedSpace.contains(s)]

    def __str__(self):
        return str(self.coordinates) + ' -> ' + str(self.boundaries)

class Bin:
    def __init__(self, boundaries):
        self.freelist = [ConvexSpace([0 for d in boundaries], list(boundaries))]
        self.spaces = []

    def canFit(self, boundariesToFit):
        return any([space.canFit(boundariesToFit) for space in self.freelist])

    # fit the given boundaries in this bin over the i-th dimension
    # (free spaces are sorted over all other boundaries first, then i-th dimension)
    # this assumes that canFit was called before
    def fit(self, boundariesToFit, i):
        freeSpace = min(filter(lambda s: s.canFit(boundariesToFit), self.freelist),
                    key=lambda s: [x for d, x in enumerate(s.coordinates) if d != i] + [s.coordinates[i]])
        assignedSpace = ConvexSpace(list(freeSpace.coordinates), list(boundariesToFit))
        # TODO: need to check if new free spaces can be joined with old ones (to form a greater convex space)
        # should only join over a given order of dimensions (order in which requests are being assigned)
        # TODO: add unit tests by checking if a configuration is possible (respects request sizes/alpha, no overlaps inside bin, all requests are inside bounds of bin)
        print('Fitting ' + str(boundariesToFit) + ' into ' + str(freeSpace))
        self.spaces.append(assignedSpace)
        self.freelist.remove(freeSpace)
        newFreeSpaces = freeSpace.minus(assignedSpace)
        for newSpace in newFreeSpaces:
            for space in self.freelist:
                space.join(newSpace, reversed(range(len(boundariesToFit))))
        print('New free spaces: ' + str([str(x) for x in newFreeSpaces]))
        self.freelist.extend(newFreeSpaces)

    def __str__(self):
        s = ''
        for space in sorted(self.spaces, key=lambda s: s.coordinates):
            s = s + str(space) + '\n'
        return s

def ffd(boundaries, listOfRequestSpaces):
    if (not listOfRequestSpaces): return []
    bins = [Bin(boundaries)]
    for space in sorted(listOfRequestSpaces, reverse=True):
        # get index of last flat dimension, this is the dimension over which the space needs
        # to be fitted into its bin
        firstNonFlatD = next( (i for i, d in enumerate(space) if d > 1), len(space)-1)
        isFit = False
        for b in bins:
            if (b.canFit(space)):
                b.fit(space, firstNonFlatD)
                isFit = True
        if (not isFit):
            newBin = Bin(boundaries)
            newBin.fit(space, firstNonFlatD)
            bins.append(newBin)
    return bins

def chooseCandidate(requestSize, total, shape_candidates, alpha):
    # get the minimum candidate (minimal size of the flattest candidates)
    # of all candidates that can hold this requestSize according to alpha
    # this is the candidate that leaves the biggest convex free space when fit into the bin,
    # according to the heuristic
    cutOffSize = min(int(requestSize+requestSize*alpha), total)
    return min([p for n, projs in shape_candidates.items()
               if n >= requestSize and n <= cutOffSize for p in projs])

def firstFitDecreasing(boundaries, requestSizes, shape_candidates, alpha):
    total = reduce(operator.mul, boundaries, 1)
    listOfRequestSpaces = map(lambda r:
                              chooseCandidate(r, total, shape_candidates, alpha), requestSizes)
    return ffd(boundaries, listOfRequestSpaces)

def randomRequestsFFD(boundaries):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)

    total = reduce(operator.mul, boundaries, 1)
    random.seed()
    # random amount of requests with random sizes (lower than half of total size)
    # TODO: there is a bug with this
    requestSizes = []
    for i in range(random.randint(1, 100)):
        requestSizes.append(random.randint(1, total/boundaries[0]))
    # get bin packing
    bins = firstFitDecreasing(boundaries, requestSizes, candidates, 1.0)
    # print the configuration
    print('Packed ' + str(len(requestSizes)) + ' requests into ' + str(len(bins)) + ' bins')
    print('Total size of packed requests is: ' + str(sum(requestSizes)) + ', total bins size: ' + str(len(bins)*total))
    print('Dimensions: ' + str(boundaries))
    for i, b in enumerate(bins):
        print('Bin ' + str(i) + ':\n' + str(b))

randomRequestsFFD([24, 24, 24])
