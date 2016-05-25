#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import math
import operator
import itertools
import random
from functools import reduce
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

    def isIntersecting(self, other):
        # check if for all dimensions, one of the other coordinates lies inside our coordinates
        for thisCoord, thisD, otherCoord, otherD in zip(self.coordinates, self.boundaries,
                                                        other.coordinates, other.boundaries):
            if ((otherCoord + otherD > thisCoord and otherCoord + otherD < thisCoord + thisD)
                or (otherCoord > thisCoord and otherCoord < thisCoord + thisD)):
                continue
            return False
        return True

    # return intersection between this space and another as
    # a list of intervals in all dimensions or an empty list (if no intersection)
    def intersectionIntervals(self, other):
        intervals = []
        for thisCoord, thisD, otherCoord, otherD in zip(self.coordinates, self.boundaries,
                                                        other.coordinates, other.boundaries):
            if (otherCoord + otherD >= thisCoord and otherCoord + otherD <= thisCoord + thisD):
                newCoord = max(thisCoord, otherCoord)
                intervals.append((newCoord, otherCoord + otherD - newCoord))
            elif (otherCoord >= thisCoord and otherCoord <= thisCoord + thisD):
                newMaxCoord = min(thisCoord + thisD, otherCoord + otherD)
                intervals.append((otherCoord, newMaxCoord - otherCoord))
            else:
                intervals = []
                break
        return intervals

    # return the intersection of this space and another as a new space
    # if no intersection possible, the new space has dimensions 0
    def intersection(self, other):
        intervals = self.intersectionIntervals(other)
        if (not intervals):
            return ConvexSpace([0 for x in self.coordinates], [0 for x in self.boundaries])
        return ConvexSpace([x[0] for x in intervals], [x[1] for x in intervals])

    # try to join two spaces over the i-th dimension
    def join(self, other, i):
        # if one space is contained in the other, return the greater one
        if (self.contains(other)):
            return [self]
        if (other.contains(self)):
            return [other]
        # check if all dimensions overlap (we can only join in that case)
        intervals = self.intersectionIntervals(other)
        #print('Intersection intervals: ' + str(intervals))
        if (not intervals):
            return []
        # if any dimension of the intervals is 0 (except for the i-th),
        # return (can not join in that case)
        if (any([d != i and s == 0 for d, (c, s) in enumerate(intervals)])):
            return []

        # take full interval of both spaces in i-th dimension
        # and all other intervals as is from intersection to create joined space
        joinedCoordinates = [x[0] for x in intervals]
        joinedCoordinates[i] = min(self.coordinates[i], other.coordinates[i])
        joinedBoundaries = [x[1] for x in intervals]
        joinedBoundaries[i] = max(self.coordinates[i] + self.boundaries[i],
                                  other.coordinates[i] + other.boundaries[i]) - joinedCoordinates[i]
        joinedSpace = ConvexSpace(joinedCoordinates, joinedBoundaries)

        # subtract the intersection from the original spaces and return all resulting spaces
        # which are not contained in the joined space already
        return [s for s in [joinedSpace] + self.minus(self.intersection(joinedSpace)) +
        other.minus(other.intersection(joinedSpace))
        if s == joinedSpace or not joinedSpace.contains(s)]

    # update this space to contain the other space, if it is adjacent or overlapping in one
    # dimension and all other dimensions overlap entirely
    # return True if the update was possible, False if not
    def joinAdjacent(self, other, doJoin = True):
        if (self.contains(other)):
            return True
        if (other.contains(self)):
            if (doJoin):
                self.coordinates = other.coordinates
                self.boundaries = other.boundaries
            return True
        intervals = self.intersectionIntervals(other)
        if (not intervals):
            return False
        # get the dimension over which we join, all others have to overlap completely
        i = -1
        for d in range(len(intervals)):
            if (intervals[d][1] != self.boundaries[d] or intervals[d][1] != other.boundaries[d]):
                if (i >= 0):
                    return False
                i = d
        if (doJoin):
            self.coordinates[i] = min(self.coordinates[i], other.coordinates[i])
            self.boundaries[i] = max(self.coordinates[i] + self.boundaries[i],
                                      other.coordinates[i] + other.boundaries[i]) - self.coordinates[i]
        return True

    def __str__(self):
        return str(self.coordinates) + ' -> ' + str(self.boundaries)

class Bin:
    def __init__(self, boundaries):
        self.freelist = [ConvexSpace([0 for d in boundaries], list(boundaries))]
        self.spaces = []

    def canFit(self, boundariesToFit):
        return any([space.canFit(boundariesToFit) for space in self.freelist])

    # choose the best join for a flat fit
    # simply try all possible joins and choose the one which leaves the biggest sizes
    # according to reverse order dimensions, this is the biggest space in the order of
    # dimensions
    # returns the list of joined spaces
    def chooseBestJoinFlat(space1, space2):
        joinedSpaces = []
        biggestSizes = max(tuple(reversed(s.boundaries)) for s in [space1, space2])
        for d in range(len(space1.boundaries)):
            spaces = space1.join(space2, d)
            if (not spaces):
                continue
            maxSize = max(tuple(reversed(s.boundaries)) for s in spaces)
            if (maxSize >= biggestSizes):
                # if a the biggest joined space has the same size than another, take it
                # (since the other could be the initial which should always be overridden if possible)
                biggestSizes = maxSize
                joinedSpaces = spaces
        return joinedSpaces


    # fit the given boundaries in this bin over the i-th dimension (boundariesToFit should
    # be flat up to that dimension), then leave the biggest possible spaces in the order
    # of dimensions
    # (free spaces are sorted over all other dimensions first, then i-th dimension)
    # this assumes that canFit was called before
    def fitFlat(self, boundariesToFit, i):
        freeSpace = min(filter(lambda s: s.canFit(boundariesToFit), self.freelist),
                    key=lambda s: [x for d, x in enumerate(s.coordinates) if d != i] + [s.coordinates[i]])
        assignedSpace = ConvexSpace(list(freeSpace.coordinates), list(boundariesToFit))
        #print('Fitting ' + str(boundariesToFit) + ' into ' + str(freeSpace))
        self.spaces.append(assignedSpace)
        self.freelist.remove(freeSpace)
        newFreeSpaces = set(freeSpace.minus(assignedSpace))

        # TODO this doesn't work. need to come up with a better idea of how to join
        # multiple spaces together (for example, for every pair, choose the dimension which
        # gives the biggest sizes according to reverse order of dimensions)
        # while it's possible to join any pair of new free spaces, replace the pair with
        # the joined spaces. join over all dimensions except first one (never needed)
        # and in reverse order
        # tryToJoin = len(newFreeSpaces) >= 2
        # while (tryToJoin):
        #     for d in reversed(range(1, len(boundariesToFit))):
        #         tryToJoin = False
        #         for idx1 in range(len(newFreeSpaces)):
        #             for idx2 in range(idx1 + 1, len(newFreeSpaces)):
        #                 fs1 = newFreeSpaces[idx1]
        #                 fs2 = newFreeSpaces[idx2]
        #                 print('Trying to join ' + str(fs1) + ' and ' + str(fs2) + ' over dimension ' + str(d))
        #                 joinedSpaces = fs1.join(fs2, d)
        #                 if (joinedSpaces):
        #                     print('Joined spaces: ' + str([str(x) for x in joinedSpaces]))
        #                     newFreeSpaces.remove(fs1)
        #                     newFreeSpaces.remove(fs2)
        #                     newFreeSpaces.extend(joinedSpaces)
        #                     tryToJoin = True
        #                     break
        #             if (tryToJoin):
        #                 break
        #         if (tryToJoin):
        #             break # start over with first dimension

        # join each pair of spaces together (according to best flat join)
        # until there is no update anymore (no pair of spaces could be joined)
        # TODO: this is not an efficient way of doing it since when a given space cannot
        # be joined to all other spaces, it doesn't need to be reconsidered anymore
        # a fast implementation would store a table with the different possible configurations
        # and do the entire loop (all joins) in one step based on that
        isUpdated = True
        while(isUpdated):
            isUpdated = False
            for fs1, fs2 in itertools.combinations(sorted(newFreeSpaces,
                                                          key=lambda s: tuple(reversed(s.boundaries)),
                                                          reverse = True), 2):
                #print('Trying to join ' + str(fs1) + ' and ' + str(fs2))
                joinedSpaces = Bin.chooseBestJoinFlat(fs1, fs2)
                if (joinedSpaces):
                    #print('Joined spaces: ' + str([str(x) for x in joinedSpaces]))
                    newFreeSpaces = newFreeSpaces.difference({fs1, fs2}).union(set(joinedSpaces))
                    isUpdated = True
                    break

        #print('New free spaces: ' + str([str(x) for x in newFreeSpaces]))
        self.freelist.extend(newFreeSpaces)

    # fit the given boundaries at the best location according to caving degree
    # i.e. always fit at a corner and try to fill all dimensions as much as possible
    def fitBest(self, boundariesToFit):
        chosenFS = None
        chosenFSCoords = None
        for space in filter(lambda s: s.canFit(boundariesToFit), self.freelist):
            # TODO
            chosenFS = space
            chosenFSCoords = list(space.coordinates)
        assignedSpace = ConvexSpace(chosenFSCoords, list(boundariesToFit))
        self.spaces.append(assignedSpace)

        # make sure that the assigned space is subtracted from all existing FS
        removedFS =[]
        addedFS = []
        for space in self.freelist:
            if (space.isIntersecting(assignedSpace)):
                removedFS.append(space)
                addedFS.extend(space.minus(assignedSpace))
        for space in removedFS:
            self.freelist.remove(space)
        self.freelist.extend(addedFS)

        # make sure that any adjacent FS are joined
        removedFS = []
        for fs1, fs2 in itertools.combinations(self.freelist, 2):
            if (s1.joinAdjacent(s2)):
                removedFS.append(s2)
        for space in removedFS:
            self.freelist.remove(space)

    # simple test to check if a bin configuration is possible according to given boundaries
    # it can be specified that free spaces must not overlap (usually, they can overlap
    # between each other, not with other spaces)
    def testPossible(self, boundaries, allowOverlappingFS = True, allowAdjFS = True):
        # check that no space goes outside of given boundaries
        allSpaces = self.freelist + self.spaces
        boundarySpace = ConvexSpace([0 for x in boundaries], boundaries)
        if (not all(boundarySpace.contains(s) for s in allSpaces)):
            return False
        # check that there is no intersection between any two spaces
        overlapSpaces = self.spaces if allowOverlappingFS else allSpaces
        for s1, s2 in itertools.combinations(overlapSpaces, 2):
            if (s1.isIntersecting(s2)):
                return False
        # if free spaces can overlap, still need to check that they don't overlap with assigned
        if (allowOverlappingFS):
            for s1, s2 in itertools.product(self.freelist, self.spaces):
                if (s1.isIntersecting(s2)):
                    return False
        if (not allowAdjFS):
            # there must not be any two adjacent free spaces
            for s1, s2 in itertools.combinations(self.freelist, 2):
                if (s1.joinAdjacent(s2, False)):
                    return False
        return True

    def __str__(self):
        s = ''
        for space in sorted(self.spaces, key=lambda s: s.coordinates):
            s = s + str(space) + '\n'
        s = s + 'Free list:\n'
        for space in sorted(self.freelist, key=lambda s: s.coordinates):
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
                b.fitFlat(space, firstNonFlatD)
                fittedBin = b
                isFit = True
        if (not isFit):
            newBin = Bin(boundaries)
            newBin.fitFlat(space, firstNonFlatD)
            fittedBin = newBin
            bins.append(newBin)
        if (not fittedBin.testPossible(boundaries, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break

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
    requestSizes = []
    for i in range(random.randint(1, 1000)):
        requestSizes.append(random.randint(1, total/boundaries[0]))

    # get bin packing
    bins = firstFitDecreasing(boundaries, requestSizes, candidates, 1.0)

    # print the configuration
    print('Packing as follows:')
    for i, b in enumerate(bins):
        print('Bin ' + str(i) + ':\n' + str(b))
    # unused space excludes the last bin's free space as it might still be filled
    unusedSpace = sum(map(lambda b:
                          sum(reduce(operator.mul, s.boundaries, 1) for s in b.freelist), bins[:-1]))
    print('Dimensions: ' + str(boundaries))
    print('Packed ' + str(len(requestSizes)) + ' requests into ' + str(len(bins)) + ' bins')
    print('Total size of packed requests is: ' + str(sum(requestSizes)) +
          ', total bins size: ' + str(len(bins)*total) + ', optimal lower bound (# bins): ' +
          str(int(math.ceil(sum(requestSizes)/total))))
    unusedPercentage = unusedSpace * 100 / (len(bins) * total)
    print('Packing unused space: ' + str(unusedSpace) + ' -> ' + str(unusedPercentage) + '%')

    # TODO: add unit test to check if request sizes/alpha are respected
    # (need to match each allocated space with request size)

    # check if bin configuration is possible for all bins
    if (not all(b.testPossible(boundaries, False) for b in bins)):
        print('Not all bin configurations are possible ! Test FAILED')


randomRequestsFFD([24, 24, 24])
