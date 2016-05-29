#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import math
import operator
import itertools
import random
import numpy as np
import uuid
from functools import reduce
from shapes import shape_candidates, char_range
from shapes import metric_compactness, metric_diameter, opt_diameter
from packingbin import Bin

def chooseCandidateFlat(requestSize, total, shape_candidates, alpha):
    # get the minimum candidate (minimal size of the flattest candidates)
    # of all candidates that can hold this requestSize according to alpha
    # this is the candidate that leaves the biggest convex free space when fit into the bin,
    # according to the heuristic
    cutOffSize = min(int(requestSize+requestSize*alpha), total)
    considered = [p for n, projs in shape_candidates.items()
               if n >= requestSize and n <= cutOffSize for p in projs]
    if (not considered):
        return requestSize
    return list(min(considered)) # return projection as list of values

def ffdFlat(dimensions, requestSizes, shape_candidates, alpha):
    total = reduce(operator.mul, dimensions.values(), 1)
    # always choose the flattest (and smallest) shape for each request size
    listOfRequestSpaces = map(lambda r:
                              chooseCandidateFlat(r, total, shape_candidates, alpha), requestSizes)
    bins = []
    # sort requested shapes by size and do a first-fit over the bins
    for space in sorted(listOfRequestSpaces, reverse=True,
                        key=lambda p: p if isinstance(p, list) else [0 for d in dimensions]):
        if (not isinstance(space, list)):
            print('Cannot allocate size ' + str(space) + ' since no shapes are available')
            continue
        # get index of last flat dimension, this is the dimension over which the space needs
        # to be fitted into its bin
        firstNonFlatD = next( (i for i, d in enumerate(space) if d > 1), len(space)-1)
        isFit = False
        for b in bins:
            if (b.canFit(space)):
                b.fitFlat(space, firstNonFlatD)
                fittedBin = b
                isFit = True
                break
        if (not isFit):
            newBin = Bin(dimensions.values())
            newBin.fitFlat(space, firstNonFlatD)
            fittedBin = newBin
            bins.append(newBin)
        if (not fittedBin.testPossible(False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break

    return bins

# return difference between second lowest metric and lowest metric for all possible shapes
# for given requestSize and alpha, or 0 if there are too few possible shapes
# TODO this could also return the ratio ?
def bestMetricsDelta(requestSize, dimensions, shape_candidates, alpha, metric):
    total = reduce(operator.mul, dimensions.values(), 1)
    cutOffSize = min(int(requestSize+requestSize*alpha), total)
    bestMetrics = sorted([metric(p, False, dimensions) for n, projs in shape_candidates.items()
    if n >= requestSize and n <= cutOffSize for p in projs])
    if (len(bestMetrics) < 2):
        return 0
    return bestMetrics[1] - bestMetrics[0]

def ffdGreatestMetricDeltaFirst(dimensions, requestSizes, shape_candidates, alpha):
    # sort the requested sizes by non-increasing size,
    # then by non-increasing delta of best metrics of possible shapes (sort must be stable)
    # TODO sorting request sizes by size is not exactly right as we want to consider the shapes
    # by non-increasing size. this still needs to be fixed
    sortedSizes = sorted(sorted(requestSizes, reverse=True), key=lambda r:
                         bestMetricsDelta(r, dimensions, shape_candidates, alpha, metric_diameter),
                         reverse=True)
    bestMetricShapes = (chooseBestMetricShape(r, dimensions, shape_candidates,
                                              alpha, metric_diameter) for r in requestSizes)
    bins = []
    fittedSpaces = []
    total = reduce(operator.mul, dimensions.values(), 1)
    for size in sortedSizes:
        # try to fit the different shapes in an existing bin, sorted by metric (best to worst)
        # and by size for the same metric
        isFit = False
        cutOffSize = min(int(size+size*alpha), total)
        spaces = sorted((p for n, projs in shape_candidates.items()
                        if n >= size and n <= cutOffSize for p in projs),
        key=lambda p: reduce(operator.mul, p, 1))
        spaces = sorted(spaces, key=lambda p:
                            metric_diameter(p, False, dimensions))
        if (not spaces):
            print('Cannot allocate size ' + str(size) + ' since no shapes are available')
            continue # ignore request sizes that cannot be placed
        for space in spaces:
            for b in bins:
                if (b.canFit(space)):
                    b.fitBest(space)
                    fittedSpaces.append(space)
                    fittedBin = b
                    isFit = True
                    break
            if (isFit):
                break
        if (not isFit):
            newBin = Bin(dimensions.values())
            # simply choose the best space in this case
            newBin.fitBest(spaces[0])
            fittedSpaces.append(spaces[0])
            bins.append(newBin)
            fittedBin = newBin
        if (not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            print('Spaces leading to this: ' + str(fittedSpaces))
            break

    return bins

def chooseBestMetricShapes(requestSize, dimensions, shape_candidates, alpha, metric):
    total = reduce(operator.mul, dimensions.values(), 1)
    cutOffSize = min(int(requestSize+requestSize*alpha), total)
    considered = [p for n, projs in shape_candidates.items()
               if n >= requestSize and n <= cutOffSize for p in projs]
    if (not considered):
        return []
    minMetric = min(metric(p, False, dimensions) for p in considered)
    # get min metric shapes, sorted by non-decreasing size
    return sorted(filter(lambda p: metric(p, False, dimensions) == minMetric, considered),
               key=lambda p: reduce(operator.mul, p, 1))

def ffdBestMetricAlways(dimensions, requestSizes, shape_candidates, alpha):
    # always pick the best metric shapes and then pack using that
    # if we choose diameter as metric, the best metric of r1 cannot have greater size
    # than the best metric of r2 if r1 < r2, since there cannot be two shapes of different sizes
    # having the best diameter for the requests
    bins = []
    for r in sorted(requestSizes, reverse=True):
        isFit = False
        bestMetricShapes = chooseBestMetricShapes(r, dimensions, shape_candidates, alpha, metric_diameter)
        if (not bestMetricShapes):
            print('Cannot allocate ' + str(r) + ' since no shapes are available')
            continue
        for shape in bestMetricShapes:
            for b in bins:
                if (b.canFit(shape)):
                    b.fitBest(shape)
                    fittedBin = b
                    isFit = True
                    break
            if (isFit):
                break
        if (not isFit):
            newBin = Bin(dimensions.values())
            newBin.fitBest(bestMetricShapes[0])
            bins.append(newBin)
            fittedBin = newBin
        if (not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break

    return bins

def ffdEachBestMetricFirst(dimensions, requestSizes, shape_candidates, alpha):
    bins = []
    total = reduce(operator.mul, dimensions.values(), 1)
    lastSize = total
    # sort requests by size
    for r in sorted(requestSizes, reverse=True):
        # try any shape which is smaller than last, from best metric to worst
        cutOffSize = min(int(r+r*alpha), total)
        # choose diameter as metric because bigger shapes should not get precedence
        # sort by size first to give precedence to smaller shapes if diameter is the same
        sortedShapes = sorted((p for n, projs in shape_candidates.items()
                               if n >= r and n <= cutOffSize for p in projs),
        key=lambda p: reduce(operator.mul, p, 1))
        sortedShapes = sorted(sortedShapes, key=lambda p: metric_diameter(p, False, dimensions))
        if (not sortedShapes):
            print('Cannot allocated ' + str(r) + ' since no shapes are available')
            continue
        isFit = False
        for shape in sortedShapes:
            if (reduce(operator.mul, shape, 1) > lastSize):
                continue # make sure that this shape is not bigger than last
            for b in bins:
                if (b.canFit(shape)):
                    b.fitBest(shape)
                    fittedBin, fittedShape = (b, shape)
                    isFit = True
                    break
            if (isFit):
                break
        if (not isFit):
            newBin = Bin(dimensions.values())
            # take the best metric shape in this case
            newBin.fitBest(sortedShapes[0])
            bins.append(newBin)
            fittedBin, fittedShape = (newBin, sortedShapes[0])
        lastSize = reduce(operator.mul, fittedShape, 1)
        if (not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break

    return bins

def chooseBestShapes(requestSize, dimensions, shape_candidates, alpha, metric):
    total = reduce(operator.mul, dimensions.values(), 1)
    cutOffSize = min(int(requestSize+requestSize*alpha), total)
    bestMetricShapes = sorted([p for n, projs in shape_candidates.items()
           if n >= requestSize and n <= cutOffSize for p in projs], key=lambda p:
            metric(p, False, dimensions))
    if (not bestMetricShapes):
        print('Cannot allocate ' + str(requestSize) + ' since no shapes are available')
        return iter([])
    bestMetricSize = reduce(operator.mul, bestMetricShapes[0], 1)
    # associate an ID with each shape for a given request size so that it can be removed
    # accordingly when all shapes are considered
    id = uuid.uuid4()
    return ((id, p) for p in
            filter(lambda p: reduce(operator.mul, p, 1) <= bestMetricSize, bestMetricShapes))

def ffdAllBestMetricFirst(dimensions, requestSizes, shape_candidates, alpha):
    bins = []
    total = reduce(operator.mul, dimensions.values(), 1)
    consideredIds = set()
    # get all best shapes, then sort all of them by non-increasing size. since sort is stable,
    # the shapes are still in the right order of metric
    for id, shape in sorted(((id, p) for l in map(lambda r:
        chooseBestShapes(r, dimensions, shape_candidates, alpha, metric_diameter),
        requestSizes) for id, p in l), reverse=True, key=lambda x: reduce(operator.mul, x[1], 1)):
        if (id in consideredIds):
            continue
        # try to fit the shape
        isFit = False
        for b in bins:
            if (b.canFit(shape)):
                b.fitBest(shape)
                fittedBin = b
                isFit = True
                break
        if (not isFit):
            newBin = Bin(dimensions.values())
            newBin.fitBest(shape)
            bins.append(newBin)
            fittedBin = newBin
        consideredIds.add(id)
        if (not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break

    return bins

def printResults(boundaries, results):
    print('Dimensions: ' + str(boundaries))
    print('Packed ' + str(results[1]) + ' requests into ' + str(results[2]) + ' bins')
    print('Total size of requests is: ' + str(results[0]) +
          ', total bins size: ' + str(results[2]*reduce(operator.mul, boundaries, 1)) +
          ', optimal lower bound (# bins): ' + str(results[3]))
    print('Packing unused space (excluding last bin): ' + str(results[4]) + ' -> ' +
          str(results[5]) + '%')

    print('Total used space: ' + str(results[6]) + ', exceeding ' +
          str(results[7]) + '% total requested space')

    print('Total convex space in bins: ' + str(results[8]) + ', exceeding ' +
          str(results[9]) + '% total requested space')
    print('Total average diameter: ' + str(results[10]) + ', total optimal average diameter: ' +
          str(results[11]))

def randomRequestsFFD(boundaries, strategy, alpha, printDetail = False):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)

    total = reduce(operator.mul, boundaries, 1)
    random.seed()
    # random amount of requests with normal distributed
    # ignore values to < 0 or > 8, then scale up to the total size
    # thus, sizes > 1/8 of the torus have a fairly small probability
    requestSizes = [int(x * (total-1) / 8) + 1 for x in
    np.random.normal(0.5, 0.5, random.randint(100, 1000)) if x >= 0 and x <= 8]

    # get flat bin packing
    bins = strategy(dimensions, requestSizes, candidates, alpha)

    totalRequestSize = sum(requestSizes)
    nRequests = len(requestSizes)
    nBins = len(bins)
    nOptBins = int(math.ceil(totalRequestSize/total))

    # unused space excludes the last bin's free space as it might still be filled
    unusedSpace = sum(map(lambda b:
                          sum(reduce(operator.mul, s.boundaries, 1) for s in b.freelist), bins[:-1]))
    unusedPercentage = unusedSpace * 100 / (len(bins) * total)

    totalUsedSpace = sum(reduce(operator.mul, s.boundaries, 1) for b in bins for s in b.spaces)
    totalUsedExceedingPercentage = (totalUsedSpace - totalRequestSize) * 100 / totalRequestSize

    def convexUsedSpace(b):
            # return space used by a convex mesh around all assigned spaces in bin b
            # get max coordinate in each dimension
            return reduce(operator.mul,
                          [max(s.coordinates[d] + s.boundaries[d] for s in b.spaces)
                          for d in range(len(b.boundaries))], 1)
    totalSubMeshSpace = sum(convexUsedSpace(b) for b in bins)
    totalSubMeshExceedingPercentage = (totalSubMeshSpace - totalRequestSize) * 100 / totalRequestSize

    averageDiameter = sum(metric_diameter(p.boundaries, False, dimensions) for b in bins
                          for p in b.spaces) / len(requestSizes)
    averageOptDiameter = sum(opt_diameter(r, False, dimensions)
                             for r in requestSizes) / len(requestSizes)

    results = [totalRequestSize, nRequests, nBins, nOptBins, unusedSpace, unusedPercentage,
    totalUsedSpace, totalUsedExceedingPercentage, totalSubMeshSpace,
    totalSubMeshExceedingPercentage, averageDiameter, averageOptDiameter]

    if (printDetail):
        # print the configuration
        print('Packing as follows:')
        for i, b in enumerate(bins):
            print('Bin ' + str(i) + ':\n' + str(b))

        printResults(boundaries, results)

    # TODO: add unit test to check if request sizes/alpha are respected
    # (need to match each allocated space with request size)

    # check if bin configuration is possible for all bins
    if ((strategy is ffdFlat and not all(b.testPossible(False, True) for b in bins)) or
        not (strategy is ffdFlat) and not all(b.testPossible(True, False) for b in bins)):
        print('Not all bin configurations are possible ! Test FAILED')

    return results

def testStrategies(boundaries):
    # test each strategy 10 times and take average measures
    for name, strategy in [('flat', ffdFlat), ('best-always', ffdBestMetricAlways),
    ('greatest-metric-delta', ffdGreatestMetricDeltaFirst),
    ('best-metric-first', ffdEachBestMetricFirst), ('all-best-metric-first', ffdAllBestMetricFirst)]:
        for alpha in [0.15, 1.0]:
            results = []
            for i in range(10):
                results.append(randomRequestsFFD(boundaries, strategy, alpha))
            npResults = np.array(results)
            avgResults = np.mean(npResults, axis=0)
            print('Average results for strategy ' + name + ' and alpha ' + str(alpha))
            printResults(boundaries, avgResults)
            print('\n')

# TODO fix bugs here
randomRequestsFFD([24, 24, 24], ffdFlat, 1.0, True)
