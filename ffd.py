#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import sys
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

def getPossibleSizes(dimensions, shape_candidates, alpha, isFull = True):
    # return a dict containing for every possible request size the list of request sizes
    # to look up in shape_candidates
    possibleSizes = {}
    total = reduce(operator.mul, dimensions.values(), 1)
    for n in reversed(range(1, total + 1)):
        cutOffSize = min(int(n+n*alpha), total)
        possibleSizes[n] = [r for r in shape_candidates if r >= n and r <= cutOffSize]
        if (isFull and not possibleSizes[n]):
            # make sure that we can have a shape for n by getting the minimal sized shapes
            possibleSizes[n] = possibleSizes[n+1]
    return possibleSizes

def chooseCandidateFlat(requestSize, possibleSizes, shape_candidates, alpha):
    # get the minimum candidate (minimal size of the flattest candidates)
    # of all candidates that can hold this requestSize according to possible sizes
    # this is the candidate that leaves the biggest convex free space when fit into the bin,
    # according to the heuristic
    if (not possibleSizes[requestSize]):
        return requestSize
    # return projection as list
    return list(min(p for n in possibleSizes[requestSize] for p in shape_candidates[n]))

def ffdFlat(dimensions, requestSizes, shape_candidates, alpha, isDebug=False):
    possibleSizes = getPossibleSizes(dimensions, shape_candidates, alpha)
    # always choose the flattest (and smallest) shape for each request size
    listOfRequestSpaces = map(lambda r:
                              chooseCandidateFlat(r, possibleSizes, shape_candidates, alpha), requestSizes)
    bins = []
    nFits = 0
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
        if (isDebug and not fittedBin.testPossible(False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break
        nFits = nFits + 1
        if (nFits % (int(len(requestSizes) / 20) + 1) == 0):
            print('.', end='', flush=True)
    print()
    return bins

# return difference between second lowest metric and lowest metric for all possible shapes
# for given requestSize and alpha, or 0 if there are too few possible shapes
# TODO most shapes have multiple rotations, hence, this will probably always return 0
# for request sizes where the best metric shapes have multiple rotations
# TODO this could also return the ratio ?
def bestMetricsDelta(requestSize, dimensions, possibleSizes, shape_candidates, metric):
    bestMetrics = sorted(metric(p, False, dimensions) for n in possibleSizes[requestSize]
                         for p in shape_candidates[n])
    if (len(bestMetrics) < 2):
        return 0
    return bestMetrics[1] - bestMetrics[0]

def ffdGreatestMetricDeltaFirst(dimensions, requestSizes, shape_candidates, alpha, isDebug=False):
    # sort the requested sizes by non-increasing size,
    # then by non-increasing delta of best metrics of possible shapes (sort must be stable)
    # TODO sorting request sizes by size is not exactly right as we want to consider the shapes
    # by non-increasing size. this still needs to be fixed
    possibleSizes = getPossibleSizes(dimensions, shape_candidates, alpha)
    sortedSizes = sorted(sorted(requestSizes, reverse=True), key=lambda r:
                         bestMetricsDelta(r, dimensions, possibleSizes, shape_candidates, metric_diameter),
                         reverse=True)
    bins = []
    fittedSpaces = []
    total = reduce(operator.mul, dimensions.values(), 1)
    nFits = 0
    for size in sortedSizes:
        # try to fit the different shapes in an existing bin, sorted by metric (best to worst)
        # and by size for the same metric
        isFit = False
        spaces = sorted((p for n in possibleSizes[size] for p in shape_candidates[n]),
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
        if (isDebug and not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            print('Spaces leading to this: ' + str(fittedSpaces))
            break
        nFits = nFits + 1
        if (nFits % (int(len(requestSizes) / 20) + 1) == 0):
            print('.', end='', flush=True)
    print()
    return bins

def chooseBestMetricShapes(requestSize, dimensions, possibleSizes, shape_candidates, metric):
    if (not possibleSizes[requestSize]):
        return []
    minMetric = min(metric(p, False, dimensions) for n in possibleSizes[requestSize]
                    for p in shape_candidates[n])
    # get min metric shapes, sorted by non-decreasing size
    considered = (p for n in possibleSizes[requestSize] for p in shape_candidates[n])
    return sorted(filter(lambda p: metric(p, False, dimensions) == minMetric, considered),
               key=lambda p: reduce(operator.mul, p, 1))

def ffdBestMetricAlways(dimensions, requestSizes, shape_candidates, alpha, isDebug=False):
    possibleSizes = getPossibleSizes(dimensions, shape_candidates, alpha)
    # always pick the best metric shapes and then pack using that
    # if we choose diameter as metric, the best metric of r1 cannot have greater size
    # than the best metric of r2 if r1 < r2, since there cannot be two shapes of different sizes
    # having the best diameter for the requests
    bins = []
    nFits = 0
    for r in sorted(requestSizes, reverse=True):
        isFit = False
        bestMetricShapes = chooseBestMetricShapes(r, dimensions, possibleSizes,
                                                  shape_candidates, metric_diameter)
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
        if (isDebug and not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break
        nFits = nFits + 1
        if (nFits % (int(len(requestSizes) / 20) + 1) == 0):
            print('.', end='', flush=True)
    print()
    return bins

def ffdEachBestMetricFirst(dimensions, requestSizes, shape_candidates, alpha, isDebug=False):
    bins = []
    possibleSizes = getPossibleSizes(dimensions, shape_candidates, alpha)
    lastSize = reduce(operator.mul, dimensions.values(), 1)
    nFits = 0
    # sort requests by size
    for r in sorted(requestSizes, reverse=True):
        # try any shape which is smaller than last, from best metric to worst
        # choose diameter as metric because bigger shapes should not get precedence
        # sort by size first to give precedence to smaller shapes if diameter is the same
        sortedShapes = sorted((p for n in possibleSizes[r] for p in shape_candidates[n]),
                              key=lambda p: reduce(operator.mul, p, 1))
        sortedShapes = sorted(sortedShapes,
                              key=lambda p: metric_diameter(p, False, dimensions))
        if (not sortedShapes):
            print('Cannot allocate ' + str(r) + ' since no shapes are available')
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
        if (isDebug and not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break
        nFits = nFits + 1
        if (nFits % (int(len(requestSizes) / 20) + 1) == 0):
            print('.', end='', flush=True)
    print()
    return bins

# seems to be same as greatest-delta atm. greatest-delta needs to change
def ffdNonStrictEachBestMetricFirst(dimensions, requestSizes, shape_candidates, alpha, isDebug=False):
    bins = []
    possibleSizes = getPossibleSizes(dimensions, shape_candidates, alpha)
    nFits = 0
    # sort requests by size
    for r in sorted(requestSizes, reverse=True):
        # try any shape which is smaller than last, from best metric to worst
        # choose diameter as metric because bigger shapes should not get precedence
        # sort by size first to give precedence to smaller shapes if diameter is the same
        sortedShapes = sorted((p for n in possibleSizes[r] for p in shape_candidates[n]),
                              key=lambda p: reduce(operator.mul, p, 1))
        sortedShapes = sorted(sortedShapes,
                              key=lambda p: metric_diameter(p, False, dimensions))
        if (not sortedShapes):
            print('Cannot allocate ' + str(r) + ' since no shapes are available')
            continue
        isFit = False
        for shape in sortedShapes:
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
        if (isDebug and not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break
        nFits = nFits + 1
        if (nFits % (int(len(requestSizes) / 20) + 1) == 0):
            print('.', end='', flush=True)
    print()
    return bins

# sort the shapes by metric (best to worst) and associate every shape with an ID
# this allows to sort all such shapes as a big batch and consider the shapes one by one
# rather than consider the request sizes one by one and then the shapes of each request size
def sortShapesByMetric(requestSize, dimensions, possibleSizes, shape_candidates, metric):
    bestMetricShapes = sorted((p for n in possibleSizes[requestSize]
                               for p in shape_candidates[n]),
    key=lambda p: metric(p, False, dimensions))
    if (not bestMetricShapes):
        print('Cannot allocate ' + str(requestSize) + ' since no shapes are available')
        return iter([])
    bestMetricSize = reduce(operator.mul, bestMetricShapes[0], 1)
    # associate an ID with each shape for a given request size so that it can be removed
    # accordingly when all shapes are considered
    id = uuid.uuid4()
    return ((id, p) for p in
            filter(lambda p: reduce(operator.mul, p, 1) <= bestMetricSize, bestMetricShapes))

def ffdAllBestMetricFirst(dimensions, requestSizes, shape_candidates, alpha, isDebug=False):
    bins = []
    possibleSizes = getPossibleSizes(dimensions, shape_candidates, alpha)
    consideredIds = set()
    nFits = 0
    # get all best shapes, then sort all of them by non-increasing size. since sort is stable,
    # the shapes are still in the right order of metric
    # TODO this is actually same as best-always since as soon as we choose the best metric shape,
    # we do not consider any other shape to make the packing better (if not isFit,
    # should check in remaining list whether we still have shapes with that ID left.
    # if yes, ignore this one and move on to next)
    for id, shape in sorted(((id, p) for l in map(lambda r:
        sortShapesByMetric(r, dimensions, possibleSizes, shape_candidates, metric_diameter),
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
        if (isDebug and not fittedBin.testPossible(True, False)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            break
        nFits = nFits + 1
        if (nFits % (int(len(requestSizes) / 20) + 1) == 0):
            print('.', end='', flush=True)
    print()
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

def getStats(dimensions, requestSizes, bins):
    total = reduce(operator.mul, dimensions.values(), 1)
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
            # get max coordinate in each dimension (assumes that spaces are packed from
            # min to max coordinate which is the case in all strategies)
            return reduce(operator.mul,
                          [max(min(s.coordinates[d] + s.boundaries[d], s.coordBoundaries[d])
                               for s in b.spaces) for d in range(len(b.boundaries))], 1)
    totalSubMeshSpace = sum(convexUsedSpace(b) for b in bins)
    totalSubMeshExceedingPercentage = (totalSubMeshSpace - totalRequestSize) * 100 / totalRequestSize

    averageDiameter = sum(metric_diameter(p.boundaries, False, dimensions) for b in bins
                          for p in b.spaces) / len(requestSizes)
    averageOptDiameter = sum(opt_diameter(r, False, dimensions)
                             for r in requestSizes) / len(requestSizes)

    return [totalRequestSize, nRequests, nBins, nOptBins, unusedSpace, unusedPercentage,
    totalUsedSpace, totalUsedExceedingPercentage, totalSubMeshSpace,
    totalSubMeshExceedingPercentage, averageDiameter, averageOptDiameter]

def performFFD(boundaries, requestSizes, strategy, alpha, printDetail = False):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)

    # get bin packing
    bins = strategy(dimensions, requestSizes, candidates, alpha)
    results = getStats(dimensions, requestSizes, bins)

    if (printDetail):
        # print the configuration
        print('Packing as follows:')
        for i, b in enumerate(bins):
            print('Bin ' + str(i) + ':\n' + str(b))

        printResults(boundaries, results)

    # TODO: add unit test to check if request sizes/alpha are respected
    # (need to match each allocated space with request size)

    # check if bin configuration is possible for all bins
    for b in bins:
        if ((strategy is ffdFlat and not b.testPossible(False, True)) or
            (not (strategy is ffdFlat) and not b.testPossible(True, False))):
            print(b)
            print('Not all bin configurations are possible ! Test FAILED')
            break

    return results

def randomRequestsFFD(boundaries, strategy, alpha, printDetail = False):
    total = reduce(operator.mul, boundaries, 1)
    random.seed()
    # random amount of requests with normal distributed
    # ignore values to < 0 or > 8, then scale up to the total size
    # thus, sizes > 1/8 of the torus have a fairly small probability
    requestSizes = [int(x * (total-1) / 8) + 1 for x in
    np.random.normal(0.5, 0.5, random.randint(100, 1000)) if x >= 0 and x <= 8]
    return performFFD(boundaries, requestSizes, strategy, alpha, printDetail)

def traceRequestsFFD(filename, boundaries, strategy, alpha, printDetail = False):
    file = open(filename)
    requestSizes = []
    for line in file.readlines():
        requestSizes.append(int(line.split()[0]))
    file.close()
    return performFFD(boundaries, requestSizes, strategy, alpha, printDetail)

def testStrategies(boundaries):
    # test each strategy 10 times and take average measures
    strategies = [('flat', ffdFlat), ('best-always', ffdBestMetricAlways),
    ('greatest-metric-delta', ffdGreatestMetricDeltaFirst),
    ('best-metric-first', ffdEachBestMetricFirst), ('all-best-metric-first', ffdAllBestMetricFirst)]
    for name, strategy in strategies:
        for alpha in [0.15, 1.0]:
            print('Testing strategy ' + name + ' with alpha ' + str(alpha))
            results = []
            for i in range(10):
                results.append(randomRequestsFFD(boundaries, strategy, alpha))
            npResults = np.array(results)
            avgResults = np.mean(npResults, axis=0)
            print('Average results for strategy ' + name + ' and alpha ' + str(alpha))
            printResults(boundaries, avgResults)
            print('\n')

def testStrategiesTrace(filename, boundaries):
    # test each strategy on the given sizes provided by trace
    strategies = [('flat', ffdFlat), ('best-always', ffdBestMetricAlways),
    ('greatest-metric-delta', ffdGreatestMetricDeltaFirst),
    ('non-strict-decreasing-best-metric-first', ffdNonStrictEachBestMetricFirst),
    ('best-metric-first', ffdEachBestMetricFirst), ('all-best-metric-first', ffdAllBestMetricFirst)]
    for name, strategy in strategies:
        for alpha in [0.15, 1.0]:
            print('Testing strategy ' + name + ' with alpha ' + str(alpha))
            results = traceRequestsFFD(filename, boundaries, strategy, alpha)
            printResults(boundaries, results)
            print('\n')

if (len(sys.argv) > 1):
    traceFile = sys.argv[1]
    testStrategiesTrace(traceFile, [24, 24, 24])
    sys.exit(0)

#randomRequestsFFD([24, 24, 24], ffdFlat, 0.15, True)
#testStrategies([24,24,24])
#randomRequestsFFD([24,24,24], ffdGreatestMetricDeltaFirst, 0.15, True)
#printResults([24,24,24],
#             traceRequestsFFD('request_sizes_scaled_5000.txt', [24,24,24],
#                              ffdNonStrictEachBestMetricFirst, 0.15))
performFFD([24,24,24],[56,34],ffdAllBestMetricFirst,1.0, True)
