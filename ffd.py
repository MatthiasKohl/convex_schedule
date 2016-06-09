#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import sys
import copy
import math
import operator
import itertools
import random
import numpy as np
import uuid
from functools import reduce
from shapes import shape_candidates, char_range, getPossibleSizes, candidates_without_rotations
from shapes import metric_compactness, metric_diameter, opt_diameter, metric_ad, opt_ad
from packingbin import Bin

def unitTestFlat(b):
    return b.testPossible(False, True)
def unitTestBest(b):
    return b.testPossible(True, False)
def placeFlat(b, shape, ID = None):
    firstNonFlatD = next( (i for i, d in enumerate(shape) if d > 1), len(shape)-1)
    b.fitFlat(shape, firstNonFlatD, ID)
def placeBest(b, shape, ID = None):
    b.fitBest(shape, ID)

def genericFirstFit(dimensions, requestSizes, shape_candidates, alpha, size_generator,
                    shape_generator, place, initialBins, unitTest, isDebug):
    possibleSizes = getPossibleSizes(dimensions, shape_candidates, alpha)
    bins = list(initialBins)
    nSizes = 0
    fittedSpaces = []
    lastShape = None
    maxFreeSpaces = 0
    for ID, size in size_generator(requestSizes, possibleSizes, shape_candidates):
        isFit = False
        for shape in shape_generator(size, lastShape, possibleSizes, shape_candidates):
            for b in bins:
                if (b.canFit(shape)):
                    isFit = True
                    place(b, shape, ID)
                    lastShape = shape
                    if (isDebug):
                        fittedSpaces.append(shape)
                        fittedBin = b
                    break
            if (isFit):
                break
        hasShapes = True
        if (not isFit):
            newBin = Bin(dimensions.values())
            # simply choose the first space in this case
            try:
                shape = next(shape_generator(size, lastShape, possibleSizes, shape_candidates))
                place(newBin, shape, ID)
                bins.append(newBin)
                lastShape = shape
                if (isDebug):
                    fittedSpaces.append(shape)
                    fittedBin = newBin
            except StopIteration:
                print('Cannot allocate size ' + str(size) + ' since no shapes are available')
                hasShapes = False
        if (isDebug and hasShapes and not unitTest(fittedBin)):
            print('Impossible configuration:\n' + str(fittedBin) + '\n')
            print('Spaces leading to this: ' + str(fittedSpaces))
            break
        nSizes = nSizes + 1
        if (isDebug and nSizes % (int(len(requestSizes) / 20) + 1) == 0):
            print('.', end='', flush=True)
        numFreeSpaces = sum(1 for b in bins for s in b.freelist)
        if (numFreeSpaces > maxFreeSpaces):
            maxFreeSpaces = numFreeSpaces
    if (isDebug):
        print()
    # print('The maximum number of total free spaces during the whole packing was: ' +
    #       str(maxFreeSpaces))
    return bins

def chooseCandidateFlat(requestSize, possibleSizes, shape_candidates):
    # get the minimum candidate (minimal size of the flattest candidates)
    # of all candidates that can hold this requestSize according to possible sizes
    # this is the candidate that leaves the biggest convex free space when fit into the bin,
    # according to the heuristic
    if (not possibleSizes[requestSize]):
        return requestSize
    # return projection as list
    return list(min(p for n in possibleSizes[requestSize] for p in shape_candidates[n]))

def ffdFlat(dimensions, requestSizes, shape_candidates, alpha,
            initialBins=[], isDebug = False):
    def size_gen(requestSizes, possibleSizes, shape_candidates):
        # sort sizes based on their assigned shape (flattest possible) with shapes being
        # sorted based on size from first to last dimension
        listOfRequestSpaces = map(lambda r:
                                  (r, chooseCandidateFlat(r[1], possibleSizes, shape_candidates)),
                                  requestSizes.items())
        return (r for r, shape in sorted(listOfRequestSpaces,
                                         key = lambda x: x[1] if isinstance(x[1], list) else
                                         [0 for d in dimensions], reverse = True))
    def shape_gen(size, lastShape, possibleSizes, shape_candidates):
        # always simply yield the flattest shape for each size
        yield chooseCandidateFlat(size, possibleSizes, shape_candidates)

    return genericFirstFit(dimensions, requestSizes, shape_candidates, alpha,
                           size_gen, shape_gen, placeFlat, initialBins, unitTestFlat, isDebug)

# return difference between second lowest metric and lowest metric for all possible shapes
# for given requestSize and alpha, or 0 if there are too few possible shapes
# TODO this could also return the ratio ?
def metricsVariance(requestSize, dimensions, possibleSizes, candidates, metric):
    bestMetrics = sorted(metric(p, False, dimensions)
                         for n in possibleSizes[requestSize] for p in candidates[n])
    if (len(bestMetrics) < 2):
         return 0
    return bestMetrics[1] - bestMetrics[0]

def ffGreatestMetricDeltaFirst(dimensions, requestSizes, shape_candidates, alpha,
                               initialBins=[], isDebug=False):
    # pre-calculate candidates without rotation to consider for sorting the sizes
    candidates_wor = candidates_without_rotations(shape_candidates)
    def size_gen(requestSizes, possibleSizes, shape_candidates):
        # sort the requested sizes by non-increasing size,
        # then by non-increasing delta of best metrics of possible shapes (sort must be stable)
        return (s for s in sorted(sorted(requestSizes.items(),
                                         reverse=True, key=lambda x: x[1]),
                                  key=lambda r: metricsVariance(r[1], dimensions, possibleSizes,
                                                                candidates_wor, metric_ad),
                                  reverse=True))
    def shape_gen(size, lastShape, possibleSizes, shape_candidates):
        # sort by metric (best to worst) and by size for same metric
        spaces = sorted((p for n in possibleSizes[size] for p in shape_candidates[n]),
                        key=lambda p: reduce(operator.mul, p, 1))
        return (s for s in sorted(spaces, key=lambda p: metric_ad(p, False, dimensions)))

    return genericFirstFit(dimensions, requestSizes, shape_candidates, alpha, size_gen,
                           shape_gen, placeBest, initialBins, unitTestBest, isDebug)

def chooseBestMetricShapes(requestSize, dimensions, possibleSizes, shape_candidates, metric):
    if (not possibleSizes[requestSize]):
        return []
    bestMetric = min(metric(p, False, dimensions) for n in possibleSizes[requestSize]
                    for p in shape_candidates[n])
    # get min metric shapes, sorted by non-decreasing size
    considered = (p for n in possibleSizes[requestSize] for p in shape_candidates[n])
    return (s for s in sorted(filter(lambda p:
                                     metric(p, False, dimensions) == bestMetric, considered),
    key=lambda p: reduce(operator.mul, p, 1)))

def ffdBestMetricAlways(dimensions, requestSizes, shape_candidates, alpha,
                        initialBins=[], isDebug=False):
    possibleSizes = getPossibleSizes(dimensions, shape_candidates, alpha)
    # always pick the best metric shapes and then pack using that
    # if we choose diameter/AD as metric, the best metric of r1 cannot have greater size
    # than the best metric of r2 if r1 < r2, since there cannot be two shapes of different sizes
    # having the best diameter/AD for the requests
    def size_gen(requestSizes, possibleSizes, shape_candidates):
        return (s for s in sorted(requestSizes.items(), reverse=True, key=lambda x: x[1]))
    def shape_gen(size, lastShape, possibleSizes, shape_candidates):
        return chooseBestMetricShapes(size, dimensions, possibleSizes, shape_candidates,
                                      metric_ad)

    return genericFirstFit(dimensions, requestSizes, shape_candidates, alpha, size_gen,
                           shape_gen, placeBest, initialBins, unitTestBest, isDebug)

def ffdEachBestMetricFirst(dimensions, requestSizes, shape_candidates, alpha,
                           initialBins=[], isDebug=False):
    def size_gen(requestSizes, possibleSizes, shape_candidates):
        return (s for s in sorted(requestSizes.items(), reverse=True, key=lambda x: x[1]))
    def shape_gen(size, lastShape, possibleSizes, shape_candidates):
        # try any shape which is smaller than last, from best metric to worst
        # choose diameter as metric because bigger shapes should not get precedence
        # sort by size first to give precedence to smaller shapes if diameter is the same
        lastSize = reduce(operator.mul,
                          (lastShape if lastShape != None else dimensions.values()), 1)
        considered = filter(lambda s: reduce(operator.mul, s, 1) <= lastSize,
                            (p for n in possibleSizes[size] for p in shape_candidates[n]))
        sortedShapes = sorted(considered, key=lambda p: reduce(operator.mul, p, 1))
        return (s for s in sorted(sortedShapes,
                                  key=lambda p: metric_ad(p, False, dimensions)))

    return genericFirstFit(dimensions, requestSizes, shape_candidates, alpha, size_gen,
                           shape_gen, placeBest, initialBins, unitTestBest, isDebug)

def ffEachBestMetricFirst(dimensions, requestSizes, shape_candidates, alpha,
                          initialBins=[], isDebug=False):
    def size_gen(requestSizes, possibleSizes, shape_candidates):
        return (s for s in sorted(requestSizes.items(), reverse=True, key=lambda x: x[1]))
    def shape_gen(size, lastShape, possibleSizes, shape_candidates):
        # try any shape, from best metric to worst
        # choose diameter/AD as metric because bigger shapes should not get precedence
        # sort by size first to give precedence to smaller shapes if diameter/AD is the same
        sortedShapes = sorted((p for n in possibleSizes[size] for p in shape_candidates[n]),
                              key=lambda p: reduce(operator.mul, p, 1))
        return (s for s in sorted(sortedShapes,
                                  key=lambda p: metric_ad(p, False, dimensions)))

    return genericFirstFit(dimensions, requestSizes, shape_candidates, alpha, size_gen,
                           shape_gen, placeBest, initialBins, unitTestBest, isDebug)

def printResults(boundaries, results):
    print('Dimensions: ' + str(boundaries))
    print('Packed ' + str(results[1]) + ' requests into ' + str(results[2]) + ' bins (' +
          str((results[2] - results[3]) * 100 / results[3]) + '% loss)')
    print('Total size of requests is: ' + str(results[0]) +
          ', total bins size: ' + str(results[2]*reduce(operator.mul, boundaries, 1)) +
          ', optimal lower bound (# bins): ' + str(results[3]))
    print('Packing unused space (excluding last bin): ' + str(results[4]) + ' -> ' +
          str(results[5]) + '%')

    print('Total used space: ' + str(results[6]) + ', exceeding ' +
          str(results[7]) + '% total requested space')

    print('Total convex space in bins: ' + str(results[8]) + ', exceeding ' +
          str(results[9]) + '% total requested space')
    print('Total average AD: ' + str(results[10]) + ', total optimal average AD: ' +
          str(results[11]) + ' (' + str((results[10]-results[11]) * 100 / results[11]) + '% loss)')

def getStats(dimensions, requestSizes, bins):
    total = reduce(operator.mul, dimensions.values(), 1)
    totalRequestSize = sum(requestSizes.values())
    nRequests = len(requestSizes)
    nBins = len(bins)
    nOptBins = int(math.ceil(totalRequestSize/total))

    totalUsedSpace = sum(reduce(operator.mul, s.boundaries, 1) for b in bins for s in b.spaces)
    totalUsedExceedingPercentage = (totalUsedSpace - totalRequestSize) * 100 / totalRequestSize

    # unused space excludes the last bin's free space as it might still be filled
    unusedSpace = totalUsedSpace - sum(reduce(operator.mul, s.boundaries, 1) for s in bins[-1].spaces)
    unusedSpace = (len(bins)-1) * total - unusedSpace
    unusedPercentage = unusedSpace * 100 / ((len(bins) - 1) * total)

    def convexUsedSpace(b):
            # return space used by a convex mesh around all assigned spaces in bin b
            # get max coordinate in each dimension (assumes that spaces are packed from
            # min to max coordinate which is the case in all strategies)
            # this is not the tightest convex space in a torus, but an approximation for it
            return reduce(operator.mul,
                          [max(min(s.coordinates[d] + s.boundaries[d], s.coordBoundaries[d])
                               for s in b.spaces) for d in range(len(b.boundaries))], 1)
    totalSubMeshSpace = sum(convexUsedSpace(b) for b in bins)
    totalSubMeshExceedingPercentage = (totalSubMeshSpace - totalRequestSize) * 100 / totalRequestSize

    averageMetric = sum(metric_ad(p.boundaries, False, dimensions) for b in bins
                          for p in b.spaces) / len(requestSizes)
    averageOptMetric = sum(opt_ad(r, False, dimensions)
                             for r in requestSizes.values()) / len(requestSizes)

    return [totalRequestSize, nRequests, nBins, nOptBins, unusedSpace, unusedPercentage,
    totalUsedSpace, totalUsedExceedingPercentage, totalSubMeshSpace,
    totalSubMeshExceedingPercentage, averageMetric, averageOptMetric]

def performFFD(boundaries, requestSizes, strategy, unitTest, alpha, printDetail = False):
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
        if (not unitTest(b)):
            print(b)
            print('Not all bin configurations are possible ! Test FAILED')
            break

    return results

def uniformRandomRequestsFFD(boundaries, strategy, unitTest, alpha, printDetail = False):
    total = reduce(operator.mul, boundaries, 1)
    random.seed()
    requestSizes = [random.randint(1, total) for i in range(random.randint(1000, 2000))]
    requestSizes = dict(enumerate(requestSizes))
    return performFFD(boundaries, requestSizes, strategy, unitTest, alpha, printDetail)

def randomRequestsFFD(boundaries, strategy, unitTest, alpha, printDetail = False):
    total = reduce(operator.mul, boundaries, 1)
    random.seed()
    # random amount of requests with normal distributed
    # ignore values to < 0 or > 8, then scale up to the total size
    # thus, sizes > 1/8 of the torus have a fairly small probability
    requestSizes = [int(x * (total-1) / 8) + 1 for x in
    np.random.normal(0.5, 0.5, random.randint(1000, 2000)) if x >= 0 and x <= 8]
    requestSizes = dict(enumerate(requestSizes))
    return performFFD(boundaries, requestSizes, strategy, unitTest, alpha, printDetail)

def traceRequestsFFD(filename, boundaries, strategy, unitTest, alpha, columnIdx = 0, printDetail = False):
    requestSizes = {}
    with open(filename) as infile:
        i = 0
        for line in infile:
            requestSizes[i] = int(line.split()[columnIdx])
            i = i + 1
    return performFFD(boundaries, requestSizes, strategy, unitTest, alpha, printDetail)

def testStrategies(boundaries):
    # test each strategy 10 times and take average measures
    strategies = [('flat', ffdFlat), ('best-always', ffdBestMetricAlways),
    ('greatest-metric-delta', ffGreatestMetricDeltaFirst),
    ('best-metric-first-strict-decreasing', ffdEachBestMetricFirst),
    ('best-metric-first', ffEachBestMetricFirst)]
    for name, strategy in strategies:
        unitTest = unitTestFlat if strategy is ffdFlat else unitTestBest
        for alpha in [0.15, 1.0]:
            print('Testing strategy ' + name + ' with alpha ' + str(alpha))
            results = []
            for i in range(10):
                results.append(randomRequestsFFD(boundaries, strategy, unitTest, alpha))
            npResults = np.array(results)
            avgResults = np.mean(npResults, axis=0)
            print('Average results for strategy ' + name + ' and alpha ' + str(alpha))
            printResults(boundaries, avgResults)
            print('\n')

def testStrategiesTrace(filename, boundaries):
    # test each strategy on the given sizes provided by trace
    strategies = [('flat', ffdFlat), ('best-always', ffdBestMetricAlways),
    ('greatest-metric-delta', ffGreatestMetricDeltaFirst),
    ('best-metric-first-strict-decreasing', ffdEachBestMetricFirst),
    ('best-metric-first', ffEachBestMetricFirst)]
    for name, strategy in strategies:
        unitTest = unitTestFlat if strategy is ffdFlat else unitTestBest
        for alpha in [0.15, 1.0]:
            print('Testing strategy ' + name + ' with alpha ' + str(alpha))
            results = traceRequestsFFD(filename, boundaries, strategy, unitTest, alpha)
            printResults(boundaries, results)
            print('\n')

if __name__ == '__main__':
    if (len(sys.argv) > 1):
        traceFile = sys.argv[1]
        testStrategiesTrace(traceFile, [24, 24, 24])
        sys.exit(0)

    # JUST TESTS
    #randomRequestsFFD([24, 24, 24], ffdFlat, 0.15, True)
    testStrategies([24,24,24])
    #randomRequestsFFD([24,24,24], ffdGreatestMetricDeltaFirst, 0.15, True)
    #printResults([24,24,24],
    #             traceRequestsFFD('request_sizes_scaled_5000.txt', [24,24,24],
    #                              ffdNonStrictEachBestMetricFirst, 0.15))
    #performFFD([24,24,24],[56,34],ffdAllBestMetricFirst,1.0, True)

# RANDOM RESULTS USING DIAMETER (best for alpha 0.15 is best-metric-first (non-strict))
# Testing strategy flat with alpha 0.15
# Average results for strategy flat and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1520.6 requests into 137.1 bins (11.2824675325% loss)
# Total size of requests is: 1696576.7, total bins size: 1895270.4, optimal lower bound (# bins): 123.2
# Packing unused space (excluding last bin): 90422.4 -> 4.77973815985%
# Total used space: 1796816.2, exceeding 5.90354108751% total requested space
# Total convex space in bins: 1839621.6, exceeding 8.43914309918% total requested space
# Total average diameter: 24.4112184039, total optimal average diameter: 22.1350806945 (10.282942903% loss)


# Testing strategy flat with alpha 1.0
# Average results for strategy flat and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1360.1 requests into 133.9 bins (22.8440366972% loss)
# Total size of requests is: 1500679.3, total bins size: 1851033.6, optimal lower bound (# bins): 109.0
# Packing unused space (excluding last bin): 7.8 -> 0.000648547254151%
# Total used space: 1842010.5, exceeding 22.7750648103% total requested space
# Total convex space in bins: 1842304.8, exceeding 22.796296621% total requested space
# Total average diameter: 24.7101412519, total optimal average diameter: 22.0568407387 (12.029376939% loss)


# Testing strategy best-always with alpha 0.15
# Average results for strategy best-always and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1211.4 requests into 108.0 bins (9.75609756098% loss)
# Total size of requests is: 1354085.8, total bins size: 1492992.0, optimal lower bound (# bins): 98.4
# Packing unused space (excluding last bin): 51315.4 -> 3.42471205067%
# Total used space: 1415055.4, exceeding 4.49403683137% total requested space
# Total convex space in bins: 1473272.8, exceeding 8.87119097307% total requested space
# Total average diameter: 22.8073748169, total optimal average diameter: 22.141566761 (3.00705032797% loss)


# Testing strategy best-always with alpha 1.0
# Average results for strategy best-always and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1309.1 requests into 116.5 bins (11.3766730402% loss)
# Total size of requests is: 1440809.8, total bins size: 1610496.0, optimal lower bound (# bins): 104.6
# Packing unused space (excluding last bin): 28547.8 -> 1.79782876857%
# Total used space: 1560654.6, exceeding 8.3189455437% total requested space
# Total convex space in bins: 1600387.2, exceeding 11.144724154% total requested space
# Total average diameter: 22.5468423198, total optimal average diameter: 22.0399633228 (2.29981778785% loss)


# Testing strategy greatest-metric-delta with alpha 0.15
# Average results for strategy greatest-metric-delta and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1213.4 requests into 105.5 bins (6.24370594159% loss)
# Total size of requests is: 1364486.0, total bins size: 1458432.0, optimal lower bound (# bins): 99.3
# Packing unused space (excluding last bin): 8334.5 -> 0.611348095604%
# Total used space: 1428576.8, exceeding 4.69266375074% total requested space
# Total convex space in bins: 1457971.2, exceeding 6.91878998458% total requested space
# Total average diameter: 23.1554850131, total optimal average diameter: 22.1867570415 (4.36624410597% loss)


# Testing strategy greatest-metric-delta with alpha 1.0
# Average results for strategy greatest-metric-delta and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1256.1 requests into 112.1 bins (9.47265625% loss)
# Total size of requests is: 1409311.0, total bins size: 1549670.4, optimal lower bound (# bins): 102.4
# Packing unused space (excluding last bin): 4321.1 -> 0.291463075941%
# Total used space: 1531224.7, exceeding 8.63664295094% total requested space
# Total convex space in bins: 1549670.4, exceeding 9.97386013879% total requested space
# Total average diameter: 22.869981998, total optimal average diameter: 22.1680402967 (3.16645807173% loss)


# Testing strategy best-metric-first-strict-decreasing with alpha 0.15
# Average results for strategy best-metric-first-strict-decreasing and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1352.3 requests into 117.7 bins (7.0% loss)
# Total size of requests is: 1513690.1, total bins size: 1627084.8, optimal lower bound (# bins): 110.0
# Packing unused space (excluding last bin): 24705.6 -> 1.53442699327%
# Total used space: 1575395.2, exceeding 4.08461647926% total requested space
# Total convex space in bins: 1622592.0, exceeding 7.28393059302% total requested space
# Total average diameter: 23.0565185126, total optimal average diameter: 22.0975234551 (4.33983047649% loss)


# Testing strategy best-metric-first-strict-decreasing with alpha 1.0
# Average results for strategy best-metric-first-strict-decreasing and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1224.9 requests into 110.6 bins (10.379241517% loss)
# Total size of requests is: 1376839.9, total bins size: 1528934.4, optimal lower bound (# bins): 100.2
# Packing unused space (excluding last bin): 14343.7 -> 0.95282570747%
# Total used space: 1487233.8, exceeding 8.02785309893% total requested space
# Total convex space in bins: 1526918.4, exceeding 10.9730063448% total requested space
# Total average diameter: 22.8993234717, total optimal average diameter: 22.1642169378 (3.31663661272% loss)


# Testing strategy best-metric-first with alpha 0.15
# Average results for strategy best-metric-first and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1251.1 requests into 108.3 bins (6.07247796278% loss)
# Total size of requests is: 1402687.3, total bins size: 1497139.2, optimal lower bound (# bins): 102.1
# Packing unused space (excluding last bin): 5654.1 -> 0.385727474062%
# Total used space: 1465450.0, exceeding 4.48223313409% total requested space
# Total convex space in bins: 1495008.0, exceeding 6.61159872674% total requested space
# Total average diameter: 23.1149526313, total optimal average diameter: 22.0884016748 (4.64746599418% loss)


# Testing strategy best-metric-first with alpha 1.0
# Average results for strategy best-metric-first and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1172.7 requests into 103.8 bins (9.84126984127% loss)
# Total size of requests is: 1298262.1, total bins size: 1434931.2, optimal lower bound (# bins): 94.5
# Packing unused space (excluding last bin): 3488.2 -> 0.248142490223%
# Total used space: 1410864.4, exceeding 8.64406950472% total requested space
# Total convex space in bins: 1433649.6, exceeding 10.4505492442% total requested space
# Total average diameter: 22.7771663937, total optimal average diameter: 21.9994457254 (3.53518301314% loss)


# RANDOM RESULTS USING AD (best for alpha 0.15 is best-metric-first (non-strict))
# Testing strategy flat with alpha 0.15
# Average results for strategy flat and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1296.4 requests into 117.8 bins (10.922787194% loss)
# Total size of requests is: 1459296.3, total bins size: 1628467.2, optimal lower bound (# bins): 106.2
# Packing unused space (excluding last bin): 75530.4 -> 4.6739512465%
# Total used space: 1545569.9, exceeding 5.90436216809% total requested space
# Total convex space in bins: 1583172.0, exceeding 8.53826695186% total requested space
# Total average AD: 11.0697690815, total optimal average AD: 9.5858411031 (15.4804149418% loss)


# Testing strategy flat with alpha 1.0
# Average results for strategy flat and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1299.8 requests into 128.9 bins (22.5285171103% loss)
# Total size of requests is: 1448509.1, total bins size: 1781913.6, optimal lower bound (# bins): 105.2
# Packing unused space (excluding last bin): 2.2 -> 0.000198929398148%
# Total used space: 1772762.0, exceeding 22.370399147% total requested space
# Total convex space in bins: 1773043.2, exceeding 22.391088617% total requested space
# Total average AD: 12.1954454263, total optimal average AD: 9.56107404452 (27.5530904737% loss)


# Testing strategy best-always with alpha 0.15
# Average results for strategy best-always and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1247.4 requests into 123.2 bins (21.6189536032% loss)
# Total size of requests is: 1394092.3, total bins size: 1703116.8, optimal lower bound (# bins): 101.3
# Packing unused space (excluding last bin): 272116.3 -> 16.2251571531%
# Total used space: 1424430.3, exceeding 2.17252101116% total requested space
# Total convex space in bins: 1554050.4, exceeding 11.7748507916% total requested space
# Total average AD: 9.73728689336, total optimal average AD: 9.56597898983 (1.7908036774% loss)


# Testing strategy best-always with alpha 1.0
# Average results for strategy best-always and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1345.7 requests into 133.2 bins (21.7550274223% loss)
# Total size of requests is: 1503794.2, total bins size: 1841356.8, optimal lower bound (# bins): 109.4
# Packing unused space (excluding last bin): 297708.4 -> 16.3803793418%
# Total used space: 1537551.6, exceeding 2.2407837047% total requested space
# Total convex space in bins: 1675237.6, exceeding 11.5973228798% total requested space
# Total average AD: 9.74279597447, total optimal average AD: 9.57343486556 (1.76907360092% loss)


# Testing strategy greatest-metric-delta with alpha 0.15
# Average results for strategy greatest-metric-delta and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1257.1 requests into 111.2 bins (9.23379174853% loss)
# Total size of requests is: 1397088.7, total bins size: 1537228.8, optimal lower bound (# bins): 101.8
# Packing unused space (excluding last bin): 92131.5 -> 6.19913654775%
# Total used space: 1437062.1, exceeding 2.85972617957% total requested space
# Total convex space in bins: 1527192.0, exceeding 9.49193669597% total requested space
# Total average AD: 9.95327040883, total optimal average AD: 9.54942847885 (4.22896439167% loss)


# Testing strategy greatest-metric-delta with alpha 1.0
# Average results for strategy greatest-metric-delta and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1189.6 requests into 106.8 bins (10.5590062112% loss)
# Total size of requests is: 1328991.2, total bins size: 1476403.2, optimal lower bound (# bins): 96.6
# Packing unused space (excluding last bin): 93406.0 -> 6.48917973632%
# Total used space: 1375754.0, exceeding 3.51772130363% total requested space
# Total convex space in bins: 1464336.0, exceeding 10.2954447486% total requested space
# Total average AD: 9.99084408754, total optimal average AD: 9.56417535456 (4.46111365765% loss)


# Testing strategy best-metric-first-strict-decreasing with alpha 0.15
# Average results for strategy best-metric-first-strict-decreasing and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1165.2 requests into 103.3 bins (9.89361702128% loss)
# Total size of requests is: 1294276.8, total bins size: 1428019.2, optimal lower bound (# bins): 94.0
# Packing unused space (excluding last bin): 94245.2 -> 6.88762537901%
# Total used space: 1325996.7, exceeding 2.44182359372% total requested space
# Total convex space in bins: 1417867.2, exceeding 9.80596563181% total requested space
# Total average AD: 9.91980045528, total optimal average AD: 9.53771176706 (4.00608340405% loss)


# Testing strategy best-metric-first-strict-decreasing with alpha 1.0
# Average results for strategy best-metric-first-strict-decreasing and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1170.5 requests into 104.3 bins (10.3703703704% loss)
# Total size of requests is: 1299761.0, total bins size: 1441843.2, optimal lower bound (# bins): 94.5
# Packing unused space (excluding last bin): 103759.9 -> 7.39138126027%
# Total used space: 1331348.6, exceeding 2.43988156117% total requested space
# Total convex space in bins: 1429869.6, exceeding 10.2151971901% total requested space
# Total average AD: 9.88961406185, total optimal average AD: 9.53360322603 (3.73427367787% loss)


# Testing strategy best-metric-first with alpha 0.15
# Average results for strategy best-metric-first and alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 1257.5 requests into 108.4 bins (7.11462450593% loss)
# Total size of requests is: 1393246.1, total bins size: 1498521.6, optimal lower bound (# bins): 101.2
# Packing unused space (excluding last bin): 59196.3 -> 4.16581030844%
# Total used space: 1432221.9, exceeding 2.8001057177% total requested space
# Total convex space in bins: 1497091.2, exceeding 7.69541800457% total requested space
# Total average AD: 10.0533148779, total optimal average AD: 9.54014325628 (5.37907668511% loss)


# Testing strategy best-metric-first with alpha 1.0
# Average results for strategy best-metric-first and alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 1347.8 requests into 115.9 bins (7.01754385965% loss)
# Total size of requests is: 1488314.7, total bins size: 1602201.6, optimal lower bound (# bins): 108.3
# Packing unused space (excluding last bin): 61427.6 -> 4.0365606703%
# Total used space: 1536160.1, exceeding 3.21609632458% total requested space
# Total convex space in bins: 1601798.4, exceeding 7.85065617432% total requested space
# Total average AD: 10.0656538615, total optimal average AD: 9.52130722455 (5.71714181789% loss)


# SCALED Curie RESULTS USING DIAMETER (best for alpha 0.15 is best-metric-first)
# Testing strategy flat with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 60 bins
# Total size of requests is: 773014, total bins size: 829440, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 819989, exceeding 6.076862773507337% total requested space
# Total convex space in bins: 820224, exceeding 6.107263257845266% total requested space
# Total average diameter: 7.1582, total optimal average diameter: 5.179844654595089


# Testing strategy flat with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 68 bins
# Total size of requests is: 773014, total bins size: 940032, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 936357, exceeding 21.130665162597314% total requested space
# Total convex space in bins: 936576, exceeding 21.158995826725% total requested space
# Total average diameter: 7.3336, total optimal average diameter: 5.179844654595089


# Testing strategy best-always with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 61 bins
# Total size of requests is: 773014, total bins size: 843264, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 15394 -> 1.8255255768063146%
# Total used space: 808384, exceeding 4.575596302266194% total requested space
# Total convex space in bins: 838080, exceeding 8.417182612475324% total requested space
# Total average diameter: 5.42, total optimal average diameter: 5.179844654595089


# Testing strategy best-always with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 61 bins
# Total size of requests is: 773014, total bins size: 843264, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 10633 -> 1.260933705221615%
# Total used space: 813428, exceeding 5.228107123544981% total requested space
# Total convex space in bins: 839808, exceeding 8.640723195181458% total requested space
# Total average diameter: 5.4146, total optimal average diameter: 5.179844654595089


# Testing strategy greatest-metric-delta with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 60 bins
# Total size of requests is: 773014, total bins size: 829440, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 810415, exceeding 4.838334105203787% total requested space
# Total convex space in bins: 829440, exceeding 7.29947969894465% total requested space
# Total average diameter: 5.5336, total optimal average diameter: 5.179844654595089


# Testing strategy greatest-metric-delta with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 61 bins
# Total size of requests is: 773014, total bins size: 843264, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 120 -> 0.014230418943533697%
# Total used space: 815185, exceeding 5.455399255382179% total requested space
# Total convex space in bins: 843264, exceeding 9.087804360593728% total requested space
# Total average diameter: 5.516, total optimal average diameter: 5.179844654595089


# Testing strategy best-metric-first-strict-decreasing with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 60 bins
# Total size of requests is: 773014, total bins size: 829440, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 807351, exceeding 4.441963534942446% total requested space
# Total convex space in bins: 829440, exceeding 7.29947969894465% total requested space
# Total average diameter: 5.6802, total optimal average diameter: 5.179844654595089


# Testing strategy best-metric-first-strict-decreasing with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 60 bins
# Total size of requests is: 773014, total bins size: 829440, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 809104, exceeding 4.668738211727084% total requested space
# Total convex space in bins: 829440, exceeding 7.29947969894465% total requested space
# Total average diameter: 5.6644, total optimal average diameter: 5.179844654595089


# Testing strategy best-metric-first with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 60 bins
# Total size of requests is: 773014, total bins size: 829440, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 810968, exceeding 4.909872266220275% total requested space
# Total convex space in bins: 829440, exceeding 7.29947969894465% total requested space
# Total average diameter: 5.5202, total optimal average diameter: 5.179844654595089


# Testing strategy best-metric-first with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 61 bins
# Total size of requests is: 773014, total bins size: 843264, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 814854, exceeding 5.412579849782798% total requested space
# Total convex space in bins: 843264, exceeding 9.087804360593728% total requested space
# Total average diameter: 5.5594, total optimal average diameter: 5.179844654595089


# SCALED Curie RESULTS USING AD (best for alpha 0.15 is best-metric-first)
# Testing strategy flat with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 60 bins (7.142857142857143% loss)
# Total size of requests is: 773014, total bins size: 829440, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 819989, exceeding 6.076862773507337% total requested space
# Total convex space in bins: 820224, exceeding 6.107263257845266% total requested space
# Total average AD: 3.258063196248265, total optimal average AD: 2.197729205525152 (48.24679892579142% loss)


# Testing strategy flat with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 68 bins (21.428571428571427% loss)
# Total size of requests is: 773014, total bins size: 940032, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 0 -> 0.0%
# Total used space: 936357, exceeding 21.130665162597314% total requested space
# Total convex space in bins: 936576, exceeding 21.158995826725% total requested space
# Total average AD: 3.4693344276095215, total optimal average AD: 2.197729205525152 (57.85995922006765% loss)


# Testing strategy best-always with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 63 bins (12.5% loss)
# Total size of requests is: 773014, total bins size: 870912, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 65124 -> 7.598286290322581%
# Total used space: 804024, exceeding 4.011570294975253% total requested space
# Total convex space in bins: 851328, exceeding 10.13099374655569% total requested space
# Total average AD: 2.2721416546416253, total optimal average AD: 2.197729205525152 (3.3858788848689154% loss)


# Testing strategy best-always with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 63 bins (12.5% loss)
# Total size of requests is: 773014, total bins size: 870912, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 65101 -> 7.595602785244922%
# Total used space: 804047, exceeding 4.0145456615274755% total requested space
# Total convex space in bins: 851328, exceeding 10.13099374655569% total requested space
# Total average AD: 2.271057051467021, total optimal average AD: 2.197729205525152 (3.3365278014016075% loss)


# Testing strategy greatest-metric-delta with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 65 bins (16.071428571428573% loss)
# Total size of requests is: 773014, total bins size: 898560, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 92688 -> 10.47634548611111%
# Total used space: 805872, exceeding 4.250634529258202% total requested space
# Total convex space in bins: 854208, exceeding 10.503561384399248% total requested space
# Total average AD: 2.3082213876863915, total optimal average AD: 2.197729205525152 (5.027561261117125% loss)


# Testing strategy greatest-metric-delta with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 62 bins (10.714285714285714% loss)
# Total size of requests is: 773014, total bins size: 857088, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 47144 -> 5.590657255616272%
# Total used space: 809944, exceeding 4.777403772764788% total requested space
# Total convex space in bins: 850752, exceeding 10.056480218986978% total requested space
# Total average AD: 2.311360613275595, total optimal average AD: 2.197729205525152 (5.170400769338211% loss)


# Testing strategy best-metric-first-strict-decreasing with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 61 bins (8.928571428571429% loss)
# Total size of requests is: 773014, total bins size: 843264, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 31838 -> 3.8384934413580245%
# Total used space: 804997, exceeding 4.137441236510593% total requested space
# Total convex space in bins: 843264, exceeding 9.087804360593728% total requested space
# Total average AD: 2.370801645021626, total optimal average AD: 2.197729205525152 (7.875057539453238% loss)


# Testing strategy best-metric-first-strict-decreasing with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 61 bins (8.928571428571429% loss)
# Total size of requests is: 773014, total bins size: 843264, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 31875 -> 3.8429542824074074%
# Total used space: 805020, exceeding 4.140416603062817% total requested space
# Total convex space in bins: 843264, exceeding 9.087804360593728% total requested space
# Total average AD: 2.369717041847023, total optimal average AD: 2.197729205525152 (7.82570645598597% loss)


# Testing strategy best-metric-first with alpha 0.15
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 61 bins (8.928571428571429% loss)
# Total size of requests is: 773014, total bins size: 843264, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 25677 -> 3.095703125%
# Total used space: 808607, exceeding 4.604444421446442% total requested space
# Total convex space in bins: 843264, exceeding 9.087804360593728% total requested space
# Total average AD: 2.3497211616161073, total optimal average AD: 2.197729205525152 (6.915863688248923% loss)


# Testing strategy best-metric-first with alpha 1.0
# Dimensions: [24, 24, 24]
# Packed 5000 requests into 62 bins (10.714285714285714% loss)
# Total size of requests is: 773014, total bins size: 857088, optimal lower bound (# bins): 56
# Packing unused space (excluding last bin): 33564 -> 3.980248178506375%
# Total used space: 812550, exceeding 5.114525739507952% total requested space
# Total convex space in bins: 857088, exceeding 10.876129022242806% total requested space
# Total average AD: 2.3480751130350614, total optimal average AD: 2.197729205525152 (6.840965990347485% loss)

