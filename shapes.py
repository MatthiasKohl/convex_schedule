#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import operator
from collections import namedtuple
from functools import reduce
from subprocess import Popen, PIPE
from math import isinf
import itertools

TMP_SHAPES = 'shapes'
TMP_EXT = '.tmp'
TMP_TORUS = 'torus'
TMP_GRID = 'grid'
TMP_SEP = '_'

def char_range(c1, numChars):
    return [chr(c) for c in range(ord(c1), ord(c1)+numChars)]

def process_candidates(s, dimensions):
    Projection = namedtuple('Projection', dimensions.keys())
    shape_candidates = {}
    for line in s.split('\n'):
        factors = line.split(',')
        if (len(factors) <= 2): continue
        n = int(factors[0])
        shape_candidates[n] = shape_candidates.get(n, []) +\
        [Projection(**{bin: int(factor) for bin, factor in zip(dimensions.keys(), factors[1:])})]
    return shape_candidates

def shape_candidates(extProc, isGrid, dimensions):
    tmp_top = TMP_TORUS
    if (isGrid): tmp_top = TMP_GRID
    tmpFileName = TMP_SHAPES + TMP_SEP + tmp_top + TMP_SEP + 'x'.join(map(str, dimensions.values())) + TMP_EXT
    result = {}
    with open(tmpFileName, 'r') as tmpFile:
        result = process_candidates(tmpFile.read(), dimensions)

    if (result):
        return result

    args = [extProc, '1' if isGrid else '0']
    args.extend(map(str, dimensions.values()))
    proc = Popen(args, stdout=PIPE)
    out, err = proc.communicate()
    exitcode = proc.returncode

    if (exitcode != 0):
        print('ERROR in external program ' + extProc + ' -> exit code: ' + str(exitcode))
        return result

    outstring = out.decode('ascii')
    with open(tmpFileName, 'w') as tmpFile:
        tmpFile.write(outstring)

    return process_candidates(outstring, dimensions)

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

# returns a deep copy of the given candidates dictionary without the shapes that are
# rotations of other shapes already in the dictionary
def candidates_without_rotations(shape_candidates):
    candidates = copy.deepcopy(shape_candidates)
    for n in candidates:
        while (True):
            hasRotations = False
            for p1, p2 in itertools.combinations(candidates[n], 2):
                if ({s for s in p1} == {s for s in p2}):
                    hasRotations = True
                    break
            if (not hasRotations):
                break
            candidates[n].remove(p1)
    return candidates

def max_potential_shape(shapes_candidates, isGrid, dimensions, minSize=1, cutOffPercentage=float('inf'), metric=lambda x, y, z : 1):
    cutOffSize = int(minSize*(1+cutOffPercentage)) if not isinf(cutOffPercentage) else len(shapes_candidates)
    shapes_metric = map(lambda x: metric(x[1], isGrid, dimensions),
                        [(n,p) for n, projection in shapes_candidates.items() if n >= minSize and n <= cutOffSize for p in projection])
    return max(shapes_metric, default=0)

def num_potential_shapes(shapes_candidates, minSize=1, cutOffPercentage=float('inf')):
    cutOffSize = int(minSize*(1+cutOffPercentage)) if not isinf(cutOffPercentage) else len(shapes_candidates)
    return sum(1 for n, proj in shapes_candidates.items() if n >= minSize and n <= cutOffSize for p in proj)

def potential_shapes(shapes_candidates, minSize=1, cutOffPercentage=float('inf'), maxShapes=-1, metric=lambda p : 1):
    cutOffSize = int(minSize*(1+cutOffPercentage)) if not isinf(cutOffPercentage) else len(shapes_candidates)
    # filter by size
    shapes = {
        n: proj
        for n, proj in shape_candidates.items()
        if n >= minSize and n <= cutOffSize
    }

    if (maxShapes >= 0):
        # sort the projections and take only the ones with smallest metric (over all dimensions)
        sortedShapes = sorted([(n, p) for n, projection in shapes.items() for p in projection],
                             key=lambda x: metric(x[1]))[:maxShapes]
        shapes = {}
        for n, p in sortedShapes:
            shapes[n] = shapes.get(n, []) + [p]

    shapes = {
        n: [
            dict(vars(p))
            for p in proj
        ]
        for n, proj in shapes.items()
        if proj
    }

    # further filtering ?
    return shapes

def metric_max_min(proj, isGrid, dimensions):
    return max(proj)/min(proj)

def metric_average_distance(proj, isGrid, dimensions):
    # TODO
    return 1

def metric_distance_from_center(proj, isGrid, dimensions):
    # TODO
    return 1

def flattened_optimal_shapes(shape, dimensions):
    # yield the possible (optimal) flattened shapes starting from the given perfect hypercube
    for n in range(2, len(dimensions) + 1):
        for l in itertools.combinations(zip(enumerate(shape), dimensions.values()), n):
            maximalSizes = sorted((x for x in l), key=lambda x: x[0][1], reverse=True)[:n-1]
            flattenFactor = reduce(operator.mul, [d / s for (i, s), d in maximalSizes], 1)
            maxIndices = [x[0][0] for x in maximalSizes]
            restSize = reduce(operator.mul, [s for i, s in enumerate(shape) if not (i in maxIndices)], 1)
            restSize = restSize / flattenFactor
            if (restSize < 1):
                continue
            newShape = [x[1] for x in maximalSizes]
            numRestSizes = len(shape) - len(maximalSizes)
            newShape.extend(restSize ** (1.0 / (numRestSizes)) for i in range(numRestSizes))
            yield newShape

# average distance for a chain or ring of size p inside a torus of size d
# according to Pahrami, exact formulas for average internode distance..., 2013
def average_distance(p, d, isGrid):
    if (isGrid or p <= d // 2):
        return (p - 1/p) / 3
    elif (p >= d):
        return (p - (p % 2)/p) / 4
    else:
        # linearly interpolate between the two values
        # this is not rigorous but only used for optimal AD calculations
        # which are debatable anyway
        v1 = (p - 1/p) / 3
        v2 = (p - (p % 2)/p) / 4
        coef = (d - p) / (d // 2 + d % 2)
        return v1 * coef + v2 * (1 - coef)

def metric_ad(proj, isGrid, dimensions):
    return sum(average_distance(p, d, isGrid) for p, d in zip(proj, dimensions.values()))

def opt_ad(size, isGrid, dimensions):
    optSide = size ** (1.0/len(dimensions))
    if isGrid:
        return len(dimensions) * average_distance(optSide, True)
    else:
        shape = [optSide for x in dimensions]
        optAD = metric_ad(shape, False, dimensions)
        for newShape in flattened_optimal_shapes(shape, dimensions):
            newAD = metric_ad(newShape, False, dimensions)
            if (newAD < optAD):
                optAD = newAD
        return optAD

def metric_diameter(proj, isGrid, dimensions):
    if isGrid:
        return sum(map(lambda x: x - 1, proj))
    else:
        return sum(min(s-1, d // 2) for s, d in zip(proj, dimensions.values()))

def opt_diameter(size, isGrid, dimensions):
    optSide = size ** (1.0/len(dimensions))
    if isGrid:
        return len(dimensions) * (optSide - 1)
    else:
        shape = [optSide for x in dimensions]
        optDiam = metric_diameter(shape, False, dimensions)
        for newShape in flattened_optimal_shapes(shape, dimensions):
            newDiam = metric_diameter(newShape, False, dimensions)
            if (newDiam < optDiam):
                optDiam = newDiam
        return optDiam

# this metric is not diameter / perfect diameter but rather shape / perfect cubic shape
def metric_compactness(proj, isGrid, dimensions):
    size = reduce(operator.mul, proj, 1.0)
    # inverse of compactness definition of paper to have smallest number = best
    return sum(proj) / (len(proj) * (size ** (1.0/ len(proj))))

# Nodes affected (NA) and links affected (LA) are not important since all projections
# are convex, thus, it simply represents the size (# of resources) of the projection
def metric_nodes_affected(proj, isGrid, dimensions):
    return reduce(operator.mul, proj, 1)
def metric_links_affected(proj, isGrid, dimensions):
    if metric_nodes_affected(proj, isGrid) <= 1:
        return 0
    # count the links dimension by dimension. for each dimension, we have (size-1) links (unless in a torus and size is full dimension),
    # then multiply by the product of the sizes of all other dimensions
    sum = 0
    for dim in proj._fields:
        dimSize = getattr(proj, dim)
        dimLinks = dimSize - 1 if isGrid or dimSize != dimensions[dim] else dimSize
        sum += dimLinks * reduce(lambda acc, x: acc * getattr(proj, x) if x != dim else acc, proj._fields, 1)
    return sum
