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
    try:
        tmpFile = open(tmpFileName, 'r')
        result = process_candidates(tmpFile.read(), dimensions)
        tmpFile.close()
    except FileNotFoundError:
        pass

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
    tmpFile = open(tmpFileName, 'w')
    tmpFile.write(outstring)
    tmpFile.close()
    return process_candidates(outstring, dimensions)

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

# what to do with projections that are bigger than half the size of a dimension ?
# because of torus topology, we could have a non-convex shape here
def metric_diameter(proj, isGrid, dimensions):
    if isGrid:
        return sum(map(lambda x: x - 1, proj))
    else:
        return sum(min(s-1, d // 2) for s, d in zip(proj, dimensions.values()))

def get_opt_diameter(shape, dimensions):
    # check if flattening gains us diameter
    currentDiam = metric_diameter(shape, False, dimensions)
    for n in range(2, len(dimensions) + 1):
        for l in itertools.combinations(zip(enumerate(shape), dimensions.values()), n):
            maximalSizes = sorted((x for x in l), key=lambda x: x[0][1], reverse=True)[:n-1]
            factor = reduce(operator.mul, [d / s for (i, s), d in maximalSizes], 1)
            maxIndices = [x[0][0] for x in maximalSizes]
            restSize = reduce(operator.mul, [s for i, s in enumerate(shape) if not (i in maxIndices)], 1)
            restSize = restSize / factor
            newShape = [x[1] for x in maximalSizes]
            numRestSizes = len(shape) - len(maximalSizes)
            newShape.extend(restSize ** (1.0 / (numRestSizes)) for i in range(numRestSizes))
            newDiam = metric_diameter(newShape, False, dimensions)
            if (newDiam < currentDiam):
                return get_opt_diameter(newShape, dimensions)
    return currentDiam

def opt_diameter(size, isGrid, dimensions):
    optSide = size ** (1.0/len(dimensions))
    if isGrid:
        return len(dimensions) * (optSide - 1)
    else:
        return get_opt_diameter([optSide for i in dimensions], dimensions)

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
