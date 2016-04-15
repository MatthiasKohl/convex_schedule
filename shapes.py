#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import operator
from collections import namedtuple
from functools import reduce
from subprocess import Popen, PIPE
from math import isinf

TMP_SHAPES = 'shapes'
TMP_EXT = '.tmp'
TMP_TORUS = 'torus'
TMP_GRID = 'grid'
TMP_SEP = '_'


def process_candidates(s, dimensions, boundaries):
    Projection = namedtuple('Projection', dimensions)
    shape_candidates = {}
    for line in s.split('\n'):
        factors = line.split(',')
        if (len(factors) <= 2): continue
        n = int(factors[0])
        shape_candidates[n] = shape_candidates.get(n, []) + [Projection(**{bin: int(factor) for bin, factor in zip(dimensions, factors[1:])})]
    return shape_candidates

def shape_candidates(extProc, isGrid, dimensions, boundaries):
    Projection = namedtuple('Projection', dimensions)
    tmp_top = TMP_TORUS
    if (isGrid): tmp_top = TMP_GRID
    tmpFileName = TMP_SHAPES + TMP_SEP + tmp_top + TMP_SEP + 'x'.join(map(str, boundaries)) + TMP_EXT
    result = {}
    try:
        tmpFile = open(tmpFileName, 'r')
        result = process_candidates(tmpFile.read(), dimensions, boundaries)
        tmpFile.close()
    except FileNotFoundError:
        pass

    if (result):
        return result

    args = [extProc, '1' if isGrid else '0']
    args.extend(map(str, boundaries))
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
    return process_candidates(outstring, dimensions, boundaries)

def max_potential_shape(shapes_candidates, minSize=1, cutOffPercentage=float('inf'), metric=lambda x : 1):
    cutOffSize = int(minSize*(1+cutOffPercentage)) if not isinf(cutOffPercentage) else len(shapes_candidates)
    shapes_metric = map(lambda x: metric(x[1]),
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

def metric_max_min(proj):
    return max(proj)/min(proj)

def metric_average_distance(proj):
    # TODO
    return 1

def metric_distance_from_center(proj):
    # TODO
    return 1

# what to do with projections that are bigger than half the size of a dimension ?
# because of torus topology, we could have a non-convex shape here
def metric_diameter(proj):
    return sum(map(lambda x: x -1, proj))

# this metric is not diameter / perfect diameter but rather shape / perfect cubic shape
def metric_compactness(proj):
    size = reduce(operator.mul, proj, 1.0)
    # inverse of compactness definition of paper to have smallest number = best
    return sum(proj) / (len(proj) * (size ** (1.0/ len(proj))))

# Nodes affected (NA) and links affected (LA) are not important since all projections
# are convex, thus, it simply represents the size (# of resources) of the projection
def metric_nodes_affected(proj):
    return reduce(operator.mul, proj, 1)
def metric_links_affected(proj):
    sum = 0
    for dim in proj._fields:
        sum += (getattr(proj, dim) - 1) * \
            reduce(lambda acc, x: acc * getattr(proj, x) if x != dim else acc, proj._fields, 1)
    return sum
