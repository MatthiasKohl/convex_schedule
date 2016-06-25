#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import operator
import copy
from collections import namedtuple
from functools import reduce
from subprocess import Popen, PIPE
from math import isinf
import itertools
import os

ABS_DIR = os.path.dirname(os.path.abspath(__file__))
GEN_PROC = ABS_DIR + '/gen_shapes'

TMP_SHAPES = 'shapes'
TMP_EXT = '.tmp'
TMP_TORUS = 'torus'
TMP_GRID = 'grid'
TMP_SEP = '_'

def char_range(start_char, n):
    return [chr(c) for c in range(ord(start_char), ord(start_char)+n)]

def process_candidates(s, dimensions):
    Projection = namedtuple('Projection', dimensions.keys())
    shape_candidates = {}
    for line in s.split('\n'):
        factors = line.split(',')
        if (len(factors) <= 2): continue
        n = int(factors[0])
        shape_candidates[n] = (shape_candidates.get(n, [])
                               + [Projection(**{bin: int(factor) for bin, factor in zip(dimensions.keys(), factors[1:])})])
    return shape_candidates

def shape_candidates(is_grid, dimensions):
    tmp_top = TMP_GRID if is_grid else TMP_TORUS
    tmp_name = (ABS_DIR + '/' + TMP_SHAPES + TMP_SEP + tmp_top
                + TMP_SEP + 'x'.join(map(str, dimensions.values())) + TMP_EXT)
    try:
        tmp_file = open(tmp_name, 'r')
    except:
        args = [GEN_PROC, '1' if is_grid else '0']
        args.extend(map(str, dimensions.values()))
        proc = Popen(args, stdout=PIPE)
        out, err = proc.communicate()
        exitcode = proc.returncode

        if (exitcode != 0):
            print('ERROR in external program ' + GEN_PROC + ' -> exit code: ' + str(exitcode))
            return result

        out_string = out.decode('ascii')
        with open(tmp_name, 'w') as tmp_file:
            tmp_file.write(out_string)

        return process_candidates(out_string, dimensions)
    else:
        with tmp_file:
            return process_candidates(tmp_file.read(), dimensions)


def get_possible_sizes(dimensions, shape_candidates, alpha, is_full = True):
    # return a dict containing for every possible request size the list of request sizes
    # to look up in shape_candidates
    possible_sizes = {}
    total = reduce(operator.mul, dimensions.values(), 1)
    for n in reversed(range(1, total + 1)):
        cutoff_size = min(int(n+n*alpha), total)
        possible_sizes[n] = [r for r in shape_candidates if r >= n and r <= cutoff_size]
        if (is_full and not possible_sizes[n]):
            # make sure that we can have a shape for n by getting the minimal sized shapes
            possible_sizes[n] = possible_sizes[n+1]
    return possible_sizes

# returns a deep copy of the given candidates dictionary without the shapes that are
# rotations of other shapes already in the dictionary
def candidates_without_rotations(shape_candidates):
    candidates = copy.deepcopy(shape_candidates)
    for n in candidates:
        while (True):
            has_rotations = False
            for p1, p2 in itertools.combinations(candidates[n], 2):
                if ({s for s in p1} == {s for s in p2}):
                    has_rotations = True
                    break
            if (not has_rotations):
                break
            candidates[n].remove(p1)
    return candidates

def max_potential_shape(shapes_candidates, is_grid, dimensions, min_size=1,
                        alpha=float('inf'), metric=lambda x, y, z : 1):
    cutoff_size = int(min_size*(1+alpha)) if not isinf(alpha) else len(shapes_candidates)
    shapes_metric = map(lambda x: metric(x[1], is_grid, dimensions),
                        [(n,p) for n, projection in shapes_candidates.items() if n >= min_size and
                        n <= cutoff_size for p in projection])
    return max(shapes_metric, default=0)

def num_potential_shapes(shapes_candidates, min_size=1, alpha=float('inf')):
    cutoff_size = int(min_size*(1+alpha)) if not isinf(alpha) else len(shapes_candidates)
    return sum(1 for n, proj in shapes_candidates.items() if n >= min_size and n <= cutoff_size for p in proj)

def potential_shapes(shapes_candidates, min_size=1, alpha=float('inf'), maxShapes=-1, metric=lambda p : 1):
    cutoff_size = int(min_size*(1+alpha)) if not isinf(alpha) else len(shapes_candidates)
    # filter by size
    shapes = {
        n: proj
        for n, proj in shape_candidates.items()
        if n >= min_size and n <= cutoff_size
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

def metric_max_min(proj, is_grid, dimensions):
    return max(proj)/min(proj)

def metric_average_distance(proj, is_grid, dimensions):
    # TODO
    return 1

def metric_distance_from_center(proj, is_grid, dimensions):
    # TODO
    return 1

def flattened_optimal_shapes(shape, dimensions):
    # yield the possible (optimal) flattened shapes starting from the given perfect hypercube
    for n in range(2, len(dimensions) + 1):
        for l in itertools.combinations(zip(enumerate(shape), dimensions.values()), n):
            max_sizes = sorted((x for x in l), key=lambda x: x[0][1], reverse=True)[:n-1]
            flatten_factor = reduce(operator.mul, [d / s for (i, s), d in max_sizes], 1)
            max_indices = [x[0][0] for x in max_sizes]
            leftover = reduce(operator.mul, [s for i, s in enumerate(shape) if not (i in max_indices)], 1)
            leftover = leftover / flatten_factor
            if (leftover < 1):
                continue
            new_shape = [x[1] for x in max_sizes]
            num_leftovers = len(shape) - len(max_sizes)
            new_shape.extend(leftover ** (1.0 / (num_leftovers)) for i in range(num_leftovers))
            yield new_shape

# average distance for a chain or ring of size p inside a torus of size d
# according to Pahrami, exact formulas for average internode distance..., 2013
def average_distance(p, d, is_grid):
    if (is_grid or p <= d // 2):
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

def metric_ad(proj, is_grid, dimensions):
    return sum(average_distance(p, d, is_grid) for p, d in zip(proj, dimensions.values()))

def opt_ad(size, is_grid, dimensions):
    opt_side = size ** (1.0/len(dimensions))
    if is_grid:
        return len(dimensions) * average_distance(opt_side, True)
    else:
        shape = [opt_side for x in dimensions]
        opt_average_distance = metric_ad(shape, False, dimensions)
        for new_shape in flattened_optimal_shapes(shape, dimensions):
            new_average_distance = metric_ad(new_shape, False, dimensions)
            if (new_average_distance < opt_average_distance):
                opt_average_distance = new_average_distance
        return opt_average_distance

def metric_diameter(proj, is_grid, dimensions):
    if is_grid:
        return sum(map(lambda x: x - 1, proj))
    else:
        return sum(min(s-1, d // 2) for s, d in zip(proj, dimensions.values()))

def opt_diameter(size, is_grid, dimensions):
    opt_side = size ** (1.0/len(dimensions))
    if is_grid:
        return len(dimensions) * (opt_side - 1)
    else:
        shape = [opt_side for x in dimensions]
        opt_diam = metric_diameter(shape, False, dimensions)
        for new_shape in flattened_optimal_shapes(shape, dimensions):
            new_diam = metric_diameter(new_shape, False, dimensions)
            if (new_diam < opt_diam):
                opt_diam = new_diam
        return opt_diam

# this metric is not diameter / perfect diameter but rather shape / perfect cubic shape
def metric_compactness(proj, is_grid, dimensions):
    size = reduce(operator.mul, proj, 1.0)
    # inverse of compactness definition of paper to have smallest number = best
    return sum(proj) / (len(proj) * (size ** (1.0/ len(proj))))

# Nodes affected (NA) and links affected (LA) are not important since all projections
# are convex, thus, it simply represents the size (# of resources) of the projection
def metric_nodes_affected(proj, is_grid, dimensions):
    return reduce(operator.mul, proj, 1)
def metric_links_affected(proj, is_grid, dimensions):
    if metric_nodes_affected(proj, is_grid) <= 1:
        return 0
    # count the links dimension by dimension. for each dimension, we have (size-1) links
    # (unless in a torus and size is full dimension),
    # then multiply by the product of the sizes of all other dimensions
    sum = 0
    for dim in proj._fields:
        dim_size = getattr(proj, dim)
        dim_links = dim_size - 1 if is_grid or dim_size != dimensions[dim] else dim_size
        sum += dim_links * reduce(lambda acc, x: acc * getattr(proj, x) if x != dim else acc, proj._fields, 1)
    return sum
