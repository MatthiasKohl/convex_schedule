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
from shapes.convex_shapes import (shape_candidates, char_range, get_possible_sizes,
                                  candidates_without_rotations, metric_compactness, metric_diameter,
                                  opt_diameter, metric_ad, opt_ad)
from packing_bin import Bin

def unit_test_flat(b):
    return b.test_possible(False, True)
def unit_test_best(b):
    return b.test_possible(True, False)
def place_flat(b, shape, ID = None):
    first_non_flat_dim = next( (i for i, d in enumerate(shape) if d > 1), len(shape)-1)
    b.fit_flat(shape, first_non_flat_dim, ID)
def place_best(b, shape, ID = None):
    b.fit_best(shape, ID)

def generic_first_fit(dimensions, request_sizes, shape_candidates, alpha, size_generator,
                    shape_generator, place, initial_bins, unit_test, is_debug):
    possible_sizes = get_possible_sizes(dimensions, shape_candidates, alpha)
    bins = list(initial_bins)
    num_sizes = 0
    fitted_spaces = []
    last_shape = None
    max_fs = 0
    for ID, size in size_generator(request_sizes, possible_sizes, shape_candidates):
        is_fit = False
        for shape in shape_generator(size, last_shape, possible_sizes, shape_candidates):
            for b in bins:
                if (b.can_fit(shape)):
                    is_fit = True
                    place(b, shape, ID)
                    last_shape = shape
                    if (is_debug):
                        fitted_spaces.append(shape)
                        fitted_bin = b
                    break
            if (is_fit):
                break
        has_shapes = True
        if (not is_fit):
            new_bin = Bin(dimensions.values())
            # simply choose the first space in this case
            try:
                shape = next(shape_generator(size, last_shape, possible_sizes, shape_candidates))
                place(new_bin, shape, ID)
                bins.append(new_bin)
                last_shape = shape
                if (is_debug):
                    fitted_spaces.append(shape)
                    fitted_bin = new_bin
            except StopIteration:
                print('Cannot allocate size ' + str(size) + ' since no shapes are available')
                has_shapes = False
        if (is_debug and has_shapes and not unit_test(fitted_bin)):
            print('Impossible configuration:\n' + str(fitted_bin) + '\n')
            print('Spaces leading to this: ' + str(fitted_spaces))
            break
        num_sizes = num_sizes + 1
        if (is_debug and num_sizes % (int(len(request_sizes) / 20) + 1) == 0):
            print('.', end='', flush=True)
        num_fs = sum(1 for b in bins for s in b.freelist)
        if (num_fs > max_fs):
            max_fs = num_fs
    if (is_debug):
        print()
    # print('The maximum number of total free spaces during the whole packing was: ' +
    #       str(max_fs))
    return bins

def choose_candidate_flat(request_size, possible_sizes, shape_candidates):
    # get the minimum candidate (minimal size of the flattest candidates)
    # of all candidates that can hold this request_size according to possible sizes
    # this is the candidate that leaves the biggest convex free space when fit into the bin,
    # according to the heuristic
    if (not possible_sizes[request_size]):
        return request_size
    # return projection as list
    return list(min(p for n in possible_sizes[request_size] for p in shape_candidates[n]))

def ffd_flat(dimensions, request_sizes, shape_candidates, alpha,
            initial_bins=[], is_debug = False):
    def size_gen(request_sizes, possible_sizes, shape_candidates):
        # sort sizes based on their assigned shape (flattest possible) with shapes being
        # sorted based on size from first to last dimension
        request_spaces = map(lambda r:
                                  (r, choose_candidate_flat(r[1], possible_sizes, shape_candidates)),
                                  request_sizes.items())
        return (r for r, shape in sorted(request_spaces,
                                         key = lambda x: x[1] if isinstance(x[1], list) else
                                         [0 for d in dimensions], reverse = True))
    def shape_gen(size, last_shape, possible_sizes, shape_candidates):
        # always simply yield the flattest shape for each size
        yield choose_candidate_flat(size, possible_sizes, shape_candidates)

    return generic_first_fit(dimensions, request_sizes, shape_candidates, alpha,
                           size_gen, shape_gen, place_flat, initial_bins, unit_test_flat, is_debug)

# return difference between second lowest metric and lowest metric for all possible shapes
# for given request_size and alpha, or 0 if there are too few possible shapes
# TODO this could also return the ratio ?
def metrics_variance(request_size, dimensions, possible_sizes, candidates, metric):
    best_metrics = sorted(metric(p, False, dimensions)
                         for n in possible_sizes[request_size] for p in candidates[n])
    if (len(best_metrics) < 2):
         return 0
    return best_metrics[1] - best_metrics[0]

def ff_greatest_metric_delta(dimensions, request_sizes, shape_candidates, alpha,
                               initial_bins=[], is_debug=False):
    # pre-calculate candidates without rotation to consider for sorting the sizes
    candidates_wor = candidates_without_rotations(shape_candidates)
    def size_gen(request_sizes, possible_sizes, shape_candidates):
        # sort the requested sizes by non-increasing size,
        # then by non-increasing delta of best metrics of possible shapes (sort must be stable)
        return (s for s in sorted(sorted(request_sizes.items(),
                                         reverse=True, key=lambda x: x[1]),
                                  key=lambda r: metrics_variance(r[1], dimensions, possible_sizes,
                                                                candidates_wor, metric_ad),
                                  reverse=True))
    def shape_gen(size, last_shape, possible_sizes, shape_candidates):
        # sort by metric (best to worst) and by size for same metric
        spaces = sorted((p for n in possible_sizes[size] for p in shape_candidates[n]),
                        key=lambda p: reduce(operator.mul, p, 1))
        return (s for s in sorted(spaces, key=lambda p: metric_ad(p, False, dimensions)))

    return generic_first_fit(dimensions, request_sizes, shape_candidates, alpha, size_gen,
                           shape_gen, place_best, initial_bins, unit_test_best, is_debug)

def choose_best_metric_shapes(request_size, dimensions, possible_sizes, shape_candidates, metric):
    if (not possible_sizes[request_size]):
        return []
    best_metric = min(metric(p, False, dimensions) for n in possible_sizes[request_size]
                    for p in shape_candidates[n])
    # get min metric shapes, sorted by non-decreasing size
    considered = (p for n in possible_sizes[request_size] for p in shape_candidates[n])
    return (s for s in sorted(filter(lambda p:
                                     metric(p, False, dimensions) == best_metric, considered),
    key=lambda p: reduce(operator.mul, p, 1)))

def ffd_best_metric_always(dimensions, request_sizes, shape_candidates, alpha,
                        initial_bins=[], is_debug=False):
    possible_sizes = get_possible_sizes(dimensions, shape_candidates, alpha)
    # always pick the best metric shapes and then pack using that
    # if we choose diameter/AD as metric, the best metric of r1 cannot have greater size
    # than the best metric of r2 if r1 < r2, since there cannot be two shapes of different sizes
    # having the best diameter/AD for the requests
    def size_gen(request_sizes, possible_sizes, shape_candidates):
        return (s for s in sorted(request_sizes.items(), reverse=True, key=lambda x: x[1]))
    def shape_gen(size, last_shape, possible_sizes, shape_candidates):
        return choose_best_metric_shapes(size, dimensions, possible_sizes, shape_candidates,
                                      metric_ad)

    return generic_first_fit(dimensions, request_sizes, shape_candidates, alpha, size_gen,
                           shape_gen, place_best, initial_bins, unit_test_best, is_debug)

def ffd_best_metric_first(dimensions, request_sizes, shape_candidates, alpha,
                           initial_bins=[], is_debug=False):
    def size_gen(request_sizes, possible_sizes, shape_candidates):
        return (s for s in sorted(request_sizes.items(), reverse=True, key=lambda x: x[1]))
    def shape_gen(size, last_shape, possible_sizes, shape_candidates):
        # try any shape which is smaller than last, from best metric to worst
        # choose diameter as metric because bigger shapes should not get precedence
        # sort by size first to give precedence to smaller shapes if diameter is the same
        last_size = reduce(operator.mul,
                          (last_shape if last_shape != None else dimensions.values()), 1)
        considered = filter(lambda s: reduce(operator.mul, s, 1) <= last_size,
                            (p for n in possible_sizes[size] for p in shape_candidates[n]))
        sorted_shapes = sorted(considered, key=lambda p: reduce(operator.mul, p, 1))
        return (s for s in sorted(sorted_shapes,
                                  key=lambda p: metric_ad(p, False, dimensions)))

    return generic_first_fit(dimensions, request_sizes, shape_candidates, alpha, size_gen,
                           shape_gen, place_best, initial_bins, unit_test_best, is_debug)

# This is the chosen strategy, yields best results over all performed tests
def ff_best_metric_first(dimensions, request_sizes, shape_candidates, alpha,
                          initial_bins=[], is_debug=False):
    def size_gen(request_sizes, possible_sizes, shape_candidates):
        return (s for s in sorted(request_sizes.items(), reverse=True, key=lambda x: x[1]))
    def shape_gen(size, last_shape, possible_sizes, shape_candidates):
        # try any shape, from best metric to worst
        # choose diameter/AD as metric because bigger shapes should not get precedence
        # sort by size first to give precedence to smaller shapes if diameter/AD is the same
        sorted_shapes = sorted((p for n in possible_sizes[size] for p in shape_candidates[n]),
                              key=lambda p: reduce(operator.mul, p, 1))
        return (s for s in sorted(sorted_shapes,
                                  key=lambda p: metric_ad(p, False, dimensions)))

    return generic_first_fit(dimensions, request_sizes, shape_candidates, alpha, size_gen,
                           shape_gen, place_best, initial_bins, unit_test_best, is_debug)

def print_restuls(boundaries, results):
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

def get_stats(dimensions, request_sizes, bins):
    total = reduce(operator.mul, dimensions.values(), 1)
    total_request_size = sum(request_sizes.values())
    num_requests = len(request_sizes)
    num_bins = len(bins)
    num_opt_bins = int(math.ceil(total_request_size/total))

    total_used_space = sum(reduce(operator.mul, s.boundaries, 1) for b in bins for s in b.spaces)
    total_used_exceeding_percentage = (total_used_space - total_request_size) * 100 / total_request_size

    # unused space excludes the last bin's free space as it might still be filled
    unused_space = total_used_space - sum(reduce(operator.mul, s.boundaries, 1) for s in bins[-1].spaces)
    unused_space = (len(bins)-1) * total - unused_space
    unused_percentage = unused_space * 100 / ((len(bins) - 1) * total)

    def convex_used_space(b):
            # return space used by a convex mesh around all assigned spaces in bin b
            # get max coordinate in each dimension (assumes that spaces are packed from
            # min to max coordinate which is the case in all strategies)
            # this is not the tightest convex space in a torus, but an approximation for it
            return reduce(operator.mul,
                          [max(min(s.coordinates[d] + s.boundaries[d], s.coordBoundaries[d])
                               for s in b.spaces) for d in range(len(b.boundaries))], 1)
    total_sub_mesh_space = sum(convex_used_space(b) for b in bins)
    total_sub_mesh_percentage = (total_sub_mesh_space - total_request_size) * 100 / total_request_size

    avg_metric = sum(metric_ad(p.boundaries, False, dimensions) for b in bins
                          for p in b.spaces) / len(request_sizes)
    avg_opt_metric = sum(opt_ad(r, False, dimensions)
                             for r in request_sizes.values()) / len(request_sizes)

    return [total_request_size, num_requests, num_bins, num_opt_bins, unused_space, unused_percentage,
    total_used_space, total_used_exceeding_percentage, total_sub_mesh_space,
    total_sub_mesh_percentage, avg_metric, avg_opt_metric]

def perform_packing(boundaries, request_sizes, strategy, unit_test, alpha, print_detail = False):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates(False, dimensions)

    # get bin packing
    bins = strategy(dimensions, request_sizes, candidates, alpha)
    results = get_stats(dimensions, request_sizes, bins)

    if (print_detail):
        # print the configuration
        print('Packing as follows:')
        for i, b in enumerate(bins):
            print('Bin ' + str(i) + ':\n' + str(b))

        print_restuls(boundaries, results)

    # TODO: add unit test to check if request sizes/alpha are respected
    # (need to match each allocated space with request size)

    # check if bin configuration is possible for all bins
    for b in bins:
        if (not unit_test(b)):
            print(b)
            print('Not all bin configurations are possible ! Test FAILED')
            break

    return results

def uniform_random_requests(boundaries, strategy, unit_test, alpha, print_detail = False):
    total = reduce(operator.mul, boundaries, 1)
    random.seed()
    request_sizes = [random.randint(1, total) for i in range(random.randint(1000, 2000))]
    request_sizes = dict(enumerate(request_sizes))
    return perform_packing(boundaries, request_sizes, strategy, unit_test, alpha, print_detail)

def normal_random_requests(boundaries, strategy, unit_test, alpha, print_detail = False):
    total = reduce(operator.mul, boundaries, 1)
    random.seed()
    # random amount of requests with normal distributed
    # ignore values to < 0 or > 8, then scale up to the total size
    # thus, sizes > 1/8 of the torus have a fairly small probability
    request_sizes = [int(x * (total-1) / 8) + 1 for x in
    np.random.normal(0.5, 0.5, random.randint(1000, 2000)) if x >= 0 and x <= 8]
    request_sizes = dict(enumerate(request_sizes))
    return perform_packing(boundaries, request_sizes, strategy, unit_test, alpha, print_detail)

def curie_trace_requests(filename, boundaries, strategy, unit_test, alpha, columnIdx = 0, print_detail = False):
    request_sizes = {}
    with open(filename) as infile:
        i = 0
        for line in infile:
            request_sizes[i] = int(line.split()[columnIdx])
            i = i + 1
    return perform_packing(boundaries, request_sizes, strategy, unit_test, alpha, print_detail)

def test_strategies(boundaries):
    # test each strategy 10 times and take average measures
    strategies = [('flat', ffd_flat), ('best-always', ffd_best_metric_always),
    ('greatest-metric-delta', ff_greatest_metric_delta),
    ('best-metric-first-strict-decreasing', ffd_best_metric_first),
    ('best-metric-first', ff_best_metric_first)]
    for name, strategy in strategies:
        unit_test = unit_test_flat if strategy is ffd_flat else unit_test_best
        for alpha in [0.15, 1.0]:
            print('Testing strategy ' + name + ' with alpha ' + str(alpha))
            results = []
            for i in range(10):
                results.append(normal_random_requests(boundaries, strategy, unit_test, alpha))
            npResults = np.array(results)
            avgResults = np.mean(npResults, axis=0)
            print('Average results for strategy ' + name + ' and alpha ' + str(alpha))
            print_restuls(boundaries, avgResults)
            print('\n')

def test_strategies_curie(filename, boundaries):
    # test each strategy on the given sizes provided by trace
    strategies = [('flat', ffd_flat), ('best-always', ffd_best_metric_always),
    ('greatest-metric-delta', ff_greatest_metric_delta),
    ('best-metric-first-strict-decreasing', ffd_best_metric_first),
    ('best-metric-first', ff_best_metric_first)]
    for name, strategy in strategies:
        unit_test = unit_test_flat if strategy is ffd_flat else unit_test_best
        for alpha in [0.15, 1.0]:
            print('Testing strategy ' + name + ' with alpha ' + str(alpha))
            results = curie_trace_requests(filename, boundaries, strategy, unit_test, alpha)
            print_restuls(boundaries, results)
            print('\n')

if __name__ == '__main__':
    if (len(sys.argv) > 1):
        traceFile = sys.argv[1]
        test_strategies_curie(traceFile, [24, 24, 24])
        sys.exit(0)

    # JUST TESTS
    #normal_random_requests([24, 24, 24], ffd_flat, 0.15, True)
    test_strategies([24,24,24])
    #normal_random_requests([24,24,24], ffdGreatestMetricDeltaFirst, 0.15, True)
    #print_restuls([24,24,24],
    #             curie_trace_requests('request_sizes_scaled_5000.txt', [24,24,24],
    #                              ffdNonStrictEachBestMetricFirst, 0.15))
    #perform_packing([24,24,24],[56,34],ffdAllBestMetricFirst,1.0, True)
