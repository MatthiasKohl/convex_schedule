#! /usr/bin/python3
# -*- encoding: utf-8 -*-

from shapes import shape_candidates, max_potential_shape, num_potential_shapes,\
metric_max_min, metric_average_distance, metric_distance_from_center,\
metric_diameter, metric_compactness, metric_nodes_affected, metric_links_affected

import operator
from collections import namedtuple
from functools import reduce
from math import isinf
import time
import matplotlib.pyplot as plt
import numpy as np

def print_max_shapes(shape_candidates, resources, boundaries, dimensions, start_time, isGrid, args):
    max_num_shapes = sorted([(n, projections) for n, projections in shape_candidates.items()],
                 key=lambda x: len(x[1]))[-3:]
    for n, projections in max_num_shapes:
        print(n)

def plot_num_shapes(shape_candidates, resources, boundaries, dimensions, start_time, isGrid, args):
    numShapes = []
    for i in resources:
            numShapes.append(num_potential_shapes(shape_candidates, i, 0))
    numCumulShapes = []
    for idx in range(len(resources)):
        numCumulShapes.append(sum(numShapes[:idx+1]))
    if (args[0]):
        plt.plot(resources, numCumulShapes, label='cumulative number of shapes')

    for alpha in args[1:]:
        numShapesAlpha = []
        for r in resources:
            cutOffSize = int(r+r*alpha)
            numShapesAlpha.append(sum(numShapes[r-1:cutOffSize-1+1]))
        print('[' + str(time.time()-start_time) + '] plotting numShapes for alpha ' + str(alpha))
        plt.plot(resources, numShapesAlpha, label='#shapes for alpha=' + str(alpha))
    plt.legend()
    plt.xlabel('Number of resources')
    plt.ylabel('Number of shapes')
    top = ' grid' if isGrid else ' torus'
    plt.title('Number of convex shapes in a ' + 'x'.join(map(str, boundaries)) + top +' for a given number of resources')
    plt.show()

def plot_max_metric(shape_candidates, resources, boundaries, dimensions, start_time, isGrid, args):
    maxMetric = []
    resources = range(1, reduce(operator.mul, boundaries, 1) + 1)
    for i in resources:
            maxMetric.append(max_potential_shape(shape_candidates, i, 0, args[0]))

    for alpha in args[2:]:
        maxMetricAlpha = []
        for idx, i in enumerate(resources):
            maxMetricAlpha.append(max(maxMetric[idx:int(idx*(1+alpha))+1]))
        print('[' + str(time.time()-start_time) + '] plotting maxMetric for alpha ' + str(alpha))
        plt.plot(resources, maxMetricAlpha, label='max metric for alpha=' + str(alpha))
    plt.legend(loc='lower right')
    plt.xlabel('Number of resources')
    plt.ylabel('Max metric of shapes')
    top = ' grid' if isGrid else ' torus'
    plt.title('Maximum possible ' + args[1] + ' metric of convex shapes in a ' + 'x'.join(map(str, boundaries)) + top + ' for a given number of resources')
    plt.show()

def plot_metric_scatter(shape_candidates, resources, boundaries, dimensions, start_time, isGrid, args):
    shapes_metrics = {
        n: {args[0](p) for p in projections}
        for n, projections in shape_candidates.items()
    }

    # TODO try to add density to this plot (color for how many shapes there are for a given point on scatter)
    for alpha in args[2:]:
        shapes_metric_scatter_alpha = {}
        for r in resources:
            cutOffSize = int(r + r*alpha)
            shapes_metric_scatter_alpha[r] = {m for n, metrics in shapes_metrics.items() if n >= r and n <= cutOffSize for m in metrics}
        scatter_metric_base = [n for n, metrics in shapes_metric_scatter_alpha.items() for m in metrics]
        scatter_metric_values = [m for n, metrics in shapes_metric_scatter_alpha.items() for m in metrics]
        plt.scatter(scatter_metric_base, scatter_metric_values, s=1, label='metric with alpha=' + str(alpha))
        #dim = len(dimensions)
        #plt.plot(resources, [dim*(r ** (1.0/dim) - 1) for r in resources])
        plt.legend(loc='upper left')
        plt.xlabel('Number of resources')
        plt.ylabel('Metric of shape')
        top = ' grid' if isGrid else ' torus'
        plt.title('Scatter plot of metric ' + args[1] + ' for convex shapes in a ' +  'x'.join(map(str, boundaries)) + top + ' for a given number of resources')
        plt.show()

def print_covering_shapes(shape_candidates, resources, boundaries, dimensions, start_time, isGrid, args):
    #TODO find set of shapes that allow to have at least 1 shape for each possible number of resources that one can ask for, given an alpha
    # + possibly find a set of shapes given an alpha and max metric that allow to allocate any possible number of resources
    shapes_metrics = {
        n: {args[0](p) for p in projections}
        for n, projections in shape_candidates.items()
    }

    total = reduce(operator.mul, boundaries, 1)
    for alpha in args[2:]:
        shapes_satisfying_alpha = {total: shape_candidates[total]}
        isCovered = False
        for minSize in reversed(resources[:-1]):
            cutOffSize = int(minSize + minSize*alpha)
            # check if list already contains a shape that can satisfy this request
            isCovered = False
            for r in range(minSize + 1, cutOffSize + 1):
                if (r in shapes_satisfying_alpha):
                    isCovered = True
                    break
            if (isCovered):
                continue
            # take the first possible shapes that cover this resource
            isCovered = False

            # for r in range(minSize, cutOffSize + 1):
            #     if (r in shape_candidates and shape_candidates[r]):
            #         shapes_satisfying_alpha[r] = shape_candidates[r]
            #         isCovered = True
            #         break

            # TODO this approach is not perfect, because it doesn't consider the whole distribution of
            # metrics vs coverage
            shapes_best_metrics = sorted([(n, m, p) for ((n, metrics), (n, projections)) in \
                                         zip(shapes_metrics.items(), shape_candidates.items()) \
                                         if n >= minSize and n <= cutOffSize for m, p in zip(metrics, projections)],
                                         key = lambda x: x[1])[:1]
            if (shapes_best_metrics):
                isCovered = True
                for n, m, p in shapes_best_metrics:
                    shapes_satisfying_alpha[n] = shapes_satisfying_alpha.get(n, []) + [p]

            if (not isCovered):
                print('Cannot satisfy request for size ' + str(minSize) + ' with shapes using alpha ' + str(alpha) + ' and metric: ' + args[1])
                break
        if (isCovered):
            print ('Number of covering shapes for alpha ' + str(alpha) + ' : ' + str(len(shapes_satisfying_alpha)) + 'using metric '+ args[1] + '\tlist of covering resources: ' + str(sorted([n for n, projections in shapes_satisfying_alpha.items()])))
            print ('Average metric (' + args[1] + '): ' + str(np.mean([args[0](p) for n, projections in shapes_satisfying_alpha.items() for p in projections])))

# MAIN function
import sys, getopt

def usage():
    print('Usage: ' + sys.argv[0] + ' [options] [list of boundaries of torus/grid]')
    print('Options:\n\
    -b --boundaries\t\trequired: boundaries of the torus/grid as space-separated list of integers\n\
    -g --grid\t\t\tif a grid should be considered and not a torus\n\
    --maxShapes\t\t\tprint the amount of resources which yield the maximum number of shapes\n\
    --numShapes <args>\t\tplot the number of shapes per resources. this requires the following arguments: [alphas] <cumulative>. alphas can be either of: list <list of alpha> specifying a list of alphas to plot or <alpha start> <alpha end> <alpha step>. <cumulative> is a boolean to plot the cumulative number of shapes\n\
    --maxMetric <args>\t\tplot the maximum metric for the given shapes per resources. this requires the following arguments: <metric> [alphas]. <metric> can be "diameter", "maxmin", "compactness". [alphas] is as in numShapes\n\
    --metricScatter <args>\tscatter plot of the given metric for the given shapes per resources. this requires the following arguments: <metric> [alphas], as in maxMetric\n\
    --coveringShapes <args>\tprint a covering set of shapes/resources given alpha and metric for the given shapes. this requires the following arguments: <metric> [alphas], as in maxMetric\n\
    -d\t\t\t\tdebug option\n\
    -h --help\t\t\tprint this and exit')

def char_range(c1, numChars):
    return [chr(c) for c in range(ord(c1), ord(c1)+numChars)]

def handle_plot(isGrid, boundaries, func, args):
    resources = range(1, reduce(operator.mul, boundaries, 1) + 1)
    dimensions = char_range('x', len(boundaries))
    start_time = time.time()
    print('[' + str(time.time()-start_time) + '] Starting shapes computation')
    candidates = shape_candidates("./primes", isGrid, dimensions, boundaries)
    print('[' + str(time.time()-start_time) + '] Computed all shapes')
    func(candidates, resources, boundaries, dimensions, start_time, isGrid, args)

def parseAlphaList(args):
    try:
        if (args[0] == 'list'):
            return [float(f) for f in args[1:]]
        else:
            return  [i for i in np.arange(float(args[0]), float(args[1]), float(args[2]))]
    except:
        print('Error specifying alphas: ' + str(args))
        return []

def parseMetric(arg):
    if (arg == 'diameter'):
        return [metric_diameter, 'diameter']
    elif (arg == 'maxmin'):
        return [metric_max_min, 'max-min']
    elif (arg == 'compactness'):
        return [metric_compactness, 'compactness']
    print('Error specifying metric: ' + str(arg))
    return []

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hdb:g", ["help", "boundaries=", "grid", "maxShapes=", "numShapes=", "maxMetric=", "metricScatter=", "coveringShapes="])
    except getopt.GetoptError:
        print('getopt error')
        usage()
        sys.exit(2)
    boundaries, fundefs, isGrid, isError = [], [], False, False
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt == '-d':
            global _debug
            _debug = 1
        elif opt in ('-g', '--grid'):
            isGrid = True
        elif opt in ("-b", "--boundaries"):
            try:
                boundaries = [int(b) for b in arg.split()]
            except:
                pass
        elif opt == '--maxShapes':
            fundefs.append((max_shapes, []))
        elif opt == '--numShapes':
            args = arg.split()
            isCumul = False
            try:
                float(args[-1])
            except ValueError:
                isCumul = True
            if (isCumul):
                alphas = parseAlphaList(args[:-1])
            else:
                alphas = parseAlphaList(args)
            if (not alphas):
                isError = True
                break
            fundefs.append((plot_num_shapes, [isCumul] + alphas))
        elif opt == '--maxMetric':
            args = arg.split()
            metric = parseMetric(args[0])
            if (not metric):
                isError = True
                break
            alphas = parseAlphaList(args[1:])
            if (not alphas):
                isError = True
                break
            fundefs.append((plot_max_metric, metric + alphas))
        elif opt == '--metricScatter':
            args = arg.split()
            metric = parseMetric(args[0])
            if (not metric): break
            alphas = parseAlphaList(args[1:])
            if (not alphas): break
            fundefs.append((plot_metric_scatter, metric + alphas))
        elif opt == '--coveringShapes':
            args = arg.split()
            metric = parseMetric(args[0])
            if (not metric):
                isError = True
                break
            alphas = parseAlphaList(args[1:])
            if (not alphas):
                isError = True
                break
            fundefs.append((print_covering_shapes, metric + alphas))
    if (isError or not boundaries):
        usage()
        sys.exit()
    for f, args in fundefs:
        handle_plot(isGrid, boundaries, f, args)

if __name__ == "__main__":
    main(sys.argv[1:])

