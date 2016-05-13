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
from scipy.stats import gaussian_kde

def print_max_shapes(shape_candidates, resources, dimensions, start_time, isGrid, args):
    max_num_shapes = sorted([(n, projections) for n, projections in shape_candidates.items()],
                 key=lambda x: len(x[1]))[-3:]
    for n, projections in max_num_shapes:
        print('#Shapes for n=' + str(n) + ' : ' + str(len(projections)))

def plot_num_shapes(shape_candidates, resources, dimensions, start_time, isGrid, args):
    # args: [isShowCumulativeNumber, alphas]
    numShapes = []
    for i in resources:
            numShapes.append(num_potential_shapes(shape_candidates, i, 0))
    numCumulShapes = []
    for idx in range(len(resources)):
        numCumulShapes.append(sum(numShapes[:idx+1]))
    if (args[0]):
        plt.plot(resources, numCumulShapes, label='cumulative number of shapes')

    total = reduce(operator.mul, dimensions.values(), 1)
    for alpha in args[1:]:
        numShapesAlpha = []
        for r in resources:
            cutOffSize = min(int(r+r*alpha), total)
            numShapesAlpha.append(sum(numShapes[r-1:cutOffSize-1+1]))
        max_num_shapes = sorted(enumerate(numShapesAlpha), key=lambda x: x[1])[-3:]
        print('[' + str(time.time()-start_time) + '] plotting numShapes for alpha ' + str(alpha) +
              ' -> max numShapes: ' + "||".join(map(lambda x: str(x[1]) + ' (resource: ' + str(x[0]) + ')', max_num_shapes)))
        plt.plot(resources, numShapesAlpha, label='#shapes for alpha=' + str(alpha))
    plt.legend()
    plt.xlabel('Number of resources')
    plt.ylabel('Number of shapes')
    top = ' grid' if isGrid else ' torus'
    plt.title('Number of convex shapes in a ' + 'x'.join(map(str, dimensions.values())) + top +' for a given number of resources')
    plt.show()

def plot_max_metric(shape_candidates, resources, dimensions, start_time, isGrid, args):
    # args: [metric function, metric name, alphas]
    maxMetric = []
    resources = range(1, reduce(operator.mul, dimensions.values(), 1) + 1)
    for i in resources:
            maxMetric.append(max_potential_shape(shape_candidates, isGrid, dimensions, i, 0, args[0]))

    for alpha in args[2:]:
        maxMetricAlpha = []
        for idx, i in enumerate(resources):
            maxMetricAlpha.append(max(maxMetric[idx:int(idx*(1+alpha))+1]))
        print('[' + str(time.time()-start_time) + '] plotting maxMetric for alpha ' + str(alpha))
        plt.plot(resources, maxMetricAlpha, label='max ' + args[1] + ' for alpha=' + str(alpha))
    plt.legend(loc='lower right')
    plt.xlabel('Number of resources')
    plt.ylabel('Max metric of shapes')
    top = ' grid' if isGrid else ' torus'
    plt.title('Maximum possible ' + args[1] + ' metric of convex shapes in a ' + 'x'.join(map(str, dimensions.values())) + top + ' for a given number of resources')
    plt.show()

def plot_metric_scatter(shape_candidates, resources, dimensions, start_time, isGrid, args):
    # args: [metric function, metric name, alphas]
    shapes_metrics = {
        n: {args[0](p, isGrid, dimensions) for p in projections}
        for n, projections in shape_candidates.items()
    }

    total = reduce(operator.mul, dimensions.values(), 1)
    for alpha in args[2:]:
        shapes_metric_scatter_alpha = {}
        for r in resources:
            cutOffSize = min(int(r + r*alpha), total)
            shapes_metric_scatter_alpha[r] = {m for n, metrics in shapes_metrics.items() if n >= r and n <= cutOffSize for m in metrics}
        scatter_metric_base = [n for n, metrics in shapes_metric_scatter_alpha.items() for m in metrics]
        scatter_metric_values = [m for n, metrics in shapes_metric_scatter_alpha.items() for m in metrics]
        plt.scatter(scatter_metric_base, scatter_metric_values, s=1, label=args[1] + ' with alpha=' + str(alpha))
        #dim = len(dimensions)
        #plt.plot(resources, [dim*(r ** (1.0/dim) - 1) for r in resources])
        plt.legend(loc='upper left')
        plt.xlabel('Number of resources')
        plt.ylabel('Metric of shape')
        top = ' grid' if isGrid else ' torus'
        plt.title('Scatter plot of metric ' + args[1] + ' for convex shapes in a ' +  'x'.join(map(str, dimensions.values())) + top + ' for a given number of resources')
        plt.show()

def plot_metric_scatter_density(shape_candidates, resources, dimensions, start_time, isGrid, args):
    # args: [metric function, metric name, nBins, alphas]
    shapes_metrics = {
        n: [args[0](p, isGrid, dimensions) for p in projections]
        for n, projections in shape_candidates.items()
    }
    total = reduce(operator.mul, dimensions.values(), 1)
    for alpha in args[3:]:
        shapes_metric_scatter_alpha = []
        for r in resources:
            cutOffSize = min(int(r + r*alpha), total)
            shapes_metric_scatter_alpha.extend([[r,m] for n, metrics in shapes_metrics.items() if n >= r and n <= cutOffSize for m in metrics])

        # densityPerResource = {}
        # for n,m in shapes_metric_scatter_alpha:
        #     densityPerResource[n] = densityPerResource.get(n, 0) + 1
        # density = [densityPerResource[n] for n,m in shapes_metric_scatter_alpha]
        #plt.scatter(metric_values[:, 0], metric_values[:, 1], c=density, s=10)

        metric_values = np.array(shapes_metric_scatter_alpha)
        bins = [args[2], len(resources)]
        hist, locy, locx = np.histogram2d(metric_values[:, 1], metric_values[:, 0], bins=bins)
        X, Y = np.meshgrid(locx, locy)
        hist = np.ma.masked_array(hist, hist < 1)
        plt.pcolormesh(X, Y, hist, vmin=1)

        print("Hist extent (first x, first y, last x, last y): (" + str(locx[0]) + ", " + str(locy[0]) + ", " + str(locx[-1]) + ", " + str(locy[-1]) + ") first/last values: " + str(hist[0][0]) + "/" + str(hist[-1][-1]))

        #plt.imshow(hist, aspect='auto', origin='lower', interpolation='nearest', extent=(locx[0], locx[-1], locy[0], locy[-1]))
        plt.colorbar()
        plt.xlabel('Number of resources')
        plt.ylabel('Metric of shape')
        top = ' grid' if isGrid else ' torus'
        plt.title('Scatter plot of metric ' + args[1] + ' for convex shapes in a ' +  'x'.join(map(str, dimensions.values())) + top + ' for a given number of resources (with alpha=' + str(alpha) + ')')
        plt.show()

# return the number of resources covered if resource r is added to shapes_satisfying_alpha, for a given alpha
def numCovered(r, alpha, shapes_satisfying_alpha, maxSize):
    coveredMinR = max(1, int(r/(1+alpha)) + 1)
    possibleCovered = r - coveredMinR + 1
    for n, projs in shapes_satisfying_alpha.items():
        coveredMinN = max(1, int(n/(1 + alpha)) + 1)
        if (r < coveredMinN or coveredMinR > n):
            continue
        if (coveredMinR >= coveredMinN and r <= n):
            return 0
        if (coveredMinR < coveredMinN):
            possibleCovered = possibleCovered - (r - coveredMinN + 1)
            r = coveredMinN - 1
        if (r > n):
            possibleCovered = possibleCovered - (n - coveredMinR + 1)
            coveredMinR = n + 1
    if (possibleCovered <= 0):
        return 0
    return possibleCovered

def print_covering_shapes(shape_candidates, resources, dimensions, start_time, isGrid, args):
    # args: [metric function, metric name, strategy, alphas]
    #TODO possibly find a set of shapes given an alpha and max metric that allow to allocate any possible number of resources
    shapes_metrics = {
        n: {args[0](p, isGrid, dimensions) for p in projections}
        for n, projections in shape_candidates.items()
    }

    total = reduce(operator.mul, dimensions.values(), 1)
    for alpha in args[3:]:
        # TODO there might be a better mechanism to check whether a resource (integer) is covered
        unCovered = -1
        shapes_satisfying_alpha = {total: shape_candidates[total]}
        if (args[2] == 'none'):
            # no strategy: simply start from the highest resource, go down to lowest and add
            # any resource which is not already covered
            # this gives the minimum-cardinality set of resources to cover any possible resource
            for minSize in reversed(resources[:-1]):
                cutOffSize = min(int(minSize + minSize*alpha), total)
                # check if list already contains a shape that can satisfy this request
                isCovered = False
                for r in range(minSize + 1, cutOffSize + 1):
                    if (r in shapes_satisfying_alpha):
                        isCovered = True
                        break
                if (isCovered):
                    continue

                # take the first possible shapes that cover this resource
                unCovered = minSize
                for r in range(minSize, cutOffSize + 1):
                    if (r in shape_candidates and shape_candidates[r]):
                        shapes_satisfying_alpha[r] = shape_candidates[r]
                        unCovered = -1
                        break
                if (unCovered >= 0):
                    break
        else:
            # This part is ignored for now, it is not directly related to the exact problem we're examining

            # Go over best metric shapes, check if all is covered. if not, include next shape
            # this makes a very important assumption about the metric: it must not depend on
            # the size of the shape, if not, high sizes get penalized

            # group by best metric, then sort by best metric
            shapes_best_metrics = [(min(args[0](p, isGrid, dimensions) for p in projections), n, projections) for n, projections in shape_candidates.items()]
            metric_values = sorted(set(map(lambda x: x[0], shapes_best_metrics)))
            shapes_best_metrics_grouped = [[x for x in shapes_best_metrics if x[0] == m] for m in metric_values]

            for projections_list in shapes_best_metrics_grouped:
                #print("size of list: " + str(len(projections_list)) + " - best metric: " + str(projections_list[0][0]))
                # choose the resource from list which covers the most space left, remove any resources that do not cover anything
                while (projections_list):
                    maxCovered = max(projections_list, key=lambda x: numCovered(x[1], alpha, shapes_satisfying_alpha, total))
                    projections_list.remove(maxCovered)
                    (m, n, projections) = maxCovered
                    if (numCovered(n, alpha, shapes_satisfying_alpha, total) <= 0):
                        continue
                    #print("n : " + str(n) + ", m : " + str(m))

                    # add shapes until all resources are covered
                    unCovered = -1
                    for r in resources:
                        isCovered = False
                        cutOffSize = min(int(r + r*alpha), total)
                        for r2 in range(r, cutOffSize + 1):
                            if (r2 in shapes_satisfying_alpha):
                                isCovered = True
                                break
                        if (not isCovered):
                            unCovered = r
                            break
                    if (unCovered < 0): # all is covered
                        break
                    # everything is not covered yet, add a request size with its shapes
                    #print(str(n) + " - unCovered: " + str(unCovered) + " - numCovered: " + str(numCovered(n, alpha, shapes_satisfying_alpha, total)))
                    shapes_satisfying_alpha[n] = projections

                if (unCovered < 0): # all is covered
                    break

            # check if anything is not covered (something might have been added in last iteration)
            unCovered = -1
            for r in resources:
                isCovered = False
                cutOffSize = min(int(r + r*alpha), total)
                for r2 in range(r, cutOffSize + 1):
                    if (r2 in shapes_satisfying_alpha):
                        isCovered = True
                        break
                if (not isCovered):
                    unCovered = r
                    break


            # # TODO this doesn't work because we need to check if we can include some shape to make the covering work completely once we know that some resource is not covered
            # # check if all resources are covered
            # # try to cover any resources which are not yet covered
            # for n in reversed(resources):
            #     cutOffSize = min(int(n + n*alpha), total)
            #     unCovered = n
            #     for r in range(n, cutOffSize + 1):
            #         if (r in shapes_satisfying_alpha):
            #             unCovered = -1
            #             break
            #     if (unCovered < 0):
            #         continue

            #     # n is not covered, try to cover it using the best possible shape
            #     covering_shapes = [(r, projections) for r, projections in shape_candidates.items() if r >= n and r <= cutOffSize]
            #     if (not covering_shapes):
            #         break
            #     shapes_satisfying_alpha[n] = min(covering_shapes, key=lambda x: min(args[0](p, isGrid, dimensions) for p in x[1]))[1]


        if (unCovered < 0):
            print ('Number of covering shapes for alpha ' + str(alpha) + ' : ' + str(len(shapes_satisfying_alpha)) + ' using metric '+ args[1] + '\nlist of covering resources (with best metric): ' + str(sorted([(n, min(args[0](p, isGrid, dimensions) for p in projections)) for n, projections in shapes_satisfying_alpha.items()], key=lambda x: x[0])))
            bestMetrics = [sorted([args[0](p, isGrid, dimensions) for p in projections])[0] for n, projections in shapes_satisfying_alpha.items()]
            print ('Average best metric (' + args[1] + '): ' + str(np.mean(bestMetrics)))
        else:
            print('Cannot satisfy request for size ' + str(unCovered) + ' with shapes using alpha ' + str(alpha) + ' and metric: ' + args[1])

def print_metrics(shape_candidates, resources, dimensions, start_time, isGrid, args):
    # args: list of resources
    numRes = len(args)
    print("Number of given resources: " + str(numRes))
    for metric, metric_name in [(metric_diameter, "diameter"), (metric_max_min, "max-min"), (metric_compactness, "compactness")]:
        metricSum = sum(map(lambda r: min(metric(p, isGrid, dimensions) for p in shape_candidates[r]), args))
        print("Average best " + metric_name + " : " + str(metricSum/numRes))

# MAIN function
import sys, getopt

def usage():
    print('Usage: ' + sys.argv[0] + ' [options]')
    print('Options:\n\
    -b --boundaries=<args>\t\trequired: boundaries of the torus/grid as space-separated list of integers\n\
    -g --grid\t\t\tif a grid should be considered and not a torus\n\
    --maxShapes\t\t\tprint the amount of resources which yield the maximum number of shapes\n\
    --numShapes=<args>\t\tplot the number of shapes per resources. this requires the following arguments: [alphas] <cumulative>. alphas can be either of: list <list of alpha> specifying a list of alphas to plot or <alpha start> <alpha end> <alpha step>. <cumulative> is a boolean to plot the cumulative number of shapes\n\
    --maxMetric=<args>\t\tplot the maximum metric for the given shapes per resources. this requires the following arguments: <metric> [alphas]. <metric> can be "diameter", "maxmin", "compactness". [alphas] is as in numShapes\n\
    --metricScatter=<args>\tscatter plot of the given metric for the given shapes per resources. this requires the following arguments: <metric> [alphas], as in maxMetric\n\
    --metricScatterDens=<args>\tdensity scatter plot of the given metric for the given shapes per resources. this requires the following arguments: <metric> [alphas], as in maxMetric\n\
    --coveringShapes=<args>\tprint a covering set of shapes/resources given alpha and metric for the given shapes. this requires the following arguments: <metric> [alphas] <strategy> - <metric> and [alphas] as in maxMetric, <strategy>: none, bestMetric\n\
    --metrics=<args>\tprint metrics about a set of resources. <args> is either a range or a list as alphas in maxShapes\n\
    -d\t\t\t\tdebug option\n\
    -h --help\t\t\tprint this and exit')

def char_range(c1, numChars):
    return [chr(c) for c in range(ord(c1), ord(c1)+numChars)]

def handle_plot(isGrid, boundaries, func, args):
    resources = range(1, reduce(operator.mul, boundaries, 1) + 1)
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    start_time = time.time()
    print('[' + str(time.time()-start_time) + '] Starting shapes computation')
    candidates = shape_candidates("./primes", isGrid, dimensions)
    print('[' + str(time.time()-start_time) + '] Computed all shapes')
    func(candidates, resources, dimensions, start_time, isGrid, args)

def parseList(args, mapping):
    try:
        if (args[0] == 'list'):
            return [mapping(f) for f in args[1:]]
        else:
            return  [i for i in np.arange(float(args[0]), float(args[1]), float(args[2]))]
    except:
        print('Error specifying list: ' + str(args))
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
        opts, args = getopt.getopt(argv, "hdb:g", ["help", "boundaries=", "grid", "maxShapes", "numShapes=", "maxMetric=", "metricScatter=", "metricScatterDens=", "coveringShapes=", "metrics="])
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
            fundefs.append((print_max_shapes, []))
        elif opt == '--numShapes':
            args = arg.split()
            isCumul = False
            try:
                float(args[-1])
            except ValueError:
                isCumul = True
            if (isCumul):
                alphas = parseList(args[:-1], float)
            else:
                alphas = parseList(args, float)
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
            alphas = parseList(args[1:], float)
            if (not alphas):
                isError = True
                break
            fundefs.append((plot_max_metric, metric + alphas))
        elif opt == '--metricScatter':
            args = arg.split()
            metric = parseMetric(args[0])
            if (not metric):
                isError = True
                break
            alphas = parseList(args[1:], float)
            if (not alphas):
                isError = True
                break
            fundefs.append((plot_metric_scatter, metric + alphas))
        elif opt == '--metricScatterDens':
            args = arg.split()
            metric = parseMetric(args[0])
            if (not metric):
                isError = True
                break
            try:
                nBins = int(args[1])
            except:
                isError = True
                break
            alphas = parseList(args[2:], float)
            if (not alphas):
                isError = True
                break
            fundefs.append((plot_metric_scatter_density, metric + [nBins] + alphas))
        elif opt == '--coveringShapes':
            args = arg.split()
            metric = parseMetric(args[0])
            if (not metric):
                isError = True
                break
            strategy = 'none'
            try:
                float(args[-1])
                alphas = parseList(args[1:], float)
            except ValueError:
                alphas = parseList(args[1:-1], float)
                strategy = args[-1]
            if (not alphas or not strategy):
                isError = True
                break
            fundefs.append((print_covering_shapes, metric + [strategy] + alphas))
        elif opt == '--metrics':
            resources = parseList(arg.split(), int)
            if (not resources):
                isError = True
                break
            fundefs.append((print_metrics, resources))
    if (isError or not boundaries):
        usage()
        sys.exit()
    for f, args in fundefs:
        handle_plot(isGrid, boundaries, f, args)

if __name__ == "__main__":
    main(sys.argv[1:])

