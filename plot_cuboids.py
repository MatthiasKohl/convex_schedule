#! /usr/bin/python3
# -*- encoding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import copy
import random
import packing_strategies
from shapes.convex_shapes import char_range, shape_candidates
from packing_space import ConvexSpace

def plot_cube_3d(fig, space, c='b'):
    ax = fig.gca(projection='3d')
    xmin, xmax = space.coordinates[0], space.coordinates[0] + space.boundaries[0]
    ymin, ymax = space.coordinates[1], space.coordinates[1] + space.boundaries[1]
    zmin, zmax = space.coordinates[2], space.coordinates[2] + space.boundaries[2]
    X, Y = np.meshgrid([xmin, xmax], [ymin, ymax])
    surf1 = ax.plot_surface(X, Y, zmin, color=c, linewidth=0, antialiased=False)
    surf2 = ax.plot_surface(X, Y, zmax, color=c, linewidth=0, antialiased=False)
    X, Z = np.meshgrid([xmin, xmax], [zmin, zmax])
    surf3 = ax.plot_surface(X, ymin, Z, color=c, linewidth=0, antialiased=False)
    surf4 = ax.plot_surface(X, ymax, Z, color=c, linewidth=0, antialiased=False)
    Y, Z = np.meshgrid([ymin, ymax], [zmin, zmax])
    surf5 = ax.plot_surface(xmin, Y, Z, color=c, linewidth=0, antialiased=False)
    surf6 = ax.plot_surface(xmax, Y, Z, color=c, linewidth=0, antialiased=False)
    ax.set_xlim(0, space.coord_boundaries[0])
    ax.set_ylim(0, space.coord_boundaries[1])
    ax.set_zlim(0, space.coord_boundaries[2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

def plot_space_3d(fig, space, c='b'):
    def cut_spaces(spaces, d):
        if (d < 0):
            return spaces
        for s in spaces:
            if (s.coordinates[d] + s.boundaries[d] > s.coord_boundaries[d]):
                s1 = copy.deepcopy(s)
                s2 = copy.deepcopy(s)
                s1.boundaries[d] = s.coord_boundaries[d] - s.coordinates[d]
                s2.boundaries[d] = s.coordinates[d] + s.boundaries[d] - s.coord_boundaries[d]
                s2.coordinates[d] = 0
                return cut_spaces([s1, s2], d-1)
        return cut_spaces(spaces, d-1)

    for s in cut_spaces([space], len(space.coord_boundaries) - 1):
        plot_cube_3d(fig, s, c)

def plot_bin(bin, fig=None):
    is_full_plot = fig == None
    if (is_full_plot):
        fig = plt.figure()
    random.seed()
    for s in bin.spaces:
        # random RGB
        plot_space_3d(fig, s, (random.random(), random.random(), random.random()))
    if (is_full_plot):
        plt.show()

def plot_sample_packing():
    boundaries = [24,24,24]
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates(False, dimensions)
    total = 24*24*24

    random.seed()
    request_sizes = [int(x * (total-1) / 8) + 1 for x in
    np.random.normal(0.5, 0.5, random.randint(100, 200)) if x >= 0 and x <= 8]
    request_sizes = dict(enumerate(request_sizes))
    bins = packing_strategies.ff_best_metric_first(dimensions, request_sizes, candidates, 0.15)

    for b in bins:
        plot_bin(b)

if __name__ == '__main__':
    plot_sample_packing()
