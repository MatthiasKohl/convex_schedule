#! /usr/bin/python3
# -*- encoding: utf-8 -*-

# Generates a schedule from the provided simplified Blue Waters trace
# Outputs some stats about the schedule on stdout and a .jed file that can be used to generate
# a Jedule plot
# Optionally, a second argument can be provided as a prefix for the Jedule .jed file


import sys
import math
import itertools
import datetime
import copy
from functools import reduce
import operator
import pytz
import jedule.schedule_common as sc
from packing_strategies import ff_best_metric_first
from shapes.convex_shapes import (char_range, shape_candidates, metric_ad, opt_ad,
                                  get_possible_sizes)
from packing_bin import Bin
from packing_space import ConvexSpace
import plot_cuboids

def test_sched(schedule, jobs, boundaries):
    scheduled_IDs = sorted(i for i, t, s in schedule)
    if (len(scheduled_IDs) > len(jobs)):
        print('Some jobs were scheduled twice: #scheduled IDs: ' + str(len(scheduled_IDs)) +
              ', #IDs: ' + str(len(jobs)))
        return False
    if (len(scheduled_IDs) < len(jobs) or scheduled_IDs != [i for i in range(len(jobs))]):
        print('Some jobs were not scheduled: There were ' +
              str(len(scheduled_IDs)) + ' scheduled jobs vs ' + str(len(jobs)) + ' requested')
        return False
    # test space and time overlaps by creating a higher dimensional bin with the allocated spaces
    cmax = max(s[1] + jobs[s[0]][1] for s in schedule)
    bin_boundaries = list(boundaries) + [cmax]
    b = Bin(bin_boundaries)
    for job_id, start_time, space in schedule:
        b.spaces.add(ConvexSpace(space.coordinates + [start_time],
                                 space.boundaries + [jobs[job_id][1]], bin_boundaries))
    b.freelist = set()
    if (not b.test_possible()):
        print('Some jobs are overlapping in the schedule')
        return False
    return True

def powers_generator(power, min, max):
    start = 0 if min == 0 else int(math.floor(math.log(min, power)))
    end = 0 if max == 0 else int(math.floor(math.log(max, power))) + 1
    return (int(power ** x) for x in range(start, end + 1))
def powers2(min, max):
    return powers_generator(2, min, max)
def inverse_powers2(min, max):
    x = max + 1
    l = []
    while (x >= min):
        l.append(x)
        x = x // 2
    return reversed(l)

# for this series, need to pay attention since the span between
# min and max must be at least M
# paper by Lee and Lee (A simple on-line bin-packing algorithm) suggests M=12
def harmonic_generator(M, min, max):
    yield min
    last = min
    k = M
    while (k >= 1):
        n = int(math.floor((max - min) / k)) + min + 1
        if (n <= last):
            return
        yield n
        last = n
        k = k - 1
def harmonic12(min, max):
    return harmonic_generator(12, min, max)

def max_bin_time(b, jobs):
    return max(jobs[b.space_IDs[s]][1] for s in b.spaces)

def pack(jobs, request_sizes, initial_bins, dimensions, candidates, alpha):
    # get bin packing of the given sizes and assign a time to each bin according to max
    # time of jobs in each bin
    bins = ff_best_metric_first(dimensions, request_sizes, candidates, alpha, initial_bins)
    # associate max time with each bin and sort by increasing time
    return bins

# schedule a list of jobs (tuple of requested size and time) onto a torus with given boundaries
# return a list of tuples with the ID (index of the job), and a location in space and time of that ID
# the given time series is used to classify the jobs in time (packing is done for each class)
# this series must be strictly increasing, the first element must be lesser or equal to the
# given minimum and last element must be strictly greater than the given max value
# the given time for each job must be integers greater or equal to 1
def schedule_strict(dimensions, candidates, alpha, jobs, time_series_generator):
    max_time = max((j[1] for j in jobs), default=0)
    min_time = min((j[1] for j in jobs), default=1)
    current_time = 0
    schedule = []
    series = time_series_generator(min_time, max_time)
    last_time_slice = next(series)
    for time_slice in series:
        request_sizes = dict((i, j[0]) for i, j in enumerate(jobs)
                          if j[1] >= last_time_slice and j[1] < time_slice)
        last_time_slice = time_slice

        if (not request_sizes):
            continue

        bins = pack(jobs, request_sizes, [], dimensions, candidates, alpha)
        print('#bins: ' + str(len(bins)) + ', slice: ' + str(last_time_slice))
        for b in sorted(bins, key=lambda b: max_bin_time(b, jobs)):
            # plot_cuboids.plot_bin(b)
            schedule.extend((b.space_IDs[s], current_time, s) for s in b.spaces)
            current_time = current_time + max_bin_time(b, jobs)

    return schedule

# same as above but allow to pack the last bin of a packing in the next packing,
# while keeping it in the same shelf in terms of time series, if this reduces the Cmax
# TODO this strategy might still be worse than strict in some cases (to be confirmed)
def schedule_last_bin(dimensions, candidates, alpha, jobs, time_series_generator):
    max_time = max((j[1] for j in jobs), default=0)
    min_time = min((j[1] for j in jobs), default=1)
    series = time_series_generator(min_time, max_time)
    last_time_slice = next(series)
    bin_size = reduce(operator.mul, dimensions.values(), 1)
    initial_bin = None
    total_bins_time = []
    for time_slice in series:
        request_sizes = dict((i, j[0]) for i, j in enumerate(jobs)
                          if j[1] >= last_time_slice and j[1] < time_slice)
        last_time_slice = time_slice

        if (not request_sizes):
            continue

        # make a reference packing without initial bins
        bins_ref = pack(jobs, request_sizes, [], dimensions, candidates, alpha)
        skip_bins = 0
        if (initial_bin != None):
            # need to deep copy the initial bins such that they are not modified if we don't choose them
            new_bins = pack(jobs, request_sizes, [copy.deepcopy(initial_bin)],
                            dimensions, candidates, alpha)
            # check if this packing is better than reference packing
            if (sum(max_bin_time(b, jobs) for b in new_bins) <
                sum(max_bin_time(b, jobs) for b in bins_ref) + max_bin_time(initial_bin, jobs)):
                # remove last bin from total bins (this is the initial bin) and replace it by new copy
                # but with old time-slice
                skip_bins = 1
                bins_ref = new_bins
                t, b = total_bins_time.pop()
                total_bins_time.append((t, new_bins[0]))

        print('#bins: ' + str(len(bins_ref) - skip_bins) + ', slice: ' + str(last_time_slice))
        # possibly skip the initial bins at the start since we already added them in the last packing
        total_bins_time.extend((time_slice, b) for b in bins_ref[skip_bins:])
        # the new initial bin is the last one from this packing
        initial_bin = bins_ref[-1]

        # bins = pack(jobs, request_sizes, [] if initial_bin == None else [initial_bin], dimensions, candidates, alpha)
        # total_bins_time.extend((time_slice, b) for b in bins[(0 if initial_bin == None else 1):])
        # initial_bin = bins[-1]

    current_time = 0
    schedule = []
    # sort the total bins by max time, then by time slice again to keep the order of time slices
    for t, b in sorted(sorted(total_bins_time, key=lambda x: max_bin_time(x[1], jobs)), key=lambda x: x[0]):
        schedule.extend((b.space_IDs[s], current_time, s) for s in b.spaces)
        current_time = current_time + max_bin_time(b, jobs)

    return schedule

def parse_datetime(dt):
    d = datetime.datetime(2000, 1, 1)
    return d.strptime(dt, '%Y-%m-%d.%H:%M:%S').replace(tzinfo=pytz.utc)

# TODO issues with comparison
# - on-line vs off-line
# - at start time, machine is full for actual schedule, for our schedule we consider it empty
# - there can be jobs which start during the time considered, but which we do not include because
# they exit after all the jobs considered have exited (these jobs can demand resources on the
# actual schedule which we don't consider in ours)
def print_stats(schedule, jobs, actual_sched, boundaries, possible_sizes):
    cmax = max(s[1] + jobs[s[0]][1] for s in schedule)
    actual_min_start = min(s[1] for s in actual_sched)
    actual_cmax = max((s[1] + s[2]) - actual_min_start for s in actual_sched)
    actual_cmax = int(actual_cmax.total_seconds())
    max_job_time = max(j[1] for j in jobs)
    jobs_work = sum(j[0] * j[1] for j in jobs)
    # lower bound for Cmax (non-convex)
    cmax_lb_non_convex = max(max_job_time, jobs_work / reduce(operator.mul, boundaries, 1))
    jobs_convex_work = sum(min(possible_sizes[j[0]]) * j[1] for j in jobs)
    cmax_lb_convex = max(max_job_time, jobs_convex_work / reduce(operator.mul, boundaries, 1))

    max_overalloc = max((reduce(operator.mul, s.boundaries, 1) - jobs[i][0]) * 100 / jobs[i][0]
                        for i, t, s in schedule)
    avg_overalloc = sum((reduce(operator.mul, s.boundaries, 1) - jobs[i][0]) * 100 / jobs[i][0]
                        for i, t, s in schedule) / len(jobs)

    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    avg_metric = sum(metric_ad(s.boundaries, False, dimensions) for i, t, s in schedule) / len(jobs)
    avg_opt_metric = sum(opt_ad(j[0], False, dimensions) for j in jobs) / len(jobs)

    print('Cmax lower bound (non-convex): ' + str(cmax_lb_non_convex) +
          ', Cmax lower bound (convex): ' + str(cmax_lb_convex) +
          ', Cmax of schedule: ' + str(cmax) +
          ' (loss of ' + str((cmax - cmax_lb_non_convex) * 100 / cmax_lb_non_convex) +
          '% compared to non-convex lower bound, ' +
          str((cmax - cmax_lb_convex) * 100 / cmax_lb_convex) +
          '% compared to convex lower bound), actual Cmax: ' + str(actual_cmax) +
          ', max overalloc: ' + str(max_overalloc) + '% / average overalloc: ' + str(avg_overalloc) +
          ', average metric: ' + str(avg_metric) + ', average optimal metric: ' + str(avg_opt_metric))

def get_resource_ids(space, boundaries):
    # coordinates here can go outside of boundaries of the torus. only when the ID is computed,
    # they are adjusted to torus boundaries
    def get_id(coordinates):
        return sum((c % boundaries[i]) *
                   reduce(operator.mul, boundaries[(i+1):], 1) for i, c in enumerate(coordinates))

    def incr_coords(coordinates, space):
        # increment coordinates except for last dimension
        if (len(coordinates) <= 1):
            return coordinates
        for d in reversed(range(len(coordinates) - 1)):
            coordinates[d] = coordinates[d] + 1
            if (coordinates[d] < space.coordinates[d] + space.boundaries[d]):
                return coordinates
            coordinates[d] = space.coordinates[d]
        return coordinates

    # make slices of coordinates aligned with last dimension
    last_boundary = space.boundaries[-1]
    resources = [(get_id(space.coordinates), last_boundary)]
    coords = incr_coords(list(space.coordinates), space)
    while (coords != space.coordinates):
        resources.append((get_id(coords), last_boundary))
        coords = incr_coords(coords, space)
    return resources

def save_jedule_output(boundaries, allocations, jobs, filename):
    total = reduce(operator.mul, boundaries, 1)
    ms = sc.MoldSchedule(total)
    for job_id, start_time, space in allocations:
        task = sc.TaskRect(job_id, 'transfer') # nodes are of transfer type
        task.set_times(start_time, start_time + jobs[job_id][1])
        task.set_procs(get_resource_ids(space, boundaries))
        ms.add_task_rect(task)
    output = ms.get_jedule_output()
    with open(filename, 'w') as outFile:
        output.dump(outFile)

def perform_schedule(filename, boundaries, time_series_generator, alpha):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates(False, dimensions)
    possible_sizes = get_possible_sizes(dimensions, candidates, alpha)

    jobs_sched = []
    actual_sched = []
    with open(filename) as infile:
        for line in infile:
            startime, endtime, size, requestedtime, walltime = line.split()
            startime, endtime = parse_datetime(startime), parse_datetime(endtime)
            size, walltime = int(size), datetime.timedelta(seconds=int(walltime))
            jobs_sched.append((size, int(walltime.total_seconds())))
            actual_sched.append((size, startime, walltime))
    sched = schedule_strict(dimensions, candidates, alpha, jobs_sched, time_series_generator)
    if (not test_sched(sched, jobs_sched, boundaries)):
        return
    print_stats(sched, jobs_sched, actual_sched, boundaries, possible_sizes)
    outFilePrefix = ''
    if (len(sys.argv) >= 3):
        outFilePrefix = sys.argv[2]
    save_jedule_output(boundaries, sched, jobs_sched, outFilePrefix + 'strict.jed')
    sched_last = schedule_last_bin(dimensions, candidates, alpha, jobs_sched, time_series_generator)
    if (not test_sched(sched_last, jobs_sched, boundaries)):
        return
    print_stats(sched_last, jobs_sched, actual_sched, boundaries, possible_sizes)
    save_jedule_output(boundaries, sched_last, jobs_sched, outFilePrefix + 'last_bin.jed')

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print('Please provide a simplified Blue Waters trace file to generate a schedule from')
        sys.exit(0)
    # alpha of 0.15 gave best results with packing
    perform_schedule(sys.argv[1], [24,24,24], inverse_powers2, 0.15)
