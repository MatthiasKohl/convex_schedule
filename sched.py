#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import sys
import math
import itertools
import datetime
import copy
from functools import reduce
import operator
import pytz
import Jedule.schedule_common as sc
from ffd import ffEachBestMetricFirst, getPossibleSizes
from shapes import char_range, shape_candidates
from packingbin import Bin
from packingspace import ConvexSpace
import plot_cuboids

def testSchedule(schedule, jobs, boundaries):
    scheduledIDs = sorted(i for i, t, s in schedule)
    if (len(scheduledIDs) > len(jobs)):
        print('Some jobs were scheduled twice: #scheduled IDs: ' + str(len(scheduledIDs)) +
              ', #IDs: ' + str(len(jobs)))
        return False
    if (len(scheduledIDs) < len(jobs) or scheduledIDs != [i for i in range(len(jobs))]):
        print('Some jobs were not scheduled: There were ' +
              str(len(scheduledIDs)) + ' scheduled jobs vs ' + str(len(jobs)) + ' requested')
        return False
    # test space and time overlaps by creating a higher dimensional bin with the allocated spaces
    cmax = max(s[1] + jobs[s[0]][1] for s in schedule)
    binBoundaries = list(boundaries) + [cmax]
    b = Bin(binBoundaries)
    for jobId, startTime, space in schedule:
        b.spaces.add(ConvexSpace(space.coordinates + [startTime],
                                 space.boundaries + [jobs[jobId][1]], binBoundaries))
    b.freelist = set()
    if (not b.testPossible()):
        print('Some jobs are overlapping in the schedule')
        return False
    return True

def powersGenerator(power, nMin, nMax):
    start = 0 if nMin == 0 else int(math.floor(math.log(nMin, power)))
    end = 0 if nMax == 0 else int(math.floor(math.log(nMax, power))) + 1
    return (int(power ** x) for x in range(start, end + 1))
def powers2(nMin, nMax):
    return powersGenerator(2, nMin, nMax)
def inversePowers2(nMin, nMax):
    x = nMax + 1
    l = []
    while (x >= nMin):
        l.append(x)
        x = x // 2
    return reversed(l)

# for this series, need to pay attention since the span between
# nMin and nMax must be at least M
# paper by Lee and Lee (A simple on-line bin-packing algorithm) suggests M=12
def harmonicGenerator(M, nMin, nMax):
    yield nMin
    last = nMin
    k = M
    while (k >= 1):
        n = int(math.floor((nMax - nMin) / k)) + nMin + 1
        if (n <= last):
            return
        yield n
        last = n
        k = k - 1
def harmonic12(nMin, nMax):
    return harmonicGenerator(12, nMin, nMax)

def max_bin_time(b, jobs):
    return max(jobs[b.spaceIDs[s]][1] for s in b.spaces)

def pack(jobs, requestSizes, initialBins, dimensions, candidates, alpha):
    # get bin packing of the given sizes and assign a time to each bin according to max
    # time of jobs in each bin
    bins = ffEachBestMetricFirst(dimensions, requestSizes, candidates, alpha, initialBins)
    # associate max time with each bin and sort by increasing time
    return bins

# schedule a list of jobs (tuple of requested size and time) onto a torus with given boundaries
# return a list of tuples with the ID (index of the job), and a location in space and time of that ID
# the given time series is used to classify the jobs in time (packing is done for each class)
# this series must be strictly increasing, the first element must be lesser or equal to the
# given minimum and last element must be strictly greater than the given max value
# the given time for each job must be integers greater or equal to 1
def schedule_strict(dimensions, candidates, alpha, jobs, time_series_generator):
    maxTime = max((j[1] for j in jobs), default=0)
    minTime = min((j[1] for j in jobs), default=1)
    currentTime = 0
    schedule = []
    series = time_series_generator(minTime, maxTime)
    lastTimeSlice = next(series)
    for timeSlice in series:
        requestSizes = dict((i, j[0]) for i, j in enumerate(jobs)
                          if j[1] >= lastTimeSlice and j[1] < timeSlice)
        lastTimeSlice = timeSlice

        if (not requestSizes):
            continue

        bins = pack(jobs, requestSizes, [], dimensions, candidates, alpha)
        print('#bins: ' + str(len(bins)) + ', slice: ' + str(lastTimeSlice))
        for b in sorted(bins, key=lambda b: max_bin_time(b, jobs)):
            plot_cuboids.plot_bin(b)
            schedule.extend((b.spaceIDs[s], currentTime, s) for s in b.spaces)
            currentTime = currentTime + max_bin_time(b, jobs)

    return schedule

# same as above but allow to pack the last bin of a packing in the next packing,
# while keeping it in the same shelf in terms of time series, if this reduces the Cmax
# TODO this strategy might still be worse than strict in some cases (to be confirmed)
def schedule_last_bin(dimensions, candidates, alpha, jobs, time_series_generator):
    maxTime = max((j[1] for j in jobs), default=0)
    minTime = min((j[1] for j in jobs), default=1)
    series = time_series_generator(minTime, maxTime)
    lastTimeSlice = next(series)
    binSize = reduce(operator.mul, dimensions.values(), 1)
    initialBin = None
    totalBinsTime = []
    for timeSlice in series:
        requestSizes = dict((i, j[0]) for i, j in enumerate(jobs)
                          if j[1] >= lastTimeSlice and j[1] < timeSlice)
        lastTimeSlice = timeSlice

        if (not requestSizes):
            continue

        # make a reference packing without initial bins
        binsRef = pack(jobs, requestSizes, [], dimensions, candidates, alpha)
        skipBins = 0
        if (initialBin != None):
            # need to deep copy the initial bins such that they are not modified if we don't choose them
            newBins = pack(jobs, requestSizes, [copy.deepcopy(initialBin)],
                            dimensions, candidates, alpha)
            # check if this packing is better than reference packing
            if (sum(max_bin_time(b, jobs) for b in newBins) <
                sum(max_bin_time(b, jobs) for b in binsRef) + max_bin_time(initialBin, jobs)):
                # remove last bin from total bins (this is the initial bin) and replace it by new copy
                # but with old time-slice
                skipBins = 1
                binsRef = newBins
                t, b = totalBinsTime.pop()
                totalBinsTime.append((t, newBins[0]))

        print('#bins: ' + str(len(binsRef) - skipBins) + ', slice: ' + str(lastTimeSlice))
        # possibly skip the initial bins at the start since we already added them in the last packing
        totalBinsTime.extend((timeSlice, b) for b in binsRef[skipBins:])
        # the new initial bin is the last one from this packing
        initialBin = binsRef[-1]

        # bins = pack(jobs, requestSizes, [] if initialBin == None else [initialBin], dimensions, candidates, alpha)
        # totalBinsTime.extend((timeSlice, b) for b in bins[(0 if initialBin == None else 1):])
        # initialBin = bins[-1]

    currentTime = 0
    schedule = []
    # sort the total bins by max time, then by time slice again to keep the order of time slices
    for t, b in sorted(sorted(totalBinsTime, key=lambda x: max_bin_time(x[1], jobs)), key=lambda x: x[0]):
        schedule.extend((b.spaceIDs[s], currentTime, s) for s in b.spaces)
        currentTime = currentTime + max_bin_time(b, jobs)

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
# - how to calculate Cmax of the actual schedule and which times should we use for our schedule ?
# at the moment, Cmax is calculated by using actual start time + requested time
# however, our schedule uses the requested time for all jobs. the on-line schedule is re-calculated
# anytime a job exits, such that the Cmax is more dependent on the actual wall-times of jobs
# whereas we are dependent on the requested wall-time (this can be much more)
def printStats(schedule, jobs, actual_sched, boundaries, possibleSizes):
    cmax = max(s[1] + jobs[s[0]][1] for s in schedule)
    actual_min_start = min(s[1] for s in actual_sched)
    actual_cmax = max((s[1] + s[2]) - actual_min_start for s in actual_sched)
    actual_cmax = int(actual_cmax.total_seconds())
    max_job_time = max(j[1] for j in jobs)
    jobs_work = sum(j[0] * j[1] for j in jobs)
    # lower bound for Cmax (non-convex)
    cmax_lb_non_convex = max(max_job_time, jobs_work / reduce(operator.mul, boundaries, 1))
    jobs_convex_work = sum(min(possibleSizes[j[0]]) * j[1] for j in jobs)
    cmax_lb_convex = max(max_job_time, jobs_convex_work / reduce(operator.mul, boundaries, 1))
    print('Cmax lower bound (non-convex): ' + str(cmax_lb_non_convex) +
          ', Cmax lower bound (convex): ' + str(cmax_lb_convex) +
          ', Cmax of schedule: ' + str(cmax) +
          ' (loss of ' + str((cmax - cmax_lb_non_convex) * 100 / cmax_lb_non_convex) +
          '% compared to non-convex lower bound, ' +
          str((cmax - cmax_lb_convex) * 100 / cmax_lb_convex) +
          '% compared to convex lower bound), actual Cmax: ' + str(actual_cmax))

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
    lastBoundary = space.boundaries[-1]
    resources = [(get_id(space.coordinates), lastBoundary)]
    coords = incr_coords(list(space.coordinates), space)
    while (coords != space.coordinates):
        resources.append((get_id(coords), lastBoundary))
        coords = incr_coords(coords, space)
    return resources

def save_jedule_output(boundaries, allocations, jobs, filename):
    total = reduce(operator.mul, boundaries, 1)
    ms = sc.MoldSchedule(total)
    for jobId, startTime, space in allocations:
        task = sc.TaskRect(jobId, 'transfer') # nodes are of transfer type
        task.set_times(startTime, startTime + jobs[jobId][1])
        task.set_procs(get_resource_ids(space, boundaries))
        ms.add_task_rect(task)
    output = ms.get_jedule_output()
    with open(filename, 'w') as outFile:
        output.dump(outFile)

def perform_schedule(filename, boundaries, time_series_generator, alpha):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)
    possibleSizes = getPossibleSizes(dimensions, candidates, alpha)

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
    if (not testSchedule(sched, jobs_sched, boundaries)):
        return
    printStats(sched, jobs_sched, actual_sched, boundaries, possibleSizes)
    # outFilePrefix = ''
    # if (len(sys.argv) >= 3):
    #     outFilePrefix = sys.argv[2]
    # save_jedule_output(boundaries, sched, jobs_sched, outFilePrefix + 'strict.jed')
    # sched_last = schedule_last_bin(dimensions, candidates, alpha, jobs_sched, time_series_generator)
    # if (not testSchedule(sched_last, jobs_sched, boundaries)):
    #     return
    # printStats(sched_last, jobs_sched, actual_sched, boundaries, possibleSizes)
    # save_jedule_output(boundaries, sched_last, jobs_sched, outFilePrefix + 'last_bin.jed')

if __name__ == '__main__':
    # alpha of 0.15 gave best results with packing
    perform_schedule(sys.argv[1], [24,24,24], inversePowers2, 0.15)

# RESULTS using requested walltimes
# bw_request_sizes_20160405_5000.txt
# Cmax lower bound non-convex: 172800, convex: 172800, actual Cmax: 302607
# strategy strict: Cmax of schedule: 441000 (loss of ~155% compared to lower bound)
# strategy last_bin: Cmax of schedule: 207900 (loss of ~20.3% compared to lower bounds)

# bw_request_sizes_20160406_5000.txt
# Cmax lower bound non-convex: 172800, convex: 172800, actual Cmax: 333918
# strategy strict: Cmax of schedule: 473880 (loss of ~174% compared to lower bound)
# strategy last_bin: Cmax of schedule: 327480 (loss of ~89.5% compared to lower bounds)

# bw_request_sizes_20160406_20000.txt
# Cmax lower bound non-convex: 387260, convex: 396851, actual Cmax: 508345
# strategy strict: Cmax of schedule: 1000050 (loss of ~158%/152% compared to lower bounds)
# strategy last_bin: Cmax of schedule: 764580 (loss of ~97.4%/92.7% compared to lower bounds)

# RESULTS using processing times (wall-time = run-time)
# bw_request_sizes_20160405_5000.txt
# Cmax lower bound non-convex: 103557, convex: 103557, actual Cmax: 130284
# strategy strict: Cmax of schedule: 235544 (loss of ~127% compared to lower bound)
# strategy last_bin: Cmax of schedule: 150883 (loss of ~45.7% compared to lower bounds)

# bw_request_sizes_20160406_5000.txt
# Cmax lower bound non-convex: 72011, convex: 72011, actual Cmax: 161359
# strategy strict: Cmax of schedule: 215555 (loss of ~199% compared to lower bound)
# strategy last_bin: Cmax of schedule: 109926 (loss of ~52.7% compared to lower bounds)

# bw_request_sizes_20160406_20000.txt
# Cmax lower bound non-convex: 173696, convex: 173696, actual Cmax: 335835
# strategy strict: Cmax of schedule: 525719 (loss of ~203% compared to lower bounds)
# strategy last_bin: Cmax of schedule: 485429 (loss of ~179% compared to lower bounds)

# RESULTS using processing times (wall-time = run-time) with inversed powers of 2
# bw_request_sizes_20160405_5000.txt
# Cmax lower bound non-convex: 103557, convex: 103557, actual Cmax: 130284
# strategy strict: Cmax of schedule: 207036 (loss of ~100% compared to lower bound)
# strategy last_bin: Cmax of schedule: 127547 (loss of ~23.2% compared to lower bounds)

# bw_request_sizes_20160406_5000.txt
# Cmax lower bound non-convex: 72011, convex: 72011, actual Cmax: 161359
# strategy strict: Cmax of schedule: 157863 (loss of ~119% compared to lower bounds)
# strategy last_bin: Cmax of schedule: 134803 (loss of ~87.2% compared to lower bounds)

# bw_request_sizes_20160406_20000.txt
# Cmax lower bound non-convex: 173696, convex: 173696, actual Cmax: 335835
# strategy strict: Cmax of schedule: 460801 (loss of ~165% compared to lower bounds)
# strategy last_bin: Cmax of schedule: 415969 (loss of ~139% compared to lower bounds)

# this strategy is not used as it does not make a lot of sense (cutOff is arbitrary etc)
# and experiments show it is performing worse than last_bin strategy
# def schedule_bin_size(boundaries, jobs, time_series_generator):
#     # get shape candidates
#     dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
#     candidates = shape_candidates("./primes", False, dimensions)

#     maxTime = max((j[1] for j in jobs), default=0)
#     minTime = min((j[1] for j in jobs), default=1)
#     currentTime = 0
#     schedule = []
#     series = time_series_generator(minTime, maxTime)
#     lastTimeSlice = next(series)
#     requestSizes = {}
#     cutOff = reduce(operator.mul, boundaries, 1) * 0.8
#     def _sched_requests(requestSizes, currentTime, schedule):
#         if (not requestSizes):
#             return currentTime, schedule
#         binsTime = pack(jobs, requestSizes, [], dimensions, candidates, 0.15)
#         for b, maxT in sorted(binsTime, key=lambda x: x[1]):
#             schedule.extend((b.spaceIDs[s], currentTime, s) for s in b.spaces)
#             currentTime = currentTime + maxT
#         return currentTime, schedule
#     for timeSlice in series:
#         sliceSizes = dict((i, j[0]) for i, j in enumerate(jobs)
#                           if j[1] >= lastTimeSlice and j[1] < timeSlice)
#         lastTimeSlice = timeSlice

#         sliceSize = sum(s for s in sliceSizes.values())
#         if (sum(s for s in requestSizes.values()) + sliceSize < cutOff):
#             requestSizes.update(sliceSizes)
#             continue
#         # schedule the last requests, then possibly the new ones
#         currentTime, schedule = _sched_requests(requestSizes, currentTime, schedule)
#         if (sliceSize < cutOff):
#             requestSizes = sliceSizes
#             continue
#         currentTime, schedule = _sched_requests(sliceSizes, currentTime, schedule)
#         # reset requestSizes
#         requestSizes = {}

#     # make sure that the last requests are scheduled
#     currentTime, schedule = _sched_requests(requestSizes, currentTime, schedule)

#     return schedule

# TODO this strategy is not smart, since we waste space by packing shorter jobs into the long bins first
# same as above but start with biggest jobs and allow to pack the smaller jobs into all bins
# of the bigger jobs. this will impact negatively on metrics like flow, since short jobs might
# get scheduled along with long jobs (in time) at the very end of the schedule
# def schedule_all_bins(boundaries, jobs, time_series_generator):
#         # get shape candidates
#     dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
#     candidates = shape_candidates("./primes", False, dimensions)

#     maxTime = max((j[1] for j in jobs), default=0)
#     minTime = min((j[1] for j in jobs), default=1)
#     # go from upper bound to lower of time slice
#     series = reversed([t for t in time_series_generator(minTime, maxTime)])
#     lastTimeSlice = next(series)
#     binSize = reduce(operator.mul, boundaries, 1)
#     initialBins = []
#     totalBinsTime = []
#     for timeSlice in series:
#         requestSizes = dict((i, j[0]) for i, j in enumerate(jobs)
#                           if j[1] < lastTimeSlice and j[1] >= timeSlice)
#         lastTimeSlice = timeSlice

#         if (not requestSizes):
#             continue

#         totalBinsTime = pack(jobs, requestSizes, initialBins, dimensions, candidates, 0.15)
#         initialBins = [b for b,maxT in totalBinsTime]

#     currentTime = 0
#     schedule = []
#     # sort the total bins by max time
#     for b, maxT in sorted(totalBinsTime, key=lambda x: x[1]):
#         schedule.extend((b.spaceIDs[s], currentTime, s) for s in b.spaces)
#         currentTime = currentTime + maxT

#     return schedule
