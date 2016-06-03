#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import sys
import math
import itertools
import datetime
from functools import reduce
import operator
import pytz
from ffd import ffEachBestMetricFirst
from shapes import char_range, shape_candidates

def powersGenerator(power, nMin, nMax):
    start = 0 if nMin == 0 else int(math.floor(math.log(nMin, power)))
    end = 0 if nMax == 0 else int(math.floor(math.log(nMax, power))) + 1
    return (int(power ** x) for x in range(start, end + 1))
def powers2(nMin, nMax):
    return powersGenerator(2, nMin, nMax)

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

def pack(jobs, requestSizes, initialBins, dimensions, candidates, alpha):
    # get bin packing of the given sizes and assign a time to each bin according to max
    # time of jobs in each bin
    bins = ffEachBestMetricFirst(dimensions, requestSizes, candidates, alpha, initialBins)
    # associate max time with each bin and sort by increasing time
    return [(b, max(jobs[b.spaceIDs[s]][1] for s in b.spaces)) for b in bins]

# schedule a list of jobs (tuple of requested size and time) onto a torus with given boundaries
# return a list of tuples with the ID (index of the job), and a location in space and time of that ID
# the given time series is used to classify the jobs in time (packing is done for each class)
# this series must be strictly increasing, the first element must be lesser or equal to the
# given minimum and last element must be strictly greater than the given max value
# the given time for each job must be integers greater or equal to 1
def schedule_strict(boundaries, jobs, time_series_generator):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)

    maxTime = max((j[1] for j in jobs), default=0)
    minTime = min((j[1] for j in jobs), default=1)
    currentTime = 0
    schedule = []
    series = time_series_generator(minTime, maxTime)
    lastTimeSlice = next(series)
    binSize = reduce(operator.mul, boundaries, 1)
    for timeSlice in series:
        requestSizes = dict((i, j[0]) for i, j in enumerate(jobs)
                          if j[1] >= lastTimeSlice and j[1] < timeSlice)
        lastTimeSlice = timeSlice

        if (not requestSizes):
            continue

        binsTime = pack(jobs, requestSizes, [], dimensions, candidates, 0.15)
        print(len(binsTime))
        for b, maxT in sorted(binsTime, key=lambda x: x[1]):
            schedule.extend((b.spaceIDs[s], currentTime, s) for s in b.spaces)
            currentTime = currentTime + maxT

    return schedule

# same as above but allow to pack the last bin of a packing in the next packing,
# while keeping it in the same shelf in terms of time series
def schedule_last_bin(boundaries, jobs, time_series_generator):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)

    maxTime = max((j[1] for j in jobs), default=0)
    minTime = min((j[1] for j in jobs), default=1)
    series = time_series_generator(minTime, maxTime)
    lastTimeSlice = next(series)
    binSize = reduce(operator.mul, boundaries, 1)
    initialBins = []
    totalBinsTime = []
    for timeSlice in series:
        requestSizes = dict((i, j[0]) for i, j in enumerate(jobs)
                          if j[1] >= lastTimeSlice and j[1] < timeSlice)
        lastTimeSlice = timeSlice

        if (not requestSizes):
            continue

        binsTime = pack(jobs, requestSizes, initialBins, dimensions, candidates, 0.15)
        # skip the initial bins at the start since we already added them in the last packing
        totalBinsTime.extend((timeSlice, b, maxT) for b, maxT in binsTime[len(initialBins):])
        # the new initial bin is the last one from this packing
        initialBins = [binsTime[-1][0]]

    currentTime = 0
    schedule = []
    # sort the total bins by max time, then by time slice again to keep the order of time slices
    for t, b, maxT in sorted(sorted(totalBinsTime, key=lambda x: x[2]), key=lambda x: x[0]):
        schedule.extend((b.spaceIDs[s], currentTime, s) for s in b.spaces)
        currentTime = currentTime + maxT

    return schedule

# TODO this strategy is not smart, since we waste space by packing shorter jobs into the long bins first
# same as above but start with biggest jobs and allow to pack the smaller jobs into all bins
# of the bigger jobs. this will impact negatively on metrics like flow, since short jobs might
# get scheduled along with long jobs (in time) at the very end of the schedule
def schedule_all_bins(boundaries, jobs, time_series_generator):
        # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)

    maxTime = max((j[1] for j in jobs), default=0)
    minTime = min((j[1] for j in jobs), default=1)
    # go from upper bound to lower of time slice
    series = reversed([t for t in time_series_generator(minTime, maxTime)])
    lastTimeSlice = next(series)
    binSize = reduce(operator.mul, boundaries, 1)
    initialBins = []
    totalBinsTime = []
    for timeSlice in series:
        requestSizes = dict((i, j[0]) for i, j in enumerate(jobs)
                          if j[1] < lastTimeSlice and j[1] >= timeSlice)
        lastTimeSlice = timeSlice

        if (not requestSizes):
            continue

        totalBinsTime = pack(jobs, requestSizes, initialBins, dimensions, candidates, 0.15)
        initialBins = [b for b,maxT in totalBinsTime]

    currentTime = 0
    schedule = []
    # sort the total bins by max time
    for b, maxT in sorted(totalBinsTime, key=lambda x: x[1]):
        schedule.extend((b.spaceIDs[s], currentTime, s) for s in b.spaces)
        currentTime = currentTime + maxT

    return schedule


def parse_datetime(dt):
    d = datetime.datetime(2000, 1, 1)
    return d.strptime(dt, '%Y-%m-%d.%H:%M:%S').replace(tzinfo=pytz.utc)

def printStats(schedule, jobs, actual_sched, boundaries):
    cmax = max(s[1] + jobs[s[0]][1] for s in schedule)
    actual_min_start = min(s[1] for s in actual_sched)
    actual_cmax = max((s[1] + s[2]) - actual_min_start for s in actual_sched)
    actual_cmax = int(actual_cmax.total_seconds())
    max_job_time = max(j[1] for j in jobs)
    jobs_volume = sum(j[0] * j[1] for j in jobs)
    # lower bound for Cmax
    cmax_lb = max(max_job_time, jobs_volume / reduce(operator.mul, boundaries, 1))
    print('Cmax lower bound: ' + str(cmax_lb) + ', Cmax of schedule: ' + str(cmax) +
          ' (loss of ' + str((cmax - cmax_lb) * 100 / cmax_lb) +
          '% compared to lower bound), actual Cmax: ' + str(actual_cmax))

def testSchedule(schedule, jobs):
    scheduledIDs = sorted(i for i, t, s in schedule)
    if (len(scheduledIDs) > len(jobs)):
        print('Some jobs were scheduled twice')
        return False
    if (len(scheduledIDs) < len(jobs) or scheduledIDs != [i for i in range(len(jobs))]):
        print('Some jobs were not scheduled')
        return False
    # TODO add test for time overlaps (space overlaps should not happen,
    # since bin packing can be tested)
    return True

def perform_schedule(filename, boundaries, time_series_generator):
    jobs_sched = []
    actual_sched = []
    with open(filename) as infile:
        for line in infile:
            startime, endtime, size, walltime = line.split()
            startime, endtime = parse_datetime(startime), parse_datetime(endtime)
            size, walltime = int(size), datetime.timedelta(seconds=int(walltime))
            jobs_sched.append((size, int(walltime.total_seconds())))
            actual_sched.append((size, startime, walltime))
    # sched = schedule_strict(boundaries, jobs_sched, time_series_generator)
    # if (not testSchedule(sched, jobs_sched)):
    #     return
    # printStats(sched, jobs_sched, actual_sched, boundaries)
    # sched_last = schedule_last_bin(boundaries, jobs_sched, time_series_generator)
    # if (not testSchedule(sched_last, jobs_sched)):
    #     return
    # printStats(sched_last, jobs_sched, actual_sched, boundaries)
    sched_all = schedule_all_bins(boundaries, jobs_sched, time_series_generator)
    if (not testSchedule(sched_all, jobs_sched)):
        return
    printStats(sched_all, jobs_sched, actual_sched, boundaries)

if __name__ == '__main__':
    perform_schedule(sys.argv[1], [24,24,24], powers2)

# RESULTS
# bw_request_sizes_20160405_5000.txt
# Cmax lower bound: 172800, actual Cmax: 302607
# strategy strict: Cmax of schedule: 441000 (loss of ~155% compared to lower bound)
# strategy last_bin: Cmax of schedule: 201660 (loss of ~16.7% compared to lower bound)
# strategy all_bins: Cmax of schedule: 203820 (loss of ~18% compared to lower bound)
# UNUSED strategy bin_size: Cmax of schedule: 213780 (loss of ~23.7% compared to lower bound)

# bw_request_sizes_20160406_5000.txt
# Cmax lower bound: 172800, actual Cmax: 333918
# strategy strict: Cmax of schedule: 473880 (loss of ~174% compared to lower bound)
# strategy last_bin: Cmax of schedule: 293880 (loss of ~70% compared to lower bound)
# strategy all_bins: Cmax of schedule: 297180 (loss of ~72% compared to lower bound)
# UNUSED strategy bin_size: Cmax of schedule: 327900 (loss of ~89.8% compared to lower bound)

# bw_request_sizes_20160406_20000.txt
# Cmax lower bound: 387260, actual Cmax: 508345
# strategy strict: Cmax of schedule: 1000050 (loss of ~158% compared to lower bound)
# strategy last_bin: Cmax of schedule: 720810 (loss of ~86.1% compared to lower bound)
# strategy all_bins: Cmax of schedule: 760020 (loss of ~96.3% compared to lower bound)

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
