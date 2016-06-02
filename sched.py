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

def pack(jobs, schedJobs, currentSchedTime, requestSizes, dimensions, candidates, alpha):
    # get bin packing of the given sizes and assign a time to each bin according to max
    # time of jobs in each bin
    bins = ffEachBestMetricFirst(dimensions, requestSizes, candidates, alpha)
    # associate max time with each bin and sort by increasing time
    binsTime = ((b, max(jobs[b.spaceIDs[s]][1] for s in b.spaces)) for b in bins)
    for b, maxT in sorted(binsTime, key=lambda x: x[1]):
        schedJobs.extend((b.spaceIDs[s], currentSchedTime, s) for s in b.spaces)
        currentSchedTime = currentSchedTime + maxT
    return schedJobs, currentSchedTime

# schedule a list of jobs (tuple of requested size and time) onto a torus with given boundaries
# return a list of tuples with the ID (index of the job), and a location in space and time of that ID
# the given time series is used to classify the jobs in time (packing is done for each class)
# this series must be strictly increasing, the first element must be lesser or equal to the
# given minimum and last element must be strictly greater than the given max value
# the given time for each job must be integers greater or equal to 1
def schedule(boundaries, jobs, time_series_generator):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)

    maxTime = max((j[1] for j in jobs), default=0)
    minTime = min((j[1] for j in jobs), default=1)
    currentSchedTime = 0
    schedJobs = []
    series = time_series_generator(minTime, maxTime)
    lastTimeSlice = next(series)
    requestSizes = {}
    binSize = reduce(operator.mul, boundaries, 1)
    for timeSlice in series:
        # TODO possibly don't just consider the jobs in this slice. instead,
        # if their volume is too small, port them over to next bin
        # volume too small could be smaller than 1 bin
        sliceSizes = dict((i, j[0]) for i, j in enumerate(jobs)
                          if j[1] >= lastTimeSlice and j[1] < timeSlice)
        lastTimeSlice = timeSlice

        # if (not sliceSizes):
        #     continue

        if (sum(s for s in requestSizes.values()) +
            sum(s for s in sliceSizes.values()) < binSize*0.8):
            requestSizes.update(sliceSizes)
            continue
        elif (requestSizes):
            # pack the last requests, then set requests to new
            schedJobs, currentSchedTime = pack(jobs, schedJobs, currentSchedTime,
                                               requestSizes, dimensions, candidates, 0.15)
        requestSizes = sliceSizes

        schedJobs, currentSchedTime = pack(jobs, schedJobs, currentSchedTime,
                                           requestSizes, dimensions, candidates, 0.15)
        # reset requestSizes
        requestSizes = {}

    return schedJobs

def parse_datetime(dt):
    d = datetime.datetime(2000, 1, 1)
    return d.strptime(dt, '%Y-%m-%d.%H:%M:%S').replace(tzinfo=pytz.utc)

def printStats(schedule, jobs, actual_sched):
    cmax = max(s[1] + jobs[s[0]][1] for s in schedule)
    actual_min_start = min(s[1] for s in actual_sched)
    actual_cmax = max((s[1] + s[2]) - actual_min_start for s in actual_sched)
    actual_cmax = int(actual_cmax.total_seconds())
    print('Cmax of schedule: ' + str(cmax) + ', actual Cmax: ' + str(actual_cmax))

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
    sched = schedule(boundaries, jobs_sched, time_series_generator)
    printStats(sched, jobs_sched, actual_sched)
    return sched

# TODO test
if __name__ == '__main__':
    perform_schedule(sys.argv[1], [24,24,24], powers2)
