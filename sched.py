#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import sys
import math
import itertools
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

    maxTime = max(j[1] for j in jobs, default=0)
    minTime = min(j[1] for j in jobs, default=1)
    currentSchedTime = 0
    schedJobs = []
    series = time_series_generator(minTime, maxTime)
    lastTimeSlice = next(series)
    for timeSlice in series:
        # TODO possibly don't just consider the jobs in this slice. instead,
        # if their volume is too small, port them over to next bin
        # volume too small could be smaller than 1 bin
        requestSizes = dict((i, j[0]) for i, j in enumerate(jobs)
                            if j[1] >= lastTimeSlice and j[1] < timeSlice)
        lastTimeSlice = timeSlice
        if (not requestSizes):
            continue

        # get bin packing of the given sizes and assign a time to each bin according to max
        # time of jobs in each bin
        bins = ffEachBestMetricFirst(dimensions, requestSizes, candidates, 0.15)
        # associate max time with each bin and sort by increasing time
        binsTime = ((b, max(jobs[b.spaceIDs[s]][1] for s in b.spaces)) for b in bins)
        for (b, maxT in sorted(binsTime, key=lambda x: x[1])):
            schedJobs.extend((b.spaceIDs[s], currentSchedTime, s) for s in b.spaces)
            currentSchedTime = currentSchedTime + maxT

    return schedJobs

