#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import sys
import itertools
from ffd import ffEachBestMetricFirst
from shapes import char_range, shape_candidates

def powersOfTwo():
    return (2 ** x for x in itertools.count())

# schedule a list of jobs (tuple of requested size and time) onto a torus with given boundaries
# return a list of tuples with the ID (index of the job), and a location in space and time of that ID
# the given time series is used to classify the jobs in time (packing is done for each class)
# this series must be strictly increasing and the first element must be strictly greater than 1
# the given time for each job must be integers greater or equal to 1
def schedule(boundaries, jobs, timeSeries):
    # get shape candidates
    dimensions = {c: boundaries[i] for i, c in enumerate(char_range('x', len(boundaries)))}
    candidates = shape_candidates("./primes", False, dimensions)

    # classify jobs by time according to powers of 2
    maxTime = max(j[1] for j in jobs, default=0)
    currentSchedTime = 0
    schedJobs = []
    series = timeSeries()
    lastTimeSlice = next(series)
    for timeSlice in series:
        if (lastTimeSlice > maxTime):
            break

        requestSizes = dict((i, j[0]) for i, j in enumerate(jobs)
                            if j[1] >= lastTimeSlice and j[1] < timeSlice)
        lastTimeSlice = timeSlice
        if (not requestSizes):
            continue

        # get bin packing of the given sizes and assign a time to each bin according to max
        # time of jobs in each bin
        bins = ffEachBestMetricFirst(dimensions, requestSizes, candidates, 0.15)
        for (b in bins):
            schedJobs.extend((b.spaceIDs[s], currentSchedTime, s) for s in b.spaces)
            maxTime = max(jobs[b.spaceIDs[s]][1] for s in b.spaces)
            currentSchedTime = currentSchedTime + maxTime

    return schedJobs

