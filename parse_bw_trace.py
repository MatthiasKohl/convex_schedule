#! /usr/bin/python3
# -*- encoding: utf-8 -*-

from trace_bluewaters.scripts.analysis.logs.torque import Event, parse_logline
import sys
import datetime
import pytz

# unused
def parse_timestamp(ts):
    date_time = ts.split()
    date = date_time[0].split('/')
    M, D, Y = ([0, 0, 0] + date)[-3:] # handle incomplete data
    Y, M, D = int(Y), int(M), int(D)
    time = date_time[1].split(':')
    h, m, s = ([0, 0, 0] + time)[-3:] # handle incomplete data
    h, m, s = int(h), int(m), int(s)
    return datetime.datetime(year=Y, month=M, day=D, hour=h, minute=m, second=s, tzinfo=pytz.utc)

def print_max_time_stamp_after(start_timestamp):
    max_job_timestamp = start_timestamp
    i = 0
    n_missing = 0
    n_jobs = 0
    with open(sys.argv[1]) as infile:
        for line in infile:
            i = i + 1
            if (i % 100000 == 0):
                print((i // 100000) % 10, end='', flush=True)
            timestamp, event_type, job_id, infodict = parse_logline(line)
            if (event_type != Event.EXIT):
                continue
            n_jobs = n_jobs + 1
            if (not ('start' in infodict) or not ('end' in infodict)):
                n_missing = n_missing + 1
                continue
            # check if the job started before start_timestamp
            job_start_time = infodict['start']
            job_end_time = infodict['end']
            if (job_start_time < start_timestamp and job_end_time > max_job_timestamp):
                max_job_timestamp = job_end_time
    print()
    print(str(n_missing) + ' jobs (out of ' + str(n_jobs) +
          ' exiting jobs) have missing start or end')
    print('The maximal date-time of a job with a start-time before ' + str(start_timestamp) +
          ' and exiting after this date-time is ' + str(max_job_timestamp))

def create_jobs_log(min_job_timestamp, total_n_jobs, min_size):
    jobs = []
    n_missing = 0
    n_jobs = 0
    n_considered = 0
    with open(sys.argv[1]) as infile:
        for line in infile:
            timestamp, event_type, job_id, infodict = parse_logline(line)
            if (event_type != Event.EXIT):
                continue
            # check if the job exited after min_job_timestamp
            job_start_time = infodict['start']
            if (job_start_time < min_job_timestamp):
                continue
            n_considered = n_considered + 1
            if (not ('Resource_List.walltime' in infodict) or not ('Resource_List.nodect' in infodict)
                or not ('resources_used.walltime' in infodict)
                or infodict['resources_used.walltime'] < datetime.timedelta(seconds=1)):
                n_missing = n_missing + 1
                print('.', end='', flush=True)
                continue
            job_size = infodict['Resource_List.nodect']
            # for blue waters, there are 2 processing nodes per network node
            # we are only interested in network nodes
            job_size = job_size // 2 + job_size % 2
            if (job_size < min_size):
                continue
            start_time_str = job_start_time.strftime('%Y-%m-%d.%H:%M:%S')
            end_time_str = infodict['end'].strftime('%Y-%m-%d.%H:%M:%S')
            jobs.append((start_time_str, end_time_str, job_size,
                         int(infodict['Resource_List.walltime'].total_seconds()),
                         int(infodict['resources_used.walltime'].total_seconds())))
            n_jobs = n_jobs + 1
            if (n_jobs >= total_n_jobs):
                break
    if (n_missing > 0):
        print()
    print(str(n_missing) + ' (out of ' + str(n_considered) +
          ' considered) jobs had missing walltime or nodect and were discarded')
    print(str(n_considered-n_missing-total_n_jobs) + ' jobs had a size smaller than ' + str(min_size))
    with open(sys.argv[2], 'w') as outfile:
        for j in jobs:
            outfile.write(' '.join(map(str, j)) + '\n')

# take all timestamps after 2016-04-01
start_timestamp = datetime.datetime(year=2016, month=4, day=1, tzinfo=pytz.utc)
# the max end time stamp for jobs starting before 2016-04-01 is 2016-04-02 23:55:53
#print_max_time_stamp_after(start_timestamp)

# it seems like jobs are bigger starting from 2016-04-05/2016-04-06, so try to consider that
min_job_timestamp = datetime.datetime(year=2016, month=4, day=6, tzinfo=pytz.utc)
# the max end time stamp for jobs starting before 2016-04-05 is 2016-04-06 20:05:30
# the max end time stamp for jobs starting before 2016-04-06 is 2016-04-06 21:01:24
#print_max_time_stamp_after(min_job_timestamp)

create_jobs_log(min_job_timestamp, int(sys.argv[3]),
                int(sys.argv[4]) if len(sys.argv) > 4 else 0)
