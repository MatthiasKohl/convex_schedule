#! /usr/bin/python3
# -*- encoding: utf-8 -*-

# taken from Raphaël Bleuse:
# bw_scripts/scripts/analysis/logs

from enum import Enum, unique
import lzma
import collections

# timestamps & timezones management
import datetime
import pytz
from pytz.exceptions import AmbiguousTimeError

#---------------------------------------

@unique
class Event(Enum):
    """
    Event types encountered in the accounting records.
    Enumeration of the different event one may encounter in the torque
    accounting logs. These fields are documented in
    http://docs.adaptivecomputing.com/torque/5-0-2/help.htm#topics/torque/9-accounting/accountingRecords.htm
    """
    ABT    = 'A'    # job has been aborted by the server
    CHKPNT = 'C'    # job has been checkpointed and held
    DEL    = 'D'    # job has been deleted
    EXIT   = 'E'    # job has exited (either successfully or unsuccessfully)
    QUEUE  = 'Q'    # job has been submitted/queued
    RERUN  = 'R'    # attempt to rerun the job has been made
    RUN    = 'S'    # attempt to start the job has been made (if the job fails
                    # to properly start, it may have multiple job start records)
    RESTRT = 'T'    # attempt to restart the job (from checkpoint) has been
                    # made (if the job fails to properly start, it may have
                    # multiple job start records)


@unique
class State(Enum):
    """
    Enumeration of the different states for a job. These constants are defined
    in include/pbs_job.h
    """
    TRANSIT  = 0
    QUEUED   = 1
    HELD     = 2
    WAITING  = 3
    RUNNING  = 4
    EXITING  = 5
    COMPLETE = 6


class Job:
    def __init__(self, id):
        self.id = id
        self.name = None

        self.queue = None
        self.events = [] # list of (datetime, event_type) associated to self

        # timestamps are updated while events are parsed
        self.ctime = None # job is created
        self.etime = None # job becomes eligible to run
        self.qtime = None # job is queued
        self.start = None # job starts to run
        self.end = None   # job finishes

        # some more interesting statistics
        self.requested = {}
        self.used = {}
        self.exit_status = None

    def update(self, timestamp, event_type, infodict):
        self.events.append((timestamp, event_type))

        if event_type is Event.QUEUE:
            self.queue = infodict['queue']
            self.qtime = timestamp
        elif event_type is Event.RUN:
            pass

    @property
    def last_event(self):
        return max(self.events, key=lambda tup: tup[0])

    @property
    def first_event(self):
        return min(self.events, key=lambda tup: tup[0])


#---------------------------------------

class MonotonicTimestampParser:
    def __init__(self, tz=pytz.utc):
        self.tz = tz
        self.previous = None


    @staticmethod
    def _is_ambiguous(timestamp, tz=pytz.utc):
        try:
            tz.dst(timestamp)
            return False
        except AmbiguousTimeError:
            return True


    def __call__(self, date_string, format):
        naive_ts = datetime.datetime.strptime(date_string, format)

        if MonotonicTimestampParser._is_ambiguous(naive_ts, self.tz):
            if self.previous is None:
                raise AmbiguousTimeError('no reference to check for warp')

            # enforce DST policy when localizing: the monotonicity of the
            # parsed sequence of timestamps allows us to choose the correct DST
            # flag
            # we need to change DST flag when the naive timestamps 'go back in
            # time': it is done with the xor
            is_dst = bool(self.previous.dst())
            is_dst ^= self.previous.replace(tzinfo=None) > naive_ts
            self.previous = self.tz.localize(naive_ts, is_dst=is_dst)
        else:
            self.previous = self.tz.localize(naive_ts)

        return self.previous


def parse_payload(payload, payload_interpreters={}):
    """
    Interpret in a pythonic way the information payload associated to an event.
    This function stores the raw string if no interpreter is defined for a
    field.
    """
    # fallback if no interpreter is defined for a given field_name
    def _identity(data):
        return data

    # actually parse values
    infodict = {}
    for chunk in payload.split():
        field_name, _, data = chunk.partition('=')
        infodict[field_name] = \
                payload_interpreters.get(field_name, _identity)(data)

    return infodict

def parse_logline(logline):
    """
    Parse a log line from the torque accounting log.
    Torque maintains accounting records for batch jobs. This function parses a
    single line, whose format is described in
    http://docs.adaptivecomputing.com/torque/5-0-2/help.htm#topics/torque/9-accounting/accountingRecords.htm
    """
    import sys

    # define field specific interpreters
    def _timestamp(data):
        return datetime.datetime.fromtimestamp(int(data), tz=pytz.utc)

    def _hms_time(data):
        split_data = data.split(':')
        h, m, s = ([0, 0, 0] + split_data)[-3:] # handle incomplete data
        h, m, s = int(h), int(m), int(s)
        return datetime.timedelta(hours=h, minutes=m, seconds=s)

    _payload_interpreters = {
        'ctime': _timestamp,
        'etime': _timestamp,
        'qtime': _timestamp,
        'start': _timestamp,
        'end': _timestamp,
        'Resource_List.walltime': _hms_time,
        'resources_used.walltime': _hms_time,
        'Exit_status': int,
        'total_execution_slots': int,
        'unique_node_count': int,
        'Resource_List.nodect': int # added to get the requested node count
    }

    # retrieve fields
    timestamp, event_type, job_id, info = logline.strip().split(';')

    # interpret fields
    event_type = Event(event_type)
    infodict = parse_payload(info, _payload_interpreters)
    # timestamp's interpretation is left in the loop as it needs context to
    # choose the correct DST flag

    return (timestamp, event_type, job_id, infodict)

