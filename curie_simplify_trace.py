#! /usr/bin/python3
# -*- encoding: utf-8 -*-

# simple script to simplify the Curie trace
# the considered jobs must have a submit time greater than 35000000, since before, the machine
# is not utilized enough
# the output is a space separated table with requested size and requested time for each job
# Arguments for the script are as follows:
# 1 - size that request sizes should be scaled to or clamped at
# 2 - 0 if request sizes should be clamped, else they will be scaled
# 3 - number of jobs to output


import sys

if (len(sys.argv) < 5):
    print('Please provide a filename to process, a clamp/scale value, ' +
          'an integer whether to clamp or scale (0/1) and a number of jobs to output')
    sys.exit(0)

filename = sys.argv[1]
try:
    total_size = int(sys.argv[2])
    is_clamp = True if int(sys.argv[3]) == 0 else False
    num_jobs = int(sys.argv[4])
except:
    print('Provided arguments are invalid or cannot open given file <' + filename + '>')
    sys.exit(2)


# SWF format is as follows (according to thesis of Emeras 28th august 2013):
# Headers are preceded by ';' and are ignored here.
# Fields are whitespace-separated as follows:
# 1. Job Number
# 2. Submit Time
# 3. Wait Time
# 4. Run Time
# 5. Number of Allocated Processors
# 6. Average CPU Time Used
# 7. Used Memory
# 8. Requested Number of Processors
# 9. Requested Time
# 10. Requested Memory
# 11. Status
# 12. User ID
# 13. Group ID
# 14. Executable (Application) Number
# 15. Queue Number
# 16. Partition Number
# 17. Preceding Job Number
# 18. Think Time from Preceding Job

# for the CEA curie trace, this is the max possible size and min submit time
max_request_size, min_submit_time = 80000, 35000000
# write to output file
out_filename = 'request_sizes_' + ('clamped_' if is_clamp else 'scaled_') + str(num_jobs) + '.txt'
print('Writing to ' + out_filename)

with open(filename) as infile:
    with open(out_filename, 'w') as outfile:
        out_sum = 0
        for line in infile:
            if (line.startswith(';')):
                continue
            fields = line.split()
            submit_time = int(fields[1])
            # only consider jobs after startup
            if (submit_time < min_submit_time):
               continue
            request_size = int(fields[7])
            if (is_clamp and request_size > total_size):
                continue
            elif (not is_clamp):
                request_size = min(int((request_size / max_request_size) * total_size) + 1, total_size)
            request_time = int(fields[8])
            outfile.write(str(request_size) + ' ' + str(request_time) + '\n')
            out_sum = out_sum + 1
            if (out_sum >= num_jobs):
                break
