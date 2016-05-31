#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import sys

if (len(sys.argv) < 5):
    print('Please provide a filename to process, a clamp/scale value, ' +
          'an integer whether to clamp or scale (0/1) and a number of jobs to output')
    sys.exit(0)

filename = sys.argv[1]
try:
    nTotalSize = int(sys.argv[2])
    isClamp = True if int(sys.argv[3]) == 0 else False
    nOutput = int(sys.argv[4])
    file = open(filename)
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

# for the CEA trace, this is the max possible size and min submit time
maxRequestSize, minSubmitTime = 80000, 35000000
# write to output file
outFileName = 'request_sizes_' + ('clamped_' if isClamp else 'scaled_') + str(nOutput) + '.txt'
outFile = open(outFileName, 'w')
print('Writing to ' + outFileName)

sumOut = 0
for line in file.readlines():
    if (line.startswith(';')):
        continue
    fields = line.split()
    submitTime = int(fields[1])
    # only consider jobs after startup
    if (submitTime < minSubmitTime):
       continue
    requestSize = int(fields[7])
    if (isClamp and requestSize > nTotalSize):
        continue
    elif (not isClamp):
        requestSize = min(int((requestSize / maxRequestSize) * nTotalSize) + 1, nTotalSize)
    requestTime = int(fields[8])
    outFile.write(str(requestSize) + ' ' + str(requestTime) + '\n')
    sumOut = sumOut + 1
    if (sumOut >= nOutput):
        break

file.close()
