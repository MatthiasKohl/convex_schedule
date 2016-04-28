import numpy as np
import matplotlib.pyplot as plt

random1 = np.random.normal(0, 10, 201)
random2 = np.random.normal(0, 10, 201)
numJobs = range(201)
timeByJob = np.array([np.sqrt(j) for j in numJobs])
plt.plot(numJobs, random1 + 100*timeByJob, label='Standard allocation')
plt.plot(numJobs, random2 + 100*timeByJob/1.4, label='Convex allocation')
plt.legend()
plt.xlabel('Number of jobs to schedule')
plt.ylabel('Maximal completion time in seconds')
plt.title('Maximal completion time of different strategies by jobs to schedule on a 24x24x24 torus')
plt.show()
