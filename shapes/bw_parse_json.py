import json
import numpy as np
from convex_shapes import (process_candidates, metric_compactness, metric_diameter,
                           metric_max_min)
from operator import mul
from functools import reduce

with open('bw_shapes.json') as json_file:
    data=json.load(json_file)

csv_string = ''
for c in data['classes']:
    for shape in c['shapes']:
        factor_strings = shape.split('x')
        line = str(reduce(mul, map(int, factor_strings))) + ',' + ','.join(factor_strings) + '\n'
        csv_string += line

dimensions = {'x': 24, 'y': 24, 'z': 24}
shapes = process_candidates(csv_string, dimensions)

#print(csv_string)
#print(shapes)

print('Number of shapes: ' + str(sum(1 for n, projs in shapes.items() for p in projs)))
for metric, name in [(metric_diameter, 'diameter'),
    (metric_compactness, 'compactness'), (metric_max_min, 'max-min')]:
    metrics = [metric(p, False, dimensions) for n, projs in shapes.items() for p in projs]
    metricSum = sum(metrics)
    metricAvg = np.mean(metrics)
    print('Average ' + name + ': ' + str(metricAvg) + ', sum: ' + str(metricSum))

