import json
import numpy as np
from shapes import process_candidates, metric_compactness, metric_diameter, metric_max_min
from operator import mul
from functools import reduce

json_file=open('blue_waters_shapes.json')
data=json.load(json_file)
json_file.close()

csvString = ''
for c in data['classes']:
    for shape in c['shapes']:
        factorStrings = shape.split('x')
        line = str(reduce(mul, map(int, factorStrings))) + ',' + ','.join(factorStrings) + '\n'
        csvString += line

dimensions = {'x': 24, 'y': 24, 'z': 24}
shapes = process_candidates(csvString, dimensions)

#print(csvString)
#print(shapes)

print('Number of shapes: ' + str(sum(1 for n, projs in shapes.items() for p in projs)))
for metric, name in [(metric_diameter, 'diameter'),
    (metric_compactness, 'compactness'), (metric_max_min, 'max-min')]:
    metrics = [metric(p, False, dimensions) for n, projs in shapes.items() for p in projs]
    metricSum = sum(metrics)
    metricAvg = np.mean(metrics)
    print('Average ' + name + ': ' + str(metricAvg) + ', sum: ' + str(metricSum))

