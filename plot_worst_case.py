#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.patches as patches

eps_space = 0.03
eps_time = 0.05

optRectangles = [((0, 0), eps_space, 2-eps_time), ((eps_space, 0), 1/2 + eps_space, 1),
((1/2 + 2*eps_space, 0), 1/2 - 2*eps_space, 1), ((eps_space, 1), 1/2 + eps_space, 1),
((1/2 + 2*eps_space, 1), 1/2 - 2*eps_space, 1)]
optText = [(0, 1 + eps_time, 'j0'), (1/4 + eps_space, 1/2, 'j1'), (1/4 + eps_space, 3/2, 'j2'),
(3/4 + eps_space, 1/2, 'j3'), (3/4 + eps_space, 3/2, 'j4')]
actualRectangles = [((0, 0), 1, 1), ((0, 1), 1, 1), ((0, 2), 1/2, 1),
((1/2, 2), 1/2, 1), ((0, 3), eps_space, 2-eps_time)]
actualText = [(0, 4 + eps_time, 'j0'), (1/2, 1/2, 'j1'), (1/2, 3/2, 'j2'), (1/4, 5/2, 'j3'), (3/4, 5/2, 'j4')]

patterns = ['/', '\\', '.', 'o', '*']
fig, axarr = plt.subplots(2, sharex=True)
axarr[0].set_ylim(0, 2)
axarr[0].set_ylabel('Time', fontweight='bold')
axarr[0].set_title('Optimal non-convex schedule', fontweight='bold')
for i, rect in enumerate(optRectangles):
    axarr[0].add_patch(patches.Rectangle(rect[0], rect[1], rect[2], hatch=patterns[i], fill=False))
    axarr[0].text(optText[i][0], optText[i][1], optText[i][2], bbox={'facecolor':'red'})

axarr[1].set_ylim(0, 5)
axarr[1].set_ylabel('TIme', fontweight='bold')
axarr[1].set_xlabel('Processors', fontweight='bold')
axarr[1].set_title('Schedule produced by proposed algorithm', fontweight='bold')
for i, rect in enumerate(actualRectangles):
    axarr[1].add_patch(patches.Rectangle(rect[0], rect[1], rect[2], hatch=patterns[i], fill=False))
    axarr[1].text(actualText[i][0], actualText[i][1], actualText[i][2], bbox={'facecolor':'red'})

plt.show()
