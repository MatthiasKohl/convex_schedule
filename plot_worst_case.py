#! /usr/bin/python3
# -*- encoding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
import numpy as np

eps_space = 0.03
eps_time = 0.05

opt_rectangles = [((0, 0), eps_space, 2-eps_time), ((eps_space, 0), 1/2 + eps_space, 1),
((1/2 + 2*eps_space, 0), 1/2 - 2*eps_space, 1), ((eps_space, 1), 1/2 + eps_space, 1),
((1/2 + 2*eps_space, 1), 1/2 - 2*eps_space, 1)]
opt_text = [(0, 1 + eps_time, 'j0'), (1/4 + eps_space, 1/2, 'j1'), (1/4 + eps_space, 3/2, 'j2'),
(3/4 + eps_space, 1/2, 'j3'), (3/4 + eps_space, 3/2, 'j4')]
actual_rectangles = [((0, 0), 1, 1), ((0, 1), 1, 1), ((0, 2), 1/2, 1),
((1/2, 2), 1/2, 1), ((0, 3), eps_space, 2-eps_time)]
actual_text = [(0, 4 + eps_time, 'j0'), (1/2, 1/2, 'j1'), (1/2, 3/2, 'j2'), (1/4, 5/2, 'j3'), (3/4, 5/2, 'j4')]

patterns = ['/', '\\', '.', 'o', '*']
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 5])
ax1 = plt.subplot(gs[0])
ax1.set_ylim(0, 2)
ax1.set_ylabel('Time unit', fontweight='bold')
ax1.set_title('Optimal non-convex schedule', fontweight='bold')
for i, rect in enumerate(opt_rectangles):
    ax1.add_patch(patches.Rectangle(rect[0], rect[1], rect[2], fill=False))
    ax1.text(opt_text[i][0], opt_text[i][1], opt_text[i][2], bbox={'facecolor':'red'})
ax1.yaxis.set_ticks(np.arange(0, 3, 1))

ax2 = plt.subplot(gs[1], sharex=ax1)
ax2.set_ylim(0, 5)
ax2.set_ylabel('TIme unit', fontweight='bold')
ax2.set_xlabel('Ratio of resources', fontweight='bold')
ax2.set_title('Schedule produced by proposed algorithm', fontweight='bold')
for i, rect in enumerate(actual_rectangles):
    ax2.add_patch(patches.Rectangle(rect[0], rect[1], rect[2], fill=False))
    ax2.text(actual_text[i][0], actual_text[i][1], actual_text[i][2], bbox={'facecolor':'red'})
ax2.yaxis.set_ticks(np.arange(0, 6, 1))

plt.tight_layout()
plt.show()
