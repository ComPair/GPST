#!/usr/bin/env python

import math
import numpy
import sys

outfilename = sys.argv[1]

background_data = [
    (0.0398, 0.0631, 0.00000001),
    (0.0631, 0.1000, 0.7200),
    (0.1000, 0.1585, 1.7000),
    (0.1585, 0.2512, 0.8400),
    (0.2512, 0.3981, 0.5500),
    (0.3981, 0.6310, 0.3700),
    (0.6310, 1.0000, 0.2800),
    (1.0000, 1.5849, 0.2000),
    (1.5849, 2.5119, 0.1000),
    (2.5119, 3.9811, 0.0800),
    (3.9811, 6.3096, 0.0600),
    (6.3096, 10.0000, 0.0600),
    (10.0000, 15.8489, 0.0400)
    ]

delta_e = 1    # need rate in 1keV steps

outfile = open(outfilename, 'w')

current_bin = 0
slope = None
scale = None
for energy in numpy.arange(100, 10000.1, delta_e):
    # find bin
    last_bin = current_bin
    while current_bin < len(background_data) and background_data[current_bin][1] < 1e-3*energy:
        current_bin += 1

    # recalculate slope and scale of power-law fit if needed
    if last_bin != current_bin or slope is None:
        slope = 0.
        cslope = 0.
        x = 0.5 * (background_data[current_bin][0] + background_data[current_bin][1])
        y = background_data[current_bin][2]
        if current_bin > 0:
            x2 = 0.5 * (background_data[current_bin-1][0] + background_data[current_bin-1][1])
            y2 = background_data[current_bin-1][2]
            slope += math.log(y/y2) / math.log(x/x2)
            cslope +=1
        if current_bin < len(background_data) - 1:
            x2 = 0.5 * (background_data[current_bin+1][0] + background_data[current_bin+1][1])
            y2 = background_data[current_bin+1][2]
            slope += math.log(y2/y) / math.log(x2/x)
            cslope +=1
        slope /= cslope
        scale = y * (1+slope)/(background_data[current_bin][1]**(1+slope) - background_data[current_bin][0]**(1+slope))
        #print 'lx : %g    ly : %g    Slope : %g     Scale : %g' % (x, y, slope, scale)

    # calculate number of events in [E, E + Delta E]
    rate = 0.
    e1 = energy*1e-3
    e2 = (energy + delta_e)*1e-3
    if slope == -1:
        rate = scale * math.log(e2/e1)
    else:
        rate = scale/(slope+1) * (e2**(slope+1) - e1**(slope+1))

    #print '%.2f %.2f %g' % (energy, energy + delta_e, rate)
    outfile.write('%g ' % rate)
