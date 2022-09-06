#!/usr/bin/env python

from math import log10
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
from numpy import *
from matplotlib.pyplot import *

rc("font", family="serif", size=18.)
rc("savefig", dpi=200)
rc("lines", linewidth=2.5, markersize=10, markeredgewidth=2.5)

##############################################################
#
#  BS3 plots; naccept = 2696226
#
##############################################################

# load the data

data_BS3 = loadtxt('shallow_water_exner_BS3_err.txt')

times_err = array([data_BS3[i][0] for i in range(len(data_BS3))])
effective_cfl_numbers = array([data_BS3[i][1] for i in range(len(data_BS3))])

# create the plot with time step as the x-axis

ax1 = gca()
time_line = ax1.plot(times_err[1:len(data_BS3)]-times_err[0:len(data_BS3)-1], label="Time step size", color="#E69F00", linestyle="-")
semilogx()

#ax1.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
xlabel("Time step")
ylabel("Time step size")
xlim(1, 10000000)
ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])

ax2 = ax1.twinx()
plot([], [])
cfl_line = ax2.plot(effective_cfl_numbers, label="Effective CFL number", color="#56B4E9", linestyle="--")
semilogx()

#ax2.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
ylabel("Effective CFL number")
ax2.set_yticks([1,2,3,4,5,6,7,8,9,10,11])

# collect the legends together

lns = time_line + cfl_line
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)

tight_layout(1)
savefig('swe_BS3_step.pdf')

clf()

# create the plot with time as the x-axis

ax1 = gca()
time_line = ax1.plot(times_err[1:len(data_BS3)], times_err[1:len(data_BS3)]-times_err[0:len(data_BS3)-1], label="Time step size", color="#E69F00", linestyle="-")
semilogx()

# ax1.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
xlabel("Time")
ylabel("Time step size")
xlim(0, 360000)
ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])

ax2 = ax1.twinx()
plot([], [])
cfl_line = ax2.plot(times_err, effective_cfl_numbers, label="Effective CFL number", color="#56B4E9", linestyle="--")
semilogx()

# ax2.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
ylabel("Effective CFL number")
ax2.set_yticks([1,2,3,4,5,6,7,8,9,10,11])

# collect the legends together

lns = time_line + cfl_line
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)

tight_layout(1)
savefig('swe_BS3_time.pdf')

##############################################################
#
#  RDPK3SpFSAL35 plots; naccept = 1374701
#
##############################################################

clf()

data_RDPK3SpFSAL35 = loadtxt('shallow_water_exner_RDPK3SpFSAL35_err.txt')

times_err = array([data_RDPK3SpFSAL35[i][0] for i in range(len(data_RDPK3SpFSAL35))])
effective_cfl_numbers = array([data_RDPK3SpFSAL35[i][1] for i in range(len(data_RDPK3SpFSAL35))])

# create the plot with time step as the x-axis

ax1 = gca()
time_line = ax1.plot(times_err[1:len(data_RDPK3SpFSAL35)]-times_err[0:len(data_RDPK3SpFSAL35)-1], label="Time step size", color="#E69F00", linestyle="-")
semilogx()

#ax1.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
xlabel("Time step")
ylabel("Time step size")
xlim(1, 10000000)
ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])

ax2 = ax1.twinx()
plot([], [])
cfl_line = ax2.plot(effective_cfl_numbers, label="Effective CFL number", color="#56B4E9", linestyle="--")
semilogx()

#ax2.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
ylabel("Effective CFL number")
ax2.set_yticks([1,2,3,4,5,6,7,8,9,10,11])

# collect the legends together

lns = time_line + cfl_line
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)

tight_layout(1)
savefig('swe_RDPK3SpFSAL35_step.pdf')

clf()

# create the plot with time as the x-axis

ax1 = gca()
time_line = ax1.plot(times_err[1:len(data_RDPK3SpFSAL35)], times_err[1:len(data_RDPK3SpFSAL35)]-times_err[0:len(data_RDPK3SpFSAL35)-1], label="Time step size", color="#E69F00", linestyle="-")
semilogx()

# ax1.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
xlabel("Time")
ylabel("Time step size")
xlim(0, 360000)
ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])

ax2 = ax1.twinx()
plot([], [])
cfl_line = ax2.plot(times_err, effective_cfl_numbers, label="Effective CFL number", color="#56B4E9", linestyle="--")
semilogx()

# ax2.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
ylabel("Effective CFL number")
ax2.set_yticks([1,2,3,4,5,6,7,8,9,10,11])

# collect the legends together

lns = time_line + cfl_line
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)

tight_layout(1)
savefig('swe_RDPK3SpFSAL35_time.pdf')

##############################################################
#
#  RDPK3SpFSAL49 plots; naccept = 715871
#
##############################################################

clf()

data_RDPK3SpFSAL49 = loadtxt('shallow_water_exner_RDPK3SpFSAL49_err.txt')

times_err = array([data_RDPK3SpFSAL49[i][0] for i in range(len(data_RDPK3SpFSAL49))])
effective_cfl_numbers = array([data_RDPK3SpFSAL49[i][1] for i in range(len(data_RDPK3SpFSAL49))])

# create the plot with time step as the x-axis

ax1 = gca()
time_line = ax1.plot(times_err[1:len(data_RDPK3SpFSAL49)]-times_err[0:len(data_RDPK3SpFSAL49)-1], label="Time step size", color="#E69F00", linestyle="-")
semilogx()

#ax1.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
xlabel("Time step")
ylabel("Time step size")
ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])

ax2 = ax1.twinx()
plot([], [])
cfl_line = ax2.plot(effective_cfl_numbers, label="Effective CFL number", color="#56B4E9", linestyle="--")
semilogx()
xlim(1, 10000000)

#ax2.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
ylabel("Effective CFL number")
ax2.set_yticks([1,2,3,4,5,6,7,8,9,10,11])

# collect the legends together

lns = time_line + cfl_line
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="lower right", fontsize="medium", fancybox=True, framealpha=0.5)

tight_layout(1)
savefig('swe_RDPK3SpFSAL49_step.pdf')

clf()

# create the plot with time as the x-axis

ax1 = gca()
time_line = ax1.plot(times_err[1:len(data_RDPK3SpFSAL49)], times_err[1:len(data_RDPK3SpFSAL49)]-times_err[0:len(data_RDPK3SpFSAL49)-1], label="Time step size", color="#E69F00", linestyle="-")
semilogx()

# ax1.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
xlabel("Time")
ylabel("Time step size")
xlim(0, 360000)
ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])

ax2 = ax1.twinx()
plot([], [])
cfl_line = ax2.plot(times_err, effective_cfl_numbers, label="Effective CFL number", color="#56B4E9", linestyle="--")
semilogx()

# ax2.legend(loc="upper right", fontsize="medium", fancybox=True, framealpha=0.5)
ylabel("Effective CFL number")
ax2.set_yticks([1,2,3,4,5,6,7,8,9,10,11])

# collect the legends together

lns = time_line + cfl_line
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="lower right", fontsize="medium", fancybox=True, framealpha=0.5)

tight_layout(1)
savefig('swe_RDPK3SpFSAL49_time.pdf')
