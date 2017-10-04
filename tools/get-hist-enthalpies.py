#!/usr/bin/env python2.7

import os, sys
from numpy import *
import pandas as pd

##################
###  For gREM  ###
##################
#
# Take walker enthalpies/energies from multiple, identical
# length trajectories, concatenate, and follow exchanges
# to get replica enthalpies/energies. Start and stop are
# file handles for each restarted simulation.
#
# Usage:
# python get-hist-enthalpies #replicas start stop
#
##################
### PARAMETERS ###
##################
NumReplica = int(sys.argv[1])
start = int(sys.argv[2])
stop = int(sys.argv[3])

# Per-run info, requires manual edits for now...
steps = 50000
exchange = 1000
thermo = 50

# Logfile info, requires manual edits for now...
header = 111
footer = 29
##################

length = steps/thermo
length_walk = steps/exchange
ratio = exchange/thermo

print "Gather enthalpies from log files, follow exchanges, save data for each replica instead of per walker"

size = ((length)*(stop-start),3*NumReplica+1)
size_walk = ((length)*(stop-start),NumReplica+1)
full_data = empty(size)
full_data[:] = NAN
# Format: step, indx1, indx2,... indxN, walk1, walk2,...walkN, rep1, rep2,...repN

# Save Exchanges
print "Input exchanges..."
for i in range(stop-start):
    wname = '../log/log.lammps' + '-' + str(start+i)
    Walk = genfromtxt(wname, skip_header=3, skip_footer=1)
    for j in range(NumReplica):
        for k in range(length):
            indx = k+(i*length)
            if (exchange == steps):
                full_data[:,j+1][indx] = int(Walk[j+1])
            else:
                full_data[:,j+1][indx] = int(Walk[:,j+1][k/ratio])

# Save Walkers
print "Input ", NumReplica, " walkers..."
for i in range(NumReplica):
    #print "\n===--- Walker", i, "---==="
    for j in range(stop-start):
        hname = '../' + str(i) + '/log.lammps.' + str(i) + '-' + str(start+j)
        Hist = genfromtxt(hname, skip_header=header, skip_footer=footer, usecols = (0,2))

        # Error checking
        if (Hist[0][0] != 0):
          print "Err. Header value is incorrect."
          sys.exit()

        if (Hist[-1][0] != steps-thermo):
          print "Err. Footer value is incorrect."
          sys.exit()

        #sys.stdout.write("Walker %d, Data index: %d\n" % (i,j+start))
        #sys.stdout.write("  Step: start=%d end=%d\n" % (Hist[0][0],Hist[-1][0]))
        sys.stdout.write("  Walker %d: Data index, %d Energy: min=%f max=%f\r" % (i,j+start,amin(Hist[:,1]),amax(Hist[:,1])))
        sys.stdout.flush()

        for k in range(length):
        # Load enthalpy data in global array
            indx = k+(j*length)
            full_data[:,0][indx] = Hist[:,0][k] + (j*steps)
            full_data[:,i+NumReplica+1][indx] = Hist[:,1][k]
print "\n"

# Save Replicas
print "Get replicas..."
for rep in range(NumReplica):
    for indx in range(len(full_data[:,0])):
        for j in range(NumReplica):
            if (full_data[:,j+1][indx] == rep):
                full_data[:,rep+(2*NumReplica)+1][indx] = full_data[:,j+NumReplica+1][indx]

#savetxt('full_replica_data.dat', full_data)
walkers = empty(size_walk, dtype=int)
print "\nSaving outputs..."
for rep in range(NumReplica):
    Rname = 'replica-' + str(rep) + '_' + str(start) + '-' + str(stop) + '.dat'
    Wname = 'walker-' + str(rep) + '_' + str(start) + '-' + str(stop) + '.dat'
    savetxt(Rname, full_data[:, (rep+(2*NumReplica)+1)])
    savetxt(Wname, full_data[:, (rep+NumReplica+1)])
    walkers[:,rep+1] = full_data[:,rep+1].astype(int)

walkers[:,0] = full_data[:,0].astype(int)
print "Saving exchanges..."
savetxt('exchange_list.dat', walkers, fmt='%d')

print "Done!"
sys.exit()
