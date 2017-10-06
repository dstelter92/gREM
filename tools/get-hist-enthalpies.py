#!/usr/bin/env python2.7

import os, sys
from numpy import *

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

# Logfile info, requires manual edits for now...
header = 139
footer = 29

# Output walkers too?
walker = 1
##################
nfiles = (stop - start) # Number of files

print "Gather enthalpies from log files, follow exchanges, save data for each replica instead of per walker"

# Save Exchanges
loop = 1.0 / (stop - start)
tmp_Walk = []
print "Read exchanges..."
for i in range(stop-start):
    wname = '../log/log.lammps' + '-' + str(start+i)
    Walk = genfromtxt(wname,dtype=int,skip_header=3,skip_footer=1)
    tmp_Walk.append(Walk)
    sys.stdout.write("  %d%% completed\r" % int(loop*(i+1)*100))
    sys.stdout.flush()
tot_Walk = array(tmp_Walk)
print "\n"

# Define useful stuff from Walk array
lwalk = len(tot_Walk[:,0])            # Number of files from exchanges
walkers = len(tot_Walk[0][0])-1       # Number of walkers/replicas
exchange = tot_Walk[0][1][0] - tot_Walk[0][0][0] # Exchange freq
steps = tot_Walk[0][-1][0] + exchange # Number of steps
nexchange = len(tot_Walk[0])          # Number of exchanges

INDX = reshape(tot_Walk,(nfiles*nexchange,NumReplica+1))

if (NumReplica != walkers):
  print "Err. Number of walkers in lammps output does not equal number of replicas specified"
  sys.exit()

# Save Walker data
loop = 1.0 / NumReplica
loop2 = (loop * 1.0) / (stop - start)
tmp_Hist = []
print "Input data for", NumReplica, "walkers..."
for i in range(NumReplica):
  walk_Hist = []
  for j in range(stop-start):
    # Read data
    hname = '../' + str(i) + '/log.lammps.' + str(i) + '-' + str(start+j)
    Hist = genfromtxt(hname,skip_header=header,skip_footer=footer,usecols=(0,2))
    walk_Hist.append(array(Hist))

    # Write % completion
    percent = (j+1)*loop2 + (i)*loop
    sys.stdout.write("  (%.1f%%) Walker %d: Data index, %d Energy: min=%.3f max=%.3f\r" % (percent*100,i,j+start,amin(Hist[:,1]),amax(Hist[:,1])))
    sys.stdout.flush()

    if (Hist[0][0] != 0):
      print "\n"
      print "Err. Header value is incorrect", Hist[0][0]
      sys.exit()
    if (Hist[-1][0] >= steps):
      print "\n"
      print "Err. Footer value is incorrect", Hist[-1][0], steps
      sys.exit()
  tmp_Hist.append(array(walk_Hist))
tot_Hist = array(tmp_Hist)
print "\n"

# Define some useful things from Hist arrays
nthermo = len(tot_Hist[0][0])             # Number of thermo outputs
thermo = steps / nthermo                  # Thermo freq
ratio = nthermo / nexchange               # Ratio of thermo / exchange data
numwalkers = len(tot_Hist)                # Number of walkers read in
lhist = len(tot_Hist[0])                  # Number of files read in
numdata = nthermo * numwalkers * lhist    # Total # of data

if ((nfiles != lwalk) or (lhist != nfiles) or (lhist != lwalk)):
  print "Err. Not all data read in correctly"
  sys.exit()

if (NumReplica != numwalkers):
  print "Err. Number of walker data read in is different than what specified"
  sys.exit()

# Save Replicas
print "Get replicas..."
f_replica_data = []
loop = 1.0 / nfiles
loop2 = (loop * 1.0) / (nexchange)
loop3 = (loop2 * 1.0) / (NumReplica)
# put all data here
replica_data = empty((NumReplica,nfiles*nthermo))
walker_data = empty((NumReplica,nfiles*nthermo))
# Start loop
for file in range(nfiles): # Loop over files
  warg = tot_Walk[file][:,1:] # walkers for each step of file
  wstep = tot_Walk[file][:,0] # walker steps for this file
  for frame in range(nexchange): # Loop over num exchanges in each file
    wlist = warg[frame]
    s = wstep[frame] + file*steps
    i = frame*ratio
    indx = frame*ratio + file*nthermo
    for walk in range(NumReplica): # Loop over each walker, assign to replica
      # Output progress
      percent = (frame)*loop2 + (file)*loop + (walk+1)*loop3
      sys.stdout.write("  %.1f%% completed\r" % (percent*100))
      sys.stdout.flush()
      this_walker = wlist[walk] # walker "walk" belongs to replica "this_walker"
      hstep = tot_Hist[walk][file][:,0][i:i+ratio] # enthalpy steps with same ex
      hvalues = tot_Hist[walk][file][:,1][i:i+ratio] # enthalpies with same ex
      replica_data[this_walker][indx:indx+ratio] = hvalues
      if (walker):
        walker_data[walk][indx:indx+ratio] = hvalues
print "\n"

print "Saving outputs..."
for rep in range(NumReplica):
    Rname = 'replica-' + str(rep) + '_' + str(start) + '-' + str(stop) + '.dat'
    Wname = 'walker-' + str(rep) + '_' + str(start) + '-' + str(stop) + '.dat'
    savetxt(Rname, replica_data[rep])
    savetxt(Wname, walker_data[rep])

print "Saving exchanges..."
savetxt('exchange_list.dat', INDX[:,1:], fmt='%d')

print "Done!"
sys.exit()
