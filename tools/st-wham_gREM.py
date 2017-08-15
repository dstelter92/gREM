#!/usr/bin/env python2.7

import os,sys
from numpy import *

########### ST-WHAM for python...##########
# Originally written in f90 by Jaegil Kim #
# Translated to Python by David Stelter   #
###########################################
#
# CITE: http://dx.doi.org/10.1063/1.3626150
# Kim, J., Keyes, T., & Straub, J. E. (2011)
#
# WARNING: This version is only for use with gREM!
# http://dx.doi.org/10.1063/1.3432176
#
# Requires list enthalpies for each replica in separate files.
# Script will automatically read into histograms, run through
# ST-WHAM machinery and calculate Ts(H) and Entropy based on
# the gREM (non-boltzmann) sampling weight.
#

## UNITS!!! IMPORTANT!!!
kb = 0.0019872041       #kcal/mol*K
#kb = 0.000086173324    #ev/K
#kb = 0.0083144621      #kj/mol*K
#kb = 1.0               #reduced/LJ


def EffTemp(lambd, H):
# Evaluates the gREM effective temperature
    w = lambd + eta*(H - H0)
    return w

def Falpha(i, j):
# Linear entropy interpolation based on Ts(H)
    Falpha = 0
    for indx in range(i+1,j):
        if (TH[indx] == TH[indx-1]):
            Falpha = Falpha + binsize/TH[indx]
        else:
            Falpha = Falpha + binsize/(TH[indx] - TH[indx-1])*log(TH[indx]/TH[indx-1])
    return Falpha


## Input parameters from inp.stwham
print "ST-WHAM for gREM\n"
ifile = open('inp.stwham', 'r')
idata = ifile.readlines()
num_lines = len(idata)

if (num_lines > 8 or num_lines < 8):
    print "Err: Invalid input: Incorrect parameters, remove any extra newlines. Should be:\n"
    print "gamma\nH0\nbinsize\nEmin\nEmax\ndata_num\ncheckLimit\nLambdas\n"
    sys.exit()

## Cast inputs...
eta = double(idata[0])
H0 = double(idata[1])
binsize = double(idata[2])
Emin = double(idata[3])
Emax = double(idata[4])
data_num = int(idata[5])
checklimit = double(idata[6]) # Cutoff for contribution from neighbor replicas
lambdas = array(map(double, idata[7].split())) # List of lambdas


## Initialize outputs...
#hout = open('histogram_stwham.dat', 'w')
tout = open('Ts_stwham.dat', 'w')
fracout = open('fract_stwham.dat', 'w')


## Calculate some constants...
nbin = int((Emax-Emin)/binsize)
nReplica = len(lambdas)
shape = (nbin, nReplica)
fullshape = (nbin, nReplica+1)


## Final checks...
if (nbin < 0):
    print "Err: Emin must be smaller than Emax.\n"
    sys.exit()
if (nReplica < 1):
    print "Err: Must supply list of lambdas.\n"
    sys.exit()


## Prepare raw enthalpy files into histograms...
print "Reading in raw data, writting to histograms...\n"
hist = zeros(shape)
edges = zeros(nbin)
print "Lambdas: "
for l in range(nReplica):
    sys.stdout.write("%s " % lambdas[l])
    sys.stdout.flush()
    data = loadtxt("./replica-%d.dat" % (l))
    #data = loadtxt("../ent_%d-%d.dat" % (lambdas[l], data_num))
    # Calculate histogram
    hist[:,l], edges = histogramdd(ravel(data), bins=nbin, range=[(Emin, Emax)])

    #for i in range(nbin):
        #hout.write("%f %f\n" % (edges[0][i], hist[i][l]))
    #hout.write("%f 0.000000\n" % Emax)
    #hout.write("\n")
print "\n"


## Assign some arrays...
TH = zeros(nbin) # Statistical Temperature
Ent = zeros(nbin) # Entropy
PDF2D_GREM = zeros(fullshape) # Enthalpy distribution for all replicas and total data
NumData = zeros(nReplica+1) # Count of data points
betaH = zeros(nbin)
betaW = zeros(nbin)
hfrac = zeros(shape) # Histogram fraction, useful for debugging

print "Run the ST-WHAM machinery...\n"

## Collect data, and normalize...
for l in range(1,nReplica+1):
    count = 0
    for i in range(nbin):
        PDF2D_GREM[i][l] = hist[i][l-1] # PDF of each replica
        count = count + PDF2D_GREM[i][l]
        PDF2D_GREM[i][0] = PDF2D_GREM[i][0] + PDF2D_GREM[i][l] # PDF of total data set

    PDF2D_GREM[:,l] = PDF2D_GREM[:,l] / count
    NumData[l] = count # Number of data in each replica
    NumData[0] = NumData[0] + count # Total number of data

PDF2D_GREM[:,0] = PDF2D_GREM[:,0] / NumData[0] # Normalized PDF


## Throw out edge data under checklimit...
bstart = None
bstop = None
for i in range(nbin):
    if (PDF2D_GREM[i][0] > checklimit):
        bstart = i + 3
        break

if (bstart == None):
    print "Err: Energy range not large enough, decrease Emin.\n"
    sys.exit()

for i in range(nbin-1,bstart,-1):
    if (PDF2D_GREM[i][0] > checklimit):
        bstop = i - 3
        break

if (bstop == None):
    print "Err: Energy range not large enough, increase Emax.\n"
    sys.exit()


## Calculate Ts(H)
for i in range(bstart,bstop):
    if (PDF2D_GREM[i+1][0] > checklimit and PDF2D_GREM[i-1][0] > checklimit):
        betaH[i] = (log(PDF2D_GREM[i+1][0] / PDF2D_GREM[i-1][0])) / (2*binsize / kb) #NOTE: UNITS important here!
    else:
        betaH[i] = 0

    for l in range(1,nReplica+1):
        # Calc B^eff_alpha
        if (PDF2D_GREM[i][0] > 0): # ensure positive, no empty bins
            w = nan_to_num(1 / EffTemp(lambdas[l-1], Emin+(i*binsize)))
            betaW[i] = betaW[i] + ((NumData[l] * PDF2D_GREM[i][l]) / (NumData[0] * PDF2D_GREM[i][0]) * w)
    TH[i] = 1 / (betaH[i] + betaW[i])
    #tout.write("%f %f %f %f\n" % (Emin+(i*binsize), TH[i], betaH[i], betaW[i]))


## Calculate histogram fraction...
for l in range(1,nReplica+1):
    for i in range(bstart,bstop):
        hfrac[i][l-1] = hist[i][l-1] / (PDF2D_GREM[i][0] * NumData[0])
        fracout.write("%f %f\n" % (Emin+(i*binsize), hfrac[i][l-1]))
    fracout.write("\n")


## Calculate entropy...
for i in range(bstart,bstop):
    Ent[i] = Falpha(bstart, i)
    tout.write("%f %f %f\n" % (Emin+(i*binsize), TH[i], Ent[i]))





sys.exit()
