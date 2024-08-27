#!/usr/bin/env python
# Run.py
#
# Run Weissinger's method code for MAE 561, programming assignment 2 and 
#   homework 3.

from __future__ import division
import scipy as sp
import sys
import Weissinger as W

# Collect patch files
#Alpha_Deg = float(sys.argv[1])
PatchFiles = sys.argv[2:]
Patches = [ W.ReadWingInput(File) for File in PatchFiles ]

# Build the 2 influence coefficient matrices
ControlMat, TrefftzMat, AllLattices = W.InfluenceCoeff(Patches)

# Calculate reference area based on area of each patch
Sref = 0
for Patch in Patches:
    Sref += Patch.Area

# Calculate effects of each freestream
Freestreams = W.ReadFsInput( sys.argv[1] )
for FS in Freestreams:
    FS.RhsMatrix = W.RhsMatrix( AllLattices, FS.AlphaW_Deg )
    FS.CalcLift(AllLattices, ControlMat, TrefftzMat, Sref)
    FS.WriteOutput( AllLattices )

# Write Integrated output for all Freestreams:
W.WriteIntegratedOutput(sys.argv[1], Freestreams)

# Build the right hand side vector:
#Freestreams[0].RhsMatrix = W.RhsMatrix(AllLattices, Alpha_Deg)

#Freestreams[0].CalcLift(AllLattices, ControlMat, TrefftzMat, Sref_M2)

### Print a bunch of output
# Set print options:
sp.set_printoptions(precision = 4, suppress = True)
#print 'ControlMat:\n', ControlMat
#print '\n\nTrefftzMat:\n', TrefftzMat
#print '\n\nRhsVec:\n', Freestreams[0].RhsMatrix
#print '\nwing_CL  =', Freestreams[0].CL
#print 'wing_CDi =', Freestreams[0].CDi
