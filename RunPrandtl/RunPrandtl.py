#!/home/cebner/anaconda/bin/python
# Prob1.py
# Function for HW2, Problem 1 for MAE 561

from __future__ import division
import PrandtlMod as PM
import Prob1Mod as P1
import scipy as sp
import sys

Wing = PM.ReadInput(sys.argv[1])

#print '\n\nY_M: ', Wing.Y_M
#print 'Thetas: ', Wing.Thetas_Rad
#print 'Chords_M: ', Wing.Chords_M
#print 'a0_pRad: ', Wing.a0_pRad
#print 'Twist_Deg: ', Wing.Twist_Deg
#print 'AlphaZL_Deg: ', Wing.AlphaZL_Deg
#print 'Alphas_Deg: ', Wing.Alphas_Deg
#print; print

Wing.IC_Mat = PM.InfluenceCoeff(Wing)

#print 'IC_Mat Shape: ', IC_Mat.shape

#Wing = PM.WriteIntegratedOutput(Wing, IC_Mat, sys.argv[1])
AlphaStall = P1.FindStall(Wing)
print 'Alpha Stall = ', AlphaStall, 'degrees'

for Alpha_Deg in Wing.Alphas_Deg:
    Wing = P1.DistributedQuantities(Wing, sys.argv[1], Alpha_Deg, True)
#PM.WriteDistributionOutput(Wing, sys.argv[1], True)
#Wing.Alphas_Deg = sp.array([AlphaStall])
Wing = PM.WriteIntegratedOutput(Wing, sys.argv[1])
