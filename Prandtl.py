# PrandtlMod.py
# 
# Module for MAE 561 Programming assignment 1 -- Fall 2014

from __future__ import division
from scipy.interpolate import interp1d
from scipy.linalg import solve
import re
import scipy as sp
import sys

################################################################################

class WingClass:
# Class to contain all input data for a wing.  Most of all this will prevent
#   the need to return a bunch of different variables from ReadInput
    def __init__(self, Span_M, Aspect, NPoints):
    # Define all wing and modeling properties
    #   The only required inputs are Span, Aspect Ratio, and number of points
        self.Span_M = float(Span_M)
        self.Aspect = float(Aspect)
        self.NPoints = int(NPoints)

        self.Thetas_Rad = sp.linspace(0, sp.pi, self.NPoints + 2)[1:-1]
        self.Y_M = (self.Span_M / -2) * sp.cos(self.Thetas_Rad)

        # Initialize empty variables
        self.Chords_M = None
        self.a0_pRad = None
        self.AlphaZL_Deg = None
        self.Twist_Deg = None
        self.Alphas_Deg = None

        self.IC_Mat = None
        self.ACoeff = None

        # Distributed quantities
        self.GammaV_M = None
        self.Cl = None
        self.ClMax = None

    ################################################################

    def SetAlphas(self, AlphaStart, AlphaEnd, NAlphas):
    # Populate the list of all angles of attack
        self.Alphas_Deg = sp.linspace(AlphaStart, AlphaEnd, NAlphas)

    ################################################################

    def SetEllipticalChords(self, Aspect, cRoot_M):
    # Populate the list of all chords in a tapered wing as a function of theta

        # Calculate root and tip chords
        if Aspect > 0:
            cRoot_M = (4 * self.Span_M) / (self.Aspect * sp.pi)
    
        # Loop over all Theta-values and calculate local chord
        self.Y_M = (self.Span_M / -2) * sp.cos(self.Thetas_Rad)
        self.Chords_M = cRoot_M * sp.sqrt(1 - 4 * \
            (self.Y_M**2 / self.Span_M**2))

    ################################################################
    
    def SetByInterp(self, InArray):
    # To be used with non-elliptical wings.  The input array will be the block
    #   of break points contained in the bottom of the input file

        # Initialize and unpack columns
        YBreak_M = InArray[:,0]
        CBreak_M = InArray[:,1]
        a0Break_pRad = InArray[:,2]
        AlphaZLBreak_Deg = InArray[:,3]
        TwistBreak_Deg = InArray[:,4]

        # Interpolation functions
        ChordFunc   = interp1d(YBreak_M, CBreak_M, kind = 'linear')
        a0Func      = interp1d(YBreak_M, a0Break_pRad, kind = 'linear')
        AlphaZLFunc = interp1d(YBreak_M, AlphaZLBreak_Deg, kind = 'linear')
        TwistFunc   = interp1d(YBreak_M, TwistBreak_Deg, kind = 'linear')

        # Interpolate data
        self.Chords_M = ChordFunc(self.Y_M)
        self.a0_pRad = a0Func(self.Y_M)
        self.AlphaZL_Deg = AlphaZLFunc(self.Y_M)
        self.Twist_Deg = TwistFunc(self.Y_M)

################################################################################

def ReadInput(InFile):
# Read the provided input file
    InData = sp.loadtxt(InFile, comments = '#')

    # Read the first line to get global wing properties
    Span_M, Aspect, Shape, NPoints, _ = InData[0,:]
    NBreak = InData[1:-1, :].shape[0]
    Wing = WingClass(Span_M, Aspect, NPoints)

    if Shape == 0:
    # Elliptical wing
        YBreak_M, CBreak_M, a0Break_pRad, AlphaZLBreak_Deg, _ = InData[1,:]
        Wing.SetEllipticalChords(Aspect, CBreak_M)

        Wing.a0_pRad = a0Break_pRad * sp.ones(Wing.Y_M.shape)
        Wing.AlphaZL_Deg = AlphaZLBreak_Deg * sp.ones(Wing.Y_M.shape)
        Wing.Twist_Deg = sp.zeros(Wing.Y_M.shape)

    elif NBreak == 1:
    # Rectangular wing with no twist and only one entry
        YBreak_M, CBreak_M, a0Break_pRad, AlphaZLBreak_Deg, _ = InData[1,:]
        cRoot_M = Wing.Span_M / Wing.Aspect

        Wing.Chords_M = cRoot_M * sp.ones(Wing.Y_M.shape)
        Wing.a0_pRad = a0Break_pRad * sp.ones(Wing.Y_M.shape)
        Wing.AlphaZL_Deg = AlphaZLBreak_Deg * sp.ones(Wing.Y_M.shape)
        Wing.Twist_Deg = sp.zeros(Wing.Y_M.shape)

    else:
    # Rectangular/Tapered/Partially Tapered wing
        Wing.SetByInterp(InData[1:NBreak + 1,:])

    # Get the Alpha values over which to run
    Wing.SetAlphas(*InData[-1, :3])

    return Wing

################################################################################

def InfluenceCoeff(Wing):
# Build an influence coefficient matrix, given a list of theta values

    # Initialize matrix
    OutMat = sp.empty((Wing.NPoints, Wing.NPoints)); OutMat.fill(sp.nan)

    for i, (Theta, a0, Chord) in \
        enumerate(zip(Wing.Thetas_Rad, Wing.a0_pRad, Wing.Chords_M)):
    # Loop over all thetas -- this will serve as looping over rows
        for j in range(Wing.NPoints):
        # Loop over number of points -- looping over columns

            Term1 = (4 * Wing.Span_M) / (a0 * Chord) * sp.sin((j + 1) * Theta)
            Term2 = ((j + 1) * sp.sin((j + 1) * Theta)) / sp.sin(Theta)

            OutMat[i,j] = Term1 + Term2

    return OutMat

################################################################################

def RhsMatrix(Wing, Alpha_Deg):
# Build the right hand side matrix

   ColumnArray = sp.radians(Alpha_Deg + Wing.Twist_Deg - Wing.AlphaZL_Deg)
   return ColumnArray.reshape(len(ColumnArray), 1)

################################################################################

def LiftCoeff(AspectRatio, AVector):
# Return CL, the wing lift coefficient given the vector of A coefficients

    return AVector[0] * sp.pi * AspectRatio

################################################################################

def Delta(AVector):
# Return delta, which is related to e, the span efficiency factor

    # Catch cases with 0 lift:
    if AVector[0] == 0:
        return 0

    delta = 0
    for n, An in enumerate(AVector[1:]):
        n += 2      # Account for n starting at 2
        delta += (n * (An / AVector[0])**2)

    return delta

################################################################################

def IndDragCoeff(AspectRatio, CL, Delta):
# Return CDi, the induced drag coefficient for the wing

    return (CL**2 * (1 + Delta)) / (sp.pi * AspectRatio)

################################################################################

def CircOverVel(Thetas_Rad, Span_M, AVector):
# Return the distribution of local circulation (gamma) 
#   divided by freestream velocity

    Out = []
    for Theta_Rad in Thetas_Rad:
        GammaV = 0
        for n, An in enumerate(AVector):
            n += 1      # Account for n starting at 1
            GammaV += An * sp.sin(n * Theta_Rad)

        GammaV *= (2 * Span_M)
        Out.append(GammaV)

    return Out

################################################################################

def LocalLiftCoeff(GammaV_M, Chords_M):
# Calculate the local lift coefficient, C_l at each theta given the chord and
#   the distribution of gamma/v_inf

#    print 'Gamma', GammaV_M
    if len(GammaV_M.shape) == 1:
        GammaV_M = GammaV_M.reshape(len(GammaV_M), 1)

    OutArray = sp.array([])
    for GammaDist in GammaV_M.T:
        OutArray = ExtendArray(OutArray, (2 / Chords_M) * GammaDist)

    return OutArray

################################################################################

def ExtendArray(Cumulative, Individual):
# All data from 'Individual' is to be stored at the end of 'Cumulative'.  This
#   will be done by appending it if they have the same number of columns.
#   Otherwise, 'Individual' will overwrite 'Cumulative'.
    
    Cumulative = sp.array(Cumulative)
    Individual = sp.array(Individual)
#    print 'Cumulative: ', Cumulative
#    print 'Individual: ', Individual

    if len(Cumulative.shape) == 1 and len(Individual.shape) == 2:
        return sp.hstack([Cumulative.reshape(len(Cumulative), 1), Individual])

    try:
        if len(Cumulative.shape) == 1:
            Cumulative = Cumulative.reshape(len(Cumulative), 1)
        return sp.hstack([Cumulative, Individual.reshape(len(Individual), 1)])
    except:
        return Individual.reshape(len(Individual), 1)

################################################################################

def WriteIntegratedOutput(Wing, InFile):
# Calculate output and print it to a file

    # Open the file for writing
    FileRoot = re.search(r'(.+).in', InFile).group(1)
    OutFile = FileRoot + '.Int.out'
    File = open(OutFile, 'w')

    # Print the first part of the file, with wing properties
    File.write('# ' + OutFile + '\n')
    File.write('#%7s %8s %8s\n' % ('NPoints', 'Span(m)', 'AR'))
    File.write('%-8i %8.3f %8.3f\n' % (Wing.NPoints, Wing.Span_M, Wing.Aspect))
    File.write('#\n#\n')

    # Print the parameters which change with angle of attack
    File.write('#%11s %10s %16s %10s\n' % \
            ('Alpha(Deg)', 'CL', 'CDi', 'Delta'))
    for Alpha_Deg in Wing.Alphas_Deg:
        RhsMat = RhsMatrix(Wing, Alpha_Deg)
        Wing.ACoeff = solve(Wing.IC_Mat, RhsMat)
#        print 'A =\n', Wing.ACoeff

        CL = LiftCoeff(Wing.Aspect, Wing.ACoeff)
        delta = Delta(Wing.ACoeff)
        CDi = IndDragCoeff(Wing.Aspect, CL, delta)

        File.write('%-12.3f %10.4f %16.10f %10.4f\n' % \
                (Alpha_Deg, CL, CDi, delta))

    # Close the file
    File.close()

    return Wing

################################################################################

def WriteDistributionOutput(Wing, InFile, WriteFlag):
# Calculate output and print it to a file

    # Open the file for writing
    FileRoot = re.search(r'(.+).in', InFile).group(1)
    GammaFile = FileRoot + '.GammaV.out'
    ClFile = FileRoot + '.Cl.out'

    for Alpha_Deg in Wing.Alphas_Deg:
        # Calculate the Gamma/V and C_l distributions for each alpha
        Wing.GammaV_M = ExtendArray(Wing.GammaV_M, \
            CircOverVel(Wing.Thetas_Rad, Wing.Span_M, Wing.ACoeff))

    Wing.Cl = LocalLiftCoeff(Wing.GammaV_M, Wing.Chords_M)

    # Create Output Arrays:
    ClOut = ExtendArray(ExtendArray(Wing.Thetas_Rad, Wing.Chords_M), Wing.Cl)
    GammaVOut = ExtendArray(ExtendArray(Wing.Thetas_Rad, Wing.Chords_M), \
        Wing.GammaV_M)

    # Write output:
    #sp.savetxt(GammaFile, ExtendArray(Wing.Thetas_Rad, Wing.GammaV_M), \
    if WriteFlag == True:
        sp.savetxt(GammaFile, GammaVOut, \
            fmt = '%10.8f', header = 'Theta(Rad)    Chord(m)     Gamma/V(m)')
    #sp.savetxt(ClFile, ExtendArray(Wing.Thetas_Rad, Wing.Cl), fmt = '%10.8f', \
        sp.savetxt(ClFile, ClOut, fmt = '%10.8f', \
            header = 'Theta(Rad)    Chord(m)     Cl')

    return Wing

################################################################################
