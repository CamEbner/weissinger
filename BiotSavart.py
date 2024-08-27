#!/usr/bin/env python
# BiotSavart.py
# Written by Cameron Ebner
#
# This python module will implement the Biot-Savart law for use in MAE 561,
#   Wing Theory.  The equations used here are presented in lecture 6 for the
#   aforementioned course and this script was written and originally used
#   for the first homework assignment.

from __future__ import division
import matplotlib.pyplot as plt
import scipy.linalg as linalg
import scipy as sp

################################################################################

# Define unit vectors
unitI = sp.array([1,0,0])
unitJ = sp.array([0,1,0])
unitK = sp.array([0,0,1])

################################################################################

def UnitVec(Vector):
# Return a unit vector in the direction of Vector
    return Vector * (1 / linalg.norm(Vector))

################################################################################

def FiniteFilament(x1, x2, Point):
# Compute the velocity induced by a vortex filament from x1 to x2 at "Point".
#   x1, x2, and Point are 3-dimensional vectors giving the positions of the 
#   filament beginning, end, and the point at which induced velocity is desired,
#   respectively.
#
# I will do my best to use the same syntax as the lecture notes in this 
#   subroutine (i.e. vectors a and b, etc.)

    a = x1 - Point
    b = x2 - Point

    # qMagnitude will be split into 3 segments:
    Seg1 = sp.cross(a,b) / sp.dot( sp.cross(a,b), sp.cross(a,b) )
    Seg2 = linalg.norm(a) + linalg.norm(b)
    Seg3 = 1 - (sp.dot(a,b) / ( linalg.norm(a) * linalg.norm(b) ) )
    

    # Calculate the direction of induced velocity
    InducedDir = UnitVec( sp.cross( x2 - x1, -1 * a ) )
    qMagnitude = abs(Seg1 * Seg2 * Seg3)

    return qMagnitude * InducedDir

################################################################################

def SemiInfiniteFilament(Start, Point, Direction = -1*unitI, \
        FlipInduced = False):
# Compute the velocity induced by a semi-infinite vortex filament at "Point".
#   Start, Direction, and Point are 3-dimensional vectors.  Start gives the 
#   position of the beginning of the vortex.  Direction is a vector telling 
#   the direction in which the vortex goes to infinity.  Point is the point 
#   at which induced velocity is desired.

    aVec = Point - Start

    Alpha_Rad = sp.arccos(sp.dot( aVec, Direction) / \
            (linalg.norm( aVec ) * linalg.norm(Direction)))

    Distance = linalg.norm( aVec ) * sp.sin(Alpha_Rad)
    InducedDir = UnitVec( sp.cross( Direction, aVec ))

    # Allow the induced velocity to be flipped.  This will be applicable for 
    #   vorticies which are left-handed about their direction vector.  Vortices
    #   originating at "S" points are left-handed, but "E" vortices are
    #   right-handed.

    qMagnitude = (sp.cos(Alpha_Rad) + 1) / Distance

    if FlipInduced:
        InducedDir *= -1.0

    return qMagnitude * InducedDir

################################################################################

def TargetPoints(Start, End, NPoints):
# Return a list of vectors at which induced downwash is desired.  These points
#   will be equally spaced between the start and end points.

    OutPoints = []
    Line = End - Start

    for Percent in sp.linspace(0, 1, NPoints):
        OutPoints.append(Start + Line * Percent)

    return sp.array(OutPoints)

################################################################################

def PlotInduced(Points, Velocity,XAxis=0, VAxis=2, File=''):
# Plot the desired component of induced velocity at a point (Default is z).
#   The first two inputs are the points and the induced velocities at those 
#   points.  The third is which component to plot (options range from 0 to 2).
#   The fourth and final input is which format in which to plot, either 'plot'
#   for a live plot or 'pdf' to plot to a file.

    Fig = plt.figure(1)
    Fig.clf()
    Ax = Fig.add_subplot(1,1,1)
    Ax.plot(sp.array(Points)[:,XAxis], sp.array(Velocity)[:,VAxis], 'ko-')
    
    AxisDict = dict(zip( [0, 1, 2], ['X', 'Y', 'Z']))
    plt.xlabel('%1s (m)'  % AxisDict[XAxis])
    plt.ylabel('Induced Velocity in the %1s-Direction (m/s)' % AxisDict[VAxis])
    plt.title('Effects of %1s-Position on Induced Velocity in the %1s-Direction' % (AxisDict[XAxis], AxisDict[VAxis]))
    plt.grid('on')

    if len(File) == 0:
        plt.show()
    else:
        plt.savefig(File, format = 'pdf')

################################################################################
#
#### This will be the main segment of problem 2 for homework 1:
#
#### Part A:
#FilamentStart = sp.array([-1, 0, 0])
#FilamentEnd = sp.array([1, 0, 0])
#
## Hand calculation verification:
#Targets = [sp.array([20, 5, 0]), sp.array([0, 5, 0])]
#for Target in Targets:
#    print 'Target point: ', Target
#    print 'V induced (m/s): ', FiniteFilament(FilamentStart, FilamentEnd,Target)
#    print
#
## Get induced velocity at points between [-20, 5, 0] and [20, 5, 0]
#Targets = TargetPoints(sp.array([-20, 5, 0]), sp.array([20, 5, 0]), 81)
#VInduced_MpS = [FiniteFilament(FilamentStart, FilamentEnd, Point) for \
#    Point in Targets]
#
## Plot Z-component of induced velocity:
#PlotInduced(Targets, VInduced_MpS, File = 'VInduced_PartA.pdf')
#
#### Part B:
## Get induced velocity at points between [0, -20, 0] and [0, 20, 0]
#Targets = TargetPoints(sp.array([0, -20, 0]), sp.array([0, 20, 0]), 80)
#VInduced_MpS = [FiniteFilament(FilamentStart, FilamentEnd, Point) for \
#    Point in Targets]
#
## Plot Z-component of induced velocity:
##PlotInduced(Targets, VInduced_MpS, XAxis = 1)
#PlotInduced(Targets, VInduced_MpS, XAxis = 1, File = 'VInduced_PartB.pdf')
