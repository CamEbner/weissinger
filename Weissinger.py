# Weissinger.py
# 
# Module for MAE 561 Programming assignment 2 -- Fall 2014
# The intention here is to write a code to implement Weissinger's method

from __future__ import division
import BiotSavart as BS
import re
from scipy.interpolate import interp1d
import scipy.linalg as linalg
import scipy as sp

################################################################################

# Define unit vectors
unitI = sp.array([1,0,0])
unitJ = sp.array([0,1,0])
unitK = sp.array([0,0,1])

################################################################################

class PatchClass:
# Class to contain all input data for a wing.  Most of all this will prevent
#   the need to return a bunch of different variables from ReadInput
    def __init__(self, NLattice, Line1, Line2):
    # Define all wing and modeling properties
    #   The only required inputs are Span, Aspect Ratio, and number of points

        self.NLattice = NLattice
        
        # Setup boundary points
        self.StartPt = Line1[:3]
        self.EndPt   = Line2[:3]

        # Calculate Span and area
        #   Span and area are calculated as they would apply to aspect ratio
        #   -- Vertical winglets would have no span and therefore no area.
        self.Span = abs( self.EndPt[1] - self.StartPt[1] )
        self.Area = self.Span * 0.5 * (Line1[3] + Line2[3])

        # Convert alphaZL and twist angle to radians
        Line1[4] = sp.radians(Line1[4]); Line1[5] = sp.radians(Line1[5])
        Line2[4] = sp.radians(Line2[4]); Line2[5] = sp.radians(Line2[5])

        # Initialize lattices
        self.Lattices = []

        for i in range(NLattice):
            self.Lattices.append(LatticeClass(self, i, Line1, Line2))

        # Distributed quantities
        self.GammaV_M = None
        self.Cl = None
        self.ClMax = None

################################################################################

def LinearInterp(StartPt, EndPt, StartVal, EndVal, InterpPt):
# Linearly interpolate from the set of points (XFrom, YFrom) to XTo

    return (linalg.norm(InterpPt - StartPt) / linalg.norm(EndPt - StartPt)) * \
        (EndVal - StartVal) + StartVal

################################################################################

class LatticeClass:
# Class to contain all necessary data to describe one lattice.  Each patch will
#   contain a list of lattices, each of which will have several parameters.

    def __init__(self, Patch, iLat, StartLine, EndLine):
    # Define all lattice properties

        # Calculate lattice points -- slide 16-14
        self.S1 = Patch.StartPt + (1 / Patch.NLattice) * \
            (Patch.EndPt - Patch.StartPt) * (iLat)
        self.E1 = Patch.StartPt + (1 / Patch.NLattice) * \
            (Patch.EndPt - Patch.StartPt) * (iLat + 1)
        self.MidPt = 0.5 * (self.S1 + self.E1)
        self.SPt = sp.array([self.MidPt[0], self.S1[1], self.S1[2]])
        self.EPt = sp.array([self.MidPt[0], self.E1[1], self.E1[2]])

        # Unit normal vector -- slide 17-8
        self.UnitNormal = BS.UnitVec( sp.cross((self.EPt - self.SPt), unitI) )
        # Dihedral angle -- slide 17-12
        self.DihedralLocal_Rad = sp.arccos(-1.0 * self.UnitNormal[2])

        # Interpolate parameters that change along the patch
        self.Chord_M = LinearInterp(Patch.StartPt, Patch.EndPt, \
            StartLine[3], EndLine[3], self.MidPt)
        self.AlphaZL_Rad = LinearInterp(Patch.StartPt, Patch.EndPt, \
            StartLine[4], EndLine[4], self.MidPt)
        self.Twist_Rad = LinearInterp(Patch.StartPt, Patch.EndPt, \
            StartLine[5], EndLine[5], self.MidPt)

        # 3c/4 point and Trefftz plane point
        self.ControlPt = sp.array([(self.MidPt[0] - self.Chord_M / 2), \
            self.MidPt[1],  self.MidPt[2] ])
        self.TrefftzPt = sp.array([ 0.0, self.MidPt[1],  self.MidPt[2] ])

        # Initialize currently unused values
        #self.AlphaLocal_Deg = None
        self.GammaV_M = None
        self.Cl = None

################################################################################

class FreestreamClass:
# Class to store all integrated data for one freestream condition.  This will
#   include patches, angles of attack, CL and CDi for each alpha, etc.

    def __init__(self, FileRoot, Alpha_Deg):
        self.AlphaW_Deg = Alpha_Deg

        self.Name = FileRoot + '.a' + str(Alpha_Deg)

        # Initialize currently unused scalars
        self.CL = sp.nan
        self.CDi = sp.nan

        # Initialize currently unused matrices
        self.RhsMatrix = None
        self.GammaV_M = None
        self.Cl = None

    ######################################################################

    def CalcLift(self, AllLattices, ControlMat, TrefftzMat, Sref_M2):
    # Once the Right-hand-side matrix is populated, calculate the circulation
    #   at each lattice and from there calculate lift and drag

        #print 'ControlMat: ', ControlMat
        self.GammaV_M = linalg.solve(ControlMat, self.RhsMatrix)
        #print 'Gamma:\n', self.GammaV_M

        ### Calculate Cl vector (slide 19-15)
        cVec = sp.array([sp.copy(Lattice.Chord_M) for Lattice in AllLattices])
        cVec = cVec.reshape( (len(cVec), 1) )
        self.Cl = (2.0 / cVec) * self.GammaV_M

        ### Calculate CL (slide 19-17)
        LatSpans_M = sp.array( \
            [sp.sqrt( (Lattice.EPt[1] - Lattice.SPt[1])**2 + (Lattice.EPt[2] \
            - Lattice.SPt[2])**2 ) for Lattice in AllLattices] )
        LatSpans_M = LatSpans_M.reshape( (len(LatSpans_M), 1) )
        Dihedrals_Rad = sp.array( [Lattice.DihedralLocal_Rad for Lattice \
            in AllLattices] )
        Dihedrals_Rad = Dihedrals_Rad.reshape( (len(Dihedrals_Rad), 1) )
        self.CL = sp.sum( self.GammaV_M * sp.cos(Dihedrals_Rad) * LatSpans_M )
        self.CL *= ( 2 / Sref_M2 )

        ### Calculate CDi (slide 19-20)
        Vnt = sp.asarray( sp.asmatrix( TrefftzMat ) * \
            sp.asmatrix( self.GammaV_M.reshape(len(self.GammaV_M), 1) ) )
        self.CDi = sp.sum( self.GammaV_M * sp.asarray(Vnt) * LatSpans_M ) / \
            Sref_M2

        # Store Gamma and Cl data with their respective lattices
        for Lattice, Gamma, Cl in zip(AllLattices, self.GammaV_M, self.Cl):
            Lattice.GammaV_M = Gamma
            Lattice.Cl = Cl

    ######################################################################

    def WriteOutput(self, AllLattices):
    # Once the above two methods have been called on this freestream condition,
    #   create output files that have all of the necessary data

    ### Data to print:
    # X, Y, Z, Chord, Gamma/V, Cl
    # Locations will be at midpoint of bound vortex, on wing quarter chord

        OutFile = self.Name + '.out'
        File = open(OutFile, 'w')

        # Print the first part of the file, with integrated properties
        File.write('# ' + OutFile + '\n')
        File.write('# %12s %8s %8s\n' % ('Alpha(Deg)', 'CL', 'CDi') )
        File.write('# %7.2f %13.4f %8.4f\n' % \
            (self.AlphaW_Deg, self.CL, self.CDi) )
        File.write('#\n#\n')

        # Print the distributed output
        File.write('# %-6s %8s %8s %10s %13s %8s\n' % ('X(m)', 'Y(m)', 'Z(m)', \
            'Chord(m)', 'Gamma/V(m)', 'Cl') )

        # Print the parameters which change with angle of attack
        for Lat in AllLattices:
            File.write('%8.3f %8.3f %8.3f %10f %13f %8f\n' % (Lat.MidPt[0], \
                Lat.MidPt[1], Lat.MidPt[2], Lat.Chord_M, Lat.GammaV_M, Lat.Cl) )

        File.close()

################################################################################

def ReadFsInput(InFile):
# Read the provided file with freestream data (angles of attack)

    FileRoot = re.search(r'(.+?)(.Alphas)?.in', InFile).group(1)
    
    # Read the file to get angles of attack
    AlphaData = sp.loadtxt(InFile)
    Alphas_Deg = sp.linspace( AlphaData[0], AlphaData[1], AlphaData[2] )

    # Create a freestream for each angle of attack
    Freestreams = [ FreestreamClass( FileRoot, Alpha ) for Alpha in Alphas_Deg ]
#    Freestreams = []
#    for Alpha_Deg in Alphas_Deg:
#        Freestreams.append( FreestreamClass( FileRoot, Alpha_Deg ) )

    return Freestreams

################################################################################

def ReadWingInput(InFile):
# Read the provided input file
    InData = sp.loadtxt(InFile, comments = '#')

    # Read the first line to get global wing properties
    NLattice = int(InData[0,0])
    Patch = PatchClass(NLattice, InData[1,:], InData[2,:])

    # Get the Alpha values over which to run
    #Patch.SetAlphas(*InData[-1, :3])

    return Patch

################################################################################

def InfluenceCoeff(Patches):
# Build an influence coefficient matrix based on all lattices
#   Info for this is given at the end of lecture 18

    # build a list of all lattices
    AllLattices = []
    for Patch in Patches:
        AllLattices.extend(Patch.Lattices)
    
    # Initialize matrices at control points and Trefftz plane
    ControlMat = sp.empty((len(AllLattices),len(AllLattices)))
    ControlMat.fill(sp.nan)
    TrefftzMat = sp.copy(ControlMat)

    # Populate matrices
    for i, Lattice in enumerate(AllLattices):
        for j, Vortex in enumerate(AllLattices):

            ### Build IC matrix for 3/4 chord control points
            # Check for colinearity of control point with bound vortex 
            #   (slide 18-6)
            Tolerance = 1e-6
            a = Vortex.SPt-Lattice.ControlPt; b = Vortex.EPt-Lattice.ControlPt
            if sp.dot( sp.cross(a,b), sp.cross(a,b) ) < Tolerance:
                qbvTerm = 0.0
            else:
                qbvTerm = sp.dot(Lattice.UnitNormal, BS.FiniteFilament( \
                    Vortex.SPt, Vortex.EPt, Lattice.ControlPt))
            qlhsTerm = sp.dot(Lattice.UnitNormal, BS.SemiInfiniteFilament( \
                Vortex.SPt, Lattice.ControlPt, FlipInduced = True))
            qrhsTerm = sp.dot(Lattice.UnitNormal, BS.SemiInfiniteFilament( \
                Vortex.EPt, Lattice.ControlPt))

            ControlMat[i,j] = (1 / (4*sp.pi)) * (qbvTerm + qlhsTerm + qrhsTerm)

            ### Build IC matrix for Trefftz plane control points
            # No bound vortex term
            sTrefftz = sp.copy(Vortex.SPt); sTrefftz[0] = 0.0
            eTrefftz = sp.copy(Vortex.EPt); eTrefftz[0] = 0.0
            qlhsTerm = 2 * sp.dot(Lattice.UnitNormal, BS.SemiInfiniteFilament( \
                sTrefftz, Lattice.TrefftzPt, FlipInduced = True))
            qrhsTerm = 2 * sp.dot(Lattice.UnitNormal, BS.SemiInfiniteFilament( \
                eTrefftz, Lattice.TrefftzPt ))

            TrefftzMat[i,j] = (1 / (4*sp.pi)) * (qlhsTerm + qrhsTerm)

    return ControlMat, TrefftzMat, AllLattices

################################################################################

def RhsMatrix(AllLattices, AlphaW_Deg):
# Given the list of lattices and an angle of attack, compute the right-hand-
#   side matrix (really a vector) for use in finding gamma values.

    # Initialize
    OutVec = sp.empty((len(AllLattices), 1)); OutVec.fill(sp.nan)
    SinAlphaW = sp.sin( sp.radians(AlphaW_Deg) )
    
    for i, Lattice in enumerate(AllLattices):
        # Calculate geometric angle of attack (slide 19-5)
        AlphaG_Rad = sp.arcsin( SinAlphaW * sp.cos(Lattice.DihedralLocal_Rad) )

        # Store the ratio of V_n / V_infinity (slide 19-9)
        OutVec[i] = -1.0 * sp.sin( AlphaG_Rad + Lattice.Twist_Rad - Lattice.AlphaZL_Rad)

    return OutVec

################################################################################

def WriteIntegratedOutput(InFile, AllFreestreams):
# Calculate output and print it to a file

    # Open the file for writing
    FileRoot = re.search(r'(.+?)(.Alphas)?.in', InFile).group(1)
    OutFile = FileRoot + '.Int.out'
    File = open(OutFile, 'w')

    # Print the first part of the file, with wing properties
    File.write('# ' + OutFile + '\n')
    File.write('#\n#\n')

    # Print the parameters which change with angle of attack
    File.write('#%11s %10s %16s\n' % ('Alpha(Deg)', 'CL', 'CDi') )
    for FS in AllFreestreams:
        File.write('%-12.3f %10.4f %16.10f\n' % (FS.AlphaW_Deg, FS.CL, FS.CDi) )

    # Close the file
    File.close()

################################################################################
