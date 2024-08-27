#!/home/cebner/anaconda/bin/python
# PlotCl.py
#
# Read output Prandtl.py files and create various plots

from __future__ import division
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import scipy as sp

# Open pdf to plot to
pp_Cl = PdfPages('ClPlots.Taper.pdf')
pp_cCl = PdfPages('cClPlots.Taper.pdf')
pp_CL = PdfPages('CLPlots.Taper.pdf')

# Collect integrated data
WIntData = sp.loadtxt('TaperWeissinger/Taper.Int.out')
WLDict = dict( zip( map(int, WIntData[:,0]), WIntData[:, 1]) )
PIntData = sp.loadtxt('RunPrandtl/Taper.Int.out', skiprows = 6)
PLDict = dict( zip( map(int, PIntData[:,0]), PIntData[:, 1]) )

# Loop over all alphas
for Alpha in range(11):
    WeissFile = 'TaperWeissinger/Taper.a%i.0.out' % int(Alpha)
    PranFile = 'RunPrandtl/Taper.a%i.0.Cl.out' % int(Alpha)

    # Load data
    WeissData = sp.loadtxt( WeissFile )
    PranData = sp.loadtxt( PranFile )

    # Average chords:
    cAvWeiss = sp.mean(WeissData[:, 3])
    cAvPran = sp.mean(PranData[:, 1])

    Pran_Y = -5 * sp.cos(PranData[:, 0])

    ### Plot vs y
    # Cl vs Y
    Fig_Cl = plt.figure(Alpha)
    Ax = Fig_Cl.add_subplot(1,1,1)
    Ax.plot(WeissData[:, 1], WeissData[:, 5], 'r-', \
        label = "Weissinger's Method" )
    Ax.plot(Pran_Y, PranData[:, 2], 'b-', label = 'Lifting Line Theory')

    Ax.set_title(r'Local Lift Coefficient Distribution For $\alpha = ' \
        + str(Alpha) + r'^\circ$')
    Ax.set_xlabel('Y (m)', fontsize = 16)
    Ax.set_ylabel(r'$C_l$', fontsize = 16)
    Ax.legend(loc = 'best')
    Ax.set_ylim([0,1.4])
    Ax.grid('on')

    Fig_Cl.savefig(pp_Cl, format = 'pdf')

    # c Cl vs Y
    Fig_cCl = plt.figure(Alpha + 100)
    Ax_cCl = Fig_cCl.add_subplot(1,1,1)
    Ax_cCl.plot(WeissData[:, 1], WeissData[:, 3] * WeissData[:, 5] / \
        (cAvWeiss * WLDict[Alpha]), 'r-', label = "Weissinger's Method" )
    Ax_cCl.plot(Pran_Y, PranData[:, 2] * PranData[:, 1] / \
        (cAvPran * PLDict[Alpha]), 'b-', label = 'Lifting Line Theory')

    Ax_cCl.set_title("Span Loading Comparison Between Weissinger's Method\n" \
        + r"and Lifting Line Theory for $\alpha =" + str(Alpha) + r'^\circ$')
    Ax_cCl.set_xlabel('Y (m)', fontsize = 16)
    Ax_cCl.set_ylabel(r'$\frac{c\ C_l}{c_{avg}\ C_L}$', fontsize = 16)
    Ax_cCl.legend(loc = 'best')
    Ax_cCl.set_ylim([0,2.0])
    Ax_cCl.grid('on')

    Fig_cCl.savefig(pp_cCl, format = 'pdf')
    
pp_Cl.close()
pp_cCl.close()

### Plot CL and CDi as functions of Alpha
# CL vs. Alpha
Fig_CL = plt.figure(1000)
Ax_CL = Fig_CL.add_subplot(1,1,1)
Ax_CL.plot(WIntData[:,0], WIntData[:,1], 'r-', \
    label = "Weissinger's Method" )
Ax_CL.plot(PIntData[:,0], PIntData[:,1], 'b-', \
    label = "Lifting Line Theory" )

Ax_CL.set_title("Lift Coefficient Comparison Between\nWeissinger's" \
    " Method and Lifting Line Theory")
Ax_CL.set_xlabel(r'$\alpha$ (Degrees)', fontsize = 16)
Ax_CL.set_ylabel(r'$C_L$', fontsize = 16)
Ax_CL.legend(loc = 'best')
#Ax_CL.set_ylim([0,1.4])
Ax_CL.grid('on')

# CDi vs. Alpha
Fig_CD = plt.figure(1001)
Ax_CD = Fig_CD.add_subplot(1,1,1)
Ax_CD.plot(WIntData[:,0], abs(WIntData[:,2]), 'r-', label = "Weissinger's Method" )
Ax_CD.plot(PIntData[:,0], PIntData[:,2], 'b-', label = "Lifting Line Theory" )

Ax_CD.set_title("Drag Coefficient Comparison Between\nWeissinger's" \
    " Method and Lifting Line Theory")
Ax_CD.set_xlabel(r'$\alpha$ (Degrees)', fontsize = 16)
Ax_CD.set_ylabel(r'$C_{D_i}$', fontsize = 16)
Ax_CD.legend(loc = 'best')
#Ax_CD.set_ylim([0,1.4])
Ax_CD.grid('on')

# CDi vs. CL
Fig_CDL = plt.figure(1002)
Ax_CDL = Fig_CDL.add_subplot(1,1,1)
Ax_CDL.set_title("Effects of Lift Coefficient on Induced Drag Coefficient")
Ax_CDL.set_xlabel(r'$C_L$', fontsize = 16)
Ax_CDL.set_ylabel(r'$C_{D_i}$', fontsize = 16)
Ax_CDL.grid('on')
Ax_CDL.plot(WIntData[:,1], abs(WIntData[:,2]), 'r-', label = "Weissinger's Method")
Ax_CDL.plot(PIntData[:,1], abs(PIntData[:,2]), 'b-', label = "Lifting Line Theory")
Theory_CL = sp.linspace(0, 1.2, 101)
Theory_CDi = Theory_CL**2 / (sp.pi * 10.0)
Ax_CDL.plot(Theory_CL, Theory_CDi, 'k--', label = "Elliptical Wing: " + \
    r'$C_{D_i} = \frac{C_L^2}{\pi\ AR}$')
Ax_CL.legend(loc = 'best')
Ax_CD.legend(loc = 'best')
Ax_CDL.legend(loc = 'best')

Fig_CL.savefig(pp_CL, format = 'pdf')
Fig_CD.savefig(pp_CL, format = 'pdf')
Fig_CDL.savefig(pp_CL, format = 'pdf')

pp_CL.close()
