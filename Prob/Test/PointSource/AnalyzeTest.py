"""
This is an analysis script to analyze a criptic test. If the
script is run with no arguments, it is produce one or more analysis
plots, labelled ProbNameN.pdf, where ProbName is the name of the test
problem and N = 0, 1, 2, .... It also prints "PASS" or "FAIL" to
report of the result is acceptable. If run with the argument
--testonly, it just reports pass/fail, and does not produce the plots.
"""

import argparse
import yt
import yt.units as u
from yt.utilities.physical_constants import c, mp
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os.path as osp
from scipy.special import erfc

# Name of the problem this script is for
prob_name = "IsoDiff"

# Check if we are supposed to produce full outputs
parser = argparse.ArgumentParser(
    description="Analysis script for the "+prob_name+" test")
parser.add_argument("--testonly",
                    help="report test results only, do not make plots",
                    default=False, action="store_true")
parser.add_argument("-p", "--plt",
                    help="number of plot file to examine (default) "
                    "is last plot file in directory",
                    default=-1, type=int)
parser.add_argument("-d", "--dir",
                    help="directory containing run to analyze",
                    default=".", type=str)
args = parser.parse_args()

# Figure out which plot file to read
if args.plt >= 0:
    pltnum = "{:05d}".format(args.plt)
else:
    pltnames = glob(osp.join(args.dir, "plt[0-9][0-9][0-9][0-9][0-9]"))
    pltnums = [int(osp.basename(p)[3:9]) for p in pltnames]
    pltnum = "{:05d}".format(np.amax(pltnums))

# Load the plot file
ds = yt.load(osp.join(args.dir, "plt"+pltnum))

# Get units used in calculation
L = ds.unit_system["length"]
T = ds.unit_system["time"]
E = ds.unit_system["energy"]

# Figure out the properties of the sources
eCR = np.array(ds.parameters["prob.source_e_min"]) * u.GeV
Tm = eCR / (mp * c**2)
pCR = mp * c * np.sqrt(Tm * (2+Tm))
lum = np.array(ds.parameters["prob.source_L"]) * E/T

# Read the diffusion coefficient formula used in the problem, and
# figure out the diffusion coefficient that applies to each of the
# sources
kappa0 = float(ds.parameters["crprop.kPar0"]) * L**2/T
kappa_idx = float(ds.parameters["crprop.kParIdx"])
kappa = kappa0 * (pCR / (mp*c))**kappa_idx

# Make profile of energy density vs radius
rmax = ds.domain_right_edge[0] / 2
nbins = int(0.2*ds.domain_dimensions[0])
dx = (ds.domain_right_edge[0] - ds.domain_left_edge[0]) / \
     ds.domain_dimensions[0]
en_flds = []
for e in eCR.to_value('GeV'):
    for _, f in ds.field_list:
        if not f.startswith('ECR'):
            continue
        else:
            T0 = float(f.split('_')[1])
            T1 = float(f.split('_')[3])
            if T0 <= e and e < T1:
                en_flds.append(f)
sp = ds.sphere("c", rmax)
prof = yt.create_profile(sp, "radius", en_flds, n_bins = nbins,
                         extrema={"radius": (2*dx,rmax)}, 
                         logs={"radius": False}, weight_field="cell_volume")
r = prof.x

# Compute analytic solution
eden_exact = []
rel_err = []
for k, l in zip(kappa, lum):
    # Note force-conversion to flat data here, because unyt doesn't
    # recognise erfc
    t = ds.current_time
    xi = np.array(r/(2*np.sqrt(k*t)))
    eden_exact.append(l / (8*np.pi * k**1.5 * t**0.5) * erfc(xi) / xi)

# Compute relative error and decide pass/fail
max_err = 0.0
for i in range(len(kappa)):
    rel_err = np.amax((prof[en_flds[i]] - eden_exact[i]) / eden_exact[i])
    if rel_err > max_err: max_err = rel_err
if rel_err < 0.1:
    print("PASS")
else:
    print("FAIL")    

# Make plots if asked to do so
if not args.testonly:

    # Make fonts look good
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    # Exact solution versus analytic solution
    fig = plt.figure(1, figsize=(3.5,2.5))
    for i in range(len(kappa)):
        plt.plot(r.to_value('pc'),
                 eden_exact[i].to_value('eV/cm**3'),
                 'C{:d}'.format(i))
        plt.plot(r.to_value('pc'),
                 prof[en_flds[i]].to_value('eV/cm**3'), 'o',
                 mec='k', mfc='C{:d}'.format(i),
                 label=str(eCR[i]) + ' (Sim)')
    plt.plot([-1,-1], [0,0], 'k', label='Exact')
    plt.xlim([0, 1.05*rmax.to_value('pc')])
    plt.yscale('log')
    plt.legend()
    plt.xlabel(r'$r$ [pc]')
    plt.ylabel(r'$U_{{\mathrm{{CR}}}}$ [eV cm$^{-3}$]')
    plt.subplots_adjust(top=0.95, right=0.95, bottom=0.18, left=0.18)

    # Save
    plt.savefig(prob_name+"1.pdf")

    # Make slice averaging together the two zones on either side of z
    # = 0 for highest energy bin
    en = -1
    slc = ds.slice('z', -dx/2)
    frb = slc.to_frb(ds.domain_right_edge[0],
                     ds.domain_dimensions[0:2])
    img = np.array(frb[en_flds[en]].to_value('eV/cm**3'))
    slc = ds.slice('z', dx/2)
    frb = slc.to_frb(ds.domain_right_edge[0],
                     ds.domain_dimensions[0:2])
    img = 0.5 * (img + np.array(frb[en_flds[en]].to_value('eV/cm**3')))
    x = np.linspace(-rmax.to_value('pc'),
                    rmax.to_value('pc'),
                    2)

    # Make plot
    fig = plt.figure(2, figsize=(3.5,2.5))
    plt.imshow(np.log10(img+1.0e-100), extent=(x[0],x[-1],x[0],x[-1]),
               vmin=-0.5, vmax=2.5)
    plt.gca().set_aspect('equal')
    plt.colorbar(label=r'$\log\,U_{\mathrm{CR}}$ [eV cm$^{-3}$]')
    plt.text(-4.2, 4, str(eCR[en]), color='w')
    plt.xlabel(r'$x$ [pc]')
    plt.ylabel(r'$y$ [pc]')
    plt.subplots_adjust(top=0.95, right=0.92, bottom=0.18, left=0.02)

    # Save
    plt.savefig(prob_name+"2.pdf")
