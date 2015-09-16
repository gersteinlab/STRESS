#!python

# This example shows how to calculate approximate low-frequency
# modes for big proteins. For a description of the techniques,
# see
#
#    K. Hinsen
#    Analysis of domain motions by approximate normal mode calculations
#    Proteins 33, 417 (1998)
#
# and
#
#    K. Hinsen, A.J. Petrescu, S. Dellerue, M.C. Bellissent-Funel, G.R. Kneller
#    Harmonicity in slow protein dynamics
#    Chem. Phys. 261, 25 (2000)
#
##
#
# Adaptation by Simon Mitternacht 2010-2011.
#


from MMTK import *
from MMTK.Proteins import Protein
from MMTK.Proteins import PeptideChain
from MMTK.ForceFields import CalphaForceField
from MMTK.FourierBasis import FourierBasis, estimateCutoff
from MMTK.NormalModes import NormalModes, VibrationalModes, SparseMatrixSubspaceNormalModes
from numpy import *
import sys
import re


# top_nm.pl only reads the first 20 modes before exiting. This results in a broken pipe error.
# Prevent calpha_modes.py from crashing by ignoring SIGPIPE signal.
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)



def modes_in_sequential_order(protein,modes):
    modes_inseq = []
    for mode in modes: #[modes.rawMode(i) for i in xrange(len(modes))]:
        tmp_mode = []
        for chain in protein:
            for residue in chain.residues():
                tmp = mode[residue.peptide.C_alpha]
                norm = 1#math.sqrt(residue.mass())
                tmp_vec = [tmp.x()/norm, tmp.y()/norm, tmp.z()/norm]
                tmp_mode.append(tmp_vec)
        modes_inseq.append(tmp_mode)
    return modes_inseq



if (len(sys.argv) != 2 and len(sys.argv) != 3):
    print "Usage: %s pdb [v|f]" % sys.argv[0]
    print "Outputs the 10 % lowest frequency normal modes,to stdout. "
    print "Second optional argument specifies if VibrationalModes or FourierBasis should be used."
    exit()


# Check CL arguments    
pdb_file = sys.argv[1]
type = "auto"
if (len(sys.argv) > 2):
    arg = sys.argv[2]
    if (arg != 'v' and arg != 'f'):
        sys.stderr.write("Argument \'%s\' invalid" % arg)
        exit()
    type = arg
        
# Construct system
universe = InfiniteUniverse(CalphaForceField(2.5))
universe.protein = Protein(pdb_file, model='calpha')

# Find a reasonable basis set size and cutoff
cutoff = None
if (type == "auto") or (type == "f"):
    nbasis = max(10, universe.numberOfAtoms()/10)
    cutoff, nbasis = estimateCutoff(universe, nbasis)

# calculate modes for first system
if (cutoff is None):
    # Do full normal mode calculation
    print "# VibrationalModes"
    modes = VibrationalModes(universe)
    print "# %d" % len(modes)
else:
    # Do subspace mode calculation with Fourier basis for large proteins
    print "# FourierBasis: cutoff = %d, nbasis = %d" % (cutoff, nbasis)
    subspace = FourierBasis(universe, cutoff)
    subspace.may_modify = 1
    modes = SparseMatrixSubspaceNormalModes(universe, subspace)

modes_seq = modes_in_sequential_order(universe.protein, modes)
N = min(len(modes_seq),3*universe.numberOfAtoms()/10+6)

for i in range(6,N):
    print "BEGIN MODE %d" % ((i+1))
    for d in modes_seq[i]:
        print "%.10f %.10f %.10f" % (d[0], d[1], d[2])
    print "END"
    print ""
    

