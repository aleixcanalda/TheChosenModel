from modeller import *
from modeller.scripts import complete_pdb
import pylab
import modeller


def DOPE_Energy(model,verbose=False):
    env = Environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
    env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

    # read model file
    #model = Model(model)
    mdl = complete_pdb(env, model)

    # Assess with DOPE:
    s = Selection(mdl)   # all atom selection
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='%s.profile'%model,
                  normalize_profile=True, smoothing_window=15)

    return '%s.profile'%model


def get_profile(profile_file):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    f = open(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)

    # Plot the template and model profiles in the same plot for comparison:
    pylab.figure(1, figsize=(10,6))
    pylab.xlabel('Alignment position')
    pylab.ylabel('DOPE per-residue score')
    pylab.plot(vals, color='red', linewidth=2, label='Model')
    pylab.legend()
    pylab.savefig('%s.png'%profile_file, dpi=65)



if __name__ == "__main__":

    file = DOPE_Energy("model.pdb")
    get_profile(file)


