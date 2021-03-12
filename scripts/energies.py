from modeller import *
from modeller.scripts import complete_pdb
import pylab
import modeller


def DOPE_Energy(model, output_directory, verbose=False):
    """ Calculates the DOPE energy of the model. """

    if verbose:
        print("Calculating the DOPE energy of the model.")
    env = Environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
    env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

    # read model file
    mdl = complete_pdb(env, model)

    print(output_directory)

    # Assess with DOPE:
    s = Selection(mdl)   # all atom selection
    file = output_directory + "/analysis/model.profile"
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file= file,
                  normalize_profile=True, smoothing_window=15)

    return file


def get_profile(profile_file, output_directory):
    """ Creates a plot for the DOPE energy of the model."""

    f = open(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10: #Read `profile_file` into a Python array, and add gaps corresponding to the alignment sequence 'seq'
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
    pylab.savefig('%s/analysis/DOPE_energy.png'%output_directory, dpi=65)


