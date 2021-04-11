from modeller import *
from modeller.scripts import complete_pdb
import pylab
import modeller


def DOPE_Energy(model, output_path,verbose=False):
    output_path.rstrip("/")
    sys.stdout = open(output_path + "/analysis/" + model + ".DOPE",'w')
    env = Environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
    env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

    # read model file
    model_path = output_path + "/structures/" + model + ".pdb"
    mdl = complete_pdb(env, model_path)

    # Assess with DOPE:
    s = Selection(mdl)   # all atom selection
    file = '%s/analysis/%s.profile'%(output_path,model)
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file= file,
                  normalize_profile=True, smoothing_window=15)

    sys.stdout.close()
    get_profile(file, model, output_path)


def get_profile(profile_file, model, output_path):
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
    pylab.figure(int(model[5]), figsize=(10,6))
    pylab.xlabel('Alignment position')
    pylab.ylabel('DOPE per-residue score')
    pylab.plot(vals, color='red', linewidth=2, label='Model')
    pylab.legend()
    pylab.savefig('%s/analysis/%s.png'%(output_path,model), dpi=65)



if __name__ == "__main__":

    file = DOPE_Energy("model.pdb")
    get_profile(file)
