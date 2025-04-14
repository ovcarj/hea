"""Get the statistics of atomic species concentrations"""

import ase.io
import ase.atoms

import numpy as np

from scipy import stats

def get_concentrations(traj_filepath):
    """
    Get statistics on the concentrations of unique atomic species in traj file

    It is assumed that all structures in the traj file have the same
        unique species as the first structure in the traj file.

    The form of the output ``chemistry`` dictionary is:

    {

        'CHEMICAL_SPECIES_1':
            {
            'concentration': concentration1,
            'uncertainty': uncertainty1
            },

        'CHEMICAL_SPECIES_2':
            {
            'concentration': concentration2,
            'uncertainty': uncertainty2
            },

            ...

    Parameters
    ----------
    traj_filepath : str
        Path to an ASE traj file

    Returns
    -------
    chemistry : dict
        Dictionary containing atomic species and concentrations

    """

    traj = ase.io.read(filename=traj_filepath, index=':')

    if isinstance(traj, ase.atoms.Atoms):
        traj = [traj]

    unique_species = traj[0].symbols.species()

    chemistry = {

            species: {'concentration': None, 'uncertainty': None,
                      'concentrations_per_struct': []}

            for species in unique_species

                }

    for atoms in traj:

        for species in unique_species:

            chemistry[species]['concentrations_per_struct'].append(
                    100.0 * float(atoms.symbols.count(species)) / float(len(atoms))
                    )

    for species in unique_species:

        chemistry[species]['concentration'] =\
                np.mean(chemistry[species]['concentrations_per_struct'])

        chemistry[species]['uncertainty'] =\
                np.std(chemistry[species]['concentrations_per_struct'], ddof=1)

    return chemistry

def report_chemistry(chemistry):
    """
    Get a formatted report on the chemistry statistics

    Parameters
    ----------
    chemistry : dict
        Dictionary containing atomic species and concentrations

    Returns
    -------
    report : str
        Report on the concentrations and standard deviations

    """

    report = '\n'

    for species in chemistry:
        report += f'{species:<2}: {np.round(chemistry[species]["concentration"], 2):<5} +/- {np.round(chemistry[species]["uncertainty"], 2)}\n'

    return report
