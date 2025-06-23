"""Builds HEA structures"""

import os

from copy import deepcopy

import numpy as np

import ase.io
import ase.build.bulk
import ase.atoms

from hea.builder.build_config import BuildConfig


class Builder: #pylint: disable=too-few-public-methods
    """
    Builds HEA structures

    """

    def __init__(self, config_filepath):
        """
        Passes ``config_filepath`` to BuildConfig and creates a unit structure

        If ``build_directory`` is None, defaults to ``{cwd}/hea_{current_datetime}``

        Parameters
        ----------
        config_filepath : str
            Path to the build configuration file

        """

        self._cfg = BuildConfig(config_filepath=config_filepath).cfg

    def build(self, build_directory=None, prefix='HEA', write_immediately=True):
        """
        Create HEA structures according to ``self._cfg`` in ``build_directory``

        If ``build_directory`` is None, defaults to ``{cwd}/hea_{current_datetime}``

        Parameters
        ----------
        build_directory : str
            Path to the directory in which to store the built structures;
                if None, defaults to CWD
        prefix : str
            Prefix for the built files; defaults to 'HEA'
        write_immediately : bool
            Immediately writes structures as they are created; does not keep structures in memory

        """

        if not build_directory:
            build_directory = os.getcwd()

        os.makedirs(build_directory, exist_ok=True)

        # Build unit template structure

        unit_template = self._build_unit_template(**self._cfg['ASE_BUILD_BULK_ARGUMENTS'])

        # Write unit template structure

        fpath_prefix = os.path.join(build_directory, prefix)

        ase.io.write(filename=f'{fpath_prefix}_template_unit.traj',
                     images=unit_template,
                     format='traj')

        # Build supercell template structure

        supercell_template = self._build_supercell_template(
                unit_template=unit_template,
                rep=self._cfg['BATCH_PARAMETERS']['supercell']
                )

        # Write supercell template structure

        ase.io.write(filename=f'{fpath_prefix}_template_supercell.traj',
                     images=supercell_template,
                     format='traj')

        # Build a batch of structures

        structures = []

        #
        # If no cell scaling, build n_structs_per_cell structures
        # from the supercell template
        #

        if not self._cfg['BATCH_PARAMETERS']['scale_cell']:

            if not write_immediately:

                structures.extend(
                    self._build_random_batch(template_structure=supercell_template,
                                    chemistry=self._cfg['CHEMISTRY'],
                                    n_structs=self._cfg['BATCH_PARAMETERS']['n_structs_per_cell'],
                                    seed=self._cfg['BATCH_PARAMETERS']['rng_seed'])
                                    )

            else:
                self._build_random_batch(template_structure=supercell_template,
                                chemistry=self._cfg['CHEMISTRY'],
                                n_structs=self._cfg['BATCH_PARAMETERS']['n_structs_per_cell'],
                                seed=self._cfg['BATCH_PARAMETERS']['rng_seed'],
                                output_filename=f'{fpath_prefix}.traj')

        #
        # If "cell scaling", build n_structs_per_cell structures
        # for each scaling value in:
        #
        # np.linspace(1.0, max_scaling, n_scaled)
        #
        # "add_edge_points" structures will be added at below (above)
        # the maximum scaling values
        #

        else:

            scaling = np.linspace(start=1.0,
                                  stop=self._cfg['BATCH_PARAMETERS']['max_scaling'],
                                  num=self._cfg['BATCH_PARAMETERS']['n_scaled'],
                                  endpoint=True)

            if self._cfg['BATCH_PARAMETERS']['add_edge_points']:

                step = scaling[1] - scaling[0]

                edge_low = np.arange(
                        start=scaling[0] -\
                                step * (self._cfg['BATCH_PARAMETERS']['add_edge_points']),
                        stop=scaling[0],
                        step=step
                        )

                edge_high = np.arange(
                        start=scaling[-1] + step,
                        stop=scaling[-1] +\
                                step * (self._cfg['BATCH_PARAMETERS']['add_edge_points'] + 1),
                        step=step
                        )

                scaling = np.concatenate((edge_low, scaling, edge_high))

            original_cell = unit_template.cell

            for scale in scaling:

                scaled_cell = original_cell * scale

                print(f'Building {self._cfg["BATCH_PARAMETERS"]["n_structs_per_cell"]} structures for cell\n\n{scaled_cell}\n\n') #pylint: disable=line-too-long

                unit_copy = deepcopy(unit_template)
                unit_copy.set_cell(cell=scaled_cell, scale_atoms=True)

                scaled_supercell_template = self._build_supercell_template(
                        unit_template=unit_copy,
                        rep=self._cfg['BATCH_PARAMETERS']['supercell']
                        )

                if not write_immediately:

                    structures.extend(
                        self._build_random_batch(template_structure=supercell_template,
                                        chemistry=self._cfg['CHEMISTRY'],
                                        n_structs=self._cfg['BATCH_PARAMETERS']['n_structs_per_cell'],
                                        seed=self._cfg['BATCH_PARAMETERS']['rng_seed'])
                                        )

                else:
                    self._build_random_batch(template_structure=supercell_template,
                                    chemistry=self._cfg['CHEMISTRY'],
                                    n_structs=self._cfg['BATCH_PARAMETERS']['n_structs_per_cell'],
                                    seed=self._cfg['BATCH_PARAMETERS']['rng_seed'],
                                    output_filename=f'{fpath_prefix}.traj')

        # Write batch .traj file

        if not write_immediately:

            ase.io.write(filename=f'{fpath_prefix}.traj',
                         images=structures,
                         format='traj')

    def _build_unit_template(self, **kwargs):
        """
        Builds a template unit structure

        Parameters
        ----------
        kwargs : dict
            Dict of parameters passed to ase.build.bulk, except for ``name``

        Returns
        -------
        unit_template : ase.atoms.Atoms
            A template unit structure

        """

        #
        # Number of atoms in primitive cell as defined by ASE
        #

        structures = {'sc': 1, 'fcc': 1, 'bcc': 1,
                  'tetragonal': 1,
                  'bct': 1,
                  'hcp': 1,
                  'rhombohedral': 1,
                  'orthorhombic': 1,
                  'mcl': 1,
                  'diamond': 1,
                  'zincblende': 2, 'rocksalt': 2, 'cesiumchloride': 2,
                  'fluorite': 3, 'wurtzite': 2}

        # Fill template with C atoms

        default_name = 'C' * structures[kwargs['crystalstructure']]

        unit_template = ase.build.bulk(name=default_name, **kwargs)

        return unit_template

    def _build_supercell_template(self, unit_template, rep):
        """
        Builds a template supercell structure

        Parameters
        ----------
        unit_template : ase.atoms.Atoms
            Template unit structure to be repeated
        rep : tuple
            How many times to repeat the structure
                in the 3 unit cell directions; e.g. (2, 3, 4)

        Returns
        -------
        supercell_template : ase.atoms.Atoms
            A template supercell structure

        """

        supercell_template = unit_template.repeat(rep=rep)

        return supercell_template

    def _build_random_batch(self,
                            template_structure,
                            chemistry,
                            n_structs,
                            seed=None,
                            output_filename=None):
        """
        Creates a batch of random HEA structures

        ``n_structs`` structures will be created.

        The ``template_structure`` is filled randomly with atoms whose
            types and probabilities are defined in the ``chemistry``
            dictionary.

            The form of the ``chemistry`` dictionary is:

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

            }

        Parameters
        ----------
        template_structure : ase.atoms.Atoms
            The template structure from which HEA will be constructed
        chemistry : dict
            Dictionary containing atomic species and concentrations
        n_structs : int
            Number of structures to create
        seed : int
            RNG seed
        output_filename : str
            Path to the output ASE .traj file.
                If not given, a list of ASE Atoms objects will be returned

        Returns
        -------
        batch : list of ase.atoms.Atoms
            Returned only if ouput_filename is None

        """

        batch = []

        positions = template_structure.get_positions()
        cell = template_structure.get_cell()
        natoms = len(template_structure)

        generator = np.random.default_rng(seed=seed)

        rng_integers = generator.integers(low=0, high=1e6, size=n_structs)

        for i in range(n_structs):

            struct_seed = rng_integers[i] + 129 * i

            random_species = self._get_random_species(
                    chemistry=chemistry, natoms=natoms, seed=struct_seed)

            atoms = ase.atoms.Atoms(symbols=random_species,
                                    positions=positions,
                                    cell=cell,
                                    pbc=True)

            if not output_filename:
                batch.append(atoms)

            else:

                with ase.io.trajectory.Trajectory(output_filename, mode='a') as traj:
                    traj.write(atoms)

        if not output_filename:
            return batch

        print(f'{n_structs} structures written to {output_filename} !')

    def _get_random_species(self, chemistry, natoms, seed=None):
        """
        Return a list of random atomic species using the ``chemistry`` dictionary

        Parameters
        ----------
        chemistry : dict
            Dictionary containing atomic species and concentrations
        natoms : int
            Desired length of the returned list of atomic species
        seed : int
            RNG seed

        Returns
        -------
        random_species_list : list
            List of atomic species

        """

        generator = np.random.default_rng(seed=seed)

        normal_chemistry = [

                {
                    'species': atomic_species,
                    'concentration': generator.normal(
                        loc=chemistry[atomic_species]['concentration'],
                        scale=chemistry[atomic_species]['uncertainty'])
                }

                for atomic_species in chemistry.keys()

                        ]

        ranged_chemistry = []

        total_concentration = 0.0

        for norm_chemistry in normal_chemistry:

            ranged_chemistry.append(
                    {
                     'species': norm_chemistry['species'],
                     'range': [total_concentration,
                               total_concentration + norm_chemistry['concentration']]
                    })

            total_concentration += norm_chemistry['concentration']

        random_numbers = generator.uniform(low=0.0,
                                           high=total_concentration,
                                           size=natoms)

        random_species_list = []

        for random_number in random_numbers:
            for range_chemistry in ranged_chemistry:

                if range_chemistry['range'][0] <= random_number < range_chemistry['range'][1]:
                    random_species_list.append(range_chemistry['species'])
                    break

        return random_species_list
