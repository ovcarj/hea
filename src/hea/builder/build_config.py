"""Handles the build.toml configuration file"""

import tomllib
import tomli_w

import ase.build.bulk


class BuildConfig:
    """
    Handles the build.toml configuration file
    """

    def __init__(self, config_filepath=None):
        """
        If ``config_filepath`` is given, reads it to self._cfg

        Parameters
        ----------
        config_filepath : str
            Path to the build configuration file

        """

        self.cfg = None

        if config_filepath:

            with open(config_filepath, mode='rb') as f:
                self.cfg = tomllib.load(f)

    @property
    def _defaults(self):
        """
        Defines default values of the build config file

        Returns
        -------
        defaults : dict
            Dictionary with default values for the build config file

        """

        defaults = {}

        ase_build_bulk_args = ase.build.bulk.__annotations__.keys()

        build_dict = {build_key: 'None'
                      for build_key in ase_build_bulk_args
                      }

        del build_dict['name']
        del build_dict['return']

        build_dict['crystalstructure'] = 'bcc'
        build_dict['a'] = 2.0

        defaults.update({'ASE_BUILD_BULK_ARGUMENTS': build_dict})

        return defaults

    def write_defaults(self, output_filepath):
        """
        Writes self._defaults to ``output_filepath``

        Parameters
        ----------
        output_filepath : str
            Path to which to write the default build.toml config file

        """

        with open(output_filepath, mode='wb') as f:
            tomli_w.dump(self._defaults, f)
