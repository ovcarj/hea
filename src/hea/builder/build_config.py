"""Handles the build.toml configuration file"""

from functools import cached_property

import tomllib
import tomli_w

import ase.build.bulk


class BuildConfig:
    """
    Handles the build.toml configuration file

    """

    def __init__(self, config_filepath=None):
        """
        Set the path to the TOML configuration file

        Parameters
        ----------
        config_filepath : str
            Path to the build configuration file

        """

        self._config_filepath = config_filepath

    @cached_property
    def cfg(self):
        """
        Read the build configuration file from ``self._config_filepath``

        Since TOML format does not suport NULL values, replace 'None' strings
            with None

        Returns
        -------
        _cfg : None | dict
            The build configuration read to a dictionary;
                if self._config_filepath is None, return None

        """

        if not self._config_filepath:
            return None

        with open(self._config_filepath, mode='rb') as f:
            cfg = tomllib.load(f)

        #
        # Replace 'None' strings with None
        #
        # For now, iterate each section individually
        #
        for key, value in cfg['ASE_BUILD_BULK_ARGUMENTS'].items():
            if value == 'None':
                cfg['ASE_BUILD_BULK_ARGUMENTS'][key] = None

        return cfg

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
