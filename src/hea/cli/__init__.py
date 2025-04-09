"""Automatically define the ``hea <command>`` CLI

The modules in this directory must be named cmd_<command>.py
    and include a function named <command>_cli()

"""

import pkgutil
import click

@click.group(name='hea')
def hea_cli():
    """
    Entry point for the hea command

    """

for loader, module_name, is_pkg in pkgutil.walk_packages(__path__):

    if module_name.startswith('cmd'):

        cmd_name = module_name.split('_')[1]
        cmd_cli = f'{cmd_name}_cli' #pylint: disable=invalid-name

        module = loader.find_spec(module_name).loader.load_module(module_name)
        command = getattr(module, cmd_cli)

        hea_cli.add_command(command, name=cmd_name)
