"""Defines the ``hea build`` CLI"""

import os

import click

from hea.builder.build_config import BuildConfig

@click.group('build', invoke_without_command=True,
    help="""*** Build HEA structures ***

    Example usage:

        hea build -c build.toml

        hea build init

""", short_help='Build HEA structures')
@click.option('-c', '--config-file', 'config_file',
              type=click.Path(dir_okay=False),
              help='Path to the TOML build configuration file')
@click.pass_context
def build_cli(ctx, config_file):
    """
    Build CLI

    """

    if ctx.invoked_subcommand is None:
        click.echo('hea build <build configuration file> not yet implemented.')

@build_cli.command('init',
    help="""Create a default build configuration file in the CWD

    Usage:
        
        hea build init

    """)
@click.option('-f', '--filename', required=False,
              default='build.toml', show_default=True,
              help='Name of the TOML build configuration file')
def init(filename):
    """
    Create a default build configuration file in CWD

    """

    cwd = os.getcwd()

    path = os.path.join(cwd, filename)

    build = BuildConfig()
    build.write_defaults(output_filepath=path)
