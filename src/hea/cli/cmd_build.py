"""Defines the ``hea build`` CLI"""

import os

import click

from hea.builder.build_config import BuildConfig
from hea.builder.builder import Builder

from hea.builder import analyze_build

@click.group('build', invoke_without_command=True,
    help="""*** Build HEA structures ***

    Example usage:

        hea build -c build.toml

        hea build init

""", short_help='Build HEA structures')
@click.option('-c', '--config-file', 'config_file',
              type=click.Path(dir_okay=False),
              help='Path to the TOML build configuration file')
@click.option('-o', '--output-directory', 'output_directory',
              type=click.Path(file_okay=False),
              required=False,
              help='Path to output build directory [default: CWD]',
              default=os.getcwd())
@click.option('-p', '--prefix', 'prefix',
              required=False,
              help='Prefix of built files',
              default='HEA',
              show_default=True)
@click.pass_context
def build_cli(ctx, config_file, output_directory, prefix):
    """
    Build CLI

    """

    if ctx.invoked_subcommand is None:

        if not config_file:
            click.echo(ctx.get_help())
            return

        builder = Builder(config_filepath=config_file)
        builder.build(build_directory=output_directory, prefix=prefix)

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

@click.argument('filename', required=True)
@build_cli.command('check',
    help="""Check concentrations in a TRAJ file

    Usage:

        hea build check HEA.traj

    """)
def check(filename):
    """
    Check concentrations in a TRAJ file

    """

    chemistry = analyze_build.get_concentrations(traj_filepath=filename)

    click.echo(analyze_build.report_chemistry(chemistry=chemistry))
