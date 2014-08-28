#!/usr/bin/env python
# encoding: utf-8
"""
run_genmod.py

Script for annotating genetic models in variant files.

Created by Måns Magnusson on 2014-01-21.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function
from __future__ import unicode_literals


import sys
import os
import argparse
import inspect
import pkg_resources
import click

from codecs import open
from tempfile import NamedTemporaryFile
from datetime import datetime

from pprint import pprint as pp

from score_mip_variants import scorer, variant_sorter
from ped_parser import parser as ped_parser
from vcf_parser import parser as vcf_parser


Version = pkg_resources.require("score_mip_variants")[0].version

def check_file_existence(infile):
    """Check is the file exists. Quit if something is wrong."""
    if not os.path.isfile(infile):
        print('The file %s does not exist!!!\nPlease check what is wrong and rerun.\nExiting...' % infile)
        sys.exit()

def get_family(family_file, family_type):
    """Return the family"""    
    my_family_parser = ped_parser.FamilyParser(family_file, family_type)
    # Stupid thing but for now when we only look at one family
    return my_family_parser.families.popitem()[1]

def add_metadata(head, command_line_string):
    """Add metadata for the information added by this script."""
    
    head.add_info('RS', '1', 'Integer', "Combined ranc score for this variant in this family.")
    head.add_version_tracking('score_mip_variants', Version, str(datetime.now()), command_line_string)
    
    return

def print_headers(head, outfile=None, silent=False):
    """Print the headers to a results file."""
    if outfile:
        with open(outfile, 'w', encoding='utf-8') as f: 
            for head_count in head.print_header():
                f.write(head_count+'\n')
    else:
        if not silent:
            for line in head.print_header():
                print(line)
    return


def print_version(ctx, param, value):
    """Callback function for printing version and exiting
    Args:
        ctx (object) : Current context
        param (object) : Click parameter(s)
        value (boolean) : Click parameter was supplied or not
    Returns:
        None:
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo('score_mip_variants version: ' + Version)
    ctx.exit()


@click.command()
@click.argument('family_file',
        nargs=1,
        type=click.Path(exists=True),
        metavar='<ped_file>'
)
@click.argument('variant_file', 
        nargs=1, 
        type=click.Path(exists=True),
        metavar='<vcf_file>'
)
@click.option('--version',
        is_flag=True,
        callback=print_version,
        expose_value=False,
        is_eager=True
)
@click.option('-s' ,'--silent', 
        is_flag=True,
        help='Do not print the variants.'
)
@click.option('-o', '--outfile', 
        type=click.Path(exists=False),
        help='Specify the path to a file where results should be stored.'
)
@click.option('--family_type', '-fam',
        type=click.Choice(['ped', 'alt', 'cmms', 'mip']),
        default='cmms',
        help='If the analysis use one of the known setups, please specify which one. Default is cmms.'
)
def variant_scorer(family_file, variant_file, family_type, silent, outfile):
    """Score the variants in a vcf file based on the MIP system."""
    
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    argument_list = [i+'='+str(values[i]) for i in values if values[i] and i != 'args' and i != 'frame' and i != 'parser']
    
    start_time_analysis = datetime.now()
    
    # Start by parsing at the pedigree file:
    
    my_family = get_family(family_file, family_type)
    
    # # Check the variants:
    
    my_vcf_parser = vcf_parser.VCFParser(variant_file)
    head = my_vcf_parser.metadata
    
    temporary_variant_file = NamedTemporaryFile(mode='w')
    
    if set(my_family.individuals.keys()) != set(my_vcf_parser.individuals):
        
        print('There must be same individuals in ped file and vcf file! Aborting...')
        print('Individuals in PED file: %s' % '\t'.join(list(my_family.individuals.keys())))
        print('Individuals in VCF file: %s' % '\t'.join(list(my_vcf_parser.individuals)))
        sys.exit()
    
    my_parser = scorer.VariantScorer(my_vcf_parser, temporary_variant_file, my_family.models_of_inheritance)
    my_parser.parse()
    temporary_variant_file.seek(0)
    
    
    # Add the new metadata to the headers:
    add_metadata(head, ','.join(argument_list))
    print_headers(head, outfile, silent)
    
    var_sorter = variant_sorter.FileSort(temporary_variant_file, outFile=outfile, silent=silent)
    var_sorter.sort()
    
    temporary_variant_file.close()
    

if __name__ == '__main__':
    variant_scorer()

