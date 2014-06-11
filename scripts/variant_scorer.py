#!/usr/bin/env python
# encoding: utf-8
"""
run_genmod.py

Script for annotating genetic models in variant files.

Created by Måns Magnusson on 2014-01-21.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
from codecs import open
from pkg_resources import require
from tempfile import NamedTemporaryFile
from datetime import datetime

from pprint import pprint as pp

from score_mip_variants import scorer, variant_sorter
from ped_parser import parser
from vcf_parser import vcf_parser


def check_file_existence(infile):
    """Check is the file exists. Quit if something is wrong."""
    if not os.path.isfile(infile):
        print('The file %s does not exist!!!\nPlease check what is wrong and rerun.\nExiting...' % infile)
        sys.exit()

def get_family(family_file):
    """Return the family"""
    family_type = 'cmms'
    
    my_family_parser = parser.FamilyParser(family_file, family_type)
    # Stupid thing but for now when we only look at one family
    return my_family_parser.families.popitem()[1]

def add_metadata(head, args):
    """Add metadata for the information added by this script."""
    
    head.add_info('RS', '1', 'Integer', "Combined ranc score for this variant in this family.")
    return

def print_headers(head, args):
    """Print the headers to a results file."""
    if args.outfile[0]:
        with open(args.outfile[0], 'w', encoding='utf-8') as f: 
            for head_count in head.print_header():
                f.write(head_count+'\n')
    else:
        if not args.silent:
            for line in head.print_header():
                print(line)
    return

def main():
    
    info_string = """Individuals that are not present in ped file will not be considered in the analysis."""
    
    parser = argparse.ArgumentParser(description="Score the variants in a vcf file based on the MIP system.")
    
    parser.add_argument('family_file',
        type=str, nargs=1,
        help='A pedigree file in .ped format.'
    )
    parser.add_argument('variant_file',
        type=str, nargs=1,
        help='A variant file. Default is vcf format.'
    )
    
    parser.add_argument('--version', 
        action="version", 
        version=require("score_mip_variants")[0].version
    )
    
    parser.add_argument('-v', '--verbose', 
        action="store_true", 
        help='Increase output verbosity.'
    )
    
    parser.add_argument('-s', '--silent', 
        action="store_true", 
        help='Do not print the variants.'
    )
    
    parser.add_argument('-o', '--outfile', 
        type=str, nargs=1, default=[None],
        help='Specify the path to a file where results should be stored.'
    )
    
    # parser.add_argument('-family', '--family_type', 
    #     type=str, nargs=1, default=['ped'], 
    #     choices=['ped', 'alt', 'cmms', 'mip'],
    #     help='If the analysis use one of the known setups, please specify which.'
    # )
    
    args = parser.parse_args()
    
    fam_file = args.family_file[0]
    var_file = args.variant_file[0]
    
    start_time_analysis = datetime.now()
    
    check_file_existence(var_file)
    check_file_existence(fam_file)
    # file_name, file_extension = os.path.splitext(var_file)
    
    # Start by parsing at the pedigree file:
    
    my_family = get_family(fam_file)
    
    # # Check the variants:
    
    my_vcf_parser = vcf_parser.VCFParser(var_file)
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
    add_metadata(head, args)
    print_headers(head, args)
    
    var_sorter = variant_sorter.FileSort(temporary_variant_file, outFile=args.outfile[0])
    var_sorter.sort()
    
    temporary_variant_file.close()
    

if __name__ == '__main__':
    main()

