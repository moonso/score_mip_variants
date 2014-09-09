#!/usr/bin/env python
# encoding: utf-8
"""
vcf_parser.py

Parse a file with variant info in vcf format.

Creates batches and put them into a queue object.
The batches are dictionary objects with overlapping features where the feature
id:s are keys and the values are dictionarys with variants.


Batch = {feature_1_id:{variant_1_id:variant_1_info, variant_2_id:
variant_2_info}, feature_2_id:... }

Created by Måns Magnusson on 2014-03-17.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function
from __future__ import unicode_literals


import sys
import os
import argparse
from datetime import datetime
from tempfile import NamedTemporaryFile

from pprint import pprint as pp

from score_mip_variants import score_model


class VariantScorer(object):
    """Creates parser objects for parsing variant files"""
    def __init__(self, variant_parser, variant_out_file, models_of_inheritance,
                 alt_dict, score_dict, value_dict, operation_dict, verbose):
        super(VariantScorer, self).__init__()
        self.variant_parser = variant_parser
        self.temp_file = variant_out_file
        self.prefered_models = models_of_inheritance
        self.alt_dict = alt_dict
        self.score_dict = score_dict
        self.value_dict = value_dict
        self.operation_dict = operation_dict
        self.verbose = verbose
        # self.verbosity = args.verbose
        # self.phased = args.phased

    def score_compounds(self, batch):
        """Score the compounds in a batch."""
        for var in batch:
            if batch[var]['info_dict'].get('Compounds', None):
                compounds = batch[var]['info_dict']['Compounds'].split(',')
                comp_list = []
                for comp in compounds:
                    comp_score = (int(batch[var].get(
                                  'Individual_rank_score', 0))
                                  + int(batch.get(comp, {}).get(
                                  'Individual_rank_score', 0)))
                    comp_list.append(comp+'>'+str(comp_score))
                    batch[var]['Compounds'] = ','.join(comp_list)

    def print_variants(self, batch, header, outfile):
        """Prints the variants to a file"""
        for variant in batch:
            # Modify the INFO field:
            info_field = batch[variant]['INFO'].split(';')
            for pos in range(len(info_field)):
                if info_field[pos].split('=')[0] == 'Compounds':
                    if info_field[pos].split('=')[-1] != '-':
                        info_field[pos] = ('Compounds=' +
                                           batch[variant]['Compounds'])
            info_field.append('RankScore=' +
                              str(batch[variant]['Individual_rank_score']))
            batch[variant]['INFO'] = ';'.join(info_field)
            print_line = [batch[variant].get(entry, '-') for entry in header]
            outfile.write('\t'.join(print_line) + '\n')

    def parse(self):
        """Start the parsing"""
        beginning = True
        batch = {}
        new_chrom = None
        current_chrom = None
        current_annotation = set([])
        nr_of_variants = 0
        nr_of_comp_cand = 0
        for variant in self.variant_parser:
            nr_of_variants += 1
            new_annotation = set(variant['vep_info'].keys())
            new_chrom = variant['CHROM']
            if beginning:
                current_annotation = new_annotation
                current_chrom = new_chrom
                beginning = False
                batch[variant['variant_id']] = variant
            else:

                if len(set.intersection(
                       new_annotation, current_annotation)) == 0:
                    score_model.score_variants(batch, self.prefered_models,
                                               self.alt_dict, self.score_dict,
                                               self.value_dict,
                                               self.operation_dict,
                                               self.verbose)
                    self.score_compounds(batch)
                    self.print_variants(batch, self.variant_parser.header,
                                        self.temp_file)

                    current_annotation = new_annotation
                    batch = {variant['variant_id']: variant}

                else:
                    batch[variant['variant_id']] = variant
                    current_annotation = current_annotation.union(new_annotation)

        score_model.score_variants(batch, self.prefered_models,
                                   self.alt_dict, self.score_dict,
                                   self.value_dict, self.operation_dict,
                                   self.verbose)
        self.score_compounds(batch)
        self.print_variants(batch, self.variant_parser.header, self.temp_file)


def main():
    from tempfile import NamedTemporaryFile
    from score_mip_variants import variant_sorter
    from vcf_parser import parser

    parser = argparse.ArgumentParser(description=
                                     "Parse different kind of pedigree files.")
    parser.add_argument('variant_file',
                        type=str, nargs=1,
                        help='A file with variant information.'
                        )
    parser.add_argument('outfile',
                        type=str, nargs=1,
                        help='Specify the path to output.'
                        )
    parser.add_argument('-phased', '--phased',
                        action="store_true",
                        help='If variant file is phased.'
                        )
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        help='Increase output verbosity.'
                        )
    parser.add_argument('-vep', '--vep',
                        action="store_true",
                        help='If variants are annotated with vep.'
                        )

    args = parser.parse_args()
    var_file = args.variant_file[0]

    my_vcf_parser = vcf_parser.VCFParser(var_file)
    temporary_variant_file = NamedTemporaryFile(mode='w')

    my_parser = VariantScorer(my_vcf_parser, temporary_variant_file)
    my_parser.parse()
    temporary_variant_file.seek(0)

    for variant_line in open(temporary_variant_file.name, 'r'):

        print(variant_line.rstrip().split('\t')[7].split(';')[-1].split('=')[-1])

    # outFile=args.outfile[0]
    var_sorter = variant_sorter.FileSort(temporary_variant_file,
                                         outFile=args.outfile[0])
    var_sorter.sort()

    temporary_variant_file.close()

if __name__ == '__main__':
    main()
