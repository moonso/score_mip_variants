#!/usr/bin/env python
# encoding: utf-8
"""
score_variant.py

Script that takes a variant as input and modify it with a score depending on its different values.

Possible names for the list of genetic models are:

AD, AD_denovo, AR, AR_denovo, AR_compound, X, X_denovo


Created by Måns Magnusson on 2013-08-14.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os

from pprint import pprint as pp

from Mip_Family_Analysis.Utils import is_number


consequence_severity = {}
# This is the rank scores for the different consequences that VEP uses:
consequence_severity['transcript_ablation'] = 5
consequence_severity['splice_donor_variant'] = 5
consequence_severity['splice_acceptor_variant'] = 5
consequence_severity['stop_gained'] = 5
consequence_severity['frameshift_variant'] = 5
consequence_severity['stop_lost'] = 5
consequence_severity['initiator_codon_variant'] = 5
consequence_severity['inframe_insertion'] = 3
consequence_severity['inframe_deletion'] = 3
consequence_severity['missense_variant'] = 3
consequence_severity['transcript_amplification'] = 3
consequence_severity['splice_region_variant'] = 3
consequence_severity['incomplete_terminal_codon_variant'] = 3
consequence_severity['synonymous_variant'] = 1
consequence_severity['stop_retained_variant'] = 1
consequence_severity['coding_sequence_variant'] = 1
consequence_severity['mature_miRNA_variant'] = 1
consequence_severity['5_prime_UTR_variant'] = 1
consequence_severity['3_prime_UTR_variant'] = 1
consequence_severity['non_coding_exon_variant'] = 1
consequence_severity['nc_transcript_variant'] = 1
consequence_severity['intron_variant'] = 1
consequence_severity['NMD_transcript_variant'] = 1
consequence_severity['upstream_gene_variant'] = 1
consequence_severity['downstream_gene_variant'] = 1
consequence_severity['TFBS_ablation'] = 1
consequence_severity['TFBS_amplification'] = 1
consequence_severity['TF_binding_site_variant'] = 1
consequence_severity['regulatory_region_variant'] = 1
consequence_severity['regulatory_region_ablation'] = 1
consequence_severity['regulatory_region_amplification'] = 1
consequence_severity['feature_elongation'] = 1
consequence_severity['feature_truncation'] = 1
consequence_severity['intergenic_variant'] = 0


def get_genetic_models(model_dict):
    """return a list with the genetic models followed"""
    models_followed = []
    for model in model_dict:
        if model_dict[model]:
            models_followed.append(model)
    return models_followed

def score_variant(variants, prefered_models = []):
    """Score a variant object according to Henriks score model. Input: A variant object and a list of genetic models."""
    
    if  prefered_models == ['NA']:
        prefered_models = []
    
    for variant_id in variants:
        variant = variants[variant_id]
        score = 0
        # Models of inheritance
        variant_models = get_genetic_models(variant.get('Inheritance_model', {}))
        
        # Predictors
        mutation_taster = variant.get('Mutation_taster', None)
        avsift = variant.get('SIFT', None)
        
        if 'Poly_phen_hdiv' in variant:
            poly_phen = variant.get('Poly_phen_hdiv', None)
        else:
            poly_phen = variant.get('Poly_phen', None)
        
        # Annotations:
        functional_annotation = variant.get('Functional_annotation', None)
        if functional_annotation:
            try:
                functional_annotation = {gene_info.split(':')[0]:gene_info.split(':')[1] for gene_info in functional_annotation.split(',')}
            except IndexError:
                functional_annotation = None
        
        # Frequency in databases:
        thousand_genomes_frequency = variant.get('1000G', None)
        dbsnp_frequency = variant.get('Dbsnp129', None)
        dbsnp_id = variant.get('Dbsnp_nonflagged', None)
        hbvdb = variant.get('HBVDB', None)
        
        # Filter
        
        filt = variant.get('GT_call_filter', None)
        
        # Conservation scores:
            # Base
        gerp_base = variant.get('GERP', None)
            # Region
        mce64way = variant.get('Phast_cons_lements', None)
        gerp_region = variant.get('GERP_elements', None)
            
            
        phylop = variant.get('Phylo_p', None)
        
        segdup = variant.get('Genomic_super_dups', None)
        
        hgmd = variant.get('HGMD', None)
        
        
        
        score += check_inheritance(variant_models, prefered_models)
        score += check_predictions(mutation_taster, avsift, poly_phen)
        score += check_functional_annotation(functional_annotation)
        score += check_frequency_score(thousand_genomes_frequency, dbsnp_frequency, hbvdb, dbsnp_id)
        score += check_filter(filt)
        score += check_region_conservation(mce64way, gerp_region)
        score += check_base_conservation(gerp_base)
        score += check_phylop_score(phylop)
        score += check_segmental_duplication(segdup)
        score += check_hgmd(hgmd)
        variant['Individual_rank_score'] = score
        
    return
    
def check_inheritance(variant_models, prefered_models):
    """Check if the models of inheritance are followed for the variant."""
    model_score = 0
    #If any of the prefered models are followed:
    for model_followed in variant_models:
        if model_followed in prefered_models:
            model_score = 3
    #Else if any model is followed
    if model_score != 3:
        if len(variant_models) > 0:
            model_score = 1
        else:
            model_score = -12
    return model_score
    
def check_predictions(mutation_taster = None, avsift = None, poly_phen = None):
    """Score the variant based on the scores from prediction databases."""
    prediction_score = 0
    if is_number.is_number(avsift):
        if float(avsift) <= 0.05:
            prediction_score += 1
    if is_number.is_number(mutation_taster):
        if float(mutation_taster) >= 0.05:
            prediction_score += 1
    if is_number.is_number(poly_phen):
        if float(poly_phen) >= 0.85:
            prediction_score += 1
    return prediction_score
    
def check_functional_annotation(functional_annotation = None):
    """Score the variant based on its functional annotation"""
    functional_annotation_score = 0
    if functional_annotation:
        for gene in functional_annotation:
            score = consequence_severity.get(functional_annotation[gene],0)
            if score > functional_annotation_score:
                functional_annotation_score = score
    return functional_annotation_score
    
def check_frequency_score(thousand_genomes_frequency = None, dbsnp_frequency = None, hbvdb_frequency = None, dbsnp_id = None):
    """Score the variant based on the frequency in population."""

    frequency_score = 0
    freq_scores = []
    
    def get_freq_score(frequency):
        """Returns a score depending on the frequency"""
        if is_number.is_number(frequency):
            if float(frequency) <= 0.005:
                return 2
            elif float(frequency) <= 0.02:
                return 1
            #If common variant:
            else:
                    return -12
        else:# If not existing in database
            return 3
    
    freq_scores.append(get_freq_score(thousand_genomes_frequency))
    freq_scores.append(get_freq_score(dbsnp_frequency))
    freq_scores.append(get_freq_score(hbvdb_frequency))
    common = False
    # If the variant if common in any database(if a score is negative) we give a low score:
    for freq_score in freq_scores:
        if freq_score < 0:
            common = True
    if common:
        frequency_score = -12
    else:
        frequency_score += sum(freq_scores) / 2
    # If variant has no ID in dbSNP it get an extra score
        if dbsnp_id == '-':
            frequency_score += 1
    return frequency_score
    
def check_filter(filt):
    """Check if variant has passed the filter process."""
    filter_score = 0
    if filt == 'PASS':
        filter_score = 3
    elif filt == 'PRES':
        filter_score = 1
    return filter_score
    
def check_region_conservation(mce64way = None, gerp_region = None):
    """Score the variant based on what annotations it has for the region conservations"""
    region_conservation_score = 0
    if mce64way != '-' and gerp_region != '-':
        region_conservation_score += 2
    elif mce64way != '-' or gerp_region != '-':
        region_conservation_score += 1
    return region_conservation_score
    
def check_base_conservation(gerp_base_score = None):
    """Score the variant based on the base level conservation."""
    base_conservation_score = 0
    if is_number.is_number(gerp_base_score):
        if float(gerp_base_score) >= 4:
            base_conservation_score += 2
        elif float(gerp_base_score) >= 2:
            base_conservation_score += 1
    return base_conservation_score
    
def check_phylop_score(phylop = None):
    """Score the variant based on the Phylop score."""
    phylop_score = 0
    if is_number.is_number(phylop):
        if float(phylop) >= 0.9984188612:
            phylop_score += 2
        elif float(phylop) >= 0.95:
            phylop_score += 1
    return phylop_score
    
def check_segmental_duplication(segdup):
    """Check if there are any annotations for segmental duplication"""
    segdup_score = 0
    if segdup != '-':
        segdup_score -= 2
    return segdup_score
    
def check_hgmd(hgmd):
    """Check if the variant have any annotation from hgmd"""
    hgmd_score = 0
    if hgmd != '-':
        hgmd_score += 1
    return hgmd_score
    
    



def main():
    pass


if __name__ == '__main__':
    main()

