#!/usr/bin/env python
# encoding: utf-8
"""
score_model.py

Script that takes a batch of variants as input and modify it with a score depending on its different values.

Possible names for the list of genetic models are:

AD, AD_denovo, AR, AR_denovo, AR_compound, X, X_denovo


Created by Måns Magnusson on 2013-08-14.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function
from __future__ import unicode_literals


import sys
import os

from pprint import pprint as pp


# from Mip_Family_Analysis.Utils import is_number
#Frågor till Henrik:
# CP = phastcons? eller ska vi köra någon av dbNSFP_phastCons100way_vertebrate, dbNSFP_phastCons46way_primate
# DB = dbsnp_noflag?
# Tidigare körde vi med Gerp base och Gerp region. Nu ser jag bara dbNSFP_GERP++_RS
# Vilken av dbNSFP_phyloP100way_vertebrate och dbNSFP_phyloP46way_primate ska vi köra med?
# Ibland 1000GMAF ibland 1000G_freq???


consequence_severity = {}
# This is the rank scores for the different consequences that VEP uses:
consequence_severity['transcript_ablation'] = 5
consequence_severity['splice_donor_variant'] = 4
consequence_severity['splice_acceptor_variant'] = 4
consequence_severity['stop_gained'] = 4
consequence_severity['frameshift_variant'] = 4
consequence_severity['stop_lost'] = 4
consequence_severity['initiator_codon_variant'] = 4
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


def is_number(number):
    """Returns true if the string is a number or False otherwise"""
    if type(number) == type(1) or type(number) == type(0.1) or type(number) == type('') or type(u''):
        try:
            float(number)
            return True
        except ValueError:
            return False
        except TypeError:
            return False
    else:
        return False


def score_variants(batch, prefered_models = []):
    """Score a variant object according to Henriks score model. Input: A variant object and a list of genetic models."""
    
    if  prefered_models == ['NA']:
        prefered_models = []
    
    for variant_id in batch:
        variant = batch[variant_id]
        score = 0
        info_dict = variant.get('info_dict', {})
        
        # Models of inheritance are annotated as GM=AR:AR_comp...
        
        variant_models = info_dict.get('GeneticModels','NA').split(',')
        # Predictors
        
        mutation_taster = info_dict.get('dbNSFP_MutationTaster_score', None)
        # print('dbNSFP_MutationTaster_score: %s' % mutation_taster)
        
        avsift = info_dict.get('Sift', None)
        if avsift:
            avsift = min([float(gene_score.split(':')[-1]) for gene_score in avsift.split(',')])
            # print('First SIFT: %s' % avsift)
        else:
            avsift = info_dict.get('dbNSFP_SIFT_score', None)
        # print('Second SIFT: %s' % avsift)
        
        poly_phen = info_dict.get('dbNSFP_Polyphen2_HVAR_score', None)
        if poly_phen:
            poly_phen = max([float(poly_score) for poly_score in poly_phen.split(',') if is_number(poly_score)])
        # print('dbNSFP_Polyphen2_HVAR_score: %s' % poly_phen)
        
        
        # Annotations:
        most_severe = info_dict.get('MostSevereConsequence', None)
        highest_consequence_score = 0
        for gene_anno in most_severe.split(','):
            new_score = consequence_severity.get(gene_anno.split(':')[-1], 0)
            if new_score > highest_consequence_score:
                highest_consequence_score = new_score
        
        #Ranges between 0-5. MAX=5
        score += highest_consequence_score
        
        # Frequency in databases:
        thousand_genomes_frequency = info_dict.get('1000GMAF', None)
        # print('Thousand g freq: %s' % thousand_genomes_frequency)
        dbsnp_frequency = info_dict.get('DbsnpMAF', None)
        # print('dbSNP freq: %s' % dbsnp_frequency)
        dbsnp129_frequency = info_dict.get('Dbsnp129MAF', None)
        # print('dbSNP 129 freq: %s' % dbsnp129_frequency)
        dbsnp_id = variant.get('ID', None)
        # print('dbSNP ID: %s' % dbsnp_id)
        esp_frequency = info_dict.get('ESPMAF', None)
        # print('ESP freq: %s' % esp_frequency)
        hbvdb = info_dict.get('BVDMAF', None)
        # print('HBVDB freq: %s' % hbvdb)
        
        # Filter
        
        filt = variant.get('FILTER', 'NOTPASS').strip()
        
        # Conservation scores:
        
        # Base
        gerp_base = info_dict.get('dbNSFP_GERP++_RS', None)
        
        # Region
        mce64way = info_dict.get('dbNSFP_phastCons46way_primate', None)
        # print('dbNSFP_phastCons46way_primate: %s' % mce64way)
        
        gerp_region = info_dict.get('dbNSFP_GERP++_NR', None)
        #Hur ska vi göra med GERP region!!!???
        phylop = info_dict.get('dbNSFP_phyloP46way_primate', None)
        
        # segdup = variant.get('Genomic_super_dups', None)
        
        # hgmd = variant.get('HGMD', None)
        
        # print('Score after most severe consequence %s' % score)
        # print('Max score after this step 5\n')
        # 
        # print('Score before inheritance: %s' % score)
        
        # Inheritance ranges between -12 - 3, MAX=3
        score += check_inheritance(variant_models, prefered_models)
        # print('Score after inheritance: %s' % score)
        # print('Max score after this step 8\n')
        
        # Predictions ranges between 0 - 3, MAX=3
        score += check_predictions(mutation_taster, avsift, poly_phen)
        # print('Predictors: %s, %s, %s' % (mutation_taster, avsift, poly_phen))
        # print('Score after predictions: %s' % score)
        # print('Max score after this step 11\n')
        # score += check_functional_annotation(functional_annotation)
        
        # Frequenciy scores ranges between -12 - 5, MAX=5
        score += check_frequency_score(thousand_genomes_frequency, dbsnp_frequency, hbvdb, esp_frequency, dbsnp_id)
        # print('Frequencies: 1000g= %s, dbsnp129= %s, hbvdb= %s, ESP= %s' % (thousand_genomes_frequency, dbsnp_frequency, hbvdb, esp_frequency))
        # print('Score after frequency: %s' % score)
        # print('Max score after this step 16\n')
        
        # Filter ranges between 0 - 3, MAX=3
        score += check_filter(filt)
        # print('Filter = %s' % filt)
        # print('Score after filter: %s' % score)
        # print('Max score after this step 19\n')
        
        # Region conservation ranges between 0 - 2, MAX=2
        score += check_region_conservation(mce64way, gerp_region)
        # print('Region conservation: MCE64= %s, GERP= %s' % (mce64way, gerp_region))
        # print('Score after region conservation: %s' % score)
        # print('Max score after this step 21\n')
        
        # Ranges between 0 - 2, MAX=2
        score += check_base_conservation(gerp_base)
        # print('Base conservation: GERP= %s' % gerp_base)
        # print('Score after base conservation: %s' % score)
        # print('Max score after this step 23\n')
        
        # Phylop ranges between 0 - 2, MAX=2
        score += check_phylop_score(phylop)
        # print('Phylop = %s' % phylop)
        # print('Score after phylop conservation: %s' % score)
        # print('Max score after this step 25\n')
        
        # Ranges between -2 - 0, MAX=0
        # score += check_segmental_duplication(segdup)
        # Ranges between -0 - 1, MAX=1
        # score += check_hgmd(hgmd)
        variant['Individual_rank_score'] = score
        
    return
    
def check_inheritance(variant_models, prefered_models):
    """Check if the models of inheritance are followed for the variant."""
    #If any of the prefered models are followed:
    for model_followed in variant_models:
        if model_followed in prefered_models:
            return 3
    #Else if any model is followed
        elif model_followed == 'NA':
            return -12
    return 1
    
def check_predictions(mutation_taster = None, avsift = None, poly_phen = None):
    """Score the variant based on the scores from prediction databases."""
    prediction_score = 0
    if is_number(avsift):
        if float(avsift) <= 0.05:
            prediction_score += 1
    if is_number(mutation_taster):
        if float(mutation_taster) >= 0.05:
            prediction_score += 1
    if is_number(poly_phen):
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
    
def check_frequency_score(thousand_genomes_frequency = None, dbsnp_frequency = None, hbvdb_frequency = None, 
                            esp_frequency = None, dbsnp_id = None):
    """Score the variant based on the frequency in population."""

    frequency_score = 0
    freq_scores = []
    
    def get_freq_score(frequency):
        """Returns a score depending on the frequency"""
        if is_number(frequency):
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
    freq_scores.append(get_freq_score(esp_frequency))
    common = False
    # If the variant if common in any database(if a score is negative) we give a low score:
    for freq_score in freq_scores:
        if freq_score < 0:
            common = True
    if common:
        frequency_score = -12
    else:
        frequency_score += sum(freq_scores) // 3
    # If variant has no ID in dbSNP it get an extra score
        if dbsnp_id in ['-','.']:
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
    if mce64way and gerp_region:
        region_conservation_score += 2
    elif mce64way or gerp_region:
        region_conservation_score += 1
    return region_conservation_score
    
def check_base_conservation(gerp_base_score = None):
    """Score the variant based on the base level conservation."""
    base_conservation_score = 0
    if is_number(gerp_base_score):
        if float(gerp_base_score) >= 4:
            base_conservation_score += 2
        elif float(gerp_base_score) >= 2:
            base_conservation_score += 1
    return base_conservation_score
    
def check_phylop_score(phylop = None):
    """Score the variant based on the Phylop score."""
    phylop_score = 0
    if is_number(phylop):
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

