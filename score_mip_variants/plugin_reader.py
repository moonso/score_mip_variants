#!/usr/bin/env python
# coding: utf-8

import sys
import configparser
from collections import defaultdict


# Import third party library
# https://github.com/mitsuhiko/logbook
from logbook import Logger, StderrHandler
log = Logger('Logbook')
log_handler = StderrHandler()


def read_config(config_file):
    """Read the config file.

    Args:
        config_file (file) : Plugin config file

    Returns:
        object: config object
    """
    config = configparser.ConfigParser()

    config.read(config_file)

    return config


def check_vcf_config(key, not_info_list, my_vcf_parser,
                     config):
    """Check that the alternatives in config is present
    within vcf INFO keys.

    Args:
        key            (string) : Plugin config alternative key
        not_info_list  (list) : List of none vcf INFO keys
        my_vcf_parser  (object) : vcf_parser object
        config         (object) : Config object

    Returns:
        None:
    """
    if key in not_info_list:

        pass
    elif key not in my_vcf_parser.metadata.info_dict:

        log.notice('Plug-in.' +
                   config['Plug-in']['name'] +
                   ': Could not find key="' + key +
                   '" in INFO fields of variant file'
                   + '\n'
                   )
        log.notice('Aborting ranking' + '\n')
        sys.exit()


def collectKeys(config_file, my_vcf_parser, verbose):
    """Collects database keys determined by plug-in and saves them in dict.

    Args:
        config_file (file) : plugin config file
        my_vcf_parser (object) : vcf_parser object

    Returns:
        dict: dict[collection][list_of_keys]
    """
    ## Create config object
    config = read_config(config_file)

    ## Create dictionnaries to store config info
    alt_dict = dict()
    score_dict = dict()
    value_dict = dict()
    operation_dict = dict()
    ## Alternatives not found in vcf INFO
    not_info_list = ["FILTER"]

    ## Collect keys for plug-in
    for section in config.items():

        key = section[0]  # Alias
        if verbose:
            log.info("Config alternative: " + str(config[key]))
            if 'version' in config[key]:
                log.info("Plugin version: " + str(config[key]['version']))

        ## Only vcf records
        if 'category' in config[section[0]]:

            ## Test for presence of vcf INFO keys
            check_vcf_config(key, not_info_list, my_vcf_parser, config)

            alt_dict[key] = dict()
            score_dict[key] = dict()
            value_dict[key] = dict()
            operation_dict[key] = dict()

            ## Save vcf_info_key and general parameters
            alt_dict[key]['category'] = config[key]['category']
            alt_dict[key]['data_type'] = config[key]['data_type']

            if 'category_aggregate' in config[section[0]]:

                alt_dict[key]['category_aggregate'] = config[key]['category_aggregate']
            if 'field_separators' in config[section[0]]:

                alt_dict[key]['field_separators'] = config[key]['field_separators']
            if 'record_aggregate' in config[section[0]]:

                alt_dict[key]['record_aggregate'] = config[key]['record_aggregate']

            for alternative in config[section[0]]:

                order = alternative.split('-')

                for term in order:

                    if term == 'score':

                        score_key = order[1]  # Alias
                        score_dict[key][score_key] = config[key][alternative]
                    if term == 'value':

                        value_key = order[1]  # Alias
                        ## Split to seperate operation from number
                        operation_list = config[key][alternative].split(':')

                        ## Check that a operation was given
                        if len(operation_list) > 1:

                            operation_dict[key][value_key] = operation_list[0]
                            value_dict[key][value_key] = operation_list[1]
                        else:

                            value_dict[key][value_key] = operation_list[0]
    # Return dicts with vcf_info_key as 1st key and alternative as 2nd key
    return alt_dict, score_dict, value_dict, operation_dict
