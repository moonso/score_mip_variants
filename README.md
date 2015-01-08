#Score Mip Variants#

Annotates variants in VCF-format with a rank score.

Score_mip_variants uses the weighted sum model (WSM) approach to score the most likely pathogenic variant.

Generally, the higher value the more likely pathogenic variant.

Score_mip_variants uses config files to define the rank model, which enables customized set-up and versioning of rank models.


###Usage###

	score_mip_variants examples/my_test.vcf examples/trio_ped.txt --plugin_file test.ini -o out.vcf


###ConfigFile###

The config file is in ini-format. INI files are simple text files with a basic structure composed of sections, properties, and values.

The sections in the ini-file are named as the VCF keys (records) that should be included in the rank model.

Each section has defined properties writen as name=value:

	data_type = data record type ('float', 'String';mandatory)

	record_aggregate = data record function ('max', 'min', 'sum';mandatory). Defined how the score from the data record should be calculated.

	category = data record category. Defined the category for this record e.g., 1000GAF and EXACAF can be grouped within the same category i.e., "allele_frequencies", since they are both allele frequencies. (mandatory)

	category_aggregate = data record category function ('max', 'min', 'sum';mandatory). Defined how each score from category data record(s) should be calculated.

	value_X = Criteria for data record for condition or string comparison. FORMAT: "Condition:Value"|"String".
	
	score_X = Performance score for X (e.g., "rare", "common") if record value fulfills criteria. FORMAT: "Integer"

	Note: "_X" is linked between value and score
	
	field_separators = data record field separator. FORMAT: "separator1_separator2_...separatorN"
	
	For Example:

	config['1000GMAF'] = {'data_type': 'float',
        	              'category': 'allele_frequency',
                	      'category_aggregate': 'min',
                	      'record_aggregate': 'max',
                	      'value-notreported': 'na:na',
                	      'value-rare': 'le:0.005',
                	      'value-intermediate': 'le:0.02',
                	      'value-common': 'gt:0.02',
                	      'score-notreported': '3',
                	      'score-rare': '2',
                	      'score-intermediate': '1',
                	      'score-common': '-12',
                	     }
