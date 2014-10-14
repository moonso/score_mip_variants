try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
    
# For making things look nice on pypi:
# try:
#     import pypandoc
#     long_description = pypandoc.convert('README.md', 'rst')
# except (IOError, ImportError):
long_description = 'Add rank scores to variants in Variant Call Format (VCF) files according to mip.'


# with open('README.txt') as file:
#     long_description = file.read()

setup(name='score_mip_variants',
    version='0.5.5',
    description='Annotate vcf variants with a rank score',
    author = 'Mans Magnusson',
    author_email = 'mans.magnusson@scilifelab.se',
    url = 'http://github.com/moonso/rank_mip_variants',
    license = 'MIT License',
    install_requires=[
        'vcf_parser', 
        'ped_parser', 
        'click',
        'configparser',
        'logbook',
    ],
    packages = [
        'score_mip_variants',
        'score_mip_variants/configs'
    ],
    package_data = {
        'score_mip_variants': ['configs/*.ini']
    },
    scripts = [
        'scripts/score_mip_variants'
    ],
    keywords = ['ranking','scoring', 'vcf', 'variants'],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 4 - Beta",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    long_description = long_description,
)