from setuptools import setup
import sys,os

with open('mtbdiff/description.txt') as f:
    long_description = f.read()

setup(
    name = 'mtbdiff',
    version = '0.1.0',
    description = 'MTB genome differences',
    long_description = long_description,
    url='https://github.com/dmnfarrell/mtbdiff',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien@gmail.com',
    packages = ['mtbdiff'],
    package_data={'mtbdiff': ['data/*.*',
                  'description.txt']
                 },
    install_requires=['numpy>=1.10',
                      'pandas>=0.24',
                      'biopython>=1.5',
                      'bcbio_gff',
                      'nucdiff',
                      'pyfaidx',
                      'future'],
    entry_points = {
        'console_scripts': [
            'mtbdiff=mtbdiff.app:main']
            },
    classifiers = ['Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.6',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords = ['bioinformatics','biology','genomics','microbiology'],
)
