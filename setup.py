#!/usr/bin/env python3
"""Setup configuration for PHAGPCR package."""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    requirements = [
        line.strip()
        for line in requirements_file.read_text().splitlines()
        if line.strip() and not line.startswith('#')
    ]
else:
    requirements = [
        'biopython>=1.81',
        'pandas>=2.0.0',
        'primer3-py>=2.0.0',
        'termcolor>=2.0.0',
    ]

setup(
    name="phagpcr",
    version="1.0.0",
    author="Nouri L. Ben Zakour",
    author_email="",  # Add email if desired
    description="qPCR primer design for phage detection in patient serum",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nbenzakour/PHAGPCR",  # Update if needed
    packages=find_packages(),
    py_modules=['phagpcr', 'bz_gene_extractor'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'phagpcr=phagpcr:main',
            'bz-gene-extractor=bz_gene_extractor:main',
        ],
    },
    include_package_data=True,
    package_data={
        '': ['bin/*', 'data/Figma_2023-03-07.png'],
    },
    keywords='qpcr primers phage detection specificity primer3 mfeprimer',
    project_urls={
        'Bug Reports': 'https://github.com/nbenzakour/PHAGPCR/issues',
        'Source': 'https://github.com/nbenzakour/PHAGPCR',
    },
)
