"""Setup script for fast_wps package."""

from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Read requirements
requirements = []
with open("requirements.txt") as f:
    requirements = [
        line.strip() for line in f if line.strip() and not line.startswith("#")
    ]

setup(
    name="fast_wps",
    version="1.0.0",
    author="JorisVermeeschLab",
    author_email="",  # Add email if available
    description="Fast Window Protection Score calculator for cell-free DNA analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JorisVermeeschLab/docker",  # Update with actual repo URL
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",  # Update if different
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
    },
    entry_points={
        "console_scripts": [
            "optwps=bin.optwps:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="bioinformatics genomics cell-free-DNA cfDNA WPS nucleosome BAM",
    project_urls={
        "Bug Reports": "https://github.com/JorisVermeeschLab/docker/issues",
        "Source": "https://github.com/JorisVermeeschLab/docker",
    },
)
