# File: setup.py
# Date: 1-Oct-2019
#
# Updates:
#
#
import re

from setuptools import find_packages
from setuptools import setup

packages = []
thisPackage = "rcsb.utils.chem"

with open("rcsb/utils/chem/__init__.py", "r") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")

setup(
    name=thisPackage,
    version=version,
    description="RCSB Python Chemical Utility Classes",
    long_description="See:  README.md",
    author="John Westbrook",
    author_email="john.westbrook@rcsb.org",
    url="https://github.com/rcsb/py-rcsb_utils_chem",
    #
    license="Apache 2.0",
    classifiers=(
        "Development Status :: 3 - Alpha",
        # 'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ),
    entry_points={"console_scripts": ["cactvs_annotate_mol=rcsb.utils.chem.cactvsAnnotateMol:main"]},
    #  The following is somewhat flakey --
    dependency_links=["https://pypi.anaconda.org/OpenEye/simple#egg=OpenEye-toolkits-2019.10.2"],
    install_requires=["mmcif >= 0.57", "rcsb.utils.io >= 0.68", "rcsb.utils.multiproc >= 0.17", "OpenEye-toolkits==2019.10.2"],
    packages=find_packages(exclude=["rcsb.mock-data", "rcsb.utils.tests-chem", "rcsb.utils.tests-*", "tests.*"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.cfg"]
    },
    #
    test_suite="rcsb.utils.tests-chem",
    tests_require=["tox"],
    #
    # Not configured ...
    extras_require={"dev": ["check-manifest"], "test": ["coverage"]},
    # Added for
    command_options={"build_sphinx": {"project": ("setup.py", thisPackage), "version": ("setup.py", version), "release": ("setup.py", version)}},
    # This setting for namespace package support -
    zip_safe=False,
)
