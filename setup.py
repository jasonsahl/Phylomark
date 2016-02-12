#!/usr/bin/env python
from distutils.core import setup

__author__ = "Jason Sahl"
__credits__ = ["Jason Sahl"]
__license__ = "GPL v3"
__version__ = "1.0"
__maintainer__ = "Jason Sahl"
__email__ = "jasonsahl@gmail.com"
__status__ = "Development"
 
long_description = """Phylomark - a method to find
representative markers from whole genome sequence data
"""

setup(name='Phylomark',
      version=__version__,
      description='Phylomark',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      packages=['phylomark'],
      long_description=long_description
)
