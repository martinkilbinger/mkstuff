# -*- coding: utf-8 -*-
  
"""INFO

Set the package information for `mkstuff'

:Author: Martin Kilbinger, <martin.kilbinger@cea.fr>

:Date: 18/04/2019
"""

# Release version
version_info = (1, 2, 0)
__version__ = '.'.join(str(c) for c in version_info)

# Package details
__author__ = 'Martin Kilbinger'
__email__ = 'martin.kilbinger@cea.fr'
__year__ = 2019
__url__ = 'https://github.com/martinkilbinger/mkstuff' 
__description__ = 'Useful but random scientific python modules and scripts'

# Dependencies
__requires__ = ['numpy', 'scipy', 'astropy']

# Default package properties
__license__ = 'MIT'
__about__ = ('{} \n\n Author: {} \n Email: {} \n Year: {} \n {} \n\n'
             ''.format(__name__, __author__, __email__, __year__,
                       __description__))
__setup_requires__ = ['pytest-runner', ]
__tests_require__ = ['pytest', 'pytest-cov', 'pytest-pep8']
