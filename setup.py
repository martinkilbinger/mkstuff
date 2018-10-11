#!/usr/bin/env python
# -*- coding: utf-8 -*-


#from setuptools import setup

from distutils.core import setup

exec(open('mkstuff/info.py').read())


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name             = __whoami__,
      author           = __author__,
      packages         = [__whoami__],
      version          = __version__,
      description      = '{}\'s stuff'.format(__author__),
      long_description = readme(),
      url              = 'https://github.com/martinkilbinger/mkstuff',
      author_email     = __email__,
      platforms        = ['posix', 'mac os'],
      license          = "GNU GPLv3",
      classifiers      = [
                          'Programming Language :: Python',
                          'Natural Language :: English',
                         ],
      scripts          = ['bin/fits2ascii.py'],
)


