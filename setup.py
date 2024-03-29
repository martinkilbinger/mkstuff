#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os

__name__ = 'mkstuff'

release_info = {}
infopath = os.path.abspath(os.path.join(os.path.dirname(__file__),
                           __name__, 'info.py'))
with open(infopath) as open_file:
    exec(open_file.read(), release_info)

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name = __name__,
    author = release_info['__author__'],
    author_email = release_info['__email__'],
    version = release_info['__version__'],
    url = release_info['__url__'],
    packages = find_packages(),
    install_requires = release_info['__requires__'],
    license = release_info['__license__'],
    description = release_info['__about__'],
    #long_description = long_description,
    long_description = release_info['__about__'],
    long_description_content_type = 'text/rst',
    setup_requires = release_info['__setup_requires__'],
    tests_require = release_info['__tests_require__'],
    #platforms = ['posix', 'mac os'],
    #classifiers= [
                  #'Programming Language :: Python',
                  #'Natural Language :: English',
                 #],
    scripts = ['bin/{}'.format(fn) for fn in ['fits2ascii.py', 'ascii2fits.py']],
)
