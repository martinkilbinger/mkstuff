# -*- coding: utf-8 -*-
  
"""UNIT TESTS FOR MKSTUFF

This module contains unit tests for the mkstuff module.

"""

from unittest import TestCase
import numpy.testing as npt
from ..mkstuff import *


class MathTestCase(TestCase):

    def setUp(self):

        self.a = np.array([4, 2], [2, 4]])

    def tearDown(self):

        self.a = None

    def test_corr_coeff(self):

        npt.assert_almost_equal(mkstuff.corr_coeff(self.a)[0, 1], 0.5,
                         err_msg='Incorrect correlation coefficient')
