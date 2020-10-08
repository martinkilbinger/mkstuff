# -*- coding: utf-8 -*-
  
"""UNIT TESTS FOR MKSTUFF

This module contains unit tests for the mkstuff module.

"""

from unittest import TestCase
import numpy.testing as npt
from ..mkstuff import *


class MathTestCase(TestCase):

    def setUp(self):

        self.a = np.array([[4, 2], [2, 4]])
        self.vector = np.ones((3))
        self.tensor = np.ones((3,2,4))
        self.zero = np.ones((4,4))
        self.zero[2,2] = 0

    def tearDown(self):

        self.a = None

    def test_corr_coeff(self):

        npt.assert_almost_equal(corr_coeff(self.a)[0, 1], 0.5,
                                err_msg='Incorrect correlation coefficient')

        npt.assert_raises(ValueError, corr_coeff, self.tensor)
        npt.assert_raises(ValueError, corr_coeff, self.vector)

        npt.assert_raises(ZeroDivisionError, corr_coeff, self.zero)
