#!/usr/bin/env python3
# cosmoEZ.py: this file is part of the pyEZmock package.
#
# pyEZmock: Python wrapper of Effective Zel'dovich approximation mock (EZmock).
#
# Github repository:
#       https://github.com/cheng-zhao/pyEZmock
#
# Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com>
#   
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import numpy as np
from scipy.special import hyp2f1
from scipy.misc import derivative

# Structure growth for flat LambdaCDM cosmology.
# ref: https://arxiv.org/abs/0906.1643
class flatLCDM:
  def __init__(self, omega_m=0.31):
    self.Om = float(omega_m)

  def delta(self, a):
    return a * hyp2f1(1./3., 1, 11./6., a**3 * (1-1./self.Om))

  def growth2z0(self, a):
    return (self.delta(a) / self.delta(1))**2

  def growthf(self, a):
    derv = derivative(self.delta, a, dx=1e-3)
    return a * derv / self.delta(a)

  def hubble(self, a):
    return 100 * np.sqrt(1 - self.Om + self.Om / a**3)
