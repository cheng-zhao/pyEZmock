#!/usr/bin/env python3
import numpy as np
from scipy.special import hyp2f1
from scipy.misc import derivative

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
