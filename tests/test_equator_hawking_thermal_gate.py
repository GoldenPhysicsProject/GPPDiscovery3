"""Thermal obstruction for natural free angular modes at Hawking temperature."""
import unittest
from math import log
from experiments.equator_hawking_thermal_gate import (
    thermal_entropy,natural_uniform_bound,h_fermion
)
class HawkingGateTests(unittest.TestCase):
    def test_fermion_entropy_zero_mode(self):
        self.assertAlmostEqual(h_fermion(0),log(2),places=13)
        self.assertGreater(h_fermion(1),0)
        self.assertLess(h_fermion(10),h_fermion(1))
    def test_natural_Hawking_entropy_uniformly_bounded(self):
        for N in (1,2,3,8,16,64,128):
            s=thermal_entropy(N)
            self.assertLessEqual(s,natural_uniform_bound()+1e-12)
        self.assertLess(natural_uniform_bound(),0.694)
    def test_extra_soft_dispersion_changes_scaling(self):
        for N in (16,32,64,128):
            s=thermal_entropy(N,soft=True)
            self.assertGreater(s,0.025*N*N)
            self.assertLess(s,0.10*N*N)
if __name__=="__main__":
    unittest.main()
