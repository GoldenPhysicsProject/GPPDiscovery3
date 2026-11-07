"""Stdlib falsifiers for principal-series circle intersection and prime phases."""
import unittest
from math import pi,log
from experiments.celestial_nodal_intersections import (
    cayley, latitude, dimension_circle, prime_phase, eta_partial
)
class CelestialNodalTests(unittest.TestCase):
    def test_circle_overlap_all_parameters(self):
        for t in (0,.2,2,14.1347251417347,30):
            s=.5+1j*t;w=cayley(s)
            self.assertAlmostEqual(abs(w),1,places=12)
            self.assertAlmostEqual(latitude(w),0,places=12)
            for d in (1,2,3,4,10):
                self.assertAlmostEqual(abs(dimension_circle(d,s)-w),0,places=12)
    def test_hemisphere_latitude(self):
        for sigma,sign in ((.25,-1),(.5,0),(.75,1)):
            z=latitude(cayley(sigma+14j))
            self.assertEqual((z>0)-(z<0),sign)
    def test_prime_intersection_is_not_first_zero(self):
        for p in (2,3,5,7,11):
            self.assertAlmostEqual(abs(prime_phase(p,0)-1),0,places=13)
        meet=2*pi/log(3/2)
        self.assertLess(abs(prime_phase(2,meet)-prime_phase(3,meet)),1e-13)
        self.assertGreater(abs(meet-14.134725141734693),1)
    def test_nonzero_finite_eta_at_all_aligned_t0(self):
        self.assertGreater(abs(eta_partial(.5,1000)),.1)
if __name__=="__main__":
    unittest.main()
