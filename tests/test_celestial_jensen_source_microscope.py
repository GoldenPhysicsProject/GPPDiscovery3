"""Exact local Cayley geometry + stdlib numerical winding falsifier."""
import unittest
import cmath
import math
from experiments.celestial_jensen_source_microscope import (
    local_cayley, local_inverse, local_latitude, fake_root_coordinate,
    jensen_lower_bound, winding_samples, jensen_samples,
)
class CelestialJensenTests(unittest.TestCase):
    def test_local_halfplane_equator_and_inverse(self):
        for t0 in (0,14.134725,70.2):
            a=.5
            for t in (-8,0,3,30):
                s=.5+1j*t
                self.assertAlmostEqual(abs(local_cayley(s,t0,a)),1,places=13)
            for s in (.75+1j*t0,1+1j*(t0+.25),2+1j*(t0-2)):
                w=local_cayley(s,t0,a)
                self.assertLess(abs(w),1)
                self.assertAlmostEqual(local_inverse(w,t0,a).real,s.real,places=10)
                self.assertAlmostEqual(local_inverse(w,t0,a).imag,s.imag,places=10)
                self.assertGreater(local_latitude(s,t0,a),0)
            self.assertGreater(abs(local_cayley(.25+1j*t0,t0,a)),1)
    def test_exact_fake_height_localization(self):
        w=fake_root_coordinate(.25,.5)
        self.assertAlmostEqual(w,-1/3,places=15)
        for t0 in (0,20,100):
            self.assertAlmostEqual(local_cayley(.75+1j*t0,t0,.5).real,w,places=14)
            self.assertAlmostEqual(local_cayley(.75+1j*t0,t0,.5).imag,0,places=14)
        self.assertAlmostEqual(jensen_lower_bound(.45,abs(w)),math.log(1.35),places=14)
    def test_poly_argument_principle_and_jensen(self):
        r=.45
        n=256
        roots=[r*cmath.exp(2j*math.pi*k/n) for k in range(n)]
        inside=lambda w: w+1/3
        outside=lambda w: w-1.5
        for f,expected in ((inside,1),(outside,0)):
            v=[f(w) for w in roots]+[f(roots[0])]
            winding,maxinc=winding_samples(v)
            self.assertAlmostEqual(winding,expected,places=12)
            self.assertLess(maxinc,.1)
        self.assertAlmostEqual(jensen_samples([inside(w) for w in roots],inside(0)),
                               math.log(1.35),places=12)
        self.assertAlmostEqual(jensen_samples([outside(w) for w in roots],outside(0)),
                               0,places=12)
    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):fake_root_coordinate(.5,.5)
        with self.assertRaises(ValueError):local_inverse(1,0,.5)
        with self.assertRaises(ValueError):jensen_lower_bound(.2,1/3)
if __name__=="__main__":
    unittest.main()
