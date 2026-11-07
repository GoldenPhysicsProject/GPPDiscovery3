"""No third-party dependencies: run with python -m unittest discover -s tests."""

import json
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from rh_lab.catalog import ROOT, load, validate
from rh_lab.controls import on_line, zero_residual
from rh_lab.runner import run
from experiments.arch_cross_superexp_null import exact_upper
from experiments.arch_ftheta_shift_gap import theta_mode_rayleigh_lower
from experiments.arch_annihilator_full_s1_probe import (
    DELTA, A0, annihilator, gaussian_q_binomial_edge,
    continuum_sewing_prime, pole_parts, prime_parts,
    prime_power_atoms, sine2_autocorrelation, folded_cross_arch,
)
from experiments.translated_arithmetic_sewing_scan import (
    convolution_coefficients, translated_prime_cross, translated_cross_arch,
)
from experiments.gamma_offsupport_laplace import (
    rational_tail_cap, laplace_offsupport,
)
from experiments.gamma_uniform_operator_bound import cross_trace_tail_cap
import math


class LabTests(unittest.TestCase):
    def test_offline_quartet_exact(self):
        self.assertEqual(zero_residual(), (0, 0))
        a, b = Fraction(1, 4), Fraction(2)
        for t in [Fraction(i, 10) for i in range(-250, 251)]:
            self.assertGreaterEqual(on_line(t, a, b), 4*a*a*b*b)

    def test_superexp_exact_witness(self):
        self.assertLess(exact_upper(6), Fraction(1, 10**18))
        self.assertLess(exact_upper(10), Fraction(1, 10**52))
        self.assertLess(exact_upper(10), exact_upper(6))

    def test_ftheta_shift_breaks_cancellation(self):
        self.assertGreater(theta_mode_rayleigh_lower(6), 10**14 * exact_upper(6))
        self.assertGreater(theta_mode_rayleigh_lower(10), 10**50 * exact_upper(10))

    def test_q_binomial_annihilator_and_zeta_pole_identity(self):
        self.assertAlmostEqual(sine2_autocorrelation(0), 3/8, places=15)
        self.assertAlmostEqual(sine2_autocorrelation(.5), 1/16, places=15)
        self.assertEqual(sine2_autocorrelation(1), 0)
        for k in (2, 3, 6, 10):
            roots, coeff = annihilator(k)
            edge = gaussian_q_binomial_edge(k)
            self.assertLess(max(abs(a-b) for a,b in zip(reversed(coeff),edge)), 4e-15)
            self.assertLess(abs(pole_parts(roots,0)["plus"]-
                                continuum_sewing_prime(roots,0)), 1e-12)
            self.assertLess(abs(pole_parts(roots,0)["plus"]+
                                pole_parts(roots,0)["minus"]), 1e-12)

    def test_prime_autocorrelation_and_squared_polynomial_grouping(self):
        k=3
        _, coeff = annihilator(k)
        for theta in (0, .25):
            prime = prime_parts(coeff,theta)
            logs,weights=prime_power_atoms(prime["prime_cutoff"])
            diag_terms,cross_terms=[],[]
            for t,w in zip(logs,weights):
                c=w if theta==0 else w*(2*math.cosh(theta*t))
                a=math.fsum(coeff[j]*coeff[h]*sine2_autocorrelation(
                    t+DELTA*(h-j)) for j in range(k) for h in range(k))
                b=math.fsum(coeff[j]*coeff[h]*sine2_autocorrelation(
                    t-(2*A0+DELTA*(j+h)+1))
                    for j in range(k) for h in range(k))
                diag_terms.append(2*c*a)
                cross_terms.append(c*b)
            self.assertAlmostEqual(prime["2D"],math.fsum(diag_terms),places=11)
            self.assertAlmostEqual(prime["X"],math.fsum(cross_terms),places=11)

    def test_translated_prime_and_arch_source_grouping(self):
        k=3
        roots,coeff=annihilator(k)
        conv=convolution_coefficients(coeff)
        for a in (0., .2, .6):
            limit=math.ceil(math.exp(2*(A0+a+DELTA*(k-1)+1)))
            logs,weights=prime_power_atoms(limit)
            grouped=translated_prime_cross(logs,weights,conv,a)
            brute=math.fsum(w*math.fsum(coeff[j]*coeff[h]*sine2_autocorrelation(
                t-2*(A0+a)-DELTA*(j+h)-1)
                for j in range(k) for h in range(k))
                for t,w in zip(logs,weights))
            self.assertAlmostEqual(grouped,brute,places=11)
            self.assertLessEqual(translated_cross_arch(roots,a),0.)
        self.assertAlmostEqual(translated_cross_arch(roots,0.),
                               folded_cross_arch(roots,0.),places=14)
        self.assertAlmostEqual(translated_prime_cross(
            *prime_power_atoms(math.ceil(math.exp(2*(A0+DELTA*(k-1)+1)))),
            conv,0.),prime_parts(coeff,0.)["X"],places=11)

    def test_offsupport_gamma_mode_tail_certificate(self):
        # Source kernel supports are length-one; d>=2 is the off-support zone.
        for theta in (Fraction(0),Fraction(1,4)):
            for d in range(2,6):
                cap=rational_tail_cap(d,theta,n_modes=60,taylor_degree=36)
                self.assertGreater(cap,0)
                self.assertLess(cap,Fraction(2,10**29))
                self.assertLess(laplace_offsupport(d,float(theta)),0.)
        with self.assertRaises(ValueError):
            laplace_offsupport(1,0.)

    def test_all_k_pure_arch_cross_trace_rational_tail(self):
        # The pure gamma Hankel block is a Gram operator. Its infinite
        # trace is bounded by a convergent geometric expansion, with
        # an EXACT rational cap on the omitted mode families.
        for theta in (Fraction(0), Fraction(1,4)):
            cap60 = cross_trace_tail_cap(theta, n_modes=60)
            cap80 = cross_trace_tail_cap(theta, n_modes=80)
            self.assertGreater(cap60, 0)
            self.assertGreater(cap80, 0)
            self.assertLess(cap80, cap60)
            self.assertLess(cap80, Fraction(1, 10**41))
        with self.assertRaises(ValueError):
            cross_trace_tail_cap(Fraction(1,3))

    def test_traversal_rejected(self):
        good = {"id": "bad", "path": "../../outside.py", "method": "exact_rational",
                "evidence": "numerical", "claim": "bad", "assumptions": ["x"],
                "falsifiers": ["x"], "limitations": ["x"],
                "timeout_s": 2, "default_run": True}
        with self.assertRaisesRegex(ValueError, "escapes repository"):
            validate(good, ROOT)

    def test_replay_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            home = Path(temp)
            folder = home / "runs"
            item = {"id": "control", "path": "rh_lab/controls.py", "method": "exact_rational",
                    "evidence": "paper_proved", "claim": "fixture", "assumptions": ["x"],
                    "falsifiers": ["x"], "limitations": ["x"],
                    "timeout_s": 15, "default_run": True}
            manifest = home / "catalog.json"
            manifest.write_text(json.dumps({"version": 1, "experiments": [item]}), encoding="utf-8")
            record = run("control", root=ROOT, catalog=manifest, output_dir=folder, replay=True)
            self.assertTrue(record["passed"])
            self.assertTrue(record["deterministic"])
            self.assertEqual(len(list(folder.glob("*.json"))), 1)
            self.assertEqual(json.loads(Path(record["artifact_path"]).read_text())["code_sha256"], record["code_sha256"])

    def test_declared_stdin_pipeline(self):
        with tempfile.TemporaryDirectory() as temp:
            home = Path(temp)
            base = {"path": "rh_lab/controls.py", "method": "exact_rational",
                    "evidence": "numerical", "claim": "test", "assumptions": ["test"],
                    "falsifiers": ["test"], "limitations": ["test"],
                    "timeout_s": 15, "default_run": False}
            manifest = home / "catalog.json"
            manifest.write_text(json.dumps({"version": 1, "experiments": [
                dict(base, id="producer"), dict(base, id="consumer", stdin_from="producer")]}),
                                encoding="utf-8")
            result = run("consumer", root=ROOT, catalog=manifest,
                         output_dir=home / "runs", replay=True)
            self.assertTrue(result["passed"])
            self.assertIsNotNone(result["first"]["input_sha256"])
            self.assertIsNotNone(result["stdin_source_sha256"])

    def test_existing_registry_when_checkout_available(self):
        # Scratch workspaces can contain only the new lab package; full GitHub
        # checkout CI verifies all registered paths.
        if not (ROOT / "experiments" / "arch_hankel_sine2_k6_exact_intervals.py").exists():
            self.skipTest("full repo checkout required for catalog validation")
        registry = load()
        self.assertIn("fake_quartet_control", registry)
        self.assertIn("sine2_cross_trace_budget", registry)


if __name__ == "__main__":
    unittest.main()
