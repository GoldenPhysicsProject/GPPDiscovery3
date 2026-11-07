"""Standard-library regression checks for the finite analytic identities."""

import unittest

from research.prime_phase_rigidity_sparse import (
    euler_correction,
    real_dilation_discontinuity,
    reflection_completion,
    upper_bound_on_line,
)


class PrimePhaseTests(unittest.TestCase):
    def test_reflection_and_real_type(self):
        speeds = {2: 0, 3: 0.25, 5: 0.6, 7: 1}
        s = 0.34 + 2.7j
        z = reflection_completion(s, speeds)
        self.assertLess(abs(z - reflection_completion(1-s, speeds)), 1e-12)
        self.assertLess(abs(z.conjugate() - reflection_completion(s.conjugate(), speeds)), 1e-12)
        self.assertGreater(abs(z), 0)

    def test_critical_line_multiplier_positive_bounded(self):
        speeds = {2: 0, 3: 0.25, 5: 0.8, 7: 1}
        upper = upper_bound_on_line(speeds)
        for t in [0, 1, 2, 14.134725141734693, 100, 100000]:
            z = reflection_completion(0.5 + 1j*t, speeds)
            self.assertLess(abs(z.imag), 2e-10 * max(1, abs(z.real)))
            self.assertGreater(z.real, 1 / upper * (1 - 1e-10))
            self.assertLess(z.real, upper * (1 + 1e-10))

    def test_normalization_and_identity_rates(self):
        self.assertLess(abs(euler_correction(0.5+14j, {2: 1, 3: 1}) - 1), 1e-14)
        self.assertLess(abs(reflection_completion(0.5+0j, {2: 0, 3: 0.1}) - 1), 1e-14)

    def test_discontinuity_witness(self):
        rows = real_dilation_discontinuity(count=12)
        self.assertGreater(len(rows), 9)
        # Select a late convergent whose original real log is tiny but
        # altered log is unbounded in size along the full sequence.
        r = rows[-1]
        self.assertLess(abs(float(r['log_n'])), 0.001)
        self.assertGreater(abs(float(r['altered_log_n'])), 100)


if __name__ == "__main__":
    unittest.main()
