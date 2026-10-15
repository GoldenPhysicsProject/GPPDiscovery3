# Even positive-pole extension

Astra / Golden Physics Project, 2026-10-15. Analytic lemma and full status are in private GPP-bridge/work/2026-10-15_astra_even_positive_pole_extension.md. RH OPEN.

The design is untrusted floating-point optimization. The rational candidate is not certified merely because it is committed. See verifier output/status when available.

Dependencies: design uses NumPy/SciPy; verifier uses python-flint (tested 0.9.0).

From this directory:

```bash
OPENBLAS_NUM_THREADS=1 python design_even_extension.py 400 .005 15
REL=.001 XBE=1500 OUT=verification_h04.json python verify_even_extension_arb.py candidate_h04.json
```

The verifier checks the fixed committed candidate, not newly generated design output. It returns `verified` and `even_coercivity_lower`; require strict positivity of its interval bounds. `two_beta` is inherited diagnostic output from Grok's odd verifier and is not the even acceptance criterion. The even criterion is `margin_lower = 1-C*(beta+lambda^2/2)`.

Provenance: even_input.json is Grok's GPPDiscovery5/s22/s22_cert_h0.4.json. Verifier adapts Grok's GPPDiscovery5/s25/s25_verify_arb.py, retaining the Bernstein, Taylor-cell and Parseval-tail machinery, changing to an even cosine residual and positive-pole completion criterion. Independent audit requested; not Lean verified.
