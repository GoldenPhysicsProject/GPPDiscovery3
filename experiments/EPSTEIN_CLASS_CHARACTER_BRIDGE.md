# The missing sector projection, in finite exact arithmetic

Daniel's compact idèle-class/Peter–Weyl proposal has a concrete
counterexample laboratory inside the classical class group of
Q(sqrt(-5)), discriminant -20, class number 2.

The integral positive forms Q0=x²+5y² and Q1=2x²+2xy+3y²
represent the two classes. Their theta/Epstein Dirichlet coefficients
r_Qj(n)/2 are nonnegative integer counts; for ideals of the maximal
order with units ±1 these are respectively the principal and
nonprincipal ideal-count coefficients at norm n (standard class-form
correspondence, not formalized here).

The C2 group algebra has idempotents
   e_+ = (1+g)/2, e_-=(1-g)/2;
its character Fourier coordinates of the ideal-class coefficient vector
are
   A_+(n) = [r_Q0(n)+r_Q1(n)]/2,
   A_-(n) = [r_Q0(n)-r_Q1(n)]/2.
They are the coefficients of the Dedekind zeta of Q(sqrt(-5))
and the nontrivial unramified class-character Hecke L-function.
The known global Euler products come from character multiplicativity.
The Python experiment verifies ONLY the finite resulting
coprime-multiplicativity checks; they do not prove Hecke theory.

**Exact counterexample**:
principal-class coefficients f(2)=f(3)=0 but f(6)=2.
Since gcd(2,3)=1, f is not multiplicative.
After character projection, A_+(2)=1,A_+(3)=2,A_+(6)=2,
and A_-(2)=-1,A_-(3)=-2,A_-(6)=2.
Both projected channels pass all coprime pairs ab<=128.

So one positive integral lattice with Poisson self-duality
can generate an Epstein zeta WITHOUT Euler product and
with off-line zeros (known theorem, Voronin for h(D)>1).
Projecting by *class characters* restores the correct
multiplicative channels. A mixing of individually legitimate
same-weight L-functions is **not** itself a single pure
Hecke eigencharacter. This explains a known off-line mechanism
without postulating an RH-strong positivity property.

**Specific implication for RHpositive:** In Daniel's compact
norm-one idèle-class K programme, do character-sector decomposition
BEFORE spectral determinant/boundary completion. Retain a single
genuinely multiplicative eigencharacter, not an arbitrary sum of
K-sector correlators. This sector projection alone certainly does NOT
prove GRH for any L-function. The missing analytic theorem remains
constructing the positive closed spectral realization and arithmetic
correspondence in that irreducible sector, rather than importing
zeros or simply asserting a polarisation.

Controls: genuine graded P1 zeta with different Tate weights;
shifted Ftheta formal fractional weights; Epstein integral-lattice
same-weight character mixture; Davenport–Heilbronn lack of Euler product.
Don't treat integer Dirichlet representation counts as the same
notion as finite-field point counts.

The theorem that Epstein zeta has off-line zeros for discriminants
with h>1 is established prior art (Voronin, cited e.g. Jutila–Srinivas);
Grok's GPP-bridge/GROK_RESEARCH_NOTES.md S99 already carried Arb
experiments for D=-20/-84. This new experiment merely makes the finite
class-character Fourier mechanism explicitly reproducible in Discovery2.
RH and BSD OPEN.
