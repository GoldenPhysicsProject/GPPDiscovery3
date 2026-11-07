# Critical n-gon rational residues and their hyperbolic shell geometry

Codex/GPT discovery track, 2026-08-25.

Two independently derived families meet at the critical numerator power.

For an n-gon, the evanescent power producing the universal finite scalar residue is

k = n-2,

so the numerator is mu^(2(n-2)) and

R_n := lim_{epsilon->0} I_n^(4-2epsilon)[mu^(2(n-2))]
     = -1/[(n-1)(n-2)],   n>=3.

For a general mu^(2k) numerator on a two-particle cut, the normalized hyperbolic radial density is

rho_k(r)=2k tanh r sech(r)^(2k).

Substituting k=n-2 gives the critical n-gon shell

boxed:

rho_n^crit(r)
 = 2(n-2) tanh r sech(r)^(2n-4),
 n>=3.

Its exact CDF and uniform coordinate are

F_n(R)=1-sech(R)^(2n-4),

U_n=sech(R)^(2n-4)=(2mu/M)^(2n-4) ~ Uniform(0,1).

The unique shell maximum occurs at

boxed:

tanh r_n^* = 1/sqrt(2n-3),

mu_n^*/M = (1/2) sqrt[2(n-2)/(2n-3)].

All transverse-mass moments are

E_n[mu^q]
 = [2(n-2)/(2(n-2)+q)] (M/2)^q.

In particular

E_n[(2mu/M)^2] = (n-2)/(n-1).

Combining this with the universal rational residue gives the exact algebraic relation

boxed:

|R_n| * E_n[(2mu/M)^2] = 1/(n-1)^2.

Thus the rational-residue hierarchy and the normalized radial shell hierarchy share a simple square law. This is a structural identity between the two scalar families; it does NOT assert that a physical YM or gravity n-gon enters with unit critical-numerator coefficient.

Examples:

n=3: R=-1/2,   E[(2mu/M)^2]=1/2,   product=1/4;
n=4: R=-1/6,   E[(2mu/M)^2]=2/3,   product=1/9;
n=5: R=-1/12,  E[(2mu/M)^2]=3/4,   product=1/16;
n=6: R=-1/20,  E[(2mu/M)^2]=4/5,   product=1/25.
