# Universal radial shell family for mu^(2k) numerator sectors

Codex/GPT discovery track, 2026-08-25.

Let a D-dimensional two-particle cut have channel mass M and transverse mass mu,

mu = M/(2 cosh r),
beta = sqrt(1-4mu^2/M^2) = tanh r.

Any numerator sector proportional to mu^(2k), k>=1, acquires from two-body phase space the universal radial factor

W_k(r) = tanh r * sech(r)^(2k)

up to the overall constant (M/2)^(2k).

## Exact normalization

Since

d/dr sech(r)^(2k) = -2k tanh r sech(r)^(2k),

we obtain

boxed:

int_0^infty W_k(r) dr = 1/(2k).

Hence

rho_k(r) = 2k tanh r sech(r)^(2k)

is a normalized radial probability density.

Its CDF is

F_k(R)=1-sech(R)^(2k).

Therefore

boxed:

U_k = sech(R)^(2k) = (2mu/M)^(2k)

is exactly Uniform(0,1) under the normalized radial law.

## Induced transverse-mass law

The normalized mu density is

boxed:

rho_{k,mu}(mu)
 = 2k (2/M)^(2k) mu^(2k-1),
0 <= mu <= M/2.

All moments close:

boxed:

E_k[mu^q]
 = [2k/(2k+q)] (M/2)^q,
q>-2k.

For k=2 this reproduces the previously derived mu^4 law

rho(mu)=64 mu^3/M^4,
E[mu^q]=4/(q+4)(M/2)^q.

## Unique shell maximum

Set x=tanh r in [0,1). Then

W_k=x(1-x^2)^k.

Differentiation gives

dW_k/dx=(1-x^2)^(k-1)[1-(2k+1)x^2].

Thus the unique interior/global maximum is

boxed:

tanh r_k^* = 1/sqrt(2k+1),

sech^2 r_k^* = 2k/(2k+1),

mu_k^* = (M/2) sqrt(2k/(2k+1)),

and

boxed:

W_k^* = (1/sqrt(2k+1)) [2k/(2k+1)]^k.

For k=2 this gives tanh r*=1/sqrt5 and mu*=M/sqrt5 exactly.

## Interpretation

This family applies to any generalized-unitarity numerator monomial mu^(2k) before insertion of the angular propagator kernel. It therefore supplies a reusable radial organization for higher-dimensional numerator sectors at higher loops or higher cuts. It does not by itself determine which powers k occur in a given renormalizable theory or topology; that remains a numerator power-counting/state-sum question.
