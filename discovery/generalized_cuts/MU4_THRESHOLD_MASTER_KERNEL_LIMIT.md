# Universal threshold limit of the massive S2 master kernel

For the fixed-radius massive-scalar cut,

\[
\mathcal J(A,B;C,D)=\int_{S^2}\frac{d\Omega}{(A+B\cdot\hat n)(C+D\cdot\hat n)},
\]

with

\[
B=-2p\,\mathbf K_1,\qquad D=-2p\,\mathbf K_2,
\qquad p=\frac M2\beta,
\qquad \beta=\sqrt{1-4\mu^2/M^2}=\tanh r.
\]

At the two-particle threshold `beta -> 0` (`r -> 0`, `mu -> M/2`), one has `B,D = O(beta)` while

\[
A\to A_0=K_1^2+M K_1^0,
\qquad
C\to C_0=K_2^2+M K_2^0.
\]

Provided `A0 C0 != 0` (the nonsingular threshold domain), dominated expansion of the affine denominators gives

\[
\frac1{(A+B\cdot\hat n)(C+D\cdot\hat n)}
=\frac1{A_0C_0}+O(\beta),
\]

and the linear angular correction integrates to zero because

\[
\int_{S^2}\hat n\,d\Omega=0.
\]

Therefore the first nonzero correction is quadratic:

\[
\boxed{
\mathcal J(r)=\frac{4\pi}{A_0C_0}+O(r^2).
}
\]

Since

\[
\mu^4=\frac{M^4}{16}+O(r^2),
\qquad
\tanh r=r+O(r^3),
\]

the isolated scalar-sector cut has the exact universal threshold asymptotic

\[
\boxed{
\int d\Pi_2\,C_s^{\rm scalar}
=
\frac{M^4}{128\pi A_0C_0}\,\Xi\,r
+O(r^3).
}
\]

Equivalently, because `beta=r+O(r^3)`,

\[
\boxed{
\int d\Pi_2\,C_s^{\rm scalar}
=
\frac{M^4}{128\pi A_0C_0}\,\Xi\,\beta
+O(\beta^3).
}
\]

Thus the threshold suppression of this D-dimensional `mu^4` sector is genuinely linear with a closed coefficient whenever the two uncut propagators remain nonsingular at threshold. This still does not supply the full D-dimensional gluon state sum, subtraction sectors, or gravity double copy.
