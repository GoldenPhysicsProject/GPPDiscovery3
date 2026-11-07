# Three-box cyclic residue hierarchy in the four-point symmetric invariant ring

Codex/GPT discovery track, 2026-08-25.

Let

\[
S_m(x,y)=\sum_{j=0}^{m} j!(m-j)!\,x^j y^{m-j}
\]

be the positive coefficient polynomial appearing in the all-order supercritical massless-box residue

\[
I_4[\mu^{2(m+2)}]\to
(-1)^{m+1}\frac{(m+1)!}{(2m+3)!}S_m(s,t).
\]

For four-point kinematics `s+t+u=0`, define the cyclic three-box sum

\[
C_m(s,t,u)=S_m(s,t)+S_m(t,u)+S_m(u,s).
\]

Because `S_m(x,y)=S_m(y,x)`, the cyclic sum is invariant under all permutations of `(s,t,u)`: every transposition turns the three unordered channel pairs `{s,t},{t,u},{u,s}` into the same set. Therefore `C_m` is a fully symmetric homogeneous polynomial of degree `m`.

On the hyperplane `e_1=s+t+u=0`, the symmetric polynomial ring reduces to

\[
\mathbb Q[e_2,e_3],\qquad e_2=st+tu+us,\quad e_3=stu.
\]

Hence there are rational/integer coefficients `c_{m,a,b}` such that

\[
\boxed{
C_m(s,t,u)=
\sum_{\substack{a,b\ge0\\2a+3b=m}}
 c_{m,a,b}\,e_2^a e_3^b
\qquad (s+t+u=0).
}
\]

This gives a finite invariant basis at every transverse-power order. The number of independent structures is the number of nonnegative solutions of `2a+3b=m`.

Exact low orders obtained by symbolic elimination of `u=-s-t` are:

\[
C_0=3,
\qquad C_1=0,
\]

\[
C_2=-7e_2,
\qquad C_3=30e_3,
\]

\[
C_4=88e_2^2,
\qquad C_5=-1092e_2e_3,
\]

\[
C_6=-2700e_2^3+3924e_3^2,
\]

\[
C_7=66096e_2^2e_3,
\]

\[
C_8=153216e_2^4-604656e_2e_3^2,
\]

\[
C_9=-6209280e_2^3e_3+2043360e_3^3.
\]

The `m=2` identity is exactly the previously used all-plus gravity collapse, since

\[
-e_2=\frac12(s^2+t^2+u^2),
\]

so

\[
C_2=\frac72(s^2+t^2+u^2).
\]

## Structural consequence

At higher `mu` powers, any three-channel box-only four-point amplitude assembled with the same channel-symmetric scalar residue family can be reduced before integration/reconstruction to a small basis of the two elementary four-point invariants `e_2,e_3`. Beginning at degree six, more than one invariant structure is available (`e_2^3` and `e_3^2`), explaining why the exceptionally simple all-plus `mu^8` polynomial does not persist as a single power at arbitrary order.

This is an algebraic organization theorem for the residue polynomials. It does not assert that a physical YM or gravity amplitude at every order is box-only, nor does it supply state-sum coefficients.
