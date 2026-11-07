# Convention-fixed helicity factor for the mu^4 MHV scalar cut

Codex/GPT discovery track, 2026-08-25.

Consider the color-ordered s-channel cut of the external MHV configuration

1^- 2^- 3^+ 4^+.

Use all-outgoing massive-scalar trees and the cyclic ordering

A_L = A_4(l1_s,1^-,2^-,l2_s),
A_R = A_4((-l2)_s,3^+,4^+,(-l1)_s).

Badger's massive-scalar tree and its parity conjugate give, in this convention,

A_L = i mu^2 <12> / ([12] [1|l1|1>),

A_R = i mu^2 [34] / (<34> <3|(-l2)|3]).

For a massive cut leg l^2=mu^2 and a massless p,

<p|l|p] = [p|l|p> = 2 l.p,

with bracket orientation chosen consistently. These are precisely the two uncut massive propagator factors in the chosen routing. Calling them D_L and D_R, the product is

C_s^scalar = A_L A_R
 = - mu^4 * (<12>[34]/([12]<34>)) /(D_L D_R).

Thus the previously symbolic external factor may be chosen as

boxed:

Xi(1^-,2^-,3^+,4^+) = - <12>[34]/([12]<34>)

for this explicit cyclic/all-outgoing convention.

The sign changes if one reverses one of the scalar momentum/orientation conventions while simultaneously redefining the corresponding propagator denominator, so the convention-independent content is the complete product above rather than Xi in isolation.

Little-group audit

Under |i> -> t_i |i>, |i] -> t_i^{-1}|i],

<12>/[12] -> t_1^2 t_2^2 <12>/[12],

[34]/<34> -> t_3^{-2} t_4^{-2} [34]/<34>.

Hence Xi has weight

t_1^2 t_2^2 t_3^{-2} t_4^{-2},

exactly the required four-point helicity weight for 1^-2^-3^+4^+.

This closes the external helicity phase of the isolated adjoint-scalar mu^4 box sector. It does NOT assemble the full D=4-2epsilon gluon state sum, nor the mu^2 triangle/bubble sectors, and therefore is not yet a pure-Yang-Mills rational amplitude.

Reference: S. D. Badger, arXiv:0806.4600, massive-scalar tree amplitudes eqs. (56)-(57), with parity applied to the all-minus scalar tree.
