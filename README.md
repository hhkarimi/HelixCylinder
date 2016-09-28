# HelixCylinder
Modeling bacterial flagella as a flexible helix as it swims through blood vessels

## LSintegralsTrapz.py
Solves for the flow-field generated by a point force (stokeslet) bounded by a cylinder [a la Liron & Shahar](http://dx.doi.org/10.1017/S0022112078001366).  The notation and equations are described in equations (4.4)--(4.19) and the appendix of that paper.  The solution is given numerically in terms of infinite summations and unbounded integrals with modified Bessel functions acting as kernels.

Since the computation time (~ 2 seconds) is too slow to be generated at run-time when incorporated with the elastic flagella model, it is generated a priori and stored, then loaded into an elastic solver (Dynamic Elastic Rods).

## LighthillSBT.hh
c++ implementation of the fluid flow-field generated in LSintegralsTrapz.py with elastic rods solver (Dynamic Elastic Rods).  This implementation is incomplete and is not intended to represent a functional program.  It is included here to demonstrate the inclusion of a complicated fluid model (LSintegralsTrapz.py) in a numerical algorithm specializing in elastic solid mechanics.
