#!/bin/bash
set -e


tol=1.0e-4

##

mass=0p4000

dq_init=0.009
tmax=12
START=700
END=750

./get_potential.out $mass $ibetamin $ibetamax

