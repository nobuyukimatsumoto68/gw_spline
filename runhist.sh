#!/bin/bash
set -e

# mass=0p4000
# mass=0p3000
# mass=0p2000
mass=0p1000

ibetamin=0
ibetamax=1000

./get_potential.o $mass $ibetamin $ibetamax

