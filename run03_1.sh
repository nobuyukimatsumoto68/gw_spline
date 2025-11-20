#!/bin/bash
set -e


tol=1.0e-4

##


mass=0p3000

dq_init=0.01
tmax=8
START=360
END=420

# dq_init=0.0007
# tmax=16
# START=750
# END=800

# dq_init=0.0007
# tmax=16
# START=780



INT=1
increment=2


maxiter=5000
tau=5.0e-5

float=$(echo "-l(${dq_init})/l(10)" | bc -l)
expn=${float%.*}
expn=$(echo "$expn+$increment" | bc)



echo "expn = $expn, dq_init = $dq_init"

for (( ibeta=$START; ibeta<=$END; ibeta+=$INT ))
do
    echo "### ibeta=$ibeta ###"

    ./hist_spline.o $ibeta $expn $tmax $mass $dq_init $maxiter $tau #

    value=$(< "./fit_params_32c_m${mass}/c_$ibeta.dat")
    echo $value
    value_=$(echo ${value} | sed 's/e/\*10\^/' | sed 's/+//')
    tol_=$(echo ${tol} | sed 's/e/\*10\^/' | sed 's/+//')

    if [ $(echo "$value_ < $tol_" | bc -l) -eq 1 ]; then
        echo "$value less than $tol. break."
        is_continue=false
    else
        echo "search fail. quit"
        exit 1
    fi

    dq_init=$(< "./fit_params_32c_m$mass/dq_init_$ibeta.dat")
    dq_init_=$(echo ${dq_init} | sed 's/e/\*10\^/' | sed 's/+//')
    float=$(echo "-l(${dq_init_})/l(10)" | bc -l)
    expn=${float%.*}
    expn=$(echo "$expn+$increment" | bc)

    echo "expn = $expn, dq_init = $dq_init"

    if [ $( echo "$expn > 15" | bc ) -eq 1 ]; then
        break
    fi
done

