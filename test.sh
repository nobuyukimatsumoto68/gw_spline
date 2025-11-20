# ./integrate.o 459 8 250. 6.8940597509908246864895887e-02 100000 1.0e-4

# make integrate.o

########### ibeta expn tmax dq_init iter_max tau
########### (iter_max for search root, tau for leapfrog)

# ./integrate.o 416 6 16. 1.0e-4 500 1.0e-5 #
# ./integrate.o 420 6 16. 1.0e-4 500 1.0e-5 #
# ./integrate.o 424 7 16. 4.0e-5 500 1.0e-5 #
# ./integrate.o 428 7 16. 4.0e-5 500 1.0e-5 #
# ./integrate.o 432 7 16. 2.0e-5 500 1.0e-5 #

# ./integrate.o 436 8 16. 1.0e-5 500 1.0e-5 #

# ./integrate.o 440 8 14. 4.0e-6 500 5.0e-5 #
# ./integrate.o 442 8 14. 1.0e-6 500 5.0e-5 #
# ./integrate.o 444 8 16. 1.0e-6 500 5.0e-5 #
# ./integrate.o 446 8 16. 1.0e-6 500 5.0e-5 #


tol=1.0e-6





##
# mass=0.4
# # dq_init=0.0000014
# # dq_init=0.025
# dq_init=0.1
# tmax=5
# START=400
# END=500
# INT=2
# increment=3







mass=0p2000
# dq_init=0.001
# dq_init=0.0001


START=400
END=500

INT=2
increment=2

# 832

tmax=12

ibeta=768
ibeta=780
dq_init=0.05

tmax=14
ibeta=784
dq_init=0.04

# ibeta=790
# dq_init=0.026






#######################################################





# mass=0p3000
# # dq_init=0.001
# # dq_init=0.0001


# START=400
# END=500

# INT=2
# increment=2




# tmax=8

# ibeta=360
# dq_init=0.01

# # ibeta=420
# # dq_init=0.001



# # tmax=8

# # ibeta=400
# # dq_init=0.003

# # ibeta=420
# # dq_init=0.001



# # tmax=12

# # ibeta=440
# # dq_init=0.0002

# # ibeta=444
# # dq_init=0.0001


# # tmax=16

# # ibeta=448
# # dq_init=0.00006


# # ibeta=450
# # dq_init=0.00001





#######################################################



# mass=0.4
# # dq_init=0.001
# dq_init=0.0001

# tmax=12
# # START=500
# # END=600

# START=600
# END=700

# INT=2
# increment=3

# # ibeta=600
# # dq_init=0.05

# # ibeta=650
# # dq_init=0.01

# # ibeta=700
# # dq_init=0.009

# # tmax=16
# # ibeta=750
# # dq_init=0.0007

# tmax=16
# ibeta=768
# dq_init=0.00001

# # tmax=16
# # ibeta=780
# # dq_init=0.0000004

# # tmax=16
# # ibeta=790
# # dq_init=0.00000002

# # tmax=16
# # ibeta=800
# # dq_init=0.0000000005

# # tmax=22
# # ibeta=810
# # dq_init=0.000000000005









increment=2


float=$(echo "-l(${dq_init})/l(10)" | bc -l)
expn=${float%.*}
expn=$(echo "$expn+$increment" | bc)


maxiter=1000
tau=5.0e-5

./a.out $ibeta $expn $tmax $mass $dq_init $maxiter $tau #




# mass=0.4
# dq_init=0.00025
# tmax=14
# START=600
# END=700
# INT=2
# increment=3





# ##
# mass=0.3
# dq_init=0.0000004
# tmax=18.
# START=450
# END=480
# INT=2


# ##
# mass=0.2
# dq_init=0.043606736254108075
# expn=13
# tmax=36.
# looprange={698..797..2}






# maxiter=250




