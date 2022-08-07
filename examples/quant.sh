#!/bin/bash
# mpirun -n 4  ./run_pselinv_linux_release_v2.0 -H H180.csc -r 2 -c 2
# 这是所有平均值小于1e-10的supernode的下标
#mpirun -n 4  ./run_pselinv_linux_release_v2.0 -H H180.csc -r 2 -c 2 -Q 21,11 20,10 21,10 21,12 7,4 26,1 26,0 7,3 17,14 20,11 17,13 22,4 22,3 23,8 22,14 23,9 21,15 21,16 22,13 19,4 23,12 23,10 19,3 12,8 12,9 26,6 22,10 23,11 22,15
# 这是所有平均值小于1e-9的supernode下标
#mpirun -n 4  ./run_pselinv_linux_release_v2.0 -H H180.csc -r 2 -c 2 -Q 21,11 20,10 21,10 21,12 7,4 26,1 26,0 7,3 17,14 20,11 17,13 22,4 22,3 23,8 22,14 23,9 21,15 21,16 22,13 19,4 23,12 23,10 19,3 12,8 12,9 26,6 22,10 23,11 22,15 22,12 26,5 22,11 24,6 26,13
# 这是所有平均值小于1e-8的supernode下标
# mpirun -n 4  ./run_pselinv_linux_release_v2.0 -H H180.csc -r 2 -c 2 -Q  21,11 20,10 21,10 21,12 7,4 26,1 26,0 7,3 17,14 20,11 17,13 22,4 22,3 23,8 22,14 23,9 21,15 21,16 22,13 19,4 23,12 23,10 19,3 12,8 12,9 26,6 22,10 23,11 22,15 22,12 26,5 22,11 24,6 26,13
# 这是所有平均值小于1e-7的supernode下标
# mpirun -n 4  ./run_pselinv_linux_release_v2.0 -H H180.csc -r 2 -c 2 -Q 21,11 20,10 21,10 21,12 7,4 26,1 26,0 7,3 17,14 20,11 17,13 22,4 22,3 23,8 22,14 23,9 21,15 21,16 22,13 19,4 23,12 23,10 19,3 12,8 12,9 26,6 22,10 23,11 22,15 22,12 26,5 22,11 24,6 26,13 26,8 25,10 26,14 26,16 26,15 26,9
# 这是所有平均值小于1e-6的supernode下标
# mpirun -n 4  ./run_pselinv_linux_release_v2.0 -H H180.csc -r 2 -c 2 -Q 21,11 20,10 21,10 21,12 7,4 26,1 26,0 7,3 17,14 20,11 17,13 22,4 22,3 23,8 22,14 23,9 21,15 21,16 22,13 19,4 23,12 23,10 19,3 12,8 12,9 26,6 22,10 23,11 22,15 22,12 26,5 22,11 24,6 26,13 26,8 25,10 26,14 26,16 26,15 26,9 26,17 26,10 23,16 24,7 23,17 25,11 24,16 24,17 19,5 23,18 22,1



# for((i=0; i<=26; i++))
# do
# mpirun -n 4 ./run_pselinv_linux_release_v2.0 -H H180.csc -Q 26,$i
# mpirun -n 4 ./run_pselinv_linux_release_v2.0 -H H180.csc -Q $i,0
# done

# mpirun -n 4  ./my_dispalySuperNode_linux_release_v2.0 -H H180.csc -r 2 -c 2 -delta 1e-6

# mpirun -n 4  ./my_dispalySuperNode_linux_release_v2.0 -H H180.csc -r 2 -c 2 -delta 1e-7

# mpirun -n 4  ./my_dispalySuperNode_linux_release_v2.0 -H H180.csc -r 2 -c 2 -delta 1e-8

# mpirun -n 4  ./my_dispalySuperNode_linux_release_v2.0 -H H180.csc -r 2 -c 2 -delta 1e-9

# mpirun -n 4  ./my_dispalySuperNode_linux_release_v2.0 -H H180.csc -r 2 -c 2 -delta 1e-10
