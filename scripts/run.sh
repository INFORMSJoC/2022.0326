#!/bin/bash
MPIJULIA="mpiexec -n 2 julia --project=@. run.jl"
JULIA="julia --project=@. run.jl"

mkdir -p logs

run() {
    $JULIA    tcuc-central decomposed/$1-2z $2 $3 $5 |& tee logs/$1-2z-tcuc-central.log
    $MPIJULIA tcuc-isf     decomposed/$1-2z $2 $3 $5 |& tee logs/$1-2z-tcuc-isf.log
    $MPIJULIA tcuc-theta   decomposed/$1-2z $2 $3 $5 |& tee logs/$1-2z-tcuc-theta.log
    $JULIA    scuc-central decomposed/$1-2z $2 $4 $5 |& tee > logs/$1-2z-scuc-central.log
    $MPIJULIA scuc-isf     decomposed/$1-2z $2 $4 $5 |& tee > logs/$1-2z-scuc-isf.log
}

case $1 in
    #                   demand  tcuc-scale scuc-scale careful?
    case300)     run $1 0.60    1.00       1.00       false ;;
    case1888rte) run $1 0.60    1.50       1.75       false ;;
    case1951rte) run $1 0.60    1.25       1.50       false ;;
    case2848rte) run $1 0.60    1.75       2.00       true  ;;
    case3012wp)  run $1 0.60    1.25       1.50       true  ;;
    case3375wp)  run $1 0.60    1.50       1.75       false ;;
    case6468rte) run $1 0.60    2.25       2.50       false ;;
    case6515rte) run $1 0.60    3.00       3.50       false ;;
    all)
        exec ./run.sh case1888rte
        exec ./run.sh case1951rte
        exec ./run.sh case2848rte
        exec ./run.sh case3012wp
        exec ./run.sh case3375wp
        exec ./run.sh case6468rte
        exec ./run.sh case6515rte
        ;;
    *) echo "Invalid case: $1" ;;
esac
