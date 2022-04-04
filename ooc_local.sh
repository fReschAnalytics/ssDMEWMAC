#! /bin/bash

export MKL_NUM_THREADS=8
export ARLS=( 100 200 500 1000 2000 )
export DELTAS=( 0.25 0.50 1.00 1.50 2.00 2.50 3.00 )
export DIMENSIONS=( 7 )
export ITERATIONS=10000
export LAMBDAS=( 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 )
export METHODS=( "ssMEWMC" )
export WAITTIMES=( 20 50 100 500 )

for d in "${DIMENSIONS[@]}"; do
    for m in "${METHODS[@]}"; do
        for l in "${LAMBDAS[@]}"; do
            for a in "${ARLS[@]}"; do
                for x in "${DELTAS[@]}"; do
                    for w in "${WAITTIMES[@]}"; do
                        echo "Running the simulation for Dimension: ${d}, Method ${m}, Lambda: ${l}, Target ARL: ${a}, Delta: ${x}, and a Pertubation at Observation: ${w}."
                          export METH="$m"
                          export DIME=$d
                          export LAMBDA=$l
                          export ARL=$a
                          export DELTA=$x
                          export WAIT_TIME=$w
                          export LOCAL_OR_GCP="LOCAL"

                          python out_of_control_simulation_local.py

                    done
                done
            done
        done
    done
done