#! /bin/bash

ARLS=(100 200 500 1000 2000)
DELTAS=(0.25 0.50 1.00 1.50 2.00 2.50 3.00)
DIMENSIONS=(2 4 6 8 10 20 50 )
ITERATIONS=10000
LAMBDAS=(0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95)
METHODS=( "MEWMA" "DMEWMA" "ssMEWMA" "ssDMEWMA" "ssMEWMC" "ssDMEWMC" )
WAITTIMES=(20 50 100 500)


for d in "${DIMENSIONS[@]}"; do
    for m in "${METHODS[@]}"; do
        for l in "${LAMBDAS[@]}"; do
            for a in "${ARLS[@]}"; do
                for x in "${DELTAS[@]}"; do
                    for w in "${WAITTIMES[@]}"; do
                        echo "Running the simulation for Dimension: ${d}, Method ${m}, Lambda: ${l}, Target ARL: ${a}, Delta: ${x}, and a Pertubation at Observation: ${w}."

                            docker run -it -v /home/rresch/Documents/Projects/ssDMEWMAC/data:/out -e SIMULATION_TYPE="OOC" -e METH="${m}" -e DIME=${d} -e ITERATIONS=${ITERATIONS} -e LAMBDA=${l} -e ARL=${a} -e DELTA=${x} -e WAIT_TIME=${w} -e LOCAL_OR_GCP="LOCAL" ssdmewmac
                    done
                done
            done
        done
    done
done
