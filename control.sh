#!/usr/bin/env bash

# This is a simple control script to run either the in_control or out_of_control simulations
# based upon the environment variable passed for SIMULATION_TYPE

if [ $SIMULATION_TYPE == "IC" ]; then
  python in_control_simulation.py
elif [ $SIMULATION_TYPE == "OOC" ]; then
  python out_of_control_simulation.py
else
  echo -e "Please enter a valid simulation type: IC - in control or OOC - out of control"
  exit 1;
fi