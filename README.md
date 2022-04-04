# ssDMEWMAC

This repository contains the code for calculating the ssMEWMAC and ssDMEWMAC control chart statistics. 

The code here is based of the 2007 hallmark paper of D.M. Hawkins and E.M. Maboudou-Tchao. 
It has taken the Fortran 95 subroutine for the "u-transformation" using Cholesky decomposition,
givens rotations, and recursive residuals to convert a stream of multivariate X vectors and return
instead the corresponding R vector (recursive residual vector) and U vector (standard normal vectors).

The output of the u-transform functions then can be passed to the function to calculate the ssMEWMAC 
and ssDMEWMAC statistics as desired. These statistics can then be plotted on control charts; functions 
to do so on seaborn charts are provided.

An example docker snippet is included below:

```
for d in 3 4 5 6 7 8 9 10 15 20 25 50 ; do
    docker run -it -v /Users/robertresch/Documents/ssDMEWMAC/data:/out -e SIMULATION_TYPE="IC" -e METH="ssDMEWMC" -e DIME=$d -e ITERATIONS=10 -e LOCAL_OR_GCP="LOCAL" ssdmewmac
done    
```

In this revision, iterations, methods,  and dimension are all now required to be passed for the simulation.