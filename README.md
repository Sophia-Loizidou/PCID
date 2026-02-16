# PCID
Permutation-based Circular IsolateDetect (PCID) is a change-point detection algorithm for multiple change-points in the mean of circular data, as proposed in Loizidou, Anastasiou and Ley (2026+).
The repository contains the following files:

* PCID.R: functions required to use the PCID algorithm 
* Simulations_code.R: functions required to run simulations (requires PCID.R)
* Simulations.R: signals and code to run the simulations as in Section 5.2 of the paper (requires Simulations_code.R)

## Using PCID
Example on how to use the method:
```r
library(circular)
set.seed(2)
x <- rvonmises(200, circular(0), 4) + c(rep(0, 100), rep(2, 100))
cpt <- PCID(x, FDR = 0.01)
```

The algorithm chooses the appropriate values of $\alpha$ (significance level of permutation tests) and $B$ (maximum number of permutations performed), according to the Type I error of the algorithm (denoted by FDR) chosen by the user.
By default, the minimum value of $\alpha$ is $\alpha = 0.001$ and the maximum value of $B$ is $B = 1000$.
The user can override this by setting 'override_default = T'.
This can slow down the computation.

```r
cpt <- PCID(x, FDR = 0.01, verbose = T, override_default = T)
```

The user can also manually choose the significance level used for the permutation tests.
In this case 'override_default = T' is required.

```r
cpt <- PCID(x, alpha = 0.002, override_default = T)
```
