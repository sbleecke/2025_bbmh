# Numerical Experiments
This directory contains code to reproduce the numerical experiments described in the manuscript.

This code is developed with Julia version 1.10.5 and Python 3.13.7. To reproduce the results, start Julia in this directory and execute the following commands in the Julia REPL to create the figures shown in the paper.

## Python code: peakon plot
The peakon plot (Figure 1) can be reproduced with the Jupyter Notebook `BBMH-traveling-waves.ipynb`.

The versions of the used software (packages) are

* Python: 3.13.7
* Numpy: 2.2.6
* Matplotlib: 3.10.1
* Ipywidgets: 8.1.5

## Julia code: remaining numerical experiments

The code is developed with Julia v1.10.5. First, you need to install
Julia. Then, you can open the Julia REPL from the console

```bash
julia --project=code
```

In the Julia REPL, you first need to load the code
(you can copy and paste the following command directly
into the Julia REPL, including the `julia>` prompt):

```julia
julia> include("code/code.jl")
```

To reproduce the numerical Petviashvili solution, execute

```julia
julia> petviashvili_plot_validation()
```


To reproduce the AP property tests with latex output in the same
order as the tables in the manuscript, execute

```julia
julia>  AP_test_in_eps(alg = AGSA342(), latex = true) # type I, GSA

julia>  AP_test_in_eps(alg = AGSA342(), latex = true, vzero_boolval = true) # type I, GSA, v=0

julia>  AP_test_in_eps(alg = ARS443(), latex = true) # type II, GSA, ARS

julia>  AP_test_in_eps(alg = ARS443(), latex = true, vzero_boolval = true) # type II, GSA, ARS, v=0

julia>  AP_test_in_eps(alg = SSP2ImEx332(), latex = true) # type I, not GSA

julia>  AP_test_in_eps(alg = BPR343(), latex = true) # type II, GSA, not ARS
```

To plot the error growth with a Petviashvili generated reference solution, execute

```julia
julia> petviashvili_error_growth(eps = 1e-3)
```

To plot the error growth with BBM as a reference solution, execute

```julia
julia> solitary_wave_error_growth(eps = 1e-1)

julia> solitary_wave_error_growth(eps = 1e-10)
```
