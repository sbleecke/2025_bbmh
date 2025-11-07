# Asymptotic-preserving and energy-conserving methods for the BBM equation and its hyperbolic approximation

Latest build of `bbmh__arXiv.pdf`:
[view](https://github.com/ranocha/2024_bbmh_dev/blob/build/bbmh__arXiv.pdf),
[download](https://github.com/ranocha/2024_bbmh_dev/raw/build/bbmh__arXiv.pdf)

The goal is to develop numerical methods for the hyperbolic approximation
of the BBM equation suggested by Gavrilyuk and Shyue with the following
properties:

- AP in time
- A semidiscretely EC methods for BBMH converges to a semidiscretely EC method for BBM
- Relaxation in time makes it fully discretely EC
- We can hopefully observe linear vs. quadratic error growth in time for both BBMH and in the limit

## Peakon Plot
The peakon plot (Figure 1) can be reproduced with the Jupyter Notebook `BBMH-traveling-waves.ipynb`.

The versions of the used software (packages) are

Python: 3.13.7
Numpy: 2.2.6
Matplotlib: 3.10.1
Ipywidgets: 8.1.5

## Using the code

The code is developed with Julia v1.10.5. First, you need to install
Julia. Then, you can open the Julia REPL from the console

```bash
julia --project=code
```

In the Julia REPL, you first need to load the code

```julia
julia> using Revise; includet("code/code.jl")
```

Using Revise ensures that changes are directly visible in the REPL
when you modify the code file. If you do not want to do that, you
can simply run

```julia
julia> include("code/code.jl")
```
instead.

Reproduce the numerical Petviashvili solution
```julia
julia>  petviashvili_plot_validation()
```


Reproducing the AP property tests with latex output

```julia
julia>  AP_test_in_eps(alg = ARS222(), latex = true)

julia>  AP_test_in_eps(alg = SSP2ImEx332(), latex = true)

julia>  AP_test_in_eps(alg = ARS443(), latex = true)

julia>  AP_test_in_eps(alg = AGSA342(), latex = true)
```

using $$v = 0$$ as initial data

```julia
julia>  AP_test_in_eps(alg = ARS222(), latex = true, vzero_boolval = true)

julia>  AP_test_in_eps(alg = SSP2ImEx332(), latex = true, vzero_boolval = true)

julia>  AP_test_in_eps(alg = ARS443(), latex = true, vzero_boolval = true)

julia>  AP_test_in_eps(alg = AGSA342(), latex = true, vzero_boolval = true)
```

Error growth with a Petviashvili generated reference solution
```julia
julia> petviashvili_error_growth(eps_ref = 1e-3)

julia> petviashvili_error_growth(eps_ref = 1e-10)
```
error growth with BBM as a reference solution
```julia
julia> solitary_wave_error_growth(eps = 1e-3)

julia> solitary_wave_error_growth(eps = 1e-10)
```

