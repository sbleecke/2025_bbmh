# Install packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load packages
using LinearAlgebra: UniformScaling, I, diag, diagind, mul!, ldiv!, lu, lu!, norm
using SparseArrays: sparse, issparse, dropzeros!

using DelimitedFiles: readdlm
using Interpolations

using Symbolics: Symbolics
using LinearSolve: LinearSolve, KLUFactorization
using OrdinaryDiffEq

using ForwardDiff: ForwardDiff
using PreallocationTools: DiffCache, get_tmp

using SummationByPartsOperators

using LaTeXStrings
using Plots: Plots, plot, plot!, scatter, scatter!, savefig

using PrettyTables: PrettyTables, pretty_table, ft_printf

#####################################################################
# Utility functions
function plot_kwargs()
    fontsizes = (
      xtickfontsize = 14, ytickfontsize = 14,
      xguidefontsize = 16, yguidefontsize = 16,
      legendfontsize = 14)

    (; linewidth = 3, gridlinewidth = 2,
       markersize = 8, markerstrokewidth = 4,
       fontsizes...)
end


function compute_eoc(Ns, errors)
    eoc = similar(errors)
    eoc[begin] = NaN # no EOC defined for the first grid
    for idx in Iterators.drop(eachindex(errors, Ns, eoc), 1)
        eoc[idx] = -( log(errors[idx] / errors[idx - 1]) /
                      log(Ns[idx] / Ns[idx - 1]) )
    end
    return eoc
end


#####################################################################
# High-level interface of the equations and IMEX ode solver

rhs_full!(du, u, parameters, t) = rhs_full!(du, u, parameters.equation, parameters, t)
rhs_stiff!(du, u, parameters, t) = rhs_stiff!(du, u, parameters.equation, parameters, t)
rhs_nonstiff!(du, u, parameters, t) = rhs_nonstiff!(du, u, parameters.equation, parameters, t)
operator(rhs_stiff!, parameters) = operator(rhs_stiff!, parameters.equation, parameters)
dot_entropy(u, v, parameters) = dot_entropy(u, v, parameters.equation, parameters)

# Coefficients
struct SP111{T}
end
SP111(T = Float64) = SP111{T}()

function coefficients(::SP111{T}) where T
    e = one(T)
    o = zero(T)

    A_stiff = [e;;]
    b_stiff = [e]
    c_stiff = [e]
    A_nonstiff = [o;;]
    b_nonstiff = [o]
    c_nonstiff = [e]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end



struct ARS111{T}

end
ARS111(T = Float64) = ARS111{T}()

function coefficients(::ARS111{T}) where T
    e = one(T)

    A_stiff = [0 0;
               0 e]
    b_stiff = [0, e]
    c_stiff = [0, e]
    A_nonstiff = [0 0;
                  e 0]
    b_nonstiff = [e, 0]
    c_nonstiff = [0, 1]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end


struct ARS222{T}

end
ARS222(T = Float64) = ARS222{T}()

function coefficients(::ARS222{T}) where T
    two = convert(T, 2)
    γ = 1 - 1 / sqrt(two)
    δ = 1 - 1 / (2 * γ)

    A_stiff = [0 0 0;
               0 γ 0;
               0 1-γ γ]
    b_stiff = [0, 1-γ, γ]
    c_stiff = [0, γ, 1]
    A_nonstiff = [0 0 0;
                  γ 0 0;
                  δ 1-δ 0]
    b_nonstiff = [δ, 1-δ, 0]
    c_nonstiff = [0, γ, 1]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end



struct ARS443{T}

end
ARS443(T = Float64) = ARS443{T}()

function coefficients(::ARS443{T}) where T
    l = one(T)

    A_stiff = [0 0 0 0 0;
               0 l/2 0 0 0;
               0 l/6 l/2 0 0;
               0 -l/2 l/2 l/2 0;
               0 3*l/2 -3*l/2 l/2 l/2]
    b_stiff = [0, 3*l/2, -3*l/2, l/2, l/2]
    c_stiff = [0, l/2, 2*l/3, l/2, l]
    A_nonstiff = [0 0 0 0 0;
                  l/2 0 0 0 0;
                  11*l/18 l/18 0 0 0;
                  5*l/6 -5*l/6 l/2 0 0;
                  l/4 7*l/4 3*l/4 -7*l/4 0]
    b_nonstiff = [l/4 7*l/4 3*l/4 -7*l/4 0]
    c_nonstiff = [0, l/2, 2*l/3, l/2, l]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end
#####################################
#type I methods following KdVH

# SSP2-IMEX(2,2,2):
# 2nd order L-stable type I ImEX-RK, explicit part - NOT FSAL, implicit part - NOT SA; NOT GSA.

struct SSP2ImEx222{T}

end
SSP2ImEx222(T = Float64) = SSP2ImEx222{T}()

function coefficients(::SSP2ImEx222{T}) where T
    γ = one(T) - 1 / sqrt(2)

    A_stiff = [γ 0;
               1 - 2 * γ γ]
    b_stiff = [1/2, 1/2]
    c_stiff = [γ, 1-γ]
    A_nonstiff = [0 0;
                  1 0]
    b_nonstiff = [1/2, 1/2]
    c_nonstiff = [0, 1]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end


# SSP2-IMEX(3,3,2):
# 2nd order L-stable type I ImEX-RK, explicit part - NOT FSAL, implicit part - SA; NOT GSA.

struct SSP2ImEx332{T} end

SSP2ImEx332(T = Float64) = SSP2ImEx332{T}()

function coefficients(::SSP2ImEx332{T}) where T

    A_stiff = [1/4 0 0;
               0 1/4 0;
               1/3 1/3 1/3]
    b_stiff = [1/3, 1/3, 1/3]
    c_stiff = [1/4, 1/4, 1]
    A_nonstiff = [0 0 0;
                  1/2 0 0;
                  1/2 1/2 0]
    b_nonstiff = [1/3, 1/3, 1/3]
    c_nonstiff = [0, 1/2, 1]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end


# AGSA(3,4,2):
# 2nd order type I ImEX-RK, explicit part - FSAL, implicit part - SA; GSA.

struct AGSA342{T}

end
AGSA342(T = Float64) = AGSA342{T}()

function coefficients(::AGSA342{T}) where T
    l = one(T)

    A_stiff = [ 168999711*l/74248304 0 0 0;
                44004295*l/24775207 202439144*l/118586105 0 0;
               -6418119*l/169001713 -748951821*l/1043823139 12015439*l/183058594 0;
               -370145222*l/355758315 l/3 0 202439144*l/118586105]

    b_stiff = [-370145222*l/355758315 l/3 0 202439144*l/118586105]
    c_stiff = [168999711*l/74248304, 10233769644823783*l/2937995298698735, -22277245178359531777915943*l/32292981880895017395247558 , 1]

    A_nonstiff = [ 0 0 0 0;
                  -139833537*l/38613965 0 0 0;
                   85870407*l/49798258 -121251843*l/1756367063 0 0;
                   l/6 l/6 2*l/3 0]

    b_nonstiff = [l/6 l/6 2*l/3 0]
    c_nonstiff = [0, -139833537*l/38613965, 144781823980515147*l/87464020145976254, 1 ]

    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end


# SSP3-IMEX(3,4,3):
# 3rd order L-stable type I ImEX-RK, explicit part - NOT FSAL, implicit part - NOT SA; NOT GSA.
struct SSP3ImEx343{T}

end
SSP3ImEx343(T = Float64) = SSP3ImEx343{T}()

function coefficients(::SSP3ImEx343{T}) where T
    α = 0.241694260788
    β = 0.0604235651970
    η = 0.12915286960590

    A_stiff = [ α 0 0 0;
               -α α 0 0;
                0 1-α α 0;
                β η 1/2 - β - η - α α]
    b_stiff = [0, 1 / 6, 1 / 6, 2 / 3]
    c_stiff = [α, 0, 1, 1/2]
    A_nonstiff = [0 0 0 0;
                  0 0 0 0;
                  0 1 0 0;
                  0 1/4 1/4 0]
    b_nonstiff = [0, 1/6, 1/6, 2/3]
    c_nonstiff = [0, 0, 1, 1/2]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end
###################################
struct BPR343{T}

end
BPR343(T = Float64) = BPR343{T}()

function coefficients(::BPR343{T}) where T
    l = one(T)

    A_stiff = [0 0 0 0 0;
               l/2 l/2 0 0 0;
               5*l/18 -l/9 l/2 0 0;
               l/2 0 0 l/2 0;
               l/4 0 3*l/4 -l/2 l/2]
    b_stiff = [l/4, 0, 3*l/4, -l/2, l/2]
    c_stiff = [0, l, 2*l/3, l, l]
    A_nonstiff = [0 0 0 0 0;
                  l 0 0 0 0;
                  4*l/9 2*l/9 0 0 0;
                  l/4 0 3*l/4 0 0;
                  l/4 0 3*l/4 0 0]
    b_nonstiff = [l/4, 0, 3*l/4, 0, 0]
    c_nonstiff = [0, l, 2*l/3, l, l]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end



"""
    LZ564(T = Float64)

H. Liu and J. Zou.
“Some new additive Runge–Kutta methods and their applications”.
In: Journal of Computational and Applied Mathematics 190.1-2 (2006), pp. 74–98.
"""
struct LZ564{T} end
LZ564(T = Float64) = LZ564{T}()

function coefficients(::LZ564{T}) where T
    l = one(T)

    A_stiff = [0 0 0 0 0 0 0;
               -l/6 l/2 0 0 0 0 0;
               l/6 -l/3 l/2 0 0 0 0;
               3*l/8 -3*l/8 0 l/2 0 0 0;
               l/8 0 3*l/8 -l/2 l/2 0 0;
               -l/2 0 3 -3 l l/2 0;
               l/6 0 0 0 2*l/3 -l/2 2*l/3]
    b_stiff = [l/6, 0, 0, 0, 2*l/3, -l/2, 2*l/3]
    c_stiff = [0, l/3, l/3, l/2, l/2, l, l]
    A_nonstiff = [0 0 0 0 0 0 0;
                  l/3 0 0 0 0 0 0;
                  l/6 l/6 0 0 0 0 0;
                  l/8 0 3*l/8 0 0 0 0;
                  l/8 0 3*l/8 0 0 0 0;
                  l/2 0 -3*l/2 0 2 0 0;
                  l/6 0 0 0 2*l/3 l/6 0]
    b_nonstiff = [l/6, 0, 0, 0, 2*l/3, l/6, 0]
    c_nonstiff = [0, l/3, l/3, l/2, l/2, l, l]
    return A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff
end

# IMEX ARK solver
# This assumes that the stiff part is linear and that the stiff solver is
# diagonally implicit.
function solve_imex(rhs_stiff!, rhs_stiff_operator, rhs_nonstiff!,
                    q0, tspan, parameters, alg;
                    dt,
                    relaxation = false,
                    callback = Returns(nothing),
                    save_everystep = false)
    A_stiff, b_stiff, c_stiff, A_nonstiff, b_nonstiff, c_nonstiff = coefficients(alg)

    s = length(b_stiff)
    @assert size(A_stiff, 1) == s && size(A_stiff, 2) == s &&
            length(b_stiff) == s && length(c_stiff) == s &&
            size(A_nonstiff, 1) == s && size(A_nonstiff, 2) == s &&
            length(b_nonstiff) == s && length(c_nonstiff) == s
    Base.require_one_based_indexing(A_stiff, b_stiff, c_stiff,
                                    A_nonstiff, b_nonstiff, c_nonstiff)

    q = copy(q0) # solution
    if save_everystep
        sol_q = [copy(q0)]
        sol_t = [first(tspan)]
    end
    y = similar(q) # stage value
    z = similar(q) # stage update value
    stage_updates = Vector{typeof(q)}(undef, s) # y^i - u^n
    t = first(tspan)
    tmp = similar(q)
    k_stiff_q = similar(q) # derivative of the previous state
    k_stiff = Vector{typeof(q)}(undef, s) # stage derivatives
    k_nonstiff = Vector{typeof(q)}(undef, s) # stage derivatives
    for i in 1:s
        k_stiff[i] = similar(q)
        k_nonstiff[i] = similar(q)
    end

    # Setup system matrix template and factorizations
    W, factorization, factorizations = let
        a = findfirst(!iszero, diag(A_stiff))
        factor = a * dt
        W = I - factor * rhs_stiff_operator

        if W isa UniformScaling
            # This happens if the stiff part is zero
            factorization = W
        else
            factorization = lu(W)
        end

        # We cache the factorizations for different factors for efficiency.
        # Since we do not use adaptive time stepping, we will only have a few
        # different factors.
        factorizations = Dict(factor => copy(factorization))
        W, factorization, factorizations
    end

    while t < last(tspan)
        dt = min(dt, last(tspan) - t)

        if t + dt >= last(tspan)
            last_step = true
        else
            last_step = false
        end

        # There are two possible formulations of a diagonally implicit RK method.
        # The "simple" one is
        #   y_i = q + h \sum_{j=1}^{i} a_{ij} f(y_j)
        # However, it can be better to use the smaller values
        #   z_i = (y_i - q) / h
        # so that the stage equations become
        #   q + h z_i = q + h \sum_{j=1}^{i} a_{ij} f(q + h z_j)
        # ⟺
        #   z_i - h a_{ii} f(q + h z_i) = \sum_{j=1}^{i-1} a_{ij} f(q + h z_j)
        # For a linear problem f(q) = T q, this becomes
        #   (I - h a_{ii} T z_i = a_{ii} T q + \sum_{j=1}^{i-1} a_{ij} T(q + h z_j)
        # We use this formulation and also avoid evaluating the stiff operator at
        # the numerical solutions (due to the large Lipschitz constant), but instead
        # rearrange the equation to obtain the required stiff RHS values as
        #   T(q + h z_i) = a_{ii}^{-1} (z_i - \sum_{j=1}^{i-1} a_{ij} f(q + h z_j))
        rhs_stiff!(k_stiff_q, q, parameters, t)

        # Compute stages
        for i in 1:s
            # RHS of linear system
            fill!(tmp, 0)
            for j in 1:(i - 1)
                @. tmp = tmp + A_stiff[i, j] * k_stiff[j] + A_nonstiff[i, j] * k_nonstiff[j]
            end
            # The right-hand side of the linear system formulated using the stages y_i
            # instead of the stage updates z_i would be
            #   @. tmp = q + dt * tmp
            # By using the stage updates z_i, we avoid the possibly different scales
            # for small dt.
            @. tmp = A_stiff[i, i] * k_stiff_q + tmp

            # Setup and solve linear system
            if iszero(rhs_stiff_operator) || iszero(A_stiff[i, i])
                copyto!(z, tmp)
            else
                factor = A_stiff[i, i] * dt

                F = let W = W, factor = factor,
                        factorization = factorization,
                        rhs_stiff_operator = rhs_stiff_operator
                    get!(factorizations, factor) do
                        fill!(W, 0)
                        W[diagind(W)] .= 1
                        @. W -= factor * rhs_stiff_operator
                        if issparse(W)
                            lu!(factorization, W)
                        else
                            factorization = lu!(W)
                        end
                        copy(factorization)
                    end
                end
                ldiv!(z, F, tmp)
            end
            # Compute new stage derivatives
            @. y = q + dt * z
            rhs_nonstiff!(k_nonstiff[i], y, parameters, t + c_nonstiff[i] * dt)
            if iszero(rhs_stiff_operator) || iszero(A_stiff[i, i])
                rhs_stiff!(k_stiff[i], y, parameters, t + c_stiff[i] * dt)
            else
                # The code below is equivalent to
                #   rhs_stiff!(k_stiff[i], y, parameters, t + c_stiff[i] * dt)
                # but avoids evaluating the stiff operator at the numerical solution.
                @. tmp = z
                for j in 1:(i-1)
                    @. tmp = tmp - A_stiff[i, j] * k_stiff[j] - A_nonstiff[i, j] * k_nonstiff[j]
                end
                @. k_stiff[i] = tmp / A_stiff[i, i]
            end

            # Store stage value in the last step
            if last_step
                stage_updates[i] = copy(z)
            end
        end

        # Update solution
        fill!(tmp, 0)
        for j in 1:s
            @. tmp += b_stiff[j] * k_stiff[j] + b_nonstiff[j] * k_nonstiff[j]
        end

        if relaxation # TODO? && (t + dt != last(tspan))
            @. y = dt * tmp # = qnew - q
            gamma = -2 * dot_entropy(q, y, parameters) / dot_entropy(y, y, parameters)
            @. q = q + gamma * y
            t += gamma * dt
        else
            @. q = q + dt * tmp
            t += dt
        end
        if save_everystep
            push!(sol_q, copy(q))
            append!(sol_t, t)
        end
        callback(q, parameters, t)

        if any(isnan, q)
            @error "NaNs in solution at time $t" q @__LINE__
            error()
        end
    end

    if save_everystep
        return (; u = sol_q,
                  t = sol_t,
                  stage_updates = stage_updates)
    else
        return (; u = (q0, q),
                  t = (first(tspan), t),
                  stage_updates = stage_updates)
    end
end



#####################################################################
# BBM discretization

struct BBM end
Base.Broadcast.broadcastable(equation::BBM) = Ref(equation)

function rhs_stiff!(du, u, equation::BBM, parameters, t)
    fill!(du, zero(eltype(du)))
    return nothing
end

operator(::typeof(rhs_stiff!), equation::BBM, parameters) = 0 * I

function rhs_nonstiff!(du, u, equation::BBM, parameters, t)
    (; D1, invImD2) = parameters
    tmp1 = get_tmp(parameters.tmp1, u)
    tmp2 = get_tmp(parameters.tmp2, u)
    one_third = one(eltype(u)) / 3

    # This semidiscretization conserves the linear and quadratic invariants
    @. tmp1 = -one_third * u^2
    mul!(tmp2, D1, tmp1)
    mul!(tmp1, D1, u)

    # There are two normalizations of the BBM equation:
    # 1. u_t - u_{txx} + u_x + u u_x = 0
    # @. tmp2 += -one_third * u * tmp1 - tmp1
    # 2. u_t - u_{txx} + u u_x = 0
    @. tmp2 += -one_third * u * tmp1

    if eltype(du) <: ForwardDiff.Dual
        # ForwardDiff does not support ldiv! for sparse matrices
        # TODO: This is not efficient for large systems
        du .= (I - Matrix(parameters.D2)) \ tmp2
    else
        ldiv!(du, invImD2, tmp2)
    end

    return nothing
end

function rhs_full!(du, u, equation::BBM, parameters, t)
    rhs_nonstiff!(du, u, equation, parameters, t)
    return nothing
end

# TODO: operator(::typeof(rhs_full!), equation::BBM, parameters) ?

function dot_entropy(u, v, equation::BBM, parameters)
    (; D1, D2, tmp1) = parameters
    mul!(tmp1, D2, v)
    half = one(eltype(u)) / 2
    @. tmp1 = half * u * (v - tmp1)
    return integrate(tmp1, D1)
end

function setup(q_func, equation::BBM, tspan, D, D2 = nothing)
    if D isa PeriodicUpwindOperators && isnothing(D2)
        D1 = D.central
        D2 = sparse(D.plus) * sparse(D.minus)
        invImD2 = lu(I - D2)
    elseif D isa FourierDerivativeOperator && isnothing(D2)
        D1 = D
        D2 = D1^2
        invImD2 = I - D2
    else
        throw(ArgumentError("Combination of operators not supported"))
    end

    x = SummationByPartsOperators.grid(D1)
    q0 = q_func(tspan[1], x, equation)
    tmp1 = similar(q0)
    tmp2 = DiffCache(similar(q0))
    parameters = (; equation, D1, D2, invImD2, tmp1, tmp2)
    return (; q0, parameters)
end


# Physical setup of a traveling wave solution with speed `c`
# and limit unity at infinity
solitary_wave_setup() = (xmin = -90.0, xmax =  90.0, c = 1.2)

function solitary_wave_solution(t, x::Number, equation::BBM)
    (; xmin, xmax, c) = solitary_wave_setup()

    A = 3 * (c - 1)
    K = 0.5 * sqrt(1 - 1 / c)
    x_t = mod(x - c * t - xmin, xmax - xmin) + xmin

    return 1 + A / cosh(K * x_t)^2
end

function solitary_wave_solution(t, x::AbstractVector, equation::BBM)
    solitary_wave_solution.(t, x, equation)
end



# TODO: development stuff - unused?
function bbm_setup(; domain_traversals = 10,
                     accuracy_order = 1)
    N = 2^8
    xmin = -90.0
    xmax =  90.0
    c = 1.2
    # upwind
    D = upwind_operators(periodic_derivative_operator;
                         derivative_order = 1, accuracy_order,
                         xmin, xmax, N)
    D1 = D.central
    D2 = sparse(D.plus) * sparse(D.minus)
    invImD2 = lu(I - D2)
    # invImD2 = cholesky(Symmetric(I - D2)) # no ldiv! CHOLMOD wrapper...

    tspan = (0.0, domain_traversals * (xmax - xmin) / c)
    x = SummationByPartsOperators.grid(D1)
    u0 = bbm_solution.(tspan[1], x)
    tmp1 = similar(u0)
    tmp2 = similar(u0)
    param = (; equation = BBM(), D1, D2, invImD2, tmp1, tmp2, bbm_solution)
    return ODEProblem(rhs_nonstiff!, u0, tspan, param)
end

function bbm_test(; tol = 1.0e-5)
    ode = bbm_setup(domain_traversals = 1.25)
    sol = solve(ode, Tsit5(); save_everystep = false, abstol = tol, reltol = tol)

    x = grid(ode.p.D1)
    fig = plot(x, sol.u[begin]; label = "u_0")
    plot!(fig, x, sol.u[end]; label = "u_end")
    return fig
end
#####################################################################

struct BBMH{T}
    c::T
    delta1::T
    delta2::T
    delta3::T
end
function BBMH(c, delta1 = 0.0, delta2 = 0.0, delta3 = 1.0)
    c, delta1, delta2, delta3 = promote(c, delta1, delta2, delta3)
    return BBMH{typeof(c)}(c, delta1, delta2, delta3)
end
Base.Broadcast.broadcastable(equation::BBMH) = (equation,)

function rhs_stiff!(dq, q, equation::BBMH, parameters, t)
    (; D1) = parameters
    N = size(D1, 2)

    du = view(dq, (0*N+1):(1*N))
    dv = view(dq, (1*N+1):(2*N))
    dw = view(dq, (2*N+1):(3*N))

    u = view(q, (0*N+1):(1*N))
    v = view(q, (1*N+1):(2*N))
    w = view(q, (2*N+1):(3*N))

    c = equation.c
    ε = 1 / sqrt(c)

    if D1 isa PeriodicUpwindOperators
        mul!(du, D1.plus, v, -(1 - equation.delta1 * ε))

        # dv .= -(1 - ε) / ε^2 * Du + 1 / ε^2 * w
        mul!(dv, D1.minus, u, -(1 - equation.delta2 * ε))
        @. dv = (dv + w) / ε^2

        # dw .= -ε^2 * (D1.central * w) - v
        mul!(dw, D1.central, w, -(1 - equation.delta3) * ε^2)
        @. dw -= v
    else
        mul!(du, D1, v)
        @. du *= -(1 - equation.delta1 * ε)

        # dv .= -(1 - ε) / ε^2 * Du + 1 / ε^2 * w
        mul!(dv, D1, u)
        @. dv = (-(1 - equation.delta2 * ε) * dv + w) / ε^2

        # dw .= -ε^2 * (D1 * w) - v
        mul!(dw, D1, w)
        @. dw = -(1 - equation.delta3) * ε^2 * dw - v
    end

    return nothing
end

function operator(::typeof(rhs_stiff!), equation::BBMH, parameters)
    D1 = parameters.D1
    c = equation.c
    ε = 1 / sqrt(c)

    if D1 isa PeriodicUpwindOperators
        D0 = sparse(D1.central)
        Dp = sparse(D1.plus)
        Dm = sparse(D1.minus)
        O = zero(D0)
        jac = [O -(1-equation.delta1*ε)*Dp O;
               -(1-equation.delta2*ε)/ε^2*Dm O 1/ε^2*I;
               O -I -((1-equation.delta3)*ε^2)*D0]
        dropzeros!(jac)
        return jac
    elseif D1 isa FourierDerivativeOperator
        D = Matrix(D1)
        O = zero(D)
        jac = [O -(1-equation.delta1*ε)*D O;
               -(1-equation.delta2*ε)/ε^2*D O 1/ε^2*I;
               O -I -((1-equation.delta3)*ε^2)*D]
        return jac
    else
        D = sparse(D1)
        O = zero(D)
        jac = [O -(1-equation.delta1*ε)*D O;
               -(1-equation.delta2*ε)/ε^2*D O 1/ε^2*I;
               O -I -((1-equation.delta3)*ε^2)*D]
        dropzeros!(jac)
        return jac
    end
end

function rhs_nonstiff!(dq, q, equation::BBMH, parameters, t)
    (; D1, tmp1) = parameters
    N = size(D1, 2)
    one_third = one(eltype(q)) / 3

    du = view(dq, (0*N+1):(1*N))
    dv = view(dq, (1*N+1):(2*N))
    dw = view(dq, (2*N+1):(3*N))

    u = view(q, (0*N+1):(1*N))
    v = view(q, (1*N+1):(2*N))
    w = view(q, (2*N+1):(3*N))

    c = equation.c
    ε = 1 / sqrt(c)
    if D1 isa PeriodicUpwindOperators
        # TODO: splitting 1
        # du_ode[:, 1] = -1 / 3 * (u .* (D1.central * u) + D1.central * (u.^2)) +
        #                 0.5 * sqrt(c) * (D1.central * u) - 0.5 * (D1.minus * v)
        # du_ode[:, 2] = -0.5 * c * (D1.plus * u) + 0.5 * sqrt(c) * (D1.central * v)
        # du_ode[:, 3] .= 0

        # splitting 3
        # du .= -1 / 3 * (u .* (D1.central * u) + D1.central * (u.^2)) - ε * (D1.minus * v)
        @. tmp1 = -one_third * u^2
        mul!(du, D1.central, tmp1)
        mul!(tmp1, D1.central, u)
        @. du += -one_third * u * tmp1
        mul!(tmp1, D1.plus, v, -equation.delta1 * ε)
        @. du += tmp1

        # dv .= -1 / ε * (D1.plus * u)
        mul!(dv, D1.minus, u, -equation.delta2 / ε)

        mul!(dw, D1.central, w, -equation.delta3 * ε^2)
    elseif D1 isa FourierDerivativeOperator
        # splitting 1
        # Du = D1 * u
        # Dv = D1 * v
        # du_ode[:, 1] = -1 / 3 * (u .* Du + D1 * (u.^2)) + 0.5 * sqrt(c) * Du - 0.5 * Dv
        # du_ode[:, 2] = -0.5 * c * Du + 0.5 * sqrt(c) * Dv
        # du_ode[:, 3] .= 0

        # splitting 3
        # du .= -1 / 3 * (u .* (D1 * u) + D1 * (u.^2)) - ε * (D1 * v)
        @. tmp1 = -one_third * u^2
        mul!(du, D1, tmp1)
        mul!(tmp1, D1, u)
        @. du += -one_third * u * tmp1
        mul!(tmp1, D1, v)
        @. du -= equation.delta1 * ε * tmp1

        # dv .= -1 / ε * (D1 * u)
        mul!(dv, D1, u)
        @. dv *= -equation.delta2 / ε

        mul!(dw, D1, w)
        @. dw *= -equation.delta3 * ε^2
    else
        # splitting 1
        # Du = D1 * u
        # Dv = D1 * v
        # du_ode[:, 1] = -1 / 3 * (u .* Du + D1 * (u.^2)) + 0.5 * sqrt(c) * Du - 0.5 * Dv
        # du_ode[:, 2] = -0.5 * c * Du + 0.5 * sqrt(c) * Dv
        # du_ode[:, 3] .= 0

        # splitting 3
        # du .= -1 / 3 * (u .* (D1 * u) + D1 * (u.^2)) - ε * (D1 * v)
        @. tmp1 = -one_third * u^2
        mul!(du, D1, tmp1)
        mul!(tmp1, D1, u)
        @. du += -one_third * u * tmp1
        mul!(tmp1, D1, v)
        @. du -= equation.delta1 * ε * tmp1

        # dv .= -1 / ε * (D1 * u)
        mul!(dv, D1, u)
        @. dv *= -equation.delta2 / ε

        mul!(dw, D1.central, w)
        @. dw *= -equation.delta3 * ε^2
    end

    return nothing
end

function rhs_full!(dq, q, equation::BBMH, parameters, t)
    # (; D1, tmp1) = parameters
    D1 = parameters.D1
    N = size(D1, 2)
    one_third = one(eltype(q)) / 3

    du = view(dq, (0*N+1):(1*N))
    dv = view(dq, (1*N+1):(2*N))
    dw = view(dq, (2*N+1):(3*N))

    u = view(q, (0*N+1):(1*N))
    v = view(q, (1*N+1):(2*N))
    w = view(q, (2*N+1):(3*N))

    # to handle auto-diff etc.
    # if eltype(tmp1) != eltype(q)
    #     tmp1 = similar(u)
    # end
    tmp1 = zero(u)

    c = equation.c

    if D1 isa PeriodicUpwindOperators
        # du = -1/3 * (u * (D1.central * u) + D1.central * (u^2))
        #      - (D1.minus * v)
        @. tmp1 = -one_third * u^2
        mul!(du, D1.central, tmp1)
        mul!(tmp1, D1.central, u)
        @. du += -one_third * u * tmp1
        mul!(tmp1, D1.minus, v)
        @. du -= tmp1

        # dv .= c * w - c * Du
        mul!(dv, D1.plus, u, -c)
        @. dv += c * w

        # dw = -v - inv(c) * (D1.central * w)
        mul!(dw, D1.central, w, -1 / c)
        @. dw -= v
    else
        # du .= -1 / 3 * (u .* (D1 * u) + D1 * (u.^2)) - (D1 * v)
        @. tmp1 = -one_third * u^2
        mul!(du, D1, tmp1)
        mul!(tmp1, D1, u)
        @. du += -one_third * u * tmp1
        mul!(tmp1, D1, v)
        @. du -= tmp1

        # dv .= c * w - c * Du
        mul!(dv, D1, u)
        @. dv = c * w - c * dv

        # dw = -v - inv(c) * (D1 * w)
        mul!(dw, D1, w)
        inv_c = inv(c)
        @. dw = -v - inv_c * dw
    end

    return nothing
end

# TODO: operator(::typeof(rhs_full!), equation::BBMH, parameters) ?


function dot_entropy(q1, q2, equation::BBMH, parameters)
    (; D1, tmp1) = parameters
    N = size(D1, 2)

    u1 = view(q1, (0*N+1):(1*N))
    v1 = view(q1, (1*N+1):(2*N))
    w1 = view(q1, (2*N+1):(3*N))

    u2 = view(q2, (0*N+1):(1*N))
    v2 = view(q2, (1*N+1):(2*N))
    w2 = view(q2, (2*N+1):(3*N))

    c = equation.c
    half = one(c) / 2
    @. tmp1 = half * (u1 * u2 + v1 * v2 / c + w1 * w2)

    return integrate(tmp1, D1)
end


function setup(q_func, equation::BBMH, tspan, D1, D2 = nothing; vzero_bool = false)
    if !isnothing(D2)
        throw(ArgumentError("Combination of operators not supported"))
    end

    x = SummationByPartsOperators.grid(D1)
    c = solitary_wave_setup().c
    q0 = q_func(tspan[1], x, equation)

    # The results are better when the discrete derivative operator is used
    # to compute the initial condition for w, in particular for low-order
    # discretizations
    if D1 isa PeriodicUpwindOperators
        u0 = view(q0, 1:length(x))
        v0 = view(q0, (1*length(x)+1):(2*length(x)))
        w0 = view(q0, (2*length(x)+1):(3*length(x)))
        mul!(w0, D1.central, u0)

        if vzero_bool
            v0 .= 0.0
        else
            tmpv = similar(u0)
            mul!(tmpv, D1.central, u0)
            mul!(v0, D1.central, tmpv)
            v0 .*= c

        end
        # TODO: v0 .= rand(length(w0))
    else
        u0 = view(q0, 1:length(x))
        v0 = view(q0, (1*length(x)+1):(2*length(x)))
        w0 = view(q0, (2*length(x)+1):(3*length(x)))
        mul!(w0, D1, u0)

        if vzero_bool
            v0 .= 0.0
        else
            tmpv = similar(u0)
            mul!(tmpv, D1.central, u0)
            mul!(v0, D1.central, tmpv)
            v0 .*= c
        end
    end

    # TODO: introduce caches for BBMH
    tmp1 = similar(u0)

    parameters = (; equation, D1, tmp1)
    return (; q0, parameters)
end


# Physical setup of a traveling wave solution with speed `c`
# and limit unity at infinity
function solitary_wave_solution(t, x::Number, equation::BBMH)
    # Use BBM solution and its derivatives as reference for now
    (; xmin, xmax, c) = solitary_wave_setup()

    A = 3 * (c - 1)
    K = 0.5 * sqrt(1 - 1 / c)
    x_t = mod(x - c * t - xmin, xmax - xmin) + xmin

    u = 1 + A / cosh(K * x_t)^2 # u
    v = 2 * A * c * K^2 * (1 - 2 * sinh(K * x_t)^2) / (cosh(K * x_t))^4 #u_tx
    w = -2 * A * K * (sinh(K * x_t) / (cosh(K * x_t))^3) # u_x

    return (u, v, w)
end

function solitary_wave_solution(t, x::AbstractVector, equation::BBMH)
    uvw = solitary_wave_solution.(t, x, equation)
    return vcat(map(uvw -> uvw[1], uvw),
                map(uvw -> uvw[2], uvw),
                map(uvw -> uvw[3], uvw))
end

#####################################################################
# BBMH reference solution
function N_rhs(phi)
    return 0.5 .* phi.^2
end
function BBMH_petviashvili(c, D; eps_sq)
    #Construct the matrix L^{-1} that will be used in the Petiashvili algorithm
    L = (c - 1) * D^0 - inv(I + c * eps_sq * (c - eps_sq) * D^2) * ((c - eps_sq) * D^2)
    #Prepare Petiashvili
    x = SummationByPartsOperators.grid(D)
    phi_0 = solitary_wave_solution.(0.0, x, BBM()) #initial_guess
    phi_current = copy(phi_0)
    res = norm(L * phi_current - N_rhs.(phi_current), Inf)
    gamma = 2.0 # factor of the Petviashvili iteration

    tol = 1e-12 #tolerance for the P-Algorithm
    iter = 0

    while res > tol
        if iter > 1e5
            @info res M_n = integrate(phi_current .* (L * phi_current), D) / integrate(phi_current .* N_rhs.(phi_current), D)
            @info "terminating since no convergence"
            break
        end
        #compute the factor M_n
        M_n = integrate(phi_current .* (L * phi_current), D) / integrate(phi_current .* N_rhs.(phi_current), D)
        phi_current .= M_n ^ gamma * (L \ N_rhs.(phi_current))
        res = norm(L * phi_current - N_rhs.(phi_current), Inf)
        iter += 1
    end
    #compute the different components of the bbmh solution
    phi_u = phi_current
    phi_w = D * phi_u + c * eps_sq * (D * (-c * phi_u + phi_u + phi_u.^2 / 2))
    #phi_v = -(inv(I + c * eps_sq * (c - eps_sq) * D^2) * ((c - eps_sq) * D^2) * phi_u)

    phi_v = (c-eps_sq) * D * phi_w
    #phi_v = -(c * D * phi_w - eps_sq * D * phi_w)
    #phi_v =  (-c * phi_u + phi_u + phi_u.^2 / 2)
    #compute the error with regards to the bbm solution
    #@info "errors" norm(phi_u - bbm_solution.(0,x), Inf) norm(phi_w - D * bbm_solution.(0,x), Inf) norm(phi_v + D * bbm_solution_t.(0,x), Inf)
    return (phi_u .+ 1, phi_v, phi_w)
  end

  function BBMH_interpolation(c, D;eps_sq = 1e-12)
    phi_u, phi_v, phi_w = BBMH_petviashvili(c, D;eps_sq)

    x = SummationByPartsOperators.grid(D)
    nodes = (x,)
    #interpolation of the petviashili
    u_itp = cubic_spline_interpolation(nodes, phi_u, extrapolation_bc = Periodic())
    v_itp = cubic_spline_interpolation(nodes, phi_v, extrapolation_bc = Periodic())
    w_itp = cubic_spline_interpolation(nodes, phi_w, extrapolation_bc = Periodic())

    return (u_itp, v_itp, w_itp)
  end

#####################################################################
# Solitary wave of the BBM equation
function plot_bbm_solitary_wave(; c = 100.0,
                                  domain_traversals = 1.25,
                                  accuracy_order = 1, N = 2^8,
                                  alg, dt, kwargs...)
    # Initialization of physical and numerical parameters
    (; xmin, xmax) = solitary_wave_setup()
    tspan = let c = solitary_wave_setup().c
        (0.0, domain_traversals * (xmax - xmin) / c)
    end

    D1 = upwind_operators(periodic_derivative_operator;
                          derivative_order = 1, accuracy_order,
                          xmin, xmax, N)


    # Setup plot
    fig_sol = plot(; xguide = L"x", yguide = L"u")


    # BBM with explicit part
    (; q0, parameters) = setup(solitary_wave_solution, BBM(), tspan, D1)
    @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters), rhs_nonstiff!,
                           q0, tspan, parameters, alg;
                           dt = dt, kwargs...)

    x = grid(parameters.D1)
    plot!(fig_sol, x, sol.u[begin]; label = L"u^0", plot_kwargs()...)
    plot!(fig_sol, x, sol.u[end]; label = "BBM", plot_kwargs()...)


    # BBMH with IMEX
    (; q0, parameters) = setup(solitary_wave_solution, BBMH(c, true), tspan, D1)
    @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters), rhs_nonstiff!,
                           q0, tspan, parameters, alg;
                           dt = dt, kwargs...)

    x = grid(parameters.D1)
    plot!(fig_sol, x, sol.u[end][1:length(x)]; label = "BBMH", plot_kwargs()...)
    plot!(fig_sol; legend = :topleft)


    return fig_sol
end
###################################################################################
#AP tests
function AP_test_in_eps(; delta1 = 0, delta2 = 0, delta3 = 1,
                                   domain_traversals = 0.13, dt = 0.01, latex = false, alg = ARS443(), N = 2^9, vzero_boolval = false, kwargs...)
    #test for different epsilon the AP property of our schemes
    #first setup the problem
    xmin, xmax, c = solitary_wave_setup()
    epss = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
    #epss = [1e-10]
    #In this case we discretize by a periodic derivative operator
    D1 = upwind_operators(periodic_derivative_operator;
                        derivative_order = 1, accuracy_order = 12,
                        xmin, xmax, N)

    #generate the BBM solution
    tspan = let c = solitary_wave_setup().c
        (0.0, domain_traversals * (xmax - xmin) / c)
    end
    errors_u = []
    errors_v = []
    errors_w = []

    u_limit, u_ini, u_limit_stage_updates = let equation = BBM()
        (; q0, parameters) = setup(solitary_wave_solution,
                                   equation, tspan, D1)
        @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters),
                                             rhs_nonstiff!,
                                             q0, tspan, parameters, alg;
                                             dt = dt, kwargs...)
        get_u(sol.u[end], equation), get_u(sol.u[begin], equation), sol.stage_updates
    end

    for eps in epss
        equation = BBMH( 1 / eps^2)
        (; q0, parameters) = setup(solitary_wave_solution,
                                   equation, tspan, D1; vzero_bool = vzero_boolval)

        @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters),
                                             rhs_nonstiff!,
                                             q0, tspan, parameters, alg;
                                             dt = dt, kwargs...)

        u = get_u(sol.u[end], equation)
        error_u = integrate(abs2, u - u_limit, parameters.D1) |> sqrt
        push!(errors_u, error_u)


        v_limit = zero(u_limit)
        a = let
            A_stiff, _ = coefficients(alg)
            if A_stiff[1,1] == 0
                Ainv = inv(A_stiff[2:end, 2:end])
            else
                Ainv = inv(A_stiff)
            end
            Ainv[end, :]
        end

        A_stiff, _ = coefficients(alg)
        if A_stiff[1,1] == 0
            for i in eachindex(a)
                v_limit .-= a[i] .* (parameters.D1.central * u_limit_stage_updates[i + 1])
            end
        else
            for i in eachindex(a)
                v_limit .-= a[i] .* (parameters.D1.central * u_limit_stage_updates[i])
            end
        end

        v = get_qi(sol.u[end], equation, 1)
        error_v = integrate(abs2, v - v_limit, parameters.D1) |> sqrt
        push!(errors_v, error_v)

        w_limit = parameters.D1.minus * u_limit
        w = get_qi(sol.u[end], equation, 2)
        error_w = integrate(abs2, w - w_limit, parameters.D1) |> sqrt
        push!(errors_w, error_w)
    end

    eoc_u =  compute_eoc(1 ./ epss.^2, errors_u)
    eoc_v = compute_eoc(1 ./ epss.^2, errors_v)
    eoc_w = compute_eoc(1 ./ epss.^2, errors_w)

    data = hcat(epss.^2, errors_u, eoc_u, errors_v, eoc_v, errors_w, eoc_w)
    header = [L"\varepsilon^2", "L2 error u", "L2 EOC u", "L2 error v", "L2 EOC v", "L2 error w", "L2 EOC w"]
    kwargs = (; header, formatters=(ft_printf("%.2e", [1, 2]),
                                        ft_printf("%.2f", [3]), ft_printf("%.2e", [4]), ft_printf("%.2f", [5]),ft_printf("%.2e", [6]), ft_printf("%.2f", [7])))
    pretty_table(data; kwargs...)
    if latex
        pretty_table(data; kwargs..., backend=Val(:latex))
    end
    return nothing
end

###################################################################################
#Petviashvili tests
function prep_data(sol)
    #prepare the solution data in such a way that the first component only contains
    # u vals 2nd v vals and 3rd w vals
    u = map(x -> x[1], sol)
    v = map(x -> x[2], sol)
    w = map(x -> x[3], sol)
    return [u,v,w]
end

function petviashvili_plot_validation()

    figdir = joinpath(@__DIR__, "PLots Petviashvili")
    isdir(figdir) || mkdir(figdir)
    N = 2^10 #grid points
    xmin, xmax, c = solitary_wave_setup()

    D = fourier_derivative_operator(xmin, xmax, N)
    x = SummationByPartsOperators.grid(D)

    sol_bbm = solitary_wave_solution.(0, x, BBM())

    epss = [0.5, 0.05]
    colors = [:green, :red]


    fig_u = plot(x, sol_bbm; color = :blue, label = "Analytical", xlabel = L"x", ylabel = L"u", xlim = (-20.0, 20.0), size = (800, 350), plot_kwargs()...)

    for i in 1:length(epss)

        u_itp, v_itp, w_itp = BBMH_interpolation(c, D; eps_sq = epss[i])

        plot!(fig_u, x, u_itp.(x); label = "Petviashvili \$\\varepsilon = $(epss[i])\$", linestyle = :dot, color = colors[i], plot_kwargs()...)


    end

    filename = joinpath(figdir,"Petviashvili_validation_u.png")
    savefig(fig_u, filename)
    @info "Saved figure" filename

    return nothing
end

function setup_pv(u0,v0,w0, equation::BBMH, D1)
    x = SummationByPartsOperators.grid(D1)
    q0 = [u0.(x); v0.(x); w0.(x)]
    tmp1 = similar(u0.(x))
    tmp2 = similar(q0)
    parameters = (;equation, D1, tmp1, tmp2, u0, v0, w0)
    return (;q0, parameters)
end

function setup_pv_bbm(u0,v0,w0, equation::BBMH, D1)
    x = SummationByPartsOperators.grid(D1)
    q0 = [u0; v0; w0]
    tmp1 = similar(u0)
    tmp2 = similar(q0)
    parameters = (;equation, D1, tmp1, tmp2)
    return (;q0, parameters)
end

function get_u(data, equation)
    if equation isa BBM
        return data
    else
        #in this case its the BBMH equation
        u = view(data, 1:Int(length(data) / 3))
        return u
    end
end

function get_qi(data, equation, index)
    n = Int(length(data) / 3)
    q_sub = view(data, (index * n + 1) : (index + 1) * n )

    return q_sub
end

function petviashvili_error_growth(; delta1 = 0, delta2 = 0, delta3 = 1,
                                     domain_traversals = 7.14, eps = 1.e-3, dt = 0.5, alg = ARS443(),  order_test = 4, N_test = 2^8, order_curves = true, legend_val =:topleft, kwargs...)
    figdir = joinpath(@__DIR__, "Plots_Error_Growth")
    isdir(figdir) || mkdir(figdir)
    #setting up the problem
    xmin, xmax, c = solitary_wave_setup()
    #setting up testvals
    N_ref = 2^12

    D_ref = fourier_derivative_operator(xmin, xmax, N_ref)
    D_test = upwind_operators(periodic_derivative_operator;
                        derivative_order = 1, accuracy_order = order_test,
                        xmin, xmax, N = N_test)

    #create the reference solution of BBMH using Petviashvili
    #solve the BBMH system with the IMEX solver
    #the initial condition of the BBMH system is given
    tspan = let c = solitary_wave_setup().c
        (0.0, domain_traversals * (xmax - xmin) / c)
    end

    # Setup plot
    fig_err = plot(; xguide = L"t", yguide = L"Error of $q$", legend = legend_val, yscale = :log10, xscale =:log10)

        #IC
        tend = tspan[2]
        u_itp, v_itp, w_itp = BBMH_interpolation(c, D_ref; eps_sq = eps^2)
        (; q0, parameters) = setup_pv(u_itp, v_itp, w_itp, BBMH(1 / eps^2, delta1, delta2, delta3), D_test)
        #Setup callback computing the error
        series_t = Vector{Float64}()
        series_error = Vector{Float64}()
        series_error_u = Vector{Float64}()
        series_error_v = Vector{Float64}()
        series_error_w = Vector{Float64}()

        callback = let series_t = series_t, series_error_u = series_error_u, series_error_v = series_error_v, series_error_w = series_error_w, series_error = series_error
            function (q, parameters, t)
                    (;tmp2, equation, u0, v0, w0) = parameters
                    x = grid(parameters.D1)
                    x_t = mod.(x .- c * t .- xmin, xmax - xmin) .+ xmin
                    D1 = parameters.D1

                    u_ref = [u0.(x_t); v0(x_t); w0.(x_t)]

                    u = q
                    N = length(x)

                    @. tmp2 = u - u_ref

                    tmp2_u = view(tmp2, 1: N)
                    tmp2_v = view(tmp2, N + 1: 2 * N)
                    tmp2_w = view(tmp2, 2 * N + 1: 3 * N)

                    err_v = integrate(abs2, tmp2_v, parameters.D1) |> sqrt
                    err_u = integrate(abs2, tmp2_u, parameters.D1) |> sqrt
                    err_w = integrate(abs2, tmp2_w, parameters.D1) |> sqrt
                    err_total = integrate(abs2, tmp2, parameters.D1) |> sqrt

                    push!(series_t, t)
                    push!(series_error_u, err_u)
                    push!(series_error_v, err_v)
                    push!(series_error_w, err_w)
                    push!(series_error, err_total)
                return nothing
            end
        end

        #Without relaxation
        @info "without relaxation"
        empty!(series_t)
        empty!(series_error)
        empty!(series_error_u)
        empty!(series_error_v)
        empty!(series_error_w)
        @show parameters.equation
        @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters),
                                rhs_nonstiff!,
                                q0, tspan, parameters, alg;
                                dt = dt, callback, relaxation = false, kwargs...)

        plot!(fig_err, series_t, series_error; label = "baseline", color = :blue, plot_kwargs()...)
        #order curves
        if order_curves
            p = 2.0
            t_start = 1e2

            # Filter time points and find corresponding errors
            idx_start = findfirst(t -> t ≥ t_start, series_t)
            t_ref = series_t[end]
            #@info series_error
            err_ref = series_error[end]

            ref_curve = [err_ref * (t / t_ref)^p for t in [series_t[idx_start], series_t[end]]]
            plot!(fig_err, [series_t[idx_start], series_t[end]], ref_curve;
                    label = L"O(t^2)", color = :gray, marker = :square)

        end


        #With relaxation
        @info "with relaxation"
        empty!(series_t)
        empty!(series_error)
        empty!(series_error_u)
        empty!(series_error_v)
        empty!(series_error_w)
        @time sol_rel = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters),
                                rhs_nonstiff!,
                                q0, tspan, parameters, alg;
                                dt = dt, callback, relaxation = true, kwargs...)

        plot!(fig_err, series_t, series_error; label = "relaxation", color = :red, plot_kwargs()...)

    #order curves
    if order_curves
        p = 1.0
        t_start = 1e2

        # Filter time points and find corresponding errors
        idx_start = findfirst(t -> t ≥ t_start, series_t)
        t_ref = series_t[end]
        err_ref = series_error[end]

        ref_curve = [err_ref * (t / t_ref)^p for t in [series_t[idx_start], series_t[end]]]
        plot!(fig_err, [series_t[idx_start], series_t[end]], ref_curve;
                label = L"O(t^1)", color = :gray, marker = :circle)

    end

    #saving the plot
    filename = joinpath(figdir,"Error_Growth__BBMH_$(alg)_$(eps)_order$(order_test)_dt($dt).pdf")
    savefig(fig_err, filename)
    @info "Saved figure" filename

    return fig_err
end

function solitary_wave_error_growth(;eps, accuracy_order = 6, alg = ARS443(), dt = 0.5, domain_traversals = 10, N = 2^8, legend_val = :topleft,
                                        delta1 = 0, delta2 = 0, delta3 = 1, kwargs...)

    figdir = joinpath(@__DIR__, "Plots_Error_Growth")
    isdir(figdir) || mkdir(figdir)
    # Initialization of physical and numerical parameters
    (; xmin, xmax) = solitary_wave_setup()
    tspan = let c = solitary_wave_setup().c
        (0.0, domain_traversals * (xmax - xmin) / c)
    end

    D1 = upwind_operators(periodic_derivative_operator;
                          derivative_order = 1, accuracy_order,
                          xmin, xmax, N)


    # Setup plot
    fig_err = plot(; xguide = L"t", yguide = L"Error of $u$", legend = legend_val)


    # Setup callback computing the error
    series_t = Vector{Float64}()
    series_error = Vector{Float64}()
    callback = let series_t = series_t, series_error = series_error
        function (q, parameters, t)
            (; tmp1, equation) = parameters

            if equation isa BBM
                u = q
            else
                u = view(q, 1:(length(q) ÷ 3))
            end
            u_ref = solitary_wave_solution(t, grid(parameters.D1), BBM())

            @. tmp1 = u - u_ref
            err = integrate(abs2, tmp1, parameters.D1) |> sqrt

            push!(series_t, t)
            push!(series_error, err)
            return nothing
        end
    end

    @info "BBM without relaxation"
    empty!(series_t)
    empty!(series_error)
    (; q0, parameters) = setup(solitary_wave_solution,
                                   BBM(), tspan, D1)

    @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters), rhs_nonstiff!,
                            q0, tspan, parameters, alg;
                            dt, callback, kwargs...)

    plot!(fig_err, series_t, series_error; label = "BBM baseline", plot_kwargs()...)

    #order curves
    p = 2.0
    t_start = 1e2

    # Filter time points and find corresponding errors
    idx_start = findfirst(t -> t ≥ t_start, series_t)
    t_ref = series_t[end]

    err_ref = series_error[end]

    ref_curve = [err_ref * (t / t_ref)^p for t in [series_t[idx_start], series_t[end]]]
    plot!(fig_err, [series_t[idx_start], series_t[end]], ref_curve;
                label = L"O(t^2)", color = :gray, marker = :square)

    savefig(fig_err, joinpath(figdir, "Error_Growth__BBM_$(alg)_$(eps)_order$(accuracy_order)_dt($dt).pdf"))

    @info "BBM with relaxation"
    empty!(series_t)
    empty!(series_error)
    (; q0, parameters) = setup(solitary_wave_solution,
                                   BBM(), tspan, D1)

    @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters), rhs_nonstiff!,
                            q0, tspan, parameters, alg;
                            dt, callback, relaxation = true, kwargs...)


    plot!(fig_err, series_t, series_error; label = "BBM relaxation", plot_kwargs()...)

    p = 1.0
    t_start = 1e2

    # Filter time points and find corresponding errors
    idx_start = findfirst(t -> t ≥ t_start, series_t)
    t_ref = series_t[end]
    err_ref = series_error[end]

    ref_curve = [err_ref * (t / t_ref)^p for t in [series_t[idx_start], series_t[end]]]
    plot!(fig_err, [series_t[idx_start], series_t[end]], ref_curve;
                label = L"O(t^1)", color = :gray, marker = :circle)

    # BBMH with IMEX
    (; q0, parameters) = setup(solitary_wave_solution,
                                BBMH(1 / eps^2, delta1, delta2, delta3),
                                tspan, D1)

    @info "BBMH without relaxation"
    empty!(series_t)
    empty!(series_error)
    @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters), rhs_nonstiff!,
                            q0, tspan, parameters, alg;
                            dt, callback, kwargs...)
    plot!(fig_err, series_t, series_error; label = "BBMH baseline", plot_kwargs()...)


    @info "BBMH with relaxation"
    empty!(series_t)
    empty!(series_error)
    @time sol = solve_imex(rhs_stiff!, operator(rhs_stiff!, parameters), rhs_nonstiff!,
                            q0, tspan, parameters, alg;
                            dt, callback, relaxation = true, kwargs...)
    plot!(fig_err, series_t, series_error; label = "BBMH relaxation", plot_kwargs()...)
    plot!(fig_err; xscale = :log10, yscale = :log10, plot_kwargs()...)

    filename = joinpath(figdir, "Error_Growth__BBM_$(alg)_$(eps)_order$(accuracy_order)_dt($dt).pdf")
    savefig(fig_err, filename)
    @info "Saved figure" filename

    return fig_err
end