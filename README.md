## KernelStat

A module for the [Julia](http://julialang.org/) language for [kernel-function based statistics](http://en.wikipedia.org/wiki/Kernel_%28statistics%29)


## Main features

This module allows Julia users to specify and evaluate a range of different kernel functions. It also includes functionality for estimating bandwidths for use with various staistical procedures, including, but not limited to, estimation of optimal block lengths for dependent bootstraps.

The module is implemented entirely in Julia and draws some inspiration from the [Distance](https://github.com/JuliaStats/Distance.jl) package. Specifically, each kernel function or bandwidth estimation method is its own type and are evaluated via the functions `evaluate` and `bandwidth` respectively, exploiting multiple dispatch.

It is worthing noting that kernel-function based statistics is a big topic, and this module is currently quite small. In fact, I only wrote this module to support a more comprehensive module for dependent bootstrapping, and so the functionality implemented thus far is all that was needed to support the package [DependentBootstrap](https://github.com/colintbowers/DependentBootstrap.jl).


## How to use KernelStat

#### Installation

This package is not yet registered, so you will not be able to get it using `Pkg.add("KernelStat")`. If you want to try it out use `Pkg.clone("https://github.com/colintbowers/KernelStat.jl")` from the Julia REPL. This should download it and install it in your `.julia` folder. Then you can load the package into the current julia session with `using KernelStat`.


#### Kernel function types and evaluation

As in the Julia Distance package, each kernel function is its own type, and all kernel function types are a subset of the abstract type `KernelFunction`. A list of kernel function types and their fields follows. Kernel functions defined on the interval [-1, 1] include:

* `KernelUniform`
* `KernelTriangular`
* `KernelEpanechnikov`
* `KernelQuartic`

Kernel functions defined on the interval [-Inf, Inf] include:

* `KernelGaussian`. Fields include:
  * p1::Float64 (scaling term almost always set to 1 / sqrt(2*pi), which is the value employed by the constructor `KernelGaussian()`)

Specific kernel functions from academic papers include:

* `KernelPR1993FlatTop`, from Politis, Romano (1993) "On a family of smoothing kernels of infinite order". Fields (defined in that paper) include:
  * m::Float64 (constrained to m > 0 and m < M)
  * M::Float64 (constrained to M > 0 and M > m)
* `KernelP2003FlatTop`, from Politis (2003) "Adaptive bandwidth choice". Note, this kernel function is equivalent to `KernelPR1993FlatTop` with `m=0.5` and `M=1`.
* `KernelPP2002Trap`, from Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences". Fields (defined in that paper) include:
  * p1::Float64 (constrained to the interval (0, 0.5])
* `KernelPP2002Smooth`, from Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences". Fields (defined in that paper) include:
  * p1::Float64 (constrained to the interval [1, Inf))

Every kernel function exhibits a constructor that requires no inputs. These constructors use sensible default values for the kernel parmeters, e.g. `KernelPR1993FlatTop() = KernelPR1993FlatTop(0.5, 1)`. 

Every `KernelFunction` exhibits the following very important function:

* evaluate{T<:Number}(x::T, kF::KernelFunction) -> returns the value of the output of the kernel function for the given input `x`.
* evaluate{T<:Number}(x::Vector{T}, kF::KernelFunction) -> a wrapper to allow a vector of inputs.

For example, if a user wishes to evaluate the uniform kernel at the number 0.3, they would use:

    evaluate(0.3, KernelUniform())

If a user wishes to evaluate the flat top kernel of Politis and Romano (1993) with parameters 1.5 and 2.5, over the integers 1, 2, 3, and 4, they would use:

    evaluate([1:4], KernelPR1993FlatTop(1.5, 2.5))

In addition to `evaluate`, every `KernelFunction` exhibits the following core methods:

* string(kF::KernelFunction) -> returns an `ASCIIString` representation of the kernel function, e.g. `string(::KernelUniform)` returns `"uniform"`
* copy(kF::KernelFunction) -> standard copy
* show(io::IO, kF::KernelFunction) -> print the name of the kernel function and list its parameter values
* show(kF::KernelFunction) -> wrapper on show(io::IO, kF::KernelFunction) that prints to `STDOUT`
* toKernel(x::ASCIIString) -> attempts to convert an `ASCIIString` to a subtype of `KernelFunction`, e.g. `toKernelFunction("uniform")` returns `KernelUniform()`. Note, this function is essentially the inverse of `string`. By construction it is *not type stable*.
* nonzeroDomain(kF::KernelFunction) -> returns a `UnitRange` that indicates the domain over which the input kernel function evaluates to a non-zero number.
* paramDomain(kF::KernelFunction, pNum::Int) -> returns a `UnitRange` indicating the valid domain of the parameter number indicated in `pNum`. Parameter numbers correspond to the order in which fields are defined for the `KernelFunction`. If a `KernelFunction` does not have any fields/parameters or if `pNum` does not correspond to a field number, then this function throws an error. 


#### Bandwidth estimation

Each method of bandwidth estimation is its own type, and they are all subtypes of the abstract type `BandwidthMethod`. Common to every bandwidth estimation type is a field called `adjustmentTerm`, which simply scales the output of bandwidth estimation by multiplication. Thus the default value to this field is `1.0` and most users will not want to alter this. Some bandwidth estimation types contain additional fields. A complete list follows

* `BandwidthWhiteNoise`. Estimate bandwidth using a white noise assumption to determine confidence bounds for autocorrelation. Bandwidth is two times the lag number of the first insignificant autocorrelation. Fields include:
  * adjustmentTerm::Float64
* `BandwidthBartlett`. Estimate bandwidth using Bartlett assumptions (ie underlying model is Moving Average (MA)) to determine confidence bounds at each step. Bandwidth is two times the lag number of the first insignificant autocorrelation. Fields include:
  * adjustmentTerm::Float64
* `BandwidthP2003`. Estimate bandwidth using the empirical rule of picking M defined in section 2.1 of Politis (2003) "Adaptive Bandwidth Choice". Fields include:
  * adjustmentTerm::Float64
  * c::Float64 (parameter from Politis (2003). Default value of 2.0 recommended in that paper)
  * K::Int (parameter from Politis (2003). Default value of 5 recommended in that paper, although the constructor `BandwidthP2003(numObs::Int)` will transform `numObs` into a better estimate for `K` than the default rule - see Politis (2003) for more detail)

Given a vector of observable data `x::Vector{T}`, where `T<:Number`, an appropriate bandwidth for `x` can be estimated using the `bandwidth` function. This function exhibits methods of the type:

    bandwidth(x, ::BandwidthMethod)

For example, to estimate the bandwidth of `x` using the Bartlett method, one would use:

    M = bandwidth(x, BandwidthBartlett())

A second example. To estimate the bandwidth of `x` using the adaptive method of Politis (2003) with default parameter value for `c`, and `K` determined optimally from the number of observations in `x`, one would use:

    M = bandwidth(x, BandwidthP2003(length(x)))

The subtypes of `BandwidthMethod` also exhibit the following core functions:

* string(b::BandwidthMethod) -> returns an `ASCIIString` representation of the bandwidth method, e.g. `string(::BandwidthP2003)` returns `"P2003"`
* copy(b::BandwidthMethod) -> standard copy
* show(io::IO, b::BandwidthMethod) -> print the name of the bandwidth estimation method and list its parameter values
* show(b::BandwidthMethod) -> wrapper on show(io::IO, b::BandwidthMethod) that prints to `STDOUT`
* toBandwidthMethod(x::ASCIIString) -> attempts to convert an `ASCIIString` to a subtype of `BandwidthMethod`, e.g. `toBandwidthMethod("P2003")` returns `BandwidthP2003()`.

This concludes the description of the KernelStat module.


