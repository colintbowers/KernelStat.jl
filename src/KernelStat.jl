module KernelStat
#-----------------------------------------------------------
#PURPOSE
#	Colin T. Bowers module for kernel functions and statistics
#NOTES
#	Currently this module only really contains a function for estimating bandwidths (needed for estimating block lengths in dependent bootstraps)
#LICENSE
#	MIT License (see github repository for more detail: https://github.com/colintbowers/KernelStat.jl.git)
#-----------------------------------------------------------


#Load any entire modules that are needed (use import ModuleName1, ModuleName2, etc)
using 	StatsBase

#Load any specific variables/functions that are needed (use import ModuleName1.FunctionName1, ModuleName2.FunctionName2, etc)
import 	Base.string,
		Base.show,
		Base.copy

#Specify the variables/functions to export (use export FunctionName1, FunctionName2, etc)
export	evaluate, #Evaluate kernel function at given input
		nonzeroDomain, #Function that returns a UnitRange indicating the range over which kernel output is non-zero
		paramDomain, #Function that returns a UnitRange indicating the feasible domain kernel parameters
		KernelFunction, #Abstract type that nests all kernel function types
		KernelUniform, #Kernel function type
		KernelTriangular, #Kernel function type
		KernelEpanechnikov, #Kernel function type
		KernelQuartic, #Kernel function type
		KernelGaussian, #Kernel function type
		KernelPR1993FlatTop, #Kernel function type
		KernelP2003FlatTop, #Kernel function type
		KernelPP2002Trap, #Kernel function type
		KernelPP2002Smooth, #Kernel function type
		bandwidth, #Function for estimating bandwidth
		BandwidthMethod, #Abstract type that nests all bandwidth method types
		BandwidthWhiteNoise, #Bandwidth type that uses Gaussian white noise assumption to get correlation confidence bounds
		BandwidthBartlett, #Bandwidth type that uses Bartlett formula to get correlation confidence bounds
		BandwidthP2003 #Bandwidth type implements the bandwidth estimation method of Politis (2003) "Adaptive Bandwidth Choice"


#******************************************************************************


#----------------------------------------------------------
#SET CONSTANTS FOR MODULE
#----------------------------------------------------------
const validBandwidthMethodString = ["whiteNoise", "bartlett", "P2003"] #Vector of valid string representations for bandwidth methods






#----------------------------------------------------------
#TYPE
#	KernelFunction (Abstract type nesting all other kernel function types)
#	KernelUniform
#	KernelTriangular
#	KernelEpanechnikov
#	KernelQuartic
#	KernelGaussian
#	KernelPR1995FlatTop
#	KernelPP2002Trap
#	KernelPP2002Smooth
#FIELDS
#	Fields store parameteres necessary to calculate each kernel function.
#BASIC METHODS (more methods below)
#	string(kF::KernelFunction)
#	copy(kF::KernelFunction)
#	show(kF::KernelFunction)
#METHODS DEFINED SEPARETELY BELOW
#	evaluate: Used to evaluate kernel functions
#	nonzeroDomain: Returns the range over which the kernel function outputs non-zero values
#	paramDomain: Returns the valid range of values that can be taken by the specified parameter of the kernel
#PURPOSE
#	These types are used to define different kernel functions, all of which can be called via the evaluate function.
#NOTES
#	Many kernel functions will have no fields, as they require no input parameters
#----------------------------------------------------------
#------- ABSTRACT TYPE DEFINITIONS ------
#Abstract type to nest all other kernel function types
abstract KernelFunction
#------- TYPE DEFINITIONS ---------------
#Common kernel functions defined on the interval [-1, 1]
type KernelUniform <: KernelFunction
end
type KernelTriangular <: KernelFunction
end
type KernelEpanechnikov <: KernelFunction
end
type KernelQuartic <: KernelFunction
end
#Common kernel functions defined on the interval [-Inf, Inf]
type KernelGaussian <: KernelFunction
	p1::Float64 #This is the scaling term, almost always set to 1 / sqrt(2*pi). We parameterize it for efficiency of computation on vector-input evaluations
	KernelGaussian(p1::Float64) = (p1 <= 0) ? error("The parameter of this kernel function must be strictly positive") : new(p1)
end
KernelGaussian() = KernelGaussian(1 / sqrt(2 * pi)) #Constructor for the most common use-case
#Politis, Romano (1993) "On a family of smoothing kernels of infinite order"
type KernelPR1993FlatTop <: KernelFunction
	m::Float64
	M::Float64
	function KernelPR1993FlatTop(m::Float64, M::Float64)
		(m <= 0 || M <= 0) && error("Parameters of this kernel function must be strictly positive")
		(m >= M) && error("First parameter of this kernel function must be stricty smaller than second parameter")
		new(m, M)
	end
end
KernelPR1993FlatTop() = KernelPR1993FlatTop(0.5, 1.0) #A popular specification of this kernel function
#Politis, Romano (2003) "Adaptive Bandwidth Choice" [NOTE, equivalent to KernelPR1993FlatTop with m = 1/2, and M = 1]
type KernelP2003FlatTop <: KernelFunction
end
#Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences"
type KernelPP2002Trap <: KernelFunction
	p1::Float64
	KernelPP2002Trap(p1::Float64) = !(0 < p1 <= 0.5) ? error("Invalid kernel function parameter") : new(p1)
end
KernelPP2002Trap() = KernelPP2002Trap(0.43) #Optimal value from Paparoditis, Politis (2002)
#Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences"
type KernelPP2002Smooth <: KernelFunction
	p1::Float64
	KernelPP2002Smooth(p1::Float64) = !(1 <= p1 < Inf) ? error("Invalid kernel function parameter") : new(p1)
end
KernelPP2002Smooth(p1::Int) = KernelPP2002Smooth(convert(Float64, p1))
KernelPP2002Smooth() = KernelPP2002Smooth(1.3) #Optimal value from Paparoditis, Politis (2002)
#------- METHODS ------
#string method gets string representation of type
string(kT::KernelUniform) = "uniform"
string(kT::KernelTriangular) = "triangular"
string(kT::KernelEpanechnikov) = "epanechnikov"
string(kT::KernelQuartic) = "quartic"
string(kT::KernelGaussian) = "gaussian"
string(kT::KernelPR1993FlatTop) = "PR1993FlatTop"
string(kT::KernelP2003FlatTop) = "P2003FlatTop"
string(kT::KernelPP2002Trap) = "PP2002Trap"
string(kT::KernelPP2002Smooth) = "PP2002Smooth"
#copy methods
copy(kT::KernelUniform) = KernelUniform()
copy(kT::KernelTriangular) = KernelTriangular()
copy(kT::KernelEpanechnikov) = KernelEpanechnikov()
copy(kT::KernelQuartic) = KernelQuartic()
copy(kT::KernelGaussian) = KernalGaussian(copy(kT.p1))
copy(kT::KernelPR1993FlatTop) = KernelPR1993FlatTop(copy(kT.m), copy(kT.M))
copy(kT::KernelP2003FlatTop) = KernelP2003FlatTop()
copy(kT::KernelPP2002Trap) = KernelPP2002Trap(copy(kT.p1))
copy(kT::KernelPP2002Smooth) = KernelPP2002Smooth(copy(kT.p1))
#show method for parameter-less kernel functions
show{T<:Union(KernelUniform, KernelTriangular, KernelEpanechnikov, KernelQuartic, KernelP2003FlatTop)}(io::IO, k::T) = println(io, "kernel function = " * string(k))
#show methods for kernel functions with parameters
function show{T<:Union(KernelGaussian, KernelPP2002Trap, KernelPP2002Smooth)}(io::IO, k::T)
	println(io, "kernel function = " * string(k))
	println(io, "    p1 = " * string(k.p1))
end
function show(io::IO, k::KernelPR1993FlatTop)
	println(io, "kernel function = " * string(k))
	println(io, "    m = " * string(k.m))
	println(io, "    M = " * string(k.M))
end
#show wrapper for STDOUT
show{T<:KernelFunction}(k::T) = show(STDOUT, k)
#toKernelFunction gets a kernel function given a string input (NOT TYPE STABLE)
function toKernelFunction(x::ASCIIString) #WARNING: FUNCTION IS NOT TYPE STABLE
	if x == "uniform"
		return(KernelUniform())
	elseif x == "triangular"
		return(KernelTriangular())
	elseif x == "epanechnikov"
		return(KernelEpanechnikov())
	elseif x == "quartic"
		return(KernelQuartic())
	elseif x == "gaussian"
		return(KernelGaussian())
	elseif x == "PR1993FlatTop"
		return(KernelPR1993FlatTop())
	elseif x == "PP2002Trap"
		return(KernelPP2002Trap())
	elseif x == "PP2002Smooth"
		return(KernelPP2002Smooth())
	else
		error("Input string not a representation of a known kernel function")
	end
end


#----------------------------------------------------------
#FUNCTION
#	evaluate
#INPUT
#	Let T<:Number. Function accepts:
#		(kT::KernelType, x::Number): Apply kernel function defined by kT to input x
#		(kT::KernelType, x::Vector{T}): Wrapper on (kT::KernelType, x::Number)
#OUTPUT
#	Output is of the same type as x for all methods
#PURPOSE
#	Obtain output of various kernel functions
#NOTES
#----------------------------------------------------------
#Common Kernel functions defined on the interval [-1, 1]
evaluate{T<:Number}(x::T, kT::KernelUniform) = indicator(x, UnitRange(-1, 1)) * 0.5
evaluate{T<:Number}(x::T, kT::KernelTriangular) = indicator(x, UnitRange(-1, 1)) * (1 - abs(x))
evaluate{T<:Number}(x::T, kT::KernelEpanechnikov) = indicator(x, UnitRange(-1, 1)) * 0.75 * (1 - x^2)
evaluate{T<:Number}(x::T, kT::KernelQuartic) = indicator(x, UnitRange(-1, 1)) * 0.9375 * (1 - x^2)^2
#Common Kernel functions defined on the interval [-Inf, Inf]
evaluate{T<:Number}(x::T, kT::KernelGaussian) = kT.p1 * exp(-0.5 * x^2)
#KernelPR1993FlatTop
function evaluate{T<:Number}(x::T, kT::KernelPR1993FlatTop)
	if abs(x) <= kT.m
		return(one(T))
	elseif abs(x) <= kT.M
		return(1 - (abs(x) - kT.m) / (kT.M - kT.m))
	else
		return(zero(T))
	end
end
#KernelP2003FlatTop (equivalent to KernelPR1993FlatTop with m=0.5 and M=1)
function evaluate{T<:Number}(x::T, kT::KernelP2003FlatTop)
	if abs(x) <= 0.5
		return(one(T))
	elseif abs(x) <= 1
		return(2 * (1 - abs(x)))
	else
		return(zero(T))
	end
end
#KernelPP2002Trap
function evaluate{T<:Number}(x::T, kT::KernelPP2002Trap)
	if x < 0
		return(zero(T))
	elseif x < kT.p1
		return(x / kT.p1)
	elseif x < 1 - kT.p1
		return(one(T))
	elseif x < 1
		return((1 - x) / kT.p1)
	else
		return(zero(T))
	end
end
#KernelPP2002Smooth
function evaluate{T<:Number}(x::T, kT::KernelPP2002Smooth)
	if 0 <= x <= 1
		return(1 - abs(2*x - 1)^kT.p1)
	else
		return(zero(T))
	end
end
#Vector input wrappers
evaluate{T<:Number}(x::Vector{T}, kT::KernelUniform) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::Vector{T}, kT::KernelTriangular) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::Vector{T}, kT::KernelEpanechnikov) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::Vector{T}, kT::KernelQuartic) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::Vector{T}, kT::KernelGaussian) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::Vector{T}, kT::KernelPR1993FlatTop) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::Vector{T}, kT::KernelP2003FlatTop) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::Vector{T}, kT::KernelPP2002Trap) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::Vector{T}, kT::KernelPP2002Smooth) = [ evaluate(x[n], kT) for n = 1:length(x) ]





#----------------------------------------------------------
#FUNCTION
#	nonzeroDomain
#	paramDomain
#INPUT
#	nonzeroRange accepts:
#		(kF::KernelFunction): kF is one of the kernel function types
#	paramDomain accepts:
#		(kF::KernelFunction, paramNum::Int): kF is a kernel function type. paramNum is an indicator of the parameter you are interested in. paramNum corresponds to the field number. Kernel function types with no fields (ie no parameters) will throw an error if input to this function.
#OUTPUT
#	Output is of both functions is a UnitRange
#PURPOSE
#	The purpose of nonzeroDomain is to return the domain (of kernel function input) over which the output of the kernel function is non-zero
#	The purpose of paramDomain is to return the domain over which the indicated parameter can be defined.
#NOTES
#	A method for each of these functions should be defined for every kernel function type
#----------------------------------------------------------
nonzeroDomain(kT::KernelUniform) = UnitRange(-1, 1)
nonzeroDomain(kT::KernelTriangular) = UnitRange(-1, 1)
nonzeroDomain(kT::KernelEpanechnikov) = UnitRange(-1, 1)
nonzeroDomain(kT::KernelQuartic) = UnitRange(-1, 1)
nonzeroDomain(kT::KernelGaussian) = UnitRange(-Inf, Inf)
nonzeroDomain(kT::KernelPR1993FlatTop) = UnitRange(-kT.M, kT.M)
nonzeroDomain(kT::KernelP2003FlatTop) = UnitRange(-1, 1)
nonzeroDomain(kT::KernelPP2002Trap) = UnitRange(0, 1)
nonzeroDomain(kT::KernelPP2002Smooth) = UnitRange(0, 1)
paramDomain(kT::KernelUniform, pNum::Int) = error("Kernel function has no parameters")
paramDomain(kT::KernelTriangular, pNum::Int) = error("Kernel function has no parameters")
paramDomain(kT::KernelEpanechnikov, pNum::Int) = error("Kernel function has no parameters")
paramDomain(kT::KernelQuartic, paramNum::Int) = error("Kernel function has no parameters")
paramDomain(kT::KernelGaussian, paramNum::Int) = (paramNum == 1) ? UnitRange(nextfloat(0.0), prevfloat(Inf)) : error("Invalid parameter number")
function paramDomain(kT::KernelPR1993FlatTop, paramNum::Int)
	if paramNum == 1
		return(UnitRange(nextfloat(0.0), Inf))
	elseif paramNum == 2
		return(UnitRange(nextfloat(0.0), Inf))
	else
		error("Invalid parameter number")
	end
end
paramDomain(kT::KernelP2003FlatTop, pNum::Int) = error("Kernel function has no parameters")
paramDomain(kT::KernelPP2002Trap) = UnitRange(nextfloat(0.0), 0.5)
paramDomain(kT::KernelPP2002Smooth) = UnitRange(1.o, prevfloat(Inf))






#----------------------------------------------------------
#TYPE
#	BandwidthWhiteNoise (Iteratively test for significance of autocorrelations using Guassian white-noise assumption for confidence interval at every lag)
#	BandwidthBartlett (Iteratively test for significance of autocorrelations using Bartlett formula at each lag k, and testing for significance of autocorrelation at lag k+1)
#	BandwidthP2003 (Politis (2003) "Adaptive Bandwidth Choice", Section 2.1, Empirical rule of picking M. One needs to specify parameters c and K_N in this rule. Constructor employs recommended values from paper of c=2.0 and K=5, or else data-driven K for the case where number of observations is supplied).
#FIELDS
#	Fields store parameteres necessary to estimate bandwidth using each method. All types contain the following field:
#		adjustmentTerm::Float64: The estimated bandwidth is multiplied by this scalar, thus most constructors set this field to 1.0. However, the field is useful as it allows empiricists to apply ad hoc rules on the go.
#	Other fields are unique to the type and described above
#METHODS
#	string(b::BandwidthMethod)
#	copy(b::BandwidthMethod)
#	show(b::BandwidthMethod)
#	toBandwidthMethod(x::ASCIIString): Convert a string representation of a bandwidth method to the appropriate type
#PURPOSE
#	These types are used to define different methods for detecting bandwidth, all of which are called via the bandwidth function.
#NOTES
#----------------------------------------------------------
#--------- ABSTRACT TYPE DEFINTIION
abstract BandwidthMethod
#--------- TYPE DEFINITIONS
#BandwidthWhiteNoise
type BandwidthWhiteNoise <: BandwidthMethod
	adjustmentTerm::Float64
	BandwidthWhiteNoise(adjustmentTerm::Float64) = (adjustmentTerm <= 0) ? error("Adjustment scalar must be strictly greater than zero") : new(adjustmentTerm)
end
BandwidthWhiteNoise() = BandwidthWhiteNoise(1.0)
#BandwidthBartlett
type BandwidthBartlett <:BandwidthMethod
	adjustmentTerm::Float64
	BandwidthBartlett(adjustmentTerm::Float64) = (adjustmentTerm <= 0) ? error("Adjustment scalar must be strictly greater than zero") : new(adjustmentTerm)
end
BandwidthBartlett() = BandwidthBartlett(1.0)
#BandwidthP2003 (Politis (2003) "Adaptive Bandwidth Choice")
type BandwidthP2003 <:BandwidthMethod
	adjustmentTerm::Float64
	c::Float64
	K::Int
	function BandwidthP2003(adjustmentTerm::Float64, c::Float64, K::Int)
		adjustmentTerm <= 0 && error("Adjustment scalar must be strictly greater than zero")
		c <= 0 && error("c must be strictly positive")
		K < 1 && error("K must be strictly positive")
		new(adjustmentTerm, c, K)
	end
end
BandwidthP2003() = BandwidthP2003(1.0, 2.0, 5)
BandwidthP2003(numObs::Int) = BandwidthP2003(1.0, 2.0, max(5, convert(Int, ceil(sqrt(log(10, numObs))))))
#-------------- METHODS---------
#string method returns string representation of BandwidthMethod
string(x::BandwidthWhiteNoise) = "whiteNoise"
string(x::BandwidthBartlett) = "bartlett"
string(x::BandwidthP2003) = "P2003"
#copy methods
copy(x::BandwidthWhiteNoise) = BandwidthWhiteNoise(copy(x.adjustmentTerm))
copy(x::BandwidthBartlett) = BandwidthBartlett(copy(x.adjustmentTerm))
copy(x::BandwidthP2003) = BandwidthP2003(copy(x.adjustmentTerm), copy(x.c), copy(x.K))
#show method for BandwidthMethod
function show{T<:Union(BandwidthWhiteNoise, BandwidthBartlett)}(io::IO, b::T)
	println(io, "bandwidth method = " * string(b))
	println(io, "    adjustment term = " * string(b.adjustmentTerm))
end
function show(io::IO, b::BandwidthP2003)
	println(io, "bandwidth method = " * string(b))
	println(io, "    adjustment term = " * string(b.adjustmentTerm))
	println(io, "    c = " * string(b.c))
	println(io, "    K = " * string(b.K))
end
#show wrapper for STDOUT
show{T<:BandwidthMethod}(b::T) = show(STDOUT, b)
#toBandwidth converts string input to BandwidthMethod
function toBandwidthMethod(x::ASCIIString, adjustmentTerm::Float64=1.0) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "whiteNoise"
		return(BandwidthWhiteNoise(adjustmentTerm))
	elseif x == "bartlett"
		return(BandwidthBartlett(adjustmentTerm))
	elseif x == "P2003"
		return(BandwidthP2003(adjustmentTerm, 2.0, 5))
	else
		error("Invalid string representation of bandwidth method type")
	end
end
function toBandwidthMethod(x::ASCIIString, adjustmentTerm::Float64, p1::Number) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "P2003"
		return(BandwidthP2003(adjustmentTerm, convert(Float64, p1), 5))
	else
		if any(x .== validBandwidthMethodString)
			error("Additional parameters unnecessary for input bandwidth")
		else
			error("Invalid string representation of bandwidth method type")
		end
	end
end
function toBandwidthMethod(x::ASCIIString, adjustmentTerm::Float64, p1::Number, p2::Number) #WARNING: FUNCTION IS NOT TYPE STABLE (ALTHOUGH PLAYS SUCH A SMALL ROLE PERFORMANCE UNLIKELY TO BE AFFECTED)
	if x == "P2003"
		return(BandwidthP2003(adjustmentTerm, convert(Float64, p1), convert(Int, p2)))
	else
		if any(x .== validBandwidthMethodString)
			error("Additional parameters unnecessary for input bandwidth")
		else
			error("Invalid string representation of bandwidth method type")
		end
	end
end





#----------------------------------------------------------
#FUNCTION
#	bandwidth
#INPUT
#	Let T<:Number. Function accepts:
#		(x::Vector{T}, method::BandwidthWhiteNoise): Estimate bandwidth using white noise assumption to determine confidence bound for autocorrelations.
#		(x::Vector{T}, method::BandwidthBartlett): Estimate bandwidth using Bartlett assumptions, e.g. underlying model is Moving Average (MA)
#		(x::Vector{T}, method::BandwidthP2003): Estimate bandwidth using the Empirical rule described in Politis (2003) "Adaptive Bandwidth Estimation"
#OUTPUT
#	Output is the tuple (M::Int, xVar::Float64, covVec::Vector{Float64})
#		M is the bandwidth estimator. Minimum possible estimate is 1, maximum possible estimate is length(x) - 1.
#		xVar is the variance of the input data x but scaled using (1/length(x)) instead of (1/(length(x)-1))
#		covVec is the vector of sample autocovariances (NOT autocorrelations) from lag 1 to lag K of input x. K is determined within the function.
#	We pass all three of these variables as output because xVar and covVec are frequently used in the calling function, and so this saves duplicating computations.
#PURPOSE
#	Obtain a bandwidth estimator for use with a kernel function applied to some function (frequently the autocorrelation function) of the input data
#NOTES
#	M can be close to length(x). Typically, this is not a good situation, so the user will probably want to control for this in the calling function and bind maximum M to a more reasonable value.
#	M has an enforced minimum value of 2. Allowing results of M < 2 typically result in kernel function evaluations that always equal 0.
#	Autocorrelations are calculated using autocor (e.g. fft methods) in blocks of 20 (i.e. if bandwidth estimate is not found using first 20 autocorrelations, the next 20 will be calculated, and so on)
#----------------------------------------------------------
function bandwidth{T<:Number}(x::Vector{T}, method::BandwidthWhiteNoise)
	length(x) < 2 && error("Input data must have at least two observations")
	mHat = 1
	corVec = Array(Float64, 0)
	append!(corVec, autocor(x, 1:min(20, length(x)-1)))
	whiteNoiseBound = sqrt(1 / length(x))
	nullRejected = false
	if abs(corVec[1]) > 1.96 * whiteNoiseBound #"1.96*" implies we are hard-coding hypothesis test at ~95% confidence level. Note, if this condition is not triggered, then we detect no autocorrelation and so leave mHat at 1
		mSt = 1
		while mSt < length(x)
			for m = mSt:length(corVec)-1
				if abs(corVec[m+1]) < 1.96 * whiteNoiseBound
					nullRejected = true
					mHat = m + 1
					break
				end
			end
			if nullRejected == true
				break
			else
				mSt = length(corVec)
				append!(corVec, autocor(x, mSt+1:min(mSt+20, length(x)))) #Grow auto-correlation vector by 20 and try again
			end
		end
		if nullRejected == false
			mHat = length(x)
		end
	end
	M = convert(Int, ceil(method.adjustmentTerm * 2 * mHat)) #"2*" is a standard rule, see e.g. Politis (2003) "Adaptive Bandwidth Choice". method.adjustmentTerm allows user to tinker with output.
	if M > length(x) - 1
		M = length(x) - 1 #Set maximum value for M
	elseif M < 2
		M = 2 #Set minimum value for M
	end
	xVar = ((length(x)-1) /  length(x)) * var(x) #Used to scale autocorrelations to autocovariances for return argument
	return(M, xVar, xVar * corVec)
end
function bandwidth{T<:Number}(x::Vector{T}, method::BandwidthBartlett)
	length(x) < 2 && error("Input data must have at least two observations")
	mHat = 1
	corVec = Array(Float64, 0)
	append!(corVec, autocor(x, 1:min(20, length(x)-1)))
	bartlettBound = sqrt(1 / length(x))
	nullRejected = false
	if abs(corVec[1]) > 1.96 * bartlettBound #"1.96*" implies we are hard-coding hypothesis test at ~95% confidence level. Note, if this condition is not triggered, then we detect no autocorrelation and so leave mHat at 1
		corSum = 0.0 #Sum of squared autocorrelations. We add lags to the overall sum iteratively within the loop
		mSt = 1
		while mSt < length(x)
			for m = mSt:length(corVec)-1
				corSum += 2*corVec[m]^2
				bartlettBound = sqrt((1 + corSum) / length(x))
				if abs(corVec[m+1]) < 1.96 * bartlettBound
					nullRejected = true
					mHat = m + 1
					break
				end
			end
			if nullRejected == true
				break
			else
				mSt = length(corVec)
				append!(corVec, autocor(x, mSt+1:min(mSt+20, length(x)))) #Grow auto-correlation vector by 20 and try again
			end
		end
		if nullRejected == false
			mHat = length(x)
		end
	end
	M = convert(Int, ceil(method.adjustmentTerm * 2 * mHat)) #"2*" is a standard rule, see e.g. Politis (2003) "Adaptive Bandwidth Choice". method.adjustmentTerm allows user to tinker with output.
	if M > length(x) - 1
		M = length(x) - 1 #Set maximum value for M
	elseif M < 2
		M = 2 #Set minimum value for M
	end
	xVar = ((length(x)-1) /  length(x)) * var(x) #Used to scale autocorrelations to autocovariances for return argument
	return(M, xVar, xVar * corVec)
end
function bandwidth{T<:Number}(x::Vector{T}, method::BandwidthP2003)
	length(x) < 2 && error("Input data must have at least two observations")
	mHat = 1
	corVec = Array(Float64, 0)
	append!(corVec, autocor(x, 1:min(20, length(x)-1))) #Add first block of autocorrelations
	corBound = method.c * sqrt(log(10, length(x)) / length(x)) #Note, use of log base 10 is deliberate and recommended in Politis (2003)
	KCounter = 0
	mHatFound = false
	for k = 1:length(x)-2
		if abs(corVec[k]) < corBound
			KCounter += 1 #Bound is satisfied so add one to counter
		else
			KCounter = 0 #Bound is not satisfied, so reset the counter
		end
		if KCounter >= method.K #We found K autocorrelations in a row that satisfy the bound, so break.
			mHat = k - method.K + 1
			mHatFound = true
			break
		end
		if k == length(corVec) #We've run out of autocorrelations to check, so we need to add another block of them to corVec
			append!(corVec, autocor(x, length(corVec)+1:min(length(corVec)+20, length(x)-1)))
		end
	end
	if mHatFound == false
		mHat = length(x) - 1
	end
	M = convert(Int, ceil(method.adjustmentTerm * 2 * mHat)) #"2*" is a standard rule, see e.g. Politis (2003) "Adaptive Bandwidth Choice". method.adjustmentTerm allows user to tinker with output.
	if M > length(x) - 1
		M = length(x) - 1 #Set maximum value for M
	elseif M < 2
		M = 2 #Set minimum value for M
	end
	xVar = ((length(x)-1) /  length(x)) * var(x) #Used to scale autocorrelations to autocovariances for return argument
	return(M, xVar, xVar * corVec)
end



#----------------------------------------------------------
#FUNCTION
#	indicator
#INPUT
#	Function accepts:
#		{T<:Number}(x::T, LB::T, UB::T): Returns one if LB <= x <= UB, zero otherwise.
#		{T<:Number}(x::T, r::UnitRange): Returns one if x lies in the range, zero otherwise
#OUTPUT
#	Output is Int of 0 or 1.
#PURPOSE
#	An indicator function
#NOTES
#----------------------------------------------------------
function indicator{T<:Number}(x::T, LB::T, UB::T)
	if LB <= x <= UB
		return(1)
	else
		return(0)
	end
end
function indicator{T<:Number}(x::T, r::UnitRange)
	if r.start <= x <= r.stop
		return(1)
	else
		return(0)
	end
end


end # module
