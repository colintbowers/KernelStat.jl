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
		activedomain, #Function that returns a UnitRange indicating the range over which kernel output is non-zero
		paramdomain, #Function that returns a UnitRange indicating the feasible domain kernel parameters
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
		KernelPR1994SB, #Kernel function for variance estimator from Politis, Romano (1994) "The Stationary Bootstrap"
		hacvariance, #Function for estimating heteroskedasticity and auto-correlation consistent variance
		HACVarianceMethod, #Abstract type for nesting all HAC variance methods
		HACVarianceBasic, #Basic HAC variance estimator
		bandwidth, #Function for estimating bandwidth
		BandwidthMethod, #Abstract type that nests all bandwidth method types
		BandwidthWhiteNoise, #Bandwidth type that uses Gaussian white noise assumption to get correlation confidence bounds
		BandwidthBartlett, #Bandwidth type that uses Bartlett formula to get correlation confidence bounds
		BandwidthP2003 #Bandwidth type implements the bandwidth estimation method of Politis (2003) "Adaptive Bandwidth Choice"


#******************************************************************************


#----------------------------------------------------------
#SET CONSTANTS FOR MODULE
#----------------------------------------------------------
#None currently needed




#----------------------------------------------------------
#TYPE
#	abstract super-type and one type for each kernel function
#PURPOSE
#	These types are used to define different kernel functions, parameterised by the fields of the type, and all of which can be called via the evaluate function.
#----------------------------------------------------------
#Abstract super type
abstract KernelFunction
#Dummy type for evaluate function that ditches the indicator function check
type KernelNoIndCheck; end
#------- TYPE DEFINITIONS ---------------
#Common kernel functions
type KernelUniform <: KernelFunction
	LB::Float64
	UB::Float64
	function KernelUniform(LB::Float64, UB::Float64)
		LB >= UB && error("lower bound greater than upper bound")
		new(LB, UB)
	end
end
KernelUniform() = KernelUniform(-1.0, 1.0)
KernelUniform(LB::Number, UB::Number) = KernelUniform(convert(Float64, LB), convert(Float64, UB))
type KernelTriangular <: KernelFunction; end #Currently only defined on [-1, 1]
type KernelEpanechnikov <: KernelFunction; end #Currently only defined on [-1, 1]
type KernelQuartic <: KernelFunction; end #Currently only defined on [-1, 1]
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
type KernelP2003FlatTop <: KernelFunction; end
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
#Politis, Romano (1994) "The Stationary Bootstrap" (kernel function for variance estimator, see eqn 7)
type KernelPR1994SB <: KernelFunction
	p1::Int #The upper bound (usually the number of observations)
	p2::Float64 #parameter of geometric distribution (1 / expected block length)
	function KernelPR1994SB(p1::Int, p2::Float64)
		p1 < 1 && error("Kernel function upper bound must be greater than zero")
		!(0 < p2 < 1) && error("Geometric distribution parameter must lie on the interval (0, 1)")
		new(p1, p2)
	end
end
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
string(kT::KernelPR1994SB) = "PR1994SB"
#copy, deepcopy methods
copy(kT::KernelUniform) = KernelUniform(copy(kT.LB), copy(kT.UB))
copy(kT::KernelTriangular) = KernelTriangular()
copy(kT::KernelEpanechnikov) = KernelEpanechnikov()
copy(kT::KernelQuartic) = KernelQuartic()
copy(kT::KernelGaussian) = KernalGaussian(copy(kT.p1))
copy(kT::KernelPR1993FlatTop) = KernelPR1993FlatTop(copy(kT.m), copy(kT.M))
copy(kT::KernelP2003FlatTop) = KernelP2003FlatTop()
copy(kT::KernelPP2002Trap) = KernelPP2002Trap(copy(kT.p1))
copy(kT::KernelPP2002Smooth) = KernelPP2002Smooth(copy(kT.p1))
copy(kT::KernelPR1994SB) = KernelPR1994SB(copy(kT.p1), copy(kT.p2))
deepcopy(kT::KernelUniform) = KernelUniform(deepcopy(kT.LB), deepcopy(kT.UB))
deepcopy(kT::KernelTriangular) = KernelTriangular()
deepcopy(kT::KernelEpanechnikov) = KernelEpanechnikov()
deepcopy(kT::KernelQuartic) = KernelQuartic()
deepcopy(kT::KernelGaussian) = KernalGaussian(deepcopy(kT.p1))
deepcopy(kT::KernelPR1993FlatTop) = KernelPR1993FlatTop(deepcopy(kT.m), deepcopy(kT.M))
deepcopy(kT::KernelP2003FlatTop) = KernelP2003FlatTop()
deepcopy(kT::KernelPP2002Trap) = KernelPP2002Trap(deepcopy(kT.p1))
deepcopy(kT::KernelPP2002Smooth) = KernelPP2002Smooth(deepcopy(kT.p1))
deepcopy(kT::KernelPR1994SB) = KernelPR1994SB(deepcopy(kT.p1), deepcopy(kT.p2))
#show method for parameter-less kernel functions
show{T<:Union(KernelUniform, KernelTriangular, KernelEpanechnikov, KernelQuartic, KernelP2003FlatTop)}(io::IO, k::T) = println(io, "kernel function = " * string(k))
#show methods for kernel functions with parameters
function show{T<:Union(KernelGaussian, KernelPP2002Trap, KernelPP2002Smooth, KernelPR1994SB)}(io::IO, k::T)
	println(io, "kernel function = " * string(k))
	println(io, "    p1 = " * string(k.p1))
end
function show{T<:Union(KernelPR1994SB)}(io::IO, k::T)
	println(io, "kernel function = " * string(k))
	println(io, "    p1 = " * string(k.p1))
	println(io, "    p2 = " * string(k.p2))
end
function show{T<:Union(KernelUniform)}(io::IO, k::T)
	println(io, "kernel function = " * string(k))
	println(io, "    lower bound = " * string(k.LB))
	println(io, "    upper bound = " * string(k.UB))
end
function show(io::IO, k::KernelPR1993FlatTop)
	println(io, "kernel function = " * string(k))
	println(io, "    m = " * string(k.m))
	println(io, "    M = " * string(k.M))
end
#show wrapper for STDOUT
show{T<:KernelFunction}(k::T) = show(STDOUT, k)
#------ EVALUATE method for evaluating kernel function at a given value ------------------------------
#Common Kernel functions
evaluate{T<:Number}(x::T, kT::KernelUniform) = indicator(x, UnitRange(kT.LB, kT.UB)) * (1 / (kT.UB - kT-LB))
evaluate{T<:Number}(x::T, kT::KernelUniform, ::KernelNoIndCheck) = 1 / (kT.UB - kT-LB)
evaluate{T<:Number}(x::T, kT::KernelTriangular) = indicator(x, UnitRange(-1, 1)) * (1 - abs(x))
evaluate{T<:Number}(x::T, kT::KernelTriangular, ::KernelNoIndCheck) = 1 - abs(x)
evaluate{T<:Number}(x::T, kT::KernelEpanechnikov) = indicator(x, UnitRange(-1, 1)) * 0.75 * (1 - x^2)
evaluate{T<:Number}(x::T, kT::KernelEpanechnikov, ::KernelNoIndCheck) = 0.75 * (1 - x^2)
evaluate{T<:Number}(x::T, kT::KernelQuartic) = indicator(x, UnitRange(-1, 1)) * 0.9375 * (1 - x^2)^2
evaluate{T<:Number}(x::T, kT::KernelQuartic, ::KernelNoIndCheck) = 0.9375 * (1 - x^2)^2
evaluate{T<:Number}(x::T, kT::KernelGaussian) = kT.p1 * exp(-0.5 * x^2)
evaluate{T<:Number}(x::T, kT::KernelGaussian, ::KernelNoIndCheck) = kT.p1 * exp(-0.5 * x^2)
#KernelPR1993FlatTop
function evaluate{T<:Number}(x::T, kT::KernelPR1993FlatTop)
	if abs(x) <= kT.m; return(one(T))
	elseif abs(x) <= kT.M; return(1 - (abs(x) - kT.m) / (kT.M - kT.m))
	else; return(zero(T))
	end
end
evaluate{T<:Number}(x::T, kT::KernelPR1993FlatTop, ::KernelNoIndCheck) = evaluate(x, kT)
#KernelP2003FlatTop (equivalent to KernelPR1993FlatTop with m=0.5 and M=1)
function evaluate{T<:Number}(x::T, kT::KernelP2003FlatTop)
	if abs(x) <= 0.5; return(one(T))
	elseif abs(x) <= 1; return(2 * (1 - abs(x)))
	else; return(zero(T))
	end
end
evaluate{T<:Number}(x::T, kT::KernelP2003FlatTop, ::KernelNoIndCheck) = evaluate(x, kT)
#KernelPP2002Trap
function evaluate{T<:Number}(x::T, kT::KernelPP2002Trap)
	if x < 0; return(zero(T))
	elseif x < kT.p1; return(x / kT.p1)
	elseif x < 1 - kT.p1; return(one(T))
	elseif x < 1; return((1 - x) / kT.p1)
	else; return(zero(T))
	end
end
evaluate{T<:Number}(x::T, kT::KernelPP2002Trap, ::KernelNoIndCheck) = evaluate(x, kT)
#KernelPP2002Smooth
function evaluate{T<:Number}(x::T, kT::KernelPP2002Smooth)
	if 0 <= x <= 1; return(1 - abs(2*x - 1)^kT.p1)
	else; return(zero(T))
	end
end
evaluate{T<:Number}(x::T, kT::KernelPP2002Smooth, ::KernelNoIndCheck) = 1 - abs(2*x - 1)^kT.p1
#KernelPR1994SB
evaluate(x::Int, kT::KernelPR1994SB) = indicator(x, UnitRange(0, kT.p1)) * ((1 - x/kT.p1)*(1 - kT.p2)^x + (x/kT.p1)*(1 - kT.p2)^(x - kT.p1))
evaluate(x::Int, kT::KernelPR1994SB, ::KernelNoIndCheck) = (1 - x/kT.p1)*(1 - kT.p2)^x + (x/kT.p1)*(1 - kT.p2)^(x - kT.p1)
#Array input wrappers
evaluate{T<:Number}(x::AbstractVector{T}, kT::KernelFunction) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::AbstractMatrix{T}, kT::KernelFunction) = [ evaluate(x[n, m], kT) for n = 1:size(x, 1), m = 1:size(x, 2) ]
#-----activedomain and paramdomain ------------------------
#The purpose of activedomain is to return the domain (of kernel function input) over which the output of the kernel function is active
#The purpose of paramdomain is to return the domain over which the indicated parameter can be defined.
activedomain(kT::KernelUniform) = UnitRange(kT.LB, kT.UB)
activedomain(kT::KernelTriangular) = UnitRange(-1, 1)
activedomain(kT::KernelEpanechnikov) = UnitRange(-1, 1)
activedomain(kT::KernelQuartic) = UnitRange(-1, 1)
activedomain(kT::KernelGaussian) = UnitRange(-Inf, Inf)
activedomain(kT::KernelPR1993FlatTop) = UnitRange(-kT.M, kT.M)
activedomain(kT::KernelP2003FlatTop) = UnitRange(-1, 1)
activedomain(kT::KernelPP2002Trap) = UnitRange(0, 1)
activedomain(kT::KernelPP2002Smooth) = UnitRange(0, 1)
activedomain(kT::KernelPR1994SB) = UnitRange(1, kT.p1)
function paramdomain(kT::KernelUniform, pNum::Int)
	pNum == 1 && return(UnitRange(-Inf, kT.UB))
	pNum == 2 && return(UnitRange(kT.LB, Inf))
	error("Invalid parameter number")
end
paramdomain(kT::KernelTriangular, pNum::Int) = error("Kernel function has no parameters")
paramdomain(kT::KernelEpanechnikov, pNum::Int) = error("Kernel function has no parameters")
paramdomain(kT::KernelQuartic, paramNum::Int) = error("Kernel function has no parameters")
paramdomain(kT::KernelGaussian, paramNum::Int) = (paramNum == 1) ? UnitRange(nextfloat(0.0), prevfloat(Inf)) : error("Invalid parameter number")
function paramdomain(kT::KernelPR1993FlatTop, paramNum::Int)
	paramNum == 1 && return(UnitRange(nextfloat(0.0), Inf))
	paramNum == 2 && return(UnitRange(nextfloat(0.0), Inf))
	error("Invalid parameter number")
end
paramdomain(kT::KernelP2003FlatTop, pNum::Int) = error("Kernel function has no parameters")
paramdomain(kT::KernelPP2002Trap, pNum::Int) = (pNum == 1) ? UnitRange(nextfloat(0.0), 0.5) : error("Invalid parameter number")
paramdomain(kT::KernelPP2002Smooth) = (pNum == 1) ? UnitRange(1, Inf) : error("Invalid parameter number")
function paramdomain(kT::KernelPR1994SB, pNum::Int)
	pNum == 1 && return(UnitRange(1, Inf))
	pNum == 2 && return(UnitRange(nextfloat(0.0), prevfloat(1.0)))
	error("Invalid parameter number")
end






#---------- TYPES FOR DIFFERENT METHODS OF BANDWIDTH ESTIMATION
#Abstract super-type
abstract BandwidthMethod
#Dummy type for using maximum possible bandwidth
type BandwidthMax <: BandwidthMethod; end
#type definitions
type BandwidthWhiteNoise <: BandwidthMethod
	adjustmentTerm::Float64
	BandwidthWhiteNoise(adjustmentTerm::Float64) = (adjustmentTerm <= 0) ? error("Adjustment scalar must be strictly greater than zero") : new(adjustmentTerm)
end
BandwidthWhiteNoise() = BandwidthWhiteNoise(1.0)
type BandwidthBartlett <:BandwidthMethod
	adjustmentTerm::Float64
	BandwidthBartlett(adjustmentTerm::Float64) = (adjustmentTerm <= 0) ? error("Adjustment scalar must be strictly greater than zero") : new(adjustmentTerm)
end
BandwidthBartlett() = BandwidthBartlett(1.0)
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
#string
string(x::BandwidthWhiteNoise) = "whiteNoise"
string(x::BandwidthBartlett) = "bartlett"
string(x::BandwidthP2003) = "P2003"
#copy, deepcopy
copy(x::BandwidthWhiteNoise) = BandwidthWhiteNoise(copy(x.adjustmentTerm))
copy(x::BandwidthBartlett) = BandwidthBartlett(copy(x.adjustmentTerm))
copy(x::BandwidthP2003) = BandwidthP2003(copy(x.adjustmentTerm), copy(x.c), copy(x.K))
deepcopy(x::BandwidthWhiteNoise) = BandwidthWhiteNoise(deepcopy(x.adjustmentTerm))
deepcopy(x::BandwidthBartlett) = BandwidthBartlett(deepcopy(x.adjustmentTerm))
deepcopy(x::BandwidthP2003) = BandwidthP2003(deepcopy(x.adjustmentTerm), deepcopy(x.c), deepcopy(x.K))
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

#---------- TYPES FOR DIFFERENT METHODS OF HAC-VARIANCE ESTIMATION
#Abstract super-type
abstract HACVarianceMethod
#type definitions
type HACVarianceBasic <: HACVarianceMethod
	kernelFunction::KernelFunction
	bandwidthMethod::BandwidthMethod
	function HACVarianceBasic(kernelFunction::KernelFunction, bandwidthMethod::BandwidthMethod)
		!(typeof(kernelFunction) <: Union(KernelUniform, KernelGaussian, KernelPR1994SB)) && error("This kernel function does not currently have a general enough form to accommodate HAC variance estimation")
		new(kernelFunction, bandwidthMethod)
	end
end
#---------- METHODS ----------------
string(::HACVarianceBasic) = "hacVarianceBasic"
copy(x::HACVarianceBasic) = HACVarianceBasic(copy(x.kernelFunction), copy(x.bandwidthMethod))
deepcopy(x::HACVarianceBasic) = HACVarianceBasic(deepcopy(x.kernelFunction), deepcopy(x.bandwidthMethod))
function show(io::IO, x::HACVarianceBasic)
	println(io, "HAC variance estimator type. Parameters are:")
	show(io, x.kernelFunction)
	show(io, x.bandwidthMethod)
end



#----------------------------------------------------------
#FUNCTION
#	hacvariance
#INPUT
#	(x::Vector{T<:Number}, method::HACVarianceMethod): Estimate variance using specified method
#OUTPUT
#	Output is a Float64 variance estimator
#NOTES
#----------------------------------------------------------
function hacvariance{T<:Number}(x::AbstractVector{T}, method::HACVarianceBasic)
	(M, xVar, xCov) = bandwidth(x, method.bandwidthMethod)
	length(xCov) < M && append!(xCov, autocov(x, length(xCov)+1:M)) #Get any additional autocovariances that we might need
	v = xVar
	noIndCheck = KernelNoIndCheck()
	for m = 1:M
		v += 2 * evaluate(m, method.kernelFunction, noIndCheck) * xCov[m]
	end
	return(v)
end








#----------------------------------------------------------
#FUNCTION
#	bandwidth
#INPUT
#	(x::Vector{T<:Number}, method::BandwidthMethod): Estimate bandwidth using specified method on input data x
#OUTPUT
#	Output is the tuple (M::Int, xVar::Float64, covVec::Vector{Float64})
#		M is the bandwidth estimator. Minimum possible estimate is 1, maximum possible estimate is length(x) - 1.
#		xVar is the variance of the input data x but scaled using (1/length(x)) instead of (1/(length(x)-1))
#		covVec is the vector of sample autocovariances (NOT autocorrelations) from lag 1 to lag K of input x. K is determined within the function.
#	We pass all three of these variables as output because xVar and covVec are frequently used in the calling function, and so this saves duplicating computations.
#NOTES
#	M can be close to length(x). Typically, this is not a good situation, so the user will probably want to control for this in the calling function and bind maximum M to a more reasonable value.
#	M has an enforced minimum value of 2. Allowing results of M < 2 typically result in kernel function evaluations that always equal 0.
#	Autocorrelations are calculated using autocor (e.g. fft methods) in blocks of 20 (i.e. if bandwidth estimate is not found using first 20 autocorrelations, the next 20 will be calculated, and so on)
#----------------------------------------------------------
function bandwidth{T<:Number}(x::AbstractVector{T}, method::BandwidthMax)
	length(x) < 2 && error("Input data must have at least two observations")
	return(length(x)-1, (length(x)-1 / length(x)) * var(x), Array(Float64, 0))
end
function bandwidth{T<:Number}(x::AbstractVector{T}, method::BandwidthWhiteNoise)
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
function bandwidth{T<:Number}(x::AbstractVector{T}, method::BandwidthBartlett)
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
function bandwidth{T<:Number}(x::AbstractVector{T}, method::BandwidthP2003)
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



#Non-exported indicator function for if a value lies in a specified range
indicator{T<:Number}(x::T, LB::T, UB::T) = (LB <= x <= UB) ? 1 : 0
indicator{T<:Number}(x::T, r::UnitRange) = (r.start <= x <= r.stop) ? 1 : 0




end # module
