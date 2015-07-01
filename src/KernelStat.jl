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
		setbandwidth!, #Function to adjust bandwidth of input kernel function
		KernelFunction, #Abstract type that nests all kernel function types
		KernelDummy, #Dummy type (useful for update! routines)
		KernelUniform, #Kernel function type
		KernelBartlett, #Kernel function type
		KernelTriangular, #Identical to KernelBartlett
		KernelEpanechnikov, #Kernel function type
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
		BandwidthDummy, #Dummy type (useful for update! functions)
		BandwidthMax, #Maximum possible bandwidth
		BandwidthWhiteNoise, #Bandwidth type that uses Gaussian white noise assumption to get correlation confidence bounds
		BandwidthBartlett, #Bandwidth type that uses Bartlett formula to get correlation confidence bounds
		BandwidthP2003 #Bandwidth type implements the bandwidth estimation method of Politis (2003) "Adaptive Bandwidth Choice"


#******************************************************************************


#----------------------------------------------------------
#SET CONSTANTS FOR MODULE
#----------------------------------------------------------
const gaussianPremult = (1 / sqrt(2*pi))::Float64



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
#General dummy type
type KernelDummy <: KernelFunction ; end
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
type KernelBartlett <: KernelFunction
	bandwidth::Float64
	KernelBartlett(bandwidth::Float64) = bandwidth <= 0 ? error("Bandwidth must be strictly positive") : new(bandwidth)
end
KernelBartlett() = KernelBartlett(1.0)
KernelBartlett(bw::Number) = KernelBartlett(convert(Float64, bw))
type KernelTriangular <: KernelFunction #Identical to KernelBartlett
	bandwidth::Float64
	KernelTriangular(bandwidth::Float64) = bandwidth <= 0 ? error("Bandwidth must be strictly positive") : new(bandwidth)
end
KernelTriangular() = KernelTriangular(1.0)
KernelTriangular(bw::Number) = KernelTriangular(convert(Float64, bw))
type KernelEpanechnikov <: KernelFunction
	bandwidth::Float64
	KernelEpanechnikov(bandwidth::Float64) = bandwidth <= 0 ? error("Bandwidth must be strictly positive") : new(bandwidth)
end
KernelEpanechnikov() = KernelEpanechnikov(1.0)
KernelEpanechnikov(bw::Number) = KernelEpanechnikov(convert(Float64, bw))
type KernelGaussian <: KernelFunction
	bandwidth::Float64
	KernelGaussian(bandwidth::Float64) = (bandwidth <= 0) ? error("Bandwidth must be strictly positive") : new(bandwidth)
end
KernelGaussian() = KernelGaussian(1.0) #Constructor for the most common use-case
KernelGaussian(bw::Number) = KernelGaussian(convert(Float64, bw))
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
	bandwidth::Float64 #The upper bound (usually the number of observations)
	p::Float64 #parameter of geometric distribution (1 / expected block length)
	function KernelPR1994SB(bandwidth::Float64, p::Float64)
		bandwidth <= 0 && error("Bandwidth must be strictly positive")
		!(0 < p < 1) && error("Geometric distribution parameter must lie on the interval (0, 1)")
		new(bandwidth, p)
	end
end
KernelPR1994SB() = KernelPR1994SB(1.0, 0.1)
KernelPR1994SB(bandwidth::Number) = KernelPR1994SB(convert(Float64, bandwidth), 0.1)
KernelPR1994SB(bandwidth::Number, p::Number) = KernelPR1994SB(convert(Float64, bandwidth), convert(Float64, p))
#------- METHODS ------
#string method gets string representation of type
string(kT::KernelUniform) = "uniform"
string(kT::KernelBartlett) = "bartlett"
string(kT::KernelTriangular) = "triangular"
string(kT::KernelEpanechnikov) = "epanechnikov"
string(kT::KernelGaussian) = "gaussian"
string(kT::KernelPR1993FlatTop) = "PR1993FlatTop"
string(kT::KernelP2003FlatTop) = "P2003FlatTop"
string(kT::KernelPP2002Trap) = "PP2002Trap"
string(kT::KernelPP2002Smooth) = "PP2002Smooth"
string(kT::KernelPR1994SB) = "PR1994SB"
#show method
function show(io::IO, kF::KernelFunction)
	println(io, "Kernel function = " * string(kF))
	fieldNames = names(kF)
	for n = 1:length(fieldNames)
		println("    field " * string(fieldNames[n]) * " = " * string(getfield(kF, n)))
	end
end
show{T<:KernelFunction}(k::T) = show(STDOUT, k)
#------ EVALUATE method for evaluating kernel function at a given value ------------------------------
#Common Kernel functions
evaluate(x::Number, kT::KernelUniform) = indicator(x, UnitRange(kT.LB, kT.UB)) * (1 / (kT.UB - kT.LB))
evaluate(x::Number, kT::KernelUniform, ::KernelNoIndCheck) = 1 / (kT.UB - kT.LB)
evaluate{T<:Union(KernelBartlett, KernelTriangular)}(x::Number, kT::T) = indicator(x / kT.bandwidth, UnitRange(-1, 1)) * (1 - abs(x / kT.bandwidth))
evaluate{T<:Union(KernelBartlett, KernelTriangular)}(x::Number, kT::T, ::KernelNoIndCheck) = 1 - abs(x / kT.bandwidth)
evaluate(x::Number, kT::KernelEpanechnikov) = indicator(x / kT.bandwidth, UnitRange(-1, 1)) * (0.75 / kT.bandwidth) * (1 - (x^2 / kT.bandwidth^2))
evaluate(x::Number, kT::KernelEpanechnikov, ::KernelNoIndCheck) = (0.75 / kT.bandwidth) * (1 - (x^2 / kT.bandwidth^2))
evaluate(x::Number, kT::KernelGaussian) = (1 / kT.bandwidth) * gaussianPremult * exp(-1 * (x^2 / (2 * kT.bandwidth^2)))
evaluate(x::Number, kT::KernelGaussian, ::KernelNoIndCheck) = evaluate(x, kT)
#Kernel functions from specific papers (note, may not satisfy standard kernel function properties)
#KernelPR1993FlatTop
function evaluate{T<:Number}(x::T, kT::KernelPR1993FlatTop)
	if abs(x) <= kT.m; return(one(T))
	elseif abs(x) <= kT.M; return(1 - (abs(x) - kT.m) / (kT.M - kT.m))
	else; return(zero(T))
	end
end
evaluate(x::Number, kT::KernelPR1993FlatTop, ::KernelNoIndCheck) = evaluate(x, kT)
#KernelP2003FlatTop (equivalent to KernelPR1993FlatTop with m=0.5 and M=1)
function evaluate{T<:Number}(x::T, kT::KernelP2003FlatTop)
	if abs(x) <= 0.5; return(one(T))
	elseif abs(x) <= 1; return(2 * (1 - abs(x)))
	else; return(zero(T))
	end
end
evaluate(x::Number, kT::KernelP2003FlatTop, ::KernelNoIndCheck) = evaluate(x, kT)
#KernelPP2002Trap
function evaluate{T<:Number}(x::T, kT::KernelPP2002Trap)
	if x < 0; return(zero(T))
	elseif x < kT.p1; return(x / kT.p1)
	elseif x < 1 - kT.p1; return(one(T))
	elseif x < 1; return((1 - x) / kT.p1)
	else; return(zero(T))
	end
end
evaluate(x::Number, kT::KernelPP2002Trap, ::KernelNoIndCheck) = evaluate(x, kT)
#KernelPP2002Smooth
function evaluate{T<:Number}(x::T, kT::KernelPP2002Smooth)
	if 0 <= x <= 1; return(1 - abs(2*x - 1)^kT.p1)
	else; return(zero(T))
	end
end
evaluate(x::Number, kT::KernelPP2002Smooth, ::KernelNoIndCheck) = 1 - abs(2*x - 1)^kT.p1
#KernelPR1994SB
evaluate(x::Number, kT::KernelPR1994SB) = indicator(x, UnitRange(0.0, kT.bandwidth)) * ((1 - x/kT.bandwidth)*(1 - kT.p)^x + (x/kT.bandwidth)*(1 - kT.p)^(kT.bandwidth - x))
evaluate(x::Number, kT::KernelPR1994SB, ::KernelNoIndCheck) = (1 - x/kT.bandwidth)*(1 - kT.p)^x + (x/kT.bandwidth)*(1 - kT.p)^(kT.bandwidth - x)
#Array input wrappers
evaluate{T<:Number}(x::AbstractVector{T}, kT::KernelFunction) = [ evaluate(x[n], kT) for n = 1:length(x) ]
evaluate{T<:Number}(x::AbstractMatrix{T}, kT::KernelFunction) = [ evaluate(x[n, m], kT) for n = 1:size(x, 1), m = 1:size(x, 2) ]
#-----activedomain and paramdomain ------------------------
#The purpose of activedomain is to return the domain (of kernel function input) over which the output of the kernel function is active
#The purpose of paramdomain is to return the domain over which the indicated parameter can be defined.
activedomain(kT::KernelUniform) = UnitRange(kT.LB, kT.UB)
activedomain(kT::KernelEpanechnikov) = UnitRange(-kT.bandwidth, kT.bandwidth)
activedomain{T<:Union(KernelBartlett, KernelTriangular)}(kT::T) = UnitRange(-kT.bandwidth, kT.bandwidth)
activedomain(kT::KernelGaussian) = UnitRange(-Inf, Inf)
activedomain(kT::KernelPR1993FlatTop) = UnitRange(-kT.M, kT.M)
activedomain(kT::KernelP2003FlatTop) = UnitRange(-1, 1)
activedomain(kT::KernelPP2002Trap) = UnitRange(0, 1)
activedomain(kT::KernelPP2002Smooth) = UnitRange(0, 1)
activedomain(kT::KernelPR1994SB) = UnitRange(0, kT.bandwidth)
function paramdomain(kT::KernelUniform, pNum::Int)
	pNum == 1 && return(UnitRange(-Inf, kT.UB))
	pNum == 2 && return(UnitRange(kT.LB, Inf))
	error("Invalid parameter number")
end
paramdomain(kT::KernelEpanechnikov, pNum::Int) = pNum == 1 ? UnitRange(nextfloat(0.0), prevfloat(Inf)) : error("Invalid parameter number")
paramdomain{T<:Union(KernelBartlett, KernelTriangular)}(kT::T, pNum::Int) = pNum == 1 ? UnitRange(nextfloat(0.0), prevfloat(Inf)) : error("Invalid parameter number")
paramdomain(kT::KernelGaussian, paramNum::Int) = pNum == 1 ? UnitRange(nextfloat(0.0), prevfloat(Inf)) : error("Invalid parameter number")
function paramdomain(kT::KernelPR1993FlatTop, paramNum::Int)
	paramNum == 1 && return(UnitRange(nextfloat(0.0), Inf))
	paramNum == 2 && return(UnitRange(nextfloat(0.0), Inf))
	error("Invalid parameter number")
end
paramdomain(kT::KernelP2003FlatTop, pNum::Int) = error("Kernel function has no parameters")
paramdomain(kT::KernelPP2002Trap, pNum::Int) = (pNum == 1) ? UnitRange(nextfloat(0.0), 0.5) : error("Invalid parameter number")
paramdomain(kT::KernelPP2002Smooth) = (pNum == 1) ? UnitRange(1, Inf) : error("Invalid parameter number")
function paramdomain(kT::KernelPR1994SB, pNum::Int)
	pNum == 1 && return(UnitRange(0, Inf))
	pNum == 2 && return(UnitRange(nextfloat(0.0), prevfloat(1.0)))
	error("Invalid parameter number")
end
function setbandwidth!(kT::KernelUniform, LB::Number, UB::Number)
	kT.LB = convert(Float64, LB)
	kT.UB = convert(Float64, UB)
end
function setbandwidth!(kT::KernelUniform, bandwidth::Number)
	kT.LB = 0.0
	kT.UB = convert(Float64, bandwidth)
end
function setbandwidth!{T<:Union(KernelEpanechnikov, KernelGaussian, KernelPR1994SB, KernelTriangular, KernelBartlett)}(kT::T, bandwidth::Number)
	kT.bandwidth = convert(Float64, bandwidth)
end
setbandwidth!(kT::KernelFunction, bandwidth::Number) = error("The input kernel function does not have a natural bandwidth measure")





#---------- TYPES FOR DIFFERENT METHODS OF BANDWIDTH ESTIMATION
#Abstract super-type
abstract BandwidthMethod
#General dummy type (useful for update! functions)
type BandwidthDummy <: BandwidthMethod; end
#Dummy type for using maximum possible bandwidth
type BandwidthMax <: BandwidthMethod; end
#type definitions
type BandwidthStockWatsonDefault
	adjustmentTerm::Float64
	BandwidthStockWatsonDefault(adjustmentTerm::Float64) = (adjustmentTerm <= 0) ? error("Adjustment scalar must be strictly greater than zero") : new(adjustmentTerm)
end
BandwidthStockWatsonDefault() = BandwidthStockWatsonDefault(1.0)
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
string(x::BandwidthMax) = "maximumBandwidth"
string(x::BandwidthStockWatsonDefault) = "stockWatsonDefault"
string(x::BandwidthWhiteNoise) = "whiteNoise"
string(x::BandwidthBartlett) = "bartlett"
string(x::BandwidthP2003) = "P2003"
#show method for BandwidthMethod
function show(io::IO, b::BandwidthMethod)
	println(io, "bandwidth method = " * string(b))
	fieldNames = names(b)
	for n = 1:length(fieldNames)
		println("    field " * string(fieldNames[n]) * " = " * string(getfield(b, n)))
	end
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
		!(typeof(kernelFunction) <: Union(KernelUniform, KernelGaussian, KernelEpanechnikov, KernelPR1994SB, KernelBartlett, KernelTriangular)) && error("This kernel function does not have appropriate properties for HAC variance estimation")
		new(kernelFunction, bandwidthMethod)
	end
end
HACVarianceBasic() = HACVarianceBasic(KernelEpanechnikov(), BandwidthP2003())
HACVarianceBasic(numObs::Int) = HACVarianceBasic(KernelEpanechnikov(), BandwidthP2003(numObs))
HACVarianceBasic{T<:Number}(x::Vector{T}) = HACVarianceBasic(KernelEpanechnikov(), BandwidthP2003(length(x)))
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
	return(length(x)-1, ((length(x)-1) / length(x)) * var(x), Array(Float64, 0))
end
function bandwidth{T<:Number}(x::AbstractVector{T}, method::BandwidthStockWatsonDefault)
	length(x) < 2 && error("Input data must have at least two observations")
	return(0.75 * length(x)^(1/3), ((length(x)-1) / length(x)) * var(x), Array(Float64, 0))
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
	setbandwidth!(method.kernelFunction, M)
	kernelAdjTerm = 1 / evaluate(0.0, method.kernelFunction) #Used to scale kernel functions that don't satisfy k(0) = 1
	v = xVar
	noIndCheck = KernelNoIndCheck() #Since we've set the bandwidth to M, we can safely ignore indicator checks
	for m = 1:M
		v += 2 * kernelAdjTerm * evaluate(m, method.kernelFunction, noIndCheck) * xCov[m]
	end
	return(v, M)
end





#Non-exported indicator function for if a value lies in a specified range
indicator{T<:Number}(x::T, LB::T, UB::T) = (LB <= x <= UB) ? 1 : 0
indicator{T<:Number}(x::T, r::UnitRange) = (r.start <= x <= r.stop) ? 1 : 0




end # module
