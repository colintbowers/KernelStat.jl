using KernelStat
using ARIMA
using DependentBootstrap


#-------------------------------
# SIMULATION FUNCTIONS
#-------------------------------
#Simulate kernel functions
const numObsBandwidthSim = 500::Int
const maxCorLag = 6::Int
const kernelList = [KernelUniform(), KernelEpanechnikov(), KernelGaussian(), KernelPR1993FlatTop(), KernelP2003FlatTop(), KernelPP2002Trap(), KernelPP2002Smooth(), KernelPR1994SB(numObsBandwidthSim)]::Vector{KernelFunction}
const bandwidthList = [BandwidthWhiteNoise(), BandwidthBartlett(), BandwidthP2003(numObsBandwidthSim)]
const hacList = [HACVarianceBasic(KernelUniform(0, maxCorLag), BandwidthP2003(numObsBandwidthSim)), HACVarianceBasic(KernelGaussian(1), BandwidthP2003(numObsBandwidthSim)), HACVarianceBasic(KernelEpanechnikov(1), BandwidthP2003(numObsBandwidthSim)), HACVarianceBasic(KernelPR1994SB(numObsBandwidthSim, 0.1), BandwidthMax())]


#------------------------------
# TEST FUNCTIONS
#------------------------------
function testevaluate()
	for k = 1:length(kernelList)
		kF = kernelList[k]
		aD = activedomain(kF)
		if aD.start == -Inf || aD.stop == Inf
			aD.start == -Inf ? (tempStart = -10) : (tempStart = aD.start)
			aD.stop == Inf ? (tempStop = 10) : (tempStop = aD.stop)
			aD = UnitRange(tempStart, tempStop)
		end
		xIn = [aD.start:(aD.stop - aD.start)/10:aD.stop]
		xOut = [ evaluate(xIn[n], kF) for n = 1:length(xIn) ]
		println("Kernel func = " * string(kF) * ". evaluate of active domain = " * string(xOut))
	end
	println("")
	println("")
	return(true)
end


function testbandwidth()
	for k = 0:maxCorLag
		k == 0 ? maCoef = Array(Float64, 0) : maCoef = [ 1 / 2^p for p = 1:k ]
		for j = 1:length(bandwidthList)
			bC = bandwidthList[j]
			mVec = Array(Int, 10)
			for i = 1:10
				x = simulate(numObsBandwidthSim, ARIMAModel(maCoef=maCoef))
				(mVec[i], xVar, xCovVec) = bandwidth(x, bC)
			end
			println("Bandwidth method = " * string(bC))
			println("    MA Order = " * string(k))
			println("    Bandwidth estimates = " * string(mVec))
			println("")
		end
	end
	println("")
	println("")
	return(true)
end


function testhacvariance()
	for k = 0:maxCorLag
		k == 0 ? maCoef = Array(Float64, 0) : maCoef = [ 1 / 2^p for p = 1:k ]
		x = simulate(numObsBandwidthSim, ARIMAModel(maCoef=maCoef))
		meanBoot = dbootstrapstatistic(x, BootstrapParam(x, statistic=mean))
		bootTrueVar = var(sqrt(numObsBandwidthSim) * meanBoot)
		xMean = [ mean(simulate(numObsBandwidthSim, ARIMAModel(maCoef=maCoef))) for i = 1:2000 ]
		simTrueVar = var(sqrt(numObsBandwidthSim) * xMean)
		for j = 1:length(hacList)
			hC = hacList[j]
			vVec = Array(Float64, 10)
			for i = 1:10
				x = simulate(numObsBandwidthSim, ARIMAModel(maCoef=maCoef))
				vVec[i] = hacvariance(x, hC)
			end
			show(hC)
			println("True MA Order = " * string(k))
			println("Simulated true variance estimate = " * string(simTrueVar))
			println("Bootstrapped true variance estimate = " * string(bootTrueVar))
			println("HAC estimates = " * string(vVec))
			println("")
		end
	end
	println("")
	println("")
	return(true)
end


testhacvariance()

function testall(KK::Int)
	testevaluate()
	testbandwidth()
	testhacvariance()
	println("Tests complete.")
end



