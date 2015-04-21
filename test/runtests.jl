using KernelStat
#Load any entire modules that are needed

#Load any specific variables/functions that are needed (use import ModuleName1.FunctionName1, ModuleName2.FunctionName2, etc)
#(none currently needed)

#Specify the variables/functions to export (use export FunctionName1, FunctionName2, etc)
# export 	kernelList, #List of kernel functions types
# 		randomKernel,
# 		randomBandwidth,
# 		testModule,
# 		testKernelFuncBasic,
# 		testBandwidthBasic,
# 		testKernelShape


#-------------------------------
# SIMULATION FUNCTIONS
#-------------------------------
#Simulate kernel functions
#WARNING: The order of kernelList is important for the kernel plotting function
const kernelList = [KernelUniform(), KernelTriangular(), KernelEpanechnikov(), KernelQuartic(), KernelGaussian(), KernelPR1993FlatTop(), KernelP2003FlatTop(), KernelPP2002Trap(), KernelPP2002Smooth()]::Vector{KernelFunction}
randomKernel(kT::KernelUniform) = KernelUniform()
randomKernel(kT::KernelTriangular) = KernelTriangular()
randomKernel(kT::KernelEpanechnikov) = KernelEpanechnikov()
randomKernel(kT::KernelQuartic) = KernelQuartic()
randomKernel(kT::KernelGaussian) = KernelGaussian(abs(randn()))
function randomKernel(kT::KernelPR1993FlatTop)
	m = abs(randn())
	M = m + abs(randn())
	return(KernelPR1993FlatTop(m, M))
end
randomKernel(kT::KernelP2003FlatTop) = KernelP2003FlatTop()
randomKernel(kT::KernelPP2002Trap) =  KernelPP2002Trap(0.5 * rand())
randomKernel(kT::KernelPP2002Smooth) =  KernelPP2002Smooth(abs(randn()) + 1.0)
#Simulate bandwidth methods
const bandwidthList = [BandwidthWhiteNoise(), BandwidthBartlett(), BandwidthP2003()]
randomBandwidth(bW::BandwidthWhiteNoise) = BandwidthWhiteNoise()
randomBandwidth(bW::BandwidthBartlett) = BandwidthBartlett()
randomBandwidth(bW::BandwidthP2003) = BandwidthP2003(1.0, 2 + 0.15 * abs(randn()), convert(Int, round(5 + 0.3 * abs(randn()))))



#------------------------------
# TEST FUNCTIONS
#------------------------------
function testModule(K::Int)
	testKernelFuncBasic(K)
	testBandwidthBasic(K)
	println("All tests passed")
end

function testKernelFuncBasic(K::Int)
	for kk = 1:K
		rK = randomKernel(kernelList[rand(1:length(kernelList))])
		kRange = nonzeroDomain(rK)
		xIn = kRange.start + (kRange.stop - kRange.start) * rand()
		x = evaluate(xIn, rK)
		if typeof(rK) == KernelGaussian
			x <= 0 && error("Invalid output from Gaussian kernel function")
		else
			if !(0 <= x <= 1)
				println(string(rK))
				error("Invalid output from kernel function")
			end
		end
	end
	println("Basic kernel function test passed")
	return(true)
end

function testBandwidthBasic(K::Int)
	for kk = 1:K
		rB = randomBandwidth(bandwidthList[rand(1:length(bandwidthList))])
		N = rand(2:200)
		x = randn(N)
		(M, cVec) = bandwidth(x, rB)
		M < 1 && error("bandwidth output too small")
		if M > N - 1
			println(M)
			println(N-1)
			error("bandwidth output too large")
		end
	end
	println("Basic bandwidth function test passed")
end

# NEED TO BRING IN GADFLY FOR THIS ONE
# function testKernelShape{T<:KernelFunction}(rK::T)
# 	kRange = nonzeroDomain(rK)
# 	if typeof(rK) != KernelGaussian
# 		xVec = [kRange.start:(kRange.stop - kRange.start) * 0.01:kRange.stop]
# 	else
# 		xVec = [-5:0.1:5]
# 	end
# 	kOut = evaluate(xVec, rK)
# 	return(plot(x=xVec, y=kOut, Geom.line, Guide.title(string(rK))))
# end

testModule(1000)
