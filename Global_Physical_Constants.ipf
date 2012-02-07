#pragma rtGlobals=1		// Use modern global access method.

Function eCharge()
	Return 1.602176487e-19
end

Function hbar()
	Return hplanck()/(2*pi)
end

Function hplanck()
	Return 6.62606896e-34
end

Function kB()
	Return 1.3806504e-23
end

Function mu0()
	Return pi*4e-7
end

Function muB()
	Return 9.27400915e-24
end

Function phi0()
	Return hplanck()/eCharge()
end