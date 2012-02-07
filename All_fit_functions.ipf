#pragma rtGlobals=1		// Use modern global access method.

//All simple fit functions that do not depend on other user defined functions

Function/S findCoefList(fitFunction)
	String fitFunction
	
	String coefList=""
	
	strswitch(fitFunction)
		Case "BesselPlusLinear":
			coefList="peakmax;peakpos;slope"
			break
	endswitch
	
	Return coefList
end

//Absolute value of the jinc function
Function AbsJinc(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = amplitude*abs(bessj(1,x*freq)/x);
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = amplitude
	//CurveFitDialog/ w[1] = frequency

	return w[0]*abs(bessj(1,x*w[1])/x);
	
End

//Interferometric signal for a low-finesse cavity with one of the mirrors attached to a simple
//harmonic oscillator.  The function gives the amplitude of reflected light at the oscillator's 
//resonant frequency.  The function also allows for a linear correction to the signal due to 
//the optical lever effect.
Function BesselPlusLinearFit(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = peakmax*abs(bessj(1,x*1.84/peakpos)/.582 +slope*x);
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = peakmax
	//CurveFitDialog/ w[1] = peakpos
	//CurveFitDialog/ w[2] = slope

	return w[0]*abs(bessj(1,x*1.84/w[1])/.582 +w[2]*x);
	
End

Function Gaussian(w,x) : FitFunc
	Wave w
	Variable x
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = offset+amplitude*(1/sqrt(2*pi*stddev^2))*exp(-(x-mean)^2/(2*stddev^2))
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ B
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = amplitude
	//CurveFitDialog/ w[1] = mean
	//CurveFitDialog/ w[2] = stddev
	//CurveFitDialog/ w[3] = offset

	Variable amplitude=w[0]
	Variable meanavg=w[1]
	Variable stddev=w[2] //standard deviation
	Variable offset=w[3]
	
	Return offset+amplitude*(1/sqrt(2*pi*stddev^2))*exp(-(x-meanavg)^2/(2*stddev^2))
end

//Magnetic field penetration depth as a function of temperature from Ginzburg Landau theory
Function GLpenetrationDepth(w,T) : FitFunc
	Wave w
	Variable T

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(T) = lambda0 / (1 - T/Tc)^2
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ B
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = lambda0
	//CurveFitDialog/ w[1] = Tc

	return w[0] / (1 - T / w[1])^2;
	
End

//Coherence length as a function of temperature from Ginzburg-Landau theory
Function GLxi(w,T) : FitFunc
	Wave w
	Variable T

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(T) = xi0 / sqrt(1 - (T/Tc))
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ B
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = xi0
	//CurveFitDialog/ w[1] = Tc

	return w[0] / sqrt(1 - (T/w[1]));
	
End

//Critical magnetic field as a function of temperature from Ginzburg Landau theory
Function HcGL(w,temp): FitFunc

	Wave w
	Variable temp
	
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(temp) = (3^.5*6.626e-34/pi/2/1.6e-19)*(1-temp/Tc)^.5/linewidth/xi0
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ temp
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = xi0
	//CurveFitDialog/ w[1] = Tc
	//CurveFitDialog/ w[2] = linewidth
	
	Variable xi0 = w[0]
	Variable Tc = w[1]
	Variable lw = w[2] //Wire linewidth (80e-9 for October 2008 cooldown)
	
	Variable hplanck=6.626e-34
	Variable echarge=1.6e-19
	
	If(temp<Tc)
		Return (3^.5)*(hplanck/(2*echarge))*(1/(lw*pi*xi0))*(1-temp/Tc)^.5
	else
		Return 0
	endif
end

Function HcGLD(w,temp): FitFunc

	Wave w
	Variable temp
	
	Variable lw = 115e-9 //Wire linewidth
	
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(temp) = something
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ temp
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = D
	//CurveFitDialog/ w[1] = Tc
	
	Variable D = w[0]
	Variable Tc = w[1]
	
	Variable hplanck=6.626e-34
	Variable hbar=hplanck/2/pi
	Variable echarge=1.6e-19
	Variable phi0=hplanck/echarge/2
	Variable kBoltzmann=1.38e-23
	
	If(temp<Tc)
		Return ((sqrt(24*kBoltzmann)*phi0)/(pi*sqrt(pi*hbar*D)*lw))*(Tc-temp)^.5
	else
		Return 0
	endif
end

Function JincPlusLinear(w,x): FitFunc

	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = gamma+alpha*bessj(1,beta*(x-epsilon))/(x-epsilon)+delta*(x-epsilon)
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = alpha
	//CurveFitDialog/ w[1] = beta
	//CurveFitDialog/ w[2] = gamma
	//CurveFitDialog/ w[3] = delta
	//CurveFitDialog/ w[4] = epsilon
	
	Return w[2]+w[0]*bessj(1,w[1]*(x-w[4]))/(w[1]*(x-w[4]))+w[3]*(x-w[4])
end

Function JincPlusLinear2(w,x): FitFunc

	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = gamma+alpha*bessj(1,beta*(x-epsilon))/(x-epsilon)+delta*(x-epsilon)
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = alpha
	//CurveFitDialog/ w[1] = beta
	//CurveFitDialog/ w[2] = gamma
	//CurveFitDialog/ w[3] = delta
	//CurveFitDialog/ w[4] = epsilon
	
	Return w[2]+w[0]*bessj(1,w[1]*(x-w[4]))/((x-w[4]))+w[3]*(x-w[4])
end

Function Lorentzian(w,frequency) : FitFunc
	Wave w
	Variable frequency

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(frequency) = offset+(Amp/Q^2)/((1-(frequency/ResFreq)^2)^2+(frequency/(ResFreq*Q))^2);
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ frequency
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = Amp
	//CurveFitDialog/ w[1] = Q
	//CurveFitDialog/ w[2] = ResFreq
	//CurveFitDialog/ w[3] = offset

	return w[3]+(w[0]/w[1]^2)/((1-(frequency/w[2])^2)^2+(frequency/(w[2]*w[1]))^2);
End

//Electron phase coherence length as a function of temperature with both the 
//electron-phonon and electron-electron contributions.
Function LphivsT(w,T) : FitFunc

	Wave w
	Variable T

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ Lphi(T) = (a0*T+a1*T^3)^(-1/2)
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ T
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = a0
	//CurveFitDialog/ w[1] = a1
	
	Variable a0 = w[0]
	Variable a1 = w[1]

	Return (a0*T^(2/3)+a1*T^3)^(-1/2)
end

Function LphivsT2(w,T) : FitFunc

	Wave w
	Variable T

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ Lphi(T) = ((Aee/D)*T^(2/3)+(Aep/D)*T^3)^(-1/2)
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ T
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = Aep
	//CurveFitDialog/ w[1] = D
	
	Variable Aep = w[0]
	Variable D = w[1]
	
	//Old values
	//Variable Rsquare=.082
	//Variable linew=80e-9
	
	Variable Rsquare=.117
	Variable linew=115e-9
	
	Variable kBoltzmann=1.3806503e-23
	Variable hplanck=6.626068e-34
	Variable echarge=1.60217646e-19
	
	//Old Wind version of the expression
	//Variable Aee=(Rsquare*echarge^2/sqrt(2)/(hplanck/2/pi))*(kBoltzmann/(hplanck/2/pi))*(sqrt(D)/linew)
	//More recent Altshuler version of the expression
	Variable Aee=(Rsquare*echarge^2/4/(hplanck/2/pi))*(kBoltzmann/(hplanck/2/pi))*(sqrt(D)/linew)
	Aee=Aee^(2/3)

	Return ((Aee/D)*T^(2/3)+(Aep/D)*T^3)^(-1/2)
end

//Temperature dependence of persistent currents according to the Riedel and von Oppen
//PRB from 1993.
Function MatsubaraIvsT(w,T) : FitFunc
	Wave w
	Variable T
	
	Variable summingTerms=20
	
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(T) = I0*sqrt((pi^6/3)*(T/Tc)^2*Sum(n*exp(-(2*pi^3*n*(T/Tc))^.5))+y_offset
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ T
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = I0
	//CurveFitDialog/ w[1] = Tc
	//CurveFitDialog/ w[2] = y_offset
	
	Variable I0=w[0],Tc=w[1],y_offset=w[2]
	Variable littleT=T/Tc
	Variable summation
	Variable n
	
	for(n=1;n<summingTerms+1;n+=1)
		summation+=n*exp(-(2*pi^3*n*littleT)^(1/2))
	endfor
	
	Variable MatsuTerm=I0*sqrt((pi^6/3)*littleT^2*summation)
	
	Return MatsuTerm+abs(y_offset)
end

Function Parabola(w,x) : FitFunc

	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0 + a0*(x-x0)^2
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = a0
	//CurveFitDialog/ w[2] = x0
	
	Variable y0 = w[0]
	Variable a0 = w[1]
	Variable x0 = w[2]

	Return y0 + a0*(x-x0)^2
end

//Describes the probability density of the amplitude of a phasor 
//which is characterized by a symmetric, two dimensional Gaussian
//centered at the origin.  For a normalized probability density, the amplitude
//should be 1 and the width is the standard deviation of the Gaussian.
//(This function just multiplies the Gaussian by 2 pi r to do the angular 
//integration).
Function RayleighDist(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = amplitude*x*exp(-x^2/(2*width^2))/width^2
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = amplitude
	//CurveFitDialog/ w[1] = width
	
	Variable amplitude=w[0]
	Variable width=w[1]

	Return amplitude*x*exp(-x^2/(2*width^2))/width^2
end

//Fit for an exponential decay in time with an offset added in quadrature.
//Supposed to represent a cantilever ringdown with noise but I am not sure
//this is the right model any more.
Function ringdownFit(w,etime) : FitFunc
	Wave w
	Variable etime

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(etime) = sqrt((amplitude*exp(-etime/tau))^2+(noise)^2);
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ etime
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = amplitude
	//CurveFitDialog/ w[1] = tau
	//CurveFitDialog/ w[2] = noise

	return sqrt((w[0]*exp(-etime / w[1]))^2+(w[2])^2)
End

//This might be a better fit for the exponential decay of an oscillator with noise.
//It assumes equal magnitudes of in-phase and out-of-phase noise added in
//quadrature with the ringdown signal.
Function ringdownFit2(w,etime) : FitFunc
	Wave w
	Variable etime

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(etime) = sqrt((amplitude*exp(-etime/tau))^2+(noise)^2);
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ etime
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = amplitude
	//CurveFitDialog/ w[1] = tau
	//CurveFitDialog/ w[2] = noise

	return sqrt((w[0]*exp(-etime / w[1]))^2+2*w[2]*w[0]*exp(-etime / w[1])+2*(w[2])^2)
End

//This was used to fit some thermometry data to check the calibration curves.
//There was a numerical error which was causing the conversion entered into Igor
//to give inaccurate values.  Not sure if  this function would ever be useful again.
Function RXChebyFit(w,R)

	Wave w
	Variable R
	
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a0+a1+a2+a3+a4+a5+a6+a7+a8
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ R
	//CurveFitDialog/ Coefficients 9
	//CurveFitDialog/ w[0] = a0
	//CurveFitDialog/ w[1] = a1
	//CurveFitDialog/ w[2] = a2
	//CurveFitDialog/ w[3] = a3
	//CurveFitDialog/ w[4] = a4
	//CurveFitDialog/ w[5] = a5
	//CurveFitDialog/ w[6] = a6
	//CurveFitDialog/ w[7] = a7
	//CurveFitDialog/ w[8] = a8
	
	Variable ZL = 2.955
	Variable ZU = 3.10855552727
	Variable order = 10
	Variable Z = log(R)
	
	Variable x = ((Z-ZL)-(ZU-Z))/(ZU-ZL)
	
	Variable i,T
	for(i=0,T=0;i<order;i+=1)
		T+= w[i] * chebyshev(i,x)
	endfor
	
	Return T
end

// Simple harmonic oscillator: amplitude versus frequency
Function SHOFit(w,frequency) : FitFunc
	Wave w
	Variable frequency

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(frequency) = Amp/sqrt((fres^2-frequency^2)^2+(fres*frequency/Q)^2)+offset
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ frequency
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = Amp
	//CurveFitDialog/ w[1] = fres
	//CurveFitDialog/ w[2] = Q
	//CurveFitDialog/ w[3] = offset

	return w[0]/sqrt((w[1]^2-frequency^2)^2+(w[1]*frequency/w[2])^2)+w[3]
End

Function SHOpsd(w,frequency) : FitFunc
	Wave w
	Variable frequency

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(frequency) = offset+(2/pi)*(xrmsSquared/(ResFreq*Q))/((1-(frequency/ResFreq)^2)^2+(frequency/(ResFreq*Q))^2);
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ frequency
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = xrmsSquared
	//CurveFitDialog/ w[1] = Q
	//CurveFitDialog/ w[2] = ResFreq
	//CurveFitDialog/ w[3] = offset

	return w[3]+(2/pi)*(w[0]/(w[1]*w[2]))/((1-(frequency/w[2])^2)^2+(frequency/(w[2]*w[1]))^2);
End