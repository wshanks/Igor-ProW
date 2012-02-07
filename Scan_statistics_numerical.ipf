#pragma rtGlobals=1		// Use modern global access method.

//Calculates the argument of complex number x+i*y
Function ArgXY(x,y)

	Variable x,y
	
	if (y == 0)
		if (x == 0)
			Return NaN
		elseif (x < 0)
			Return pi
		endif
	endif
	
	Return 2*atan(y/(sqrt(x^2+y^2)+x))
end

//Autocorrelate() creates a wave containing the autocorrelation of targetWave.
//Autocorrelate() returns the name of the created wave (which is the name of 
//targetWave with "_AC" appended to it).
//
//The linear correlation, targetWave_AC[j], between targetWave 
//and itself shifted by j is calculated with weighting given by the overlap between 
//targetWave and its shifted version.  That is, 
//targetWave_AC[j]=Sum(targetWave[i]*targetWave[i+j])/(numpnts(targetWave)-j) 
//with i running from 0 to numpnts(targetWave)-j-1.
//targetWave_AC is given the same scaling as targetWave.  A note is attached 
//recording the location of targetWave.
Function/S Autocorrelate(targetWave)

	Wave targetWave
	
	String ACwaveName=nameofwave(targetWave)+"_AC"
	
	Duplicate/O targetWave, $ACwaveName
	Wave ACwave=$ACwaveName
	
	ACwave=0
	Variable i,j
	for(j=0;j<numpnts(targetWave);j+=1)
		for(i=0;i<numpnts(targetWave)-j;i+=1)
			ACwave[j]+=targetWave[i]*targetWave[i+j]/(numpnts(targetWave)-j)
		endfor
	endfor
	
	SetScale/P x 0,DimDelta(targetWave,0),ACwave
	
	String NoteStr="GENERATEDBY:Autocorrelate();SRCWAVE:"
	NoteStr+=GetWavesDataFolder(targetWave,2)+";"
	
	Note ACwave, NoteStr
	
	Return ACwaveName
end

//Saves sections of all the waves in the current data folder into subfolders.
//This function was created to take data files consisting of many scans and
//break them up into waves for each scan.  It assumes that the first point in 
//each wave is the first point of the scan (it doesn't search for a turning point).
//Also sets the wave scaling of the wave sections to match the scan 
//(assumes the scan points are equally spaced).
//
//The subfolder are created with names determined by currentVal+endNameTag
//where currentVal is startVal+stepVal*(the number of the scan).  All the scan
//waves saved to the subfolder also have this tag appended to the end of their
//names.  (eg if taking scans while stepping an RF power from 0 dBm, to 
//2 dBm to 4 dBm etc, startVal=0,stepVal=2,endNameTag="dBm").  They also 
//have a "p" or "n" added to the beginning of this end tag depending on whether 
//the scan was positive or negative.
//
//scanWaveTag -- partial wave name used to identify the scan wave from 
//which the scan step and direction are extracted.
//
//The section lengths are singleScanLength.
//
//This function could be improved to cut off the incomplete scan at the 
//beginning of the data and to find the scan length itself.
Function breakUpSeriesStep(scanWaveTag,singleScanLength,startVal,stepVal,endNameTag)

	Variable singleScanLength,startVal,stepVal
	String endNameTag,ScanWaveTag
	
	Variable totWaves = CountObjects("",1)
	
	String tempWaveName=GetIndexedObjName("",1,0)
	Wave tempWave=$tempWaveName
	Variable totScans=floor(numpnts(tempWave)/singleScanLength)
	
	String origFolder = GetDataFolder(1);
	String tempScanTag
	String tempWaveTag, tempScanWaveName,tempTestWaveTag
	Variable i,j
	Variable scanDelta,scanStart
	for(i=0;i<totScans;i+=1)
		//Trim big scan to appropriate data range
		tempScanTag=num2str(startVal+stepVal*i)+endNameTag
		DupWaves2Folder(tempScanTag)
		SetDataFolder $tempScanTag
		For(j=0;j<totWaves;j+=1)
			tempWaveName=GetIndexedObjName("",1,j)
			trimWave($tempWaveName,i*singleScanLength,(i+1)*singleScanLength-1)
		endfor
		
		//Find the scan wave
		For(j=0;j<totWaves;j+=1)
			tempWaveName=GetIndexedObjName("",1,j)
			tempTestWaveTag=tempWaveName[strlen(tempWaveName)-strlen(ScanWaveTag),strlen(tempWaveName)-1]
			if (stringmatch(ScanWaveTag,tempTestWaveTag))
				tempScanWaveName=tempWaveName
			endif
		endfor
		
		//Reverse negative scans
		Wave tempScanWave=$tempScanWaveName
		if (FindScanDirection(tempScanWave)<0)
			For(j=0;j<totWaves;j+=1)
				tempWaveName=GetIndexedObjName("",1,j)
				Reverse/P $tempWaveName
			endfor
		endif
		
		//Set wave scaling
		scanStart=tempScanWave[0]
		//Find scan step by finding the first point in the scan wave different from the 
		//zeroth point and then looking at the difference in value between the two points.
		scanDelta=tempScanWave[findSecLength(tempScanWave,0)]-tempScanWave[0]
		For(j=0;j<totWaves;j+=1)
			tempWaveName=GetIndexedObjName("",1,j)
			SetScale/P x, scanStart, scanDelta, $tempWaveName
		endfor
		
		//Add tag to end of waves
		//Mark positive or negative scans
		if((startVal+stepVal*i)<0)
			tempWaveTag="n"
		else
			tempWaveTag="p"
		endif
		tempWaveTag+=tempScanTag+endNameTag
		AddTagAllWaves(tempWaveTag)
		
		SetDataFolder $origFolder
	endfor
end

//This function is for scans that don't double up at the turnaround points.
//
//Saves sections of all the waves in the current data folder into subfolders.
//This function was created to take data files consisting of many scans and
//break them up into waves for each scan.  It assumes that the first point in 
//each wave is the first point of the scan (it doesn't search for a turning point).
//Also sets the wave scaling of the wave sections to match the scan 
//(assumes the scan points are equally spaced).
//
//The subfolder are created with names determined by currentVal+endNameTag
//where currentVal is startVal+stepVal*(the number of the scan).  All the scan
//waves saved to the subfolder also have this tag appended to the end of their
//names.  (eg if taking scans while stepping an RF power from 0 dBm, to 
//2 dBm to 4 dBm etc, startVal=0,stepVal=2,endNameTag="dBm").  They also 
//have a "p" or "n" added to the beginning of this end tag depending on whether 
//the scan was positive or negative.
//
//scanWaveTag -- partial wave name used to identify the scan wave from 
//which the scan step and direction are extracted.
//
//The section lengths are singleScanLength.
//
//This function could be improved to cut off the incomplete scan at the 
//beginning of the data and to find the scan length itself.
Function breakUpSeriesStep2(,ScanWaveTag,singleScanLength,startVal,stepVal,endNameTag)

	//This function is for scans that don't double up at the turnaround points.
	Variable singleScanLength,startVal,stepVal
	String endNameTag,ScanWaveTag
	
	Variable totWaves = CountObjects("",1)
	
	String tempWaveName=GetIndexedObjName("",1,0)
	Wave tempWave=$tempWaveName
	//Different from when the scan doubles up at the turnaround points
	Variable totScans=floor((numpnts(tempWave)-1)/(singleScanLength-1))
	
	String origFolder = GetDataFolder(1);
	String tempFolderName
	String tempWaveTag, tempScanWaveName,tempTestWaveTag
	Variable i,j
	Variable scanDelta,scanStart
	for(i=0;i<totScans;i+=1)
		//Trim big scan to appropriate data range
		tempFolderName=num2str(startVal+stepVal*i)+endNameTag
		DupWaves2Folder(tempFolderName)
		SetDataFolder $tempFolderName
		For(j=0;j<totWaves;j+=1)
			tempWaveName=GetIndexedObjName("",1,j)
			//Different from case with doubled up turnaround points
			trimWave($tempWaveName,i*(singleScanLength-1),(i+1)*(singleScanLength-1))
		endfor
		
		//Find the scan wave
		For(j=0;j<totWaves;j+=1)
			tempWaveName=GetIndexedObjName("",1,j)
			tempTestWaveTag=tempWaveName[strlen(tempWaveName)-strlen(ScanWaveTag),strlen(tempWaveName)-1]
			if (stringmatch(ScanWaveTag,tempTestWaveTag))
				tempScanWaveName=tempWaveName
			endif
		endfor
		
		//Reverse negative scans
		Wave tempScanWave=$tempScanWaveName
		if (FindScanDirection(tempScanWave)<0)
			For(j=0;j<totWaves;j+=1)
				tempWaveName=GetIndexedObjName("",1,j)
				Reverse/P $tempWaveName
			endfor
		endif
		
		//Set wave scaling
		scanStart=tempScanWave[0]
		scanDelta=tempScanWave[findSecLength(tempScanWave,0)]-tempScanWave[0]
		For(j=0;j<totWaves;j+=1)
			tempWaveName=GetIndexedObjName("",1,j)
			SetScale/P x, scanStart, scanDelta, $tempWaveName
		endfor
		
		//Add tag to end of waves
		if((startVal+stepVal*i)<0)
			tempWaveTag="n"
		else
			tempWaveTag="p"
		endif
		tempWaveTag+=num2str(abs(startVal+stepVal*i))+endNameTag
		AddTagAllWaves(tempWaveTag)
		
		SetDataFolder $origFolder
	endfor
end

//Calculates the moment of targetWave specified by moment
//ie the average of targetWave^moment
Function calcMomentData(targetWave,moment)

	Wave targetWave
	Variable moment
	
	Variable momentVal,i
	for(i=0,momentVal=0;i<numpnts(targetWave);i+=1)
		momentVal+=targetWave[i]^moment
	endfor
	momentVal/=numpnts(targetWave)
	
	Return momentVal
end

//Treats targetWave as an unnormalized distribution function and calculates
//the moment of the distribution of order moment.  Uses the scaling of 
//targetWave to sum up x^moment*targetWave(x) over the whole wave
//and then normalizes by the sum of targetWave(x) over the whole wave.
Function calcMomentDFunc(targetWave,moment)

	Wave targetWave
	Variable moment
	
	Variable i,sum1,norm1
	for(i=0,sum1=0,norm1=0;i<numpnts(targetWave);i+=1)
		sum1+=((leftx(targetWave)+i*deltax(targetWave))^moment)*targetWave[i]
		norm1+=targetWave[i]
	endfor
	
	Return sum1/norm1
end
//calcMoment() on only a subset of targetWave
Function calcMomentDFuncSubset(targetWave,moment,firstPoint,lastPoint)

	Wave targetWave
	Variable moment,firstPoint,lastPoint
	
	Variable i,sum1,norm1
	for(i=firstPoint,sum1=0,norm1=0;i<lastPoint+1;i+=1)
		sum1+=((leftx(targetWave)+i*deltax(targetWave))^moment)*targetWave[i]
		norm1+=targetWave[i]
	endfor
	
	Return sum1/norm1
end

//Calculates the moments as estimated from the variance assuming a Gaussian 
//distribution centered at zero.
Function calcMomentGaussian(targetWave,moment)

	Wave targetWave
	Variable moment
	
	if (mod(moment,2)!=0)
		Return 0
	endif
	
	Return calcMomentGaussianFV(calcMomentData(targetWave,2),moment)
end

//Calculates the moments as estimated from the variance assuming a Gaussian 
//distribution centered at zero with input variance.  (FV="from variance")
Function calcMomentGaussianFV(variance,moment)

	Variable variance,moment
	
	if (mod(moment,2)!=0)
		Return 0
	endif
	
	Return variance^(moment/2)*doublefactorial(moment-1)
end

//Calculates the normalized standard error of the moment for a normal/Gaussian distribution
//Assume indSamp is the number of independent samples
Function calcMomentGaussianSENorm(moment,indSamp)

	Variable moment,indSamp
	
	Return sqrt(1/indSamp)*sqrt(doublefactorial(2*moment-1)/doublefactorial(moment-1)^2-1)
end

//Creates an output wave with the first totMoments moments of targetWave
//as calculated using the raw second moment and the formula for moments 
//of a Gaussian distribution
Function/S calcMomentGWave(targetWave,totMoments)

	Wave targetWave
	Variable totMoments
	
	String momWaveName=NameOfWave(targetWave)+"_GM"
	Make/D/O/N=(totMoments+1) $momWaveName
	Wave momWave=$momWaveName
	
	momWave[0]=0
	Variable i
	for(i=1;i<totMoments+1;i+=1)
		momWave[i]=calcMomentGaussian(targetWave,i)
	endfor
	
	Note momWave,"Gaussian moments of "+GetWavesDataFolder(targetWave,2)
	Note momWave,"(Moments calculated using raw 2nd moment and assuming a Gaussian distribution"
	
	Return momWaveName
end

//Creates an output wave with the first totMoments raw moments of targetWave
Function/S calcMomentWave(targetWave,totMoments)

	Wave targetWave
	Variable totMoments
	
	String momWaveName=NameOfWave(targetWave)+"_M"
	Make/D/O/N=(totMoments+1) $momWaveName
	Wave momWave=$momWaveName
	
	momWave[0]=0
	Variable i
	for(i=1;i<totMoments+1;i+=1)
		momWave[i]=calcMomentData(targetWave,i)
	endfor
	
	Note momWave,"Raw moments of "+GetWavesDataFolder(targetWave,2)
	
	Return momWaveName
end

//Creates an output wave with the square root of the variance of the 
//first totMoments raw moments of targetWave
Function/S calcMomentVarWave(targetWave,totMoments)

	Wave targetWave
	Variable totMoments
	
	String momUWaveName=NameOfWave(targetWave)+"_MU"
	Make/D/O/N=(totMoments+1) $momUWaveName
	Wave momUWave=$momUWaveName
	
	momUWave[0]=0
	Variable i
	for(i=1;i<totMoments+1;i+=1)
		momUWave[i]=sqrt((1/numpnts(targetWave))*(calcMomentData(targetWave,2*i)-calcMomentData(targetWave,i)^2))
	endfor
	
	String momUNote="Variance of the raw moments of "
	momUNote+=GetWavesDataFolder(targetWave,2)
	Note momUWave,momUNote
	
	Return momUWaveName
end

//Creates an output wave with the square root of the variance of the 
//first totMoments raw moments of targetWave.  The variance is 
//scaled by correlationLength/deltax(targetWave) in order to account 
//for the finite correlation length between the data points (ie the variance 
//is scaled by the effective number of measurements (data range / correlation 
//length) rather than actual number of points in targetWave).  correlationLength
//should be in the same units as the point spacing of targetWave.
Function/S calcMomentVarCorWave(targetWave,totMoments,correlationLength)

	Wave targetWave
	Variable totMoments,correlationLength
	
	String momUName=calcMomentVarWave(targetWave,totMoments)
	Wave momU=$momUName
	
	momU*=sqrt(correlationLength/deltax(targetWave))
	
	String momUnote
	sprintf momUnote,"Correlation length %g",correlationLength
	Note momU,momUnote
	
	Return momUName
end

//Creates an output wave with the square root of the variance of the 
//first totMoments raw moments of targetWave.  In calculating the 
//higher order moments to find the variances, only the second order moment 
//of the data is used (it is taken to the appropriate power and multiplied by the 
//appropriate coefficient).  The variance is 
//scaled by correlationLength/deltax(targetWave) in order to account 
//for the finite correlation length between the data points (ie the variance 
//is scaled by the effective number of measurements (data range / correlation 
//length) rather than actual number of points in targetWave).  correlationLength
//should be in the same units as the point spacing of targetWave.
Function/S calcMomentGauVarCorWave(targetWave,totMoments,correlationLength)

	Wave targetWave
	Variable totMoments,correlationLength
	
	String momUWaveName=NameOfWave(targetWave)+"_GMU"
	Make/D/O/N=(totMoments+1) $momUWaveName
	Wave momUWave=$momUWaveName
	
	momUWave[0]=0
	Variable i
	for(i=1;i<totMoments+1;i+=1)
		momUWave[i]=sqrt((1/numpnts(targetWave))*(calcMomentGaussian(targetWave,2*i)-calcMomentGaussian(targetWave,i)^2))
	endfor
	
	momUWave*=sqrt(correlationLength/deltax(targetWave))
	
	String momUNote="Variance of the raw moments of "
	momUNote+=GetWavesDataFolder(targetWave,2)
	momUNote+=" assuming a normal distribution (higher moments calculated as multiples of "
	momUNote+=" powers of the variance)"
	Note momUWave,momUNote
	
	sprintf momUnote,"Correlation length %g",correlationLength
	Note momUWave,momUnote
	
	Return momUWaveName
end

//Calculates moment of Rayleigh distribution with the same first moment as targetWave
Function calcMomentRayleigh(targetWave,moment)

	Wave targetWave
	Variable moment
	
	Variable sig=sqrt(2/pi)*calcMomentData(targetWave,1)
	
	Return calcMomentRayleighSig(sig,moment)
end

//Calculates moment of Rayleigh distribution with the same first moment as targetWave
Function calcMomentRayleighSig(sigma,moment)

	Variable sigma, moment
	
	Return sigma^moment*2^(moment/2)*gamma(1+moment/2)
end

//Calculates the variance in the estimated moment for a Rayleigh distribution with mean m1
//from samples number of independent measurements
Function calcMomentRayleighVar(m1,moment,samples)

	Variable m1,moment,samples
	
	Variable sigma=sqrt(2/pi)*m1
	
	Return (1/samples)*(calcMomentRayleighSig(sigma,2*moment)-calcMomentRayleighSig(sigma,moment)^2)
end

//Analyzes some statistics of the values of targetWave.  Specifically the function
//calculates the first four moments and the skewness and kurtosis (excess) 
//coefficients.  Also, calculates the Kolmogorov-Smirnov test statistic assuming 
//the distribution set by distribution ID and paramWave (see KSTestStatistic() ).
//
//Returns all of this information as a keyword indexed string (see code for keywords)
//along with the standard deviation.
//Also appends the string to targetWave's note if notationBit == 1.
Function/S checkGaussianity(targetWave,paramWave,distributionID,notationBit)

	Wave targetWave,paramWave
	Variable distributionID,notationBit
	
	Variable m1=calcMomentData(targetWave,1)
	Variable m2=calcMomentData(targetWave,2)
	Variable m3=calcMomentData(targetWave,3)
	Variable m4=calcMomentData(targetWave,4)
	Variable skewness=m3/m2^1.5
	Variable excessCoefficient=m4/m2^2-3
	Variable KStstat=KSTestStatistic(targetWave,distributionID,paramWave)
	
	//Output string
	String gaussanityStr=""
	sprintf gaussanityStr, gaussanityStr+"avg: %g;",m1
	sprintf gaussanityStr, gaussanityStr+"stdDev: %g;",sqrt(m2)
	sprintf gaussanityStr, gaussanityStr+"Moment2: %g;",m2
	sprintf gaussanityStr, gaussanityStr+"Moment3: %g;",m3
	sprintf gaussanityStr, gaussanityStr+"Moment4: %g;",m4
	sprintf gaussanityStr, gaussanityStr+"Skewness: %g;",skewness
	sprintf gaussanityStr, gaussanityStr+"Kurtosis: %g;",excessCoefficient
	sprintf gaussanityStr, gaussanityStr+"KSTestStatistic: %g;",KStstat
	
	if (notationBit==1)
		Note targetWave, gaussanityStr
	endif
	Return gaussanityStr
end

Function checkPtExists(targetWave,value,tolerance)
	//Search sequentially through targetWave and return the index of the first point
	//within tolerance of value
	//Return -1 if there is no such point.
	Wave targetWave
	Variable value,tolerance
	
	Variable i,vIndex=-1
	for(i=0;i<numpnts(targetWave);i+=1)
		if (abs(targetWave[i]-value)<tolerance)
			vIndex=i
			break
		endif
	endfor
	Return vIndex
end

// Determines how many loops are in a frequency wave from FFTLoop (where the frequency data
// is repeated over and over for each loop in the file)
Function CountLoops(loopcountwave)

	Wave loopcountwave;
	
	Variable i,psdlength;

	for(i=0;i<numpnts(loopcountwave);i+=1)
		psdlength = i + 1;
		if (loopcountwave[i] > loopcountwave[i+1])
			break;
		endif
	endfor
	
	Variable loops = (numpnts(loopcountwave) / psdlength);
	Return loops;

end

//Returns the number of oscillations in wiggleWave about the value midpoint
//with a tolerance of threshold.
//
//The function steps through the wave and adds one to its count each time 
//the wave passes reaches a point above midpoint-threshold following a point 
//below midpoint-threshold
Function CountWiggles(wiggleWave,midpoint,threshold)

	Wave wiggleWave
	Variable midpoint,threshold
	
	Variable wiggleCount=0,i,updown=1
	//updown keeps track of whether the function is looking for the wave
	//to go up or down.  1=wave has crossed the low threshold is and should
	//go up soon.  -1=wave has crossed the upper threshold and should go 
	//down soon.
	
	//step through wave
	for(i=0;i<numpnts(wiggleWave);i+=1)
		switch (updown)
			case 1: 
				if (wiggleWave[i]>threshold+midpoint)
					//crossing upper threshold
					wiggleCount+=1
					updown=-1
				endif
				break;
			case -1:
				if (wiggleWave[i]<midpoint-threshold)
					//crossing lower threshold
					updown=1
				endif
				break;
		endswitch
	endfor
	
	return wiggleCount
end

//Must be included by ringdown loading procedure.
// Create a histogram of the data in dataWave called histName.  Make enough bins to have
// bins with width of have the standard deviation of the data covering the entire range of the data.
Function CreateHist(dataWave, histName)

	Wave dataWave;
	String histName;
	
	Variable numBins;
	
	WaveStats/C=1/W dataWave;
	Wave M_WaveStats;
	
	numBins = ceil(2 * (M_WaveStats[12] - M_WaveStats[10]) / M_WaveStats[4] + .0001);
	
	Make/N=1 w1;	
	Histogram/B={M_WaveStats[10], M_WaveStats[4] / 2,numBins} dataWave, w1;

	Duplicate/O w1 $histName;
	KillWaves w1;

end

//Crosscorrelate() creates a wave containing the cross correlation of targetWave1
//and targetWave2.
//Crosscorrelate() returns the name of the created wave (which is the name of 
//targetWave1 with "_CC" appended to it).
//
//The linear correlation, targetWave1_CC[j], between targetWave1
//and targetWave2 lagged by j is calculated with weighting given by the overlap 
//between targetWave1 and targetWave2 for lag j.  That is, 
//targetWave1_CC[j]=Sum(targetWave1[i]*targetWave2[i+j])/(imax-imin-1) 
//with i running from imin=max(0,-j) to 
//imax=min(numpnts(targetWave1)-1,numpnts(targetWave1)-1-j)
//targetWave1_CC is given the same scaling as targetWave1.  A note is attached 
//recording the full path of targetWave1.
Function/S Crosscorrelate(targetWave1,targetWave2)

	Wave targetWave1,targetWave2
	
	String CCwaveName=nameofwave(targetWave1)+"_CC"
	
	Make/O/N=(2*numpnts(targetWave1)-1) $CCwaveName
	Wave CCwave=$CCwaveName
	
	CCwave=0
	Variable i,j,imin,imax
	for(j=-numpnts(targetWave1)+1;j<numpnts(targetWave1);j+=1)
		imin=max(0,-j)
		imax=min(numpnts(targetWave1)-1,numpnts(targetWave1)-1-j)
		for(i=imin;i<imax+1;i+=1)
			CCwave[j+numpnts(targetWave1)-1]+=targetWave1[i]*targetWave2[i+j]/(imax-imin+1)
		endfor
	endfor
	
	Variable CCmin=DimDelta(targetWave1,0)*(1-numpnts(targetWave1))
	SetScale/P x CCmin,DimDelta(targetWave1,0),CCwave
	
	String NoteStr="GENERATEDBY:Autocorrelate();SRCWAVE:"
	NoteStr+=GetWavesDataFolder(targetWave1,2)+";"
	
	Note CCwave, NoteStr
	
	Return CCwaveName
end

//Calculates the double factorial integerN!!
Function doublefactorial(integerN)

	Variable integerN
	
	Variable dfact=1,i
	for(i=integerN;i>0;i-=2)
		dfact*=i
	endfor
	
	Return dfact
end

Function FindFitEndPoint(xWave,lowxValue,highxValue)

	Wave xWave
	Variable lowxValue,highxValue
	
	Variable temp
	if (highxValue < lowxValue)
		temp = lowxValue
		lowxValue = highxValue
		highxValue = temp
	endif
	
	Variable lowPoint = FindValPointOrLess(xwave,lowxValue)
	Variable highPoint = FindValPointOrMore(xwave,highxValue)
	
	Return max(lowPoint,highPoint)
end

Function FindFitStartPoint(xWave,lowxValue,highxValue)

	Wave xWave
	Variable lowxValue,highxValue
	
	Variable temp
	if (highxValue < lowxValue)
		temp = lowxValue
		lowxValue = highxValue
		highxValue = temp
	endif
	
	Variable lowPoint = FindValPointOrLess(xwave,lowxValue)
	Variable highPoint = FindValPointOrMore(xwave,highxValue)
	
	Return min(lowPoint,highPoint)
end

Function findMedian(opWave)

	Wave opWave;
	
	Duplicate/O opWave,tempopWave;
	Sort tempopWave,tempopWave;
	
	Variable midpoint,median;
	if ( mod(numpnts(tempopWave),2) == 1)
		midpoint = (numpnts(tempopWave) + 1)/2;
		median = tempopWave[midpoint];
	else
		midpoint = numpnts(tempopWave) / 2;
		median = (tempopWave[midpoint]+tempopWave[midpoint+1])/2;
	endif
	
	KillWaves tempopWave;
	Return median;
end

Function FindNaNIndex(targetWave)

	Wave targetWave
	
	Variable i
	for(i=0;i<numpnts(targetWave);i+=1)
		if (numtype(targetWave[i])==2)
			Return i
		endif
	endfor
	
	Return nan
end

Function FindNearestValue(targetWave,value)
	//Returns the *index* of the entry of targetWave with value closest to value
	Wave targetWave
	Variable value
	
	Variable i
	Variable bestIndex=0,bestProximity=abs(targetWave[0]-value)
	for(i=0;i<numpnts(targetWave);i+=1)
		if (abs(targetWave[i]-value)<bestProximity)
			bestProximity=abs(targetWave[i]-value)
			bestIndex=i
		endif
	endfor

	Return bestIndex
end

Function FindPeakLocation(srchWave)

	Wave srchWave;
	
	Variable wLength = numpnts(srchWave);
	Variable maxval = srchWave[0],peakindex=0;
	Variable i;
	for(i=0;i<wLength;i+=1)
		if (srchWave[i] > maxval)
			maxval = srchWave[i];
			peakindex = i;
		endif
	endfor
	
	Return peakindex;
end

Function FindScanDirection(multiScanWave)

	Wave multiScanWave;
	
	Variable i,direction;
	Variable totPoints = numpnts(multiScanWave);
	for(i=1;i<totPoints;i+=1)
		if ((multiscanWave[i]-multiscanWave[i-1])!=0)
			direction = sign(multiScanWave[i]-multiScanWave[i-1]);
		endif
	endfor
	
	Return direction;
end

Function FindScanLength(multiScanWave)

	Wave multiScanWave;
	
	Variable direction = FindScanDirection(multiScanWave);
	
	Variable i;
	Variable totPoints = numpnts(multiScanWave);
	for(i=1;i<totPoints;i+=1)
		if ((multiScanWave[i]-multiScanWave[i-1])*direction < 0)
			break;
		endif
	endfor
	
	Return i;
end

Function findSecLength(xwave,firstPnt)

	Wave xwave;
	Variable firstPnt;
	
	Variable pntCnt=1,i;
	Variable maxpnt = numpnts(xwave) - 1;
	
	if (firstPnt > maxpnt)
		Return 0;
	elseif (firstPnt == maxpnt)
		Return 1;
	endif
	
	for (pntCnt=1;(firstPnt+pntCnt<maxpnt);pntCnt+=1)
		if (xwave[firstPnt+pntCnt-1]!=xwave[firstPnt+pntCnt])
			Return pntCnt;
		endif
	endfor
	
	Return pntCnt;
end

Function FindValPointOrLess(xwave,value)
	// Finds the point number of xwave where xwave == value
	// if that point does not exist it finds the nearest point with a value less than input value
	Wave xwave
	Variable value

	Variable direction
	if (xwave[1] > xwave[0])
		direction = 1
	else
		direction = -1
	endif
	
	Variable i,otherend
	if (direction == 1)
		i = numpnts(xwave) - 1
		otherend=0
	else
		i = 0
		otherend = numpnts(xwave) - 1
	endif
	
	for(;xwave[i]>value;i-=direction)
		if (i==otherend)
			break
		endif
		if (abs(xwave[i] - value) < 1e-10)
			break
		endif
	endfor
	
	Return i
end
	
Function FindValPointOrMore(xwave,value)
	// Finds the point number of xwave where xwave == value.
	// If that point does not exist, it finds the nearest point with a value more than input value
	Wave xwave
	Variable value

	Variable direction
	if (xwave[1] > xwave[0])
		direction = 1
	else
		direction = -1
	endif
	
	Variable i,otherend
	if (direction == -1)
		i = numpnts(xwave) - 1
		otherend = 0
	else
		i = 0
		otherend = numpnts(xwave) - 1
	endif
	
	for(;xwave[i]<value;i+=direction)
		if (i == otherend)
			break
		endif

		if (abs(xwave[i] - value) < 1e-10)
			break
		endif
	endfor
	
	Return i
end

Function FindValueIndex(targetWave,value)

	Wave targetWave
	Variable value
	
	Variable i
	for(i=0;i<numpnts(targetWave);i+=1)
		if (targetWave[i]==value)
			Return i
		endif
	endfor
	
	Return nan
end

Function FindWaveDirection(targetWave,targetPoint)

	Wave targetWave
	Variable targetPoint
	targetPoint=floor(targetPoint)
	
	If(numpnts(targetWave)<2)
		Return 0
	endif
	
	Variable firstPoint,secondPoint
	if(targetPoint==numpnts(targetWave)-1)
		targetPoint-=1
	endif
	
	if(targetWave[targetPoint+1]>targetWave[targetPoint])
		Return 1
	elseif(targetWave[targetPoint+1]<targetWave[targetPoint])
		Return -1
	else
		Return 0
	endif

end

Function findWaveInsertionPt(targetWave,value)
	// Find the index in targetWave where value would go if it were to be inserted
	// (assumes targetWave is sorted in ascending order).
	Wave targetWave
	Variable value
	
	Variable i
	if (numpnts(targetWave)==0)
		Return 0
	elseif (value>targetWave[numpnts(targetWave)-1])
		Return numpnts(targetWave)
	else
		for(i=0;value>targetWave[i];i+=1)
		endfor
		Return i
	endif
end

Function firstNumber(targetWave)
	Wave targetWave

	Variable i,targetVal

	for(i=0;i<numpnts(targetWave);i+=1)
		targetVal=targetWave[i]
		if (targetVal!=nan)
			return targetVal
		endif
	endfor
	Return nan
end

//Returns the value of the cumulative distribution function for the Gaussian 
//distribution with average avg and standard deviation stdDev at point xVal.
//
//avg and stdDev are stored in the 0th and 1st entries of paramWave 
//respectively.
Function GaussianCDF(xVal,paramWave)
	
	Variable xVal
	Wave paramWave
	
	Variable avg=paramWave[0],stdDev=paramWave[1]
	
	Return .5*(1+erf((xVal-avg)/(sqrt(2)*stdDev)))
end

//This function creates a 2D histogram from the coordinate pairs formed by 
//the data in xWave and yWave.  The row dimension of the output histogram 
//wave corresponds to the xWave axis and the column dimension to the 
//yWave axis.  The data in the columns is reversed because Igor likes to 
//plot 2D data with the column number increasing from top to bottom of the 
//graph (whereas usually graph axes are plotted with increasing values from 
//bottom to top).
//
//Returns a string with the name of the output wave.  The output wave scaling
//is changed to match the center point of each bin.
//
//binStart -- the lowest value of the first bin (bins are defined in ascending order)
//binSize -- the width of each bin
//binTotal -- the total number of bins to use
Function/S Histogram2D(xWave,yWave,binStartX,binStartY,binSizeX,binSizeY,binTotalX,binTotalY)
	
	Wave xWave,yWave
	Variable binStartX,binStartY,binSizeX,binSizeY,binTotalX,binTotalY
	
	//Create the output histogram
	String Hist2DWaveName=NameOfWave(xWave)+"_H2D"
	Make/O/D/N=(binTotalX,binTotalY) $Hist2DWaveName
	Wave Hist2DWave=$Hist2DWaveName
	SetScale/P x, binStartX+binSizeX/2,binSizeX,Hist2DWave
	//As mentioned in the notes above.  The y axis is reversed for easier plotting in Igor.
	SetScale/P y, binStartY+binSizeY*(binTotalY-.5),-binSizeY,Hist2DWave
	
	Make/O/D/N=(binTotalY) tempYHist
	Variable i,j
	for(i=0;i<binTotalX;i+=1)
		Make/O/D/N=0 tempYSubset
		//Loop through the xWave looking for values that fall into the current bin.
		//Record the yWave values associated with these xWave points into the 
		//tempYSubset wave.
		for(j=0;j<numpnts(xWave);j+=1)
			if((xWave[j]>=binStartX+i*binSizeX)&&(xWave[j]<=binStartX+(i+1)*binSizeX))
				AppendPoint2Wave(tempYSubset,yWave[j])
			endif
		endfor
		//Make a histogram of the tempYSubset wave found in the loop above 
		//and save it as a row in the output wave (for the appropriate x bin).
		if (numpnts(tempYSubset)>0)
			Histogram/B={binStartY,binSizeY,binTotalY} tempYSubset, tempYHist
		else
			Make/O/D/N=(binTotalY) tempYSubset
			tempYHist=0
		endif
		Reverse tempYHist //Need to reverse the columns so that Igor plots the images intuitively
		Hist2DWave[i][]=tempYHist[q]
		KillWaves tempYSubset
	endfor
	KillWaves tempYHist
	
	Return Hist2DWaveName
end

//Gives the argument of the first bessel function for the amplitude corresponding to bessVal
//Valid only below the first zero of the Bessel function.  (Coefficients were found by fitting 
//the first Bessel function to a polynomial over this range).
Function InvertBessel1(bessVal)

	Variable bessVal
	
	Make/D/N=20 Bessel1CoefficientWave
	Bessel1CoefficientWave={-0.002117290037,3.251295895,-125.5703328,5063.940571,-106426.9817,1314316.756,-10080126.84,48440683.95,-137859556.2,175968734.7,127903202,-584110637.6,-139941281.5,1716780639,485578382.2,-5332602515,-1207112382,1.93751381e+10,-2.364567245e+10,9260978305}
	
	Variable bessArg= poly(Bessel1CoefficientWave,bessVal)
	KillWaves Bessel1CoefficientWave
	Return bessArg
end

Function jinc(x)

	Variable x
	
	Return 2*BesselJ(1,x)/x
	
end

//Calculatives the Kolmogorov-Smirnov test statistic (often denoted as M in the literature) 
//for the data contained in dataWave assuming a distribution identified by distributionID 
//and paramWave.
//distributionID: 0=Gaussian,1=Rayleigh
//paramWave: usually {mean,standard deviation} see documentation for the distribution
//function you want to use for more.
//KS test follows description in Handbook of parametric and nonparametric statistical 
//procedures by David Sheskin
Function KSTestStatistic(dataWave,distributionID,paramWave)

	Wave dataWave,paramWave
	Variable distributionID
	
	Switch(distributionID)
		default:
		Case 0:
			FUNCREF GaussianCDF cdfunc=GaussianCDF
			break
		Case 1:
			FUNCREF GaussianCDF cdfunc=RayleighCDF
			break
	endswitch
	
	//Make new wave with the data sorted in ascending order
	Duplicate/O dataWave,KStempData
	Sort dataWave,KStempData
	
	Variable KStstat=0 //holds the current test statistic value
	Variable i,j,Sxi,FXi
	for(i=0;i<2;i+=1)
	//Loop through the data twice comparing the theoretical cumulative distribution function
	//to the data's cumulative proportion, comparing the theoretical cumulative distribution
	//function to the cumulative proportion of the ith and (i-1)th data points (reason for the 
	//two loops).
		for(j=0;j<numpnts(KStempData);j+=1)
			SXi=(i+j)/numpnts(KStempData) //S(X_(i+j-1)) in the usual notation (cumulative proportion)
			FXi=cdfunc(KStempData[j],paramWave) //theoretical cumulative proportion of data point i
			if (abs(SXi-FXi)>KStstat)
				KStstat=abs(SXi-FXi)
			endif
		endfor
	endfor
	
	KillWaves KStempData
	Return KStstat
end

Function lastNumber(targetWave)
	Wave targetWave

	Variable i,targetVal

	for(i=numpnts(targetWave)-1;i>=0;i-=1)
		targetVal=targetWave[i]
		if (targetVal!=nan)
			return targetVal
		endif
	endfor
	Return nan
end

//Creates a wave with the same number of entries as phaseWave but with all values
//integer multiples of 2 pi.  The procedure tracks successive entries of phaseWave to
//create an output that can be added to phaseWave in order to create a continous phase
//function by correcting for jumps across the branch cut at phase = +/- pi.
//
//Returns the name of the output wave.
Function/S MakePhaseJumps(phaseWave)

	Wave phaseWave
	
	//Make output wave
	String outWaveName=NameOfWave(phaseWave)+"_J"
	
	Duplicate/O phaseWave, $outWaveName
	Wave outWave=$outWaveName
	
	outWave=0
	Variable i,windingNumber=0
	for(i=1;i<numpnts(outWave);i+=1)
		//Check if consecutive points in the phase wave jump between a value greater than 3pi/4 
		//and less than -3pi/4.  If so, change the winding number appropriately.
		if ((abs(phaseWave[i-1])>3*pi/4)&&(abs(phaseWave[i])>3*pi/4)&&(sign(phaseWave[i-1])!=sign(phaseWave[i])))
			windingNumber=windingNumber-sign(phaseWave[i]-phaseWave[i-1])
		endif
		outWave[i]=windingNumber*2*pi
	endfor
	
	Return outWaveName
end

//This function takes in a wave containing the moments of a distribution and returns
//the CumulantIndex order cumulant.  The moment indices are assumed to be equal
//to the wave indices (the zeroth entry is ignored).
Function Moments2CumulantValue(MomWave,CumulantIndex)

	Wave MomWave
	Variable CumulantIndex
	
	Variable cumulant
	Switch(CumulantIndex)
		Case 1:
			Return MomWave[1]
		Case 2:
			Return MomWave[2]-MomWave[1]^2
		Case 3:
			Return MomWave[3]-3*MomWave[1]*MomWave[2]+2*MomWave[1]^3
		Case 4:
			cumulant=MomWave[4]-4*MomWave[1]*MomWave[3]-3*MomWave[2]^2
			cumulant+=12*MomWave[1]^2*MomWave[2]-6*MomWave[1]^4
			Return cumulant
		Case 5:
			cumulant=MomWave[5]-5*MomWave[1]*MomWave[4]-10*MomWave[2]*MomWave[3]
			cumulant+=20*MomWave[1]^2*MomWave[3]+30*MomWave[1]*MomWave[2]^2
			cumulant+=-60*MomWave[1]^3*MomWave[2]+24*MomWave[1]^5
			Return cumulant
		Case 6:
			cumulant=MomWave[6]-6*MomWave[1]*MomWave[5]+30*MomWave[1]^2*MomWave[4]
			cumulant+=-10*MomWave[3]^2+120*MomWave[1]*MomWave[2]*MomWave[3]
			cumulant+=-120*MomWave[1]^3*MomWave[3]+30*MomWave[2]^3
			cumulant+=-270*MomWave[1]^2*MomWave[2]^2+360*MomWave[1]^4*MomWave[2]
			cumulant+=-120*MomWave[1]^6
			Return cumulant
		Case 7:
			cumulant=MomWave[7]-7*MomWave[1]*MomWave[6]-21*MomWave[2]*MomWave[5]
			cumulant+=42*MomWave[1]^2*MomWave[5]-35*MomWave[3]*MomWave[4]
			cumulant+=210*MomWave[1]*MomWave[2]*MomWave[4]-210*MomWave[1]^3*MomWave[4]
			cumulant+=140*MomWave[1]*MomWave[3]^2+210*MomWave[2]^2*MomWave[3]
			cumulant+=-1260*MomWave[1]^2*MomWave[2]*MomWave[3]+840*MomWave[1]^4*MomWave[3]
			cumulant+=-630*MomWave[1]*MomWave[2]^3+2520*MomWave[1]^3*MomWave[2]^2
			cumulant+=-2520*MomWave[1]^5*MomWave[2]+720*MomWave[1]^7
			Return cumulant
		default:
			Return NaN
		endswitch
end

//This function takes in a wave containing each of the moments of a wave
//and returns a wave with the cumulants.  The order of the moment should be
//the index of the wave.  The 0th entry is ignored.  The cumulant wave 
//is named CumulantName.  If CumulantName is the empty string, the 
//cumulant wave is the moment wave name plus "_C"
Function Moments2CumulantWave(MomWave,CumulantName)

	Wave MomWave
	String CumulantName
	
	if (strlen(CumulantName)==0)
		CumulantName=NameOfWave(MomWave)+"_C"
	endif
	
	Duplicate/O MomWave,$CumulantName
	Wave CumulantWave=$CumulantName
	
	Variable i
	for(i=1;i<numpnts(MomWave);i+=1)
		CumulantWave[i]=Moments2CumulantValue(MomWave,i)
	endfor
	
	
end

//Creates a wave (named by signWaveName) with the same number of 
//points as phaseWave that has all entries of 1 or -1
//The 1 entries occur for the entries of phaseWave with values between 
//Pos2NegPhase - 180 and Pos2NegPhase (modulo 360, inclusive) and 
//the -1 values occur for the other half of the circle.
//Assumes all entries in phaseWave are >-180 and <=180.
Function Phase2Sign(phaseWave,Pos2NegPhase,signWaveName)

	Wave phaseWave
	Variable Pos2NegPhase
	String signWaveName
	
	Duplicate/O phaseWave,$signWaveName
	Wave signWave=$signWaveName
	
	//Determine whether the positive or negative values are the ones 
	//with the contiguous interval between -180 and +180.
	Variable lowPosPhase,highPosPhase,signVariable
	if (Pos2NegPhase>0)
		highPosPhase=Pos2NegPhase
		lowPosPhase=Pos2NegPhase-180
		signVariable=1
	else
		lowPosPhase=Pos2NegPhase
		highPosPhase=Pos2NegPhase+180
		signVariable=-1
	endif
	
	//Loop through the phase wave and save 1 for values within the contiguous
	//phase region and -1 outside of it.  Save signVariable for the boundary points.
	Variable i
	for(i=0;i<numpnts(phaseWave);i+=1)
		if ((phaseWave[i]>lowPosPhase) && (phaseWave[i]<highPosPhase))
			signWave[i]=1
		elseif ((phaseWave[i]==lowPosPhase) || (phaseWave[i]==highPosPhase))
			signWave[i]=signVariable
		else
			signWave[i]=-1
		endif
	endfor
	//Negate the output wave if the contiguous phase region corresponds to the 
	//negative values.
	signWave*=signVariable
end

//Returns the projection of the complex number cmplxNum onto the unit vector
//in the complex plane with angle PhasorAngle relative to the x axis.
Function ProjectCmplx(cmplxNum,PhasorAngle)

	Variable/C cmplxNum
	Variable PhasorAngle
	
	Variable projection=cos(PhasorAngle)*real(cmplxNum)+sin(PhasorAngle)*imag(cmplxNum)
	
	Return projection
end

//The value of the function f(x)=y0+A*x^pow for the value x and 
//coefWave={y0,A,pow}.
Function Power(coefWave,x)

	Wave coefWave
	Variable x
	
	Variable y0=coefWave[0]
	Variable A=coefWave[1]
	Variable pow=coefWave[2]
	
	Return y0+A*x^pow
end

//Returns the value of the cumulative distribution function for the Rayleigh
//distribution with standard deviation stdDev at point xVal
//
//stdDev is stored in the 1st entry of paramWave (0th entry is not used).
Function RayleighCDF(xVal,paramWave)

	Variable xVal
	Wave paramWave
	
	Variable stdDev=paramWave[1]
	
	Return 1-exp(-xVal^2/(2*stdDev^2))
end

//This function calls Histogram2D() to create a square 2D histogram from the sets of 
//coordinate pairs in xWave, yWave.  The histogram's bins run from -dMax to +dMax 
//where dMax is the maximum value of any point in either xWave or yWave.  The bin 
//widths are given by 2*dMax/numBinsLinear.
//
//The name of the 2D output wave is returned.  The wave is normalized so that 
//integration over the histogram gives 1.
Function/S SquareHist2D(xWave,yWave,numBinsLinear)
	
	Wave xWave,yWave
	Variable numBinsLinear
	
	//Find the largest data point magnitude to set the range for the histogram to 
	//cover all data.
	Wavestats/Q xWave
	Variable dataRange=max(abs(V_max),abs(V_min))
	Wavestats/Q yWave
	dataRange=max(dataRange,max(abs(V_max),abs(V_min)))
	Variable binSize=2*dataRange/numBinsLinear
	String HistName=Histogram2D(xWave,yWave,-dataRange,-dataRange,binSize,binSize,numBinsLinear,numBinsLinear)
	
	//Normalize the data by the weight of the output histogram and the area covered
	//by the histogram.
	Wave HistWave=$HistName
	Variable normalization=sum(HistWave)*4*dataRange^2
	HistWave=HistWave/normalization
	
	Return HistName
end

//Computes the weighted average of dataWave with weights given by
//weightWave.  The entries of weightWave should be proportional to 
//the standard deviations of the entries in dataWave
Function weightedAverage(dataWave,weightWave)

	Wave dataWave,weightWave
	
	Variable i,wAvg,normWeight
	for(i=0,wAvg=0,normWeight=0;i<numpnts(dataWave);i+=1)
		wAvg+=dataWave[i]/weightWave[i]^2
		normWeight+=1/weightWave[i]^2
	endfor
	wAvg/=normWeight
	
	Return wAvg
end

//Computes the standard deviation of the weighted average calculated by 
//weightedAverage().  The entries of weightWave must be the standard 
//deviations of the values being averaged (can't just be proportional as in 
//weightedAverage()).
Function weightedAverageStdDev(weightWave)

	Wave weightWave
	
	Variable i,normWeight
	for(i=0,normWeight=0;i<numpnts(weightWave);i+=1)
		normWeight+=1/weightWave[i]^2
	endfor
	
	Return sqrt(1/normWeight)
end