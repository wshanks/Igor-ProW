#pragma rtGlobals=1		// Use modern global access method.

Function fitPoly(yWave,xWave,pOrd,outCoefName)
	//Fit yWave vs xWave to a polynomial of order pOrd.  Create wave
	//named outCoefName for the fit coefficients.
	Wave yWave,xWave
	Variable pOrd
	String outCoefname
	
	Make/D/N=(pOrd+1) $outCoefName
	Wave coefWave=$outCoefName
	
	Switch (pOrd)
		case 0:
			WaveStats/M=1 yWave;
			coefWave = V_avg;
			break;
		case 1:
			CurveFit/Q/M=2/W=0 line, yWave /X=xWave /D=coefWave	
			break;
		default:
			CurveFit/Q/M=2/W=0 poly (pOrd+1), yWave /X=xWave /D=coefWave
			break;
	endswitch
end

Function MakeBSWave(fullBwave,fullFwave,BGwaveName,BSwaveName)

	//Need to start in Background subtraction folder
	Wave fullBwave,fullFwave //Need to have "::" at beginning of wave names
	String BGwaveName,BSwaveName
	
	String origFolder = GetDataFolder(1)
	SetDataFolder ::
	Dup1Wave2Folder(fullBwave,origFolder)
	Dup1Wave2Folder(fullFwave,origFolder)
	SetDataFolder $origFolder
	
	MakePolyWave(fullBwave,BGwaveName)
	Wave BGwave = $BGwaveName
	
	Duplicate/O fullFwave, $BSwaveName
	Wave BSwave = $BSwaveName
	BSwave = fullFwave - BGwave
	
	SetDataFolder ::
	Duplicate/O BSwave $BSwaveName
	
end

Function MakePolyWave(xWave,OutWaveName)

	Wave xWave;
	String OutWaveName;
	
	Duplicate/O xWave $OutWaveName;
	Wave outWave = $OutWaveName;
	
	Wave W_coef;
	
	Variable i,j;
	
	for(i=0;i<numpnts(xWave);i+=1)
		outWave[i]=0;
		for(j=0;j<numpnts(W_coef);j+=1)
			outWave[i]+=W_coef[j]*xWave[i]^j;
		endfor
		if (xWave[i] == 0)
			outWave[i] = W_coef[0];
		endif
	endfor
	
end

Function MakePolyWave2(xWave,coefWave,OutWaveName)

	Wave xWave,coefWave;
	String OutWaveName;
	
	Duplicate/O xWave $OutWaveName;
	Wave outWave = $OutWaveName;
	
	Variable i,j;
	
	for(i=0;i<numpnts(xWave);i+=1)
		outWave[i]=0;
		for(j=0;j<numpnts(coefWave);j+=1)
			outWave[i]+=coefWave[j]*xWave[i]^j;
		endfor
		if (xWave[i] == 0)
			outWave[i] = coefWave[0];
		endif
	endfor
	
end

Function/S SubPoly(yWave,xWave,pOrder,[suppressPlots])

// This function takes in a frequncy wave and xWave, subtracts off a 2nd order polynomial and computes the variance
// FreqWave is the input wave and xWave is the corresponding xwave (time or B usually). 
// The subtracted slope is a polynomial fit to FreqWave 

	Wave yWave, xWave;
	Variable pOrder,suppressPlots

	//Create a subfolder to hold all the waves generated for the smooth
	//background fitting and subtraction
	String origFolder = GetDataFolder(1)
	String ywName = NameofWave(yWave)
	String folderName
	sprintf folderName "%s_%s%d",ywName,"Poly",pOrder
	NewDataFolder /O $folderName
	SetDataFolder $folderName
	
	//Copy yWave into subfolder
	Duplicate /O yWave $ywName

	//If there is a separate xWave, copy it into subfolder
	if (waveMatch(yWave,xWave)!=1)
		String xwName=NameOfWave(xWave)
		Duplicate /O xWave $xwName
	endif
	
	//Create a wave to store the fit
	String FitWaveName;
	sprintf FitWaveName "fit_%s_%s%d",ywName,"Poly",pOrder;
	Duplicate/O yWave $FitWaveName
	Wave PolyFitWave=$FitWaveName
	
	//Create a wave to store the output wave containing yWave minus the polynomial fit
	String ySlopeSubWaveName
	sprintf ySlopeSubWaveName "%s_%s%d",ywName,"P",pOrder
	Duplicate/O yWave $ySlopeSubWaveName
	Wave ySlopeSubWave = $ySlopeSubWaveName
	
	//Make the x-axis start at 0.  For some reason Igor gives better fits when this is the case
	Duplicate yWave yWave4Fit;Duplicate xWave xWave4Fit
	if (WaveMatch(yWave,xWave))
		SetScale/P x, 0, 1, yWave4Fit
	else
		Variable xWaveOffset=xWave[0]
		xWave4Fit-=xWaveOffset
	endif
	
	//Fit the yWave with the appropriate x wave and polynomial order
	Variable V_avg;
	Switch (pOrder)
		case 0:
			WaveStats/M=1 yWave4Fit
			PolyFitWave = V_avg;
			break
		case 1:
			if(WaveMatch(yWave,xWave))
				CurveFit/Q/M=2/W=0 line, yWave4Fit /D=PolyFitWave
			else
				CurveFit/Q/M=2/W=0 line, yWave4Fit /X=xWave4Fit /D=PolyFitWave
			endif
			break
		default:
			if(WaveMatch(yWave,xWave))
				CurveFit/Q/M=2/W=0 poly (pOrder+1), yWave4Fit /D=PolyFitWave
			else
				CurveFit/Q/M=2/W=0 poly (pOrder+1), yWave4Fit /X=xWave4Fit /D=PolyFitWave
			endif
			break
	endswitch
	KillWaves yWave4Fit,xWave4Fit
	
	ySlopeSubWave=yWave-PolyFitWave
	
	if(suppressPlots==0)
		if(WaveMatch(yWave,xWave))
			Display $ywName,PolyFitWave
		else
			Display $ywName,PolyFitWave vs $xwName
		endif
		DoGMacro(11)
	
		if(WaveMatch(yWave,xWave))
			Display ySlopeSubWave
		else
			Display ySlopeSubWave vs $xwName
		endif
		DoGMacro(0)
		MoveWindow 450,0,750,400;
	endif

	//Copy subtracted wave up a folder	
	SetDataFolder $origFolder;
	Duplicate /O ySlopeSubWave $ySlopeSubWaveName
	
	Return ySlopeSubWaveName
end

//For loop which runs subpoly() on each chunksize chunk of ywave vs xwave
//for polynomial order pord.  The function also creates a copy of ywave with
//each of the polynomial fits subtracted from the corresponding chunk.
//Returns the name of the wave with the fits subtracted.
Function/S SubPolyChunks(ywave,xwave,pord,chunksize)

	Wave ywave,xwave
	Variable chunksize,pord
	
	//If this does not equal zero, plots of individual subPoly calls will be suppressed.
	Variable noPlots=1
	
	String tempName=NameOfWave(ywave)+"_pch"+num2str(pord)
	Duplicate/O ywave, $tempName
	Wave subYwave=$tempName
	
	Variable i
	for(i=0;i<floor(numpnts(ywave)/chunksize);i+=1)
		tempName=NameOfWave(ywave)+num2str(i)
		Duplicate/O/R=(i*chunksize,(i+1)*chunksize-1) ywave,$tempName
		Wave tempY=$tempName
		
		tempName=NameOfWave(xwave)+num2str(i)
		Duplicate/O/R=(i*chunksize,(i+1)*chunksize-1) xwave,$tempName
		Wave tempX=$tempName
		
		tempName=Subpoly(tempY,tempX,pord,suppressPlots=noPlots)
		Wave tempSub=$tempName
		subYWave[i*chunkSize,(i+1)*chunkSize-1]=tempSub[p-i*chunkSize]
	endfor
	
	Return NameOfWave(subYwave)
end

//Fits a polynomial of order pOrder to yWave vs. xWave.  mWave is a mask wave for 
//yWave.  It is expected to have the same number of points as yWave and to have 
//0 and 1 for entries.  Only the points in yWave for which the corresponding point in 
//mWave !=0 are used in the fit.
//Calling this function creates a new folder containing duplicates of yWave and xWave 
//as well as a wave containing the polynomial fit with the same number of points and 
//scaling as yWave, another wave with difference between the yWave and the 
//polynomial fit, and the coefficient and sigma waves for the polynomial fit.
//
//The function also displays yWave vs xWave together with the polynomial fit in one 
//window and the difference wave in another.
//
//Returns the name of the difference wave.
Function/S SubPolyMask(yWave,xWave,mWave,pOrder)

	Wave yWave, xWave,mWave
	Variable pOrder

	//Create a subfolder to hold all the waves generated for the smooth
	//background fitting and subtraction
	String origFolder = GetDataFolder(1)
	String ywName = NameofWave(yWave)
	String folderName
	sprintf folderName "%s_%s%d",ywName,"Poly",pOrder
	NewDataFolder /O $folderName
	SetDataFolder $folderName
	
	//Copy yWave into subfolder
	Duplicate /O yWave $ywName

	//If there is a separate xWave, copy it into subfolder
	if (waveMatch(yWave,xWave)!=1)
		String xwName=NameOfWave(xWave)
		Duplicate /O xWave $xwName
	endif
	
	//Create a wave to store the fit
	String FitWaveName;
	sprintf FitWaveName "f_%s_%s%d",ywName,"P",pOrder;
	
	//Create a wave to store the output wave containing yWave minus the polynomial fit
	String ySlopeSubWaveName
	sprintf ySlopeSubWaveName "%s_%s%d",ywName,"P",pOrder
	Duplicate/O yWave $ySlopeSubWaveName
	Wave ySlopeSubWave = $ySlopeSubWaveName
	
	//Make the x-axis start at 0.  For some reason Igor gives better fits when this is the case
	Duplicate/O xWave xWave4Fit
	if (WaveMatch(yWave,xWave))
		//Making a separate wave that has values that match the scaling of yWave.
		//This just allows for the yWave to be fit versus the xWave4Fit in the code below
		//so that there doesn't need to be separate code for the fit vs. scaling case from 
		//the fit vs. xWave case.
		xWave4Fit=p*DimDelta(yWave,0)
	else
		Variable xWaveOffset=xWave[0]
		xWave4Fit-=xWaveOffset
	endif
	
	//Fit the yWave with the appropriate x wave and polynomial order
	Variable V_avg;
	Switch (pOrder)
		case 0:
			//Not really a fit in Igor.  Just subtract the average.
			String maskedyWave4FitName=ApplyFilterMask(yWave4Fit,mWave)
			Wave maskedyWave4Fit=$maskedyWave4FitName
			WaveStats/M=1 maskedyWave4Fit
			Make/O/N=1 w_coef
			w_coef=V_avg
			KillWaves maskedyWave4Fit
			break;
		case 1:
			CurveFit/Q/M=2/W=0 line, yWave /X=xWave4Fit /M=mWave
			break
		default:
			CurveFit/Q/M=2/W=0 poly (pOrder+1), yWave /X=xWave4Fit /M=mWave
			break
	endswitch
	//Make output fit wave.  Making it explicitly here rather than using Igor's /D option in CurveFit
	//ensures that the polynomial fit is calculated even for the masked out points.
	MakePolyWave(xWave4Fit,FitWaveName)	
	Wave PolyFitWave=$FitWaveName
	SetScale/P x, leftx(yWave),deltax(yWave)
	KillWaves xWave4Fit
	
	ySlopeSubWave=yWave-PolyFitWave

	//Display the input wave and the polynomial fit	
	if(WaveMatch(yWave,xWave))
		Display $ywName,PolyFitWave
	else
		Display $ywName,PolyFitWave vs $xwName
	endif
	DoGMacro(11)
	
	//Display the input wave after subtraction of the polynomial fit
	if(WaveMatch(yWave,xWave))
		Display ySlopeSubWave
	else
		Display ySlopeSubWave vs $xwName
	endif
	DoGMacro(0)
	MoveWindow 450,0,750,400;

	//Copy subtracted wave up a folder	
	SetDataFolder $origFolder;
	Duplicate /O ySlopeSubWave $ySlopeSubWaveName
	
	Return ySlopeSubWaveName
end

Function RemoveSlope()

	Wave yWave = CsrWaveRef(A);
	Wave xWave = CsrXWaveRef(A);
	
	Variable m = (yWave[pCsr(B)] - yWave[pCsr(A)]) / (xWave[pCsr(B)] - xWave[pCsr(A)]);
	Variable b = (yWave[pCsr(A)] * xWave[pCsr(B)] - yWave[pCsr(B)] * xWave[pCsr(A)])/ (xWave[pCsr(B)] - xWave[pCsr(A)]);

	String flatWaveName = NameofWave(yWave) + "_fl";
	Duplicate/O yWave $flatWaveName;
	Wave flatWave = $flatWaveName;
	flatWave = yWave - (m * xWave);

end

//Fits a polynomial (of order pOrd) to the wave with WaveTag in its name with the points outside of 
//TailA and TailB and the points between FeatA and FeatB masked out.  Subtracts
//the polynomial fit from the original wave (see SubPolyMask()) and records the value of the 
//wave with the polynomial subtracted at point RecPt to the wave RecWave.
//
//This function was used to subtract a polynomial background from an FFT and record 
//the amplitude of one point in the FFT.
Function RmFeatWaveBG(WaveTag,TailA,TailB,FeatA,FeatB,pOrd,RecPt,RecWave)
	
	Wave RecWave
	String WaveTag
	Variable TailA,TailB,FeatA,FeatB,pOrd,RecPt
	
	String targetWaveName = FindWaveWithTag(waveTag)
	Wave targetWave=$targetWaveName
	
	String MaskWaveName
	MaskWaveName = MaskTails(targetWave,TailA,TailB)
	MaskWavePoints(targetWave,FeatA,FeatB)
	Wave mWave = $MaskWaveName
	
	String SubWaveName = SubPolyMask(targetWave,targetWave,mWave,pOrd)
	Wave SubWave = $SubWaveName
	AppendPoint2Wave(RecWave,SubWave[RecPt])
end

//OLD FUNCTIONS
//
//

//At some point, I'd like to make this function more general by adding some sort of graph style macro
//selector instead of always using frequency and magnetic field to label the axes.
Function SubPolyOrd(FreqWave,xWave,pOrder)

// This function takes in a frequncy wave and xWave, subtracts off a 2nd order polynomial and computes the variance
// FreqWave is the input wave and xWave is the corresponding xwave (time or B usually). 
// The subtracted slope is a polynomial fit to FreqWave 

	Wave FreqWave, xWave;
	Variable pOrder;
	String FreqSlopeSubWaveName,xSlopeSubWaveName;
	Wave xSecWave,FreqSecWave,PolyFitWave, FreqSlopeSubWave;

	String FreqwName = NameofWave(FreqWave);
	String folderName;
	sprintf folderName "%s_%s%d",freqwName,"Poly",pOrder;
	NewDataFolder /O $folderName;
	
	String origFolder = GetDataFolder(1);
	
	String FreqSecWaveName = folderName+":FreqSecWave";
	String xSecWavename = folderName+":xSecWave";
	
	Duplicate /O FreqWave FreqSecWave; // The Frequency Wave in the range I want to fit
	Duplicate /O xWave xSecWave; // The Time wave in the range I want to fit
	SetDataFolder $folderName;
	Duplicate /O FreqSecWave $FreqwName;
	Duplicate FreqSecWave PolyFitWave;						// The Line Fit in the range I want to fit
	Variable V_avg;
	
	Switch (pOrder)
		case 0:
			WaveStats/M=1 FreqSecWave;
			PolyFitWave = V_avg;
			break;
		case 1:
			CurveFit/Q/M=2/W=0 line, FreqSecWave /X=xSecWave /D=PolyFitWave;		
			break;
		default:
			CurveFit/Q/M=2/W=0 poly (pOrder+1), FreqSecWave /X=xSecWave /D=PolyFitWave;
			break;
	endswitch
	
	Duplicate PolyFitWave FreqSlopeSubWave;
	FreqSlopeSubWave = FreqWave-PolyFitWave;	// The Frequency Wave with slope subtracted in the range I want to fit
	

	sprintf FreqSlopeSubWaveName "%s_%s%d",freqwName,"P",pOrder;
	String FitWaveName;
	sprintf FitWaveName "fit_%s_%s%d",freqwName,"Poly",pOrder;
	String xwName = NameofWave(xWave);
	sprintf xSlopeSubWaveName "%s_%s%d",xwName,"P",pOrder;
	
	Duplicate /O FreqSlopeSubWave $FreqSlopeSubWaveName;	// Copying the slope-subtracted wave to the input freq wave name
	Duplicate /O xSecWave $xSlopeSubWaveName;		//copying the section of the time wave to the input wave name
	Duplicate /O PolyFitWave $FitWaveName;
	
	Display $FreqwName,$FitWaveName vs xWave;
	ModifyGraph mode ($FreqwName)=3,mode ($FitWaveName)=0,marker ($FreqwName)=19,msize ($FreqwName)=2;
	ModifyGraph tick=2,mirror=1,fStyle=1,fSize=16,font="Arial";
	ModifyGraph rgb ($FreqwName)=(0,0,65535),rgb ($FitWaveName) = (65535,0,0);
	
	Display $FreqSlopeSubWaveName vs $xSlopeSubWaveName;
	Label bottom "V\Bdrive\M (mV)"; Label left "Frequency (µHz)";ModifyGraph prescaleExp(left)=6;
	ModifyGraph mode=4,mirror=1,marker($FreqSlopeSubWaveName)=19,msize($FreqSlopeSubWaveName)=1;
	ModifyGraph rgb($FreqSlopeSubWaveName)=(65535,0,0),tick=2;
	ModifyGraph prescaleExp(left)=6, fStyle=1,fSize=16,font="Arial";
	MoveWindow 500,0,900,400;

	//Copy subtracted wave up a folder	
	SetDataFolder $origFolder;
	Duplicate /O FreqSlopeSubWave $FreqSlopeSubWaveName

	SetDataFolder $folderName	
	KillWaves FreqSlopeSubWave, PolyFitWave,FreqSecWave,xSecWave;
	SetDataFolder $origFolder;
end