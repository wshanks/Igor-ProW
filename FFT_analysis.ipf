#pragma rtGlobals=1		// Use modern global access method.

//You can change the graph macro used by FFTindexD() here
Function FFTGraphMacro()
	Return 10
end

Function AvgFFT(sweepList,Bwave,FFTNamesWave)

	Wave/T sweepList, FFTNamesWave;
	Wave Bwave;
	
	Variable totalsweeps = numpnts(sweepList);
	Variable i,j;
	String currentSweep;
	
	for(i=0;i<totalsweeps;i+=1)
		currentSweep = sweepList[i];
		FreqFFTALL($currentSweep,BWave);
	endfor
	
	String sweepListName=NameofWave(sweepList);
	String windowName = sweepListName + "_AvgFFT";
	DoWindow/K $windowName;
	Display/N=$windowName;
	DoWindow/T $windowName, "Avg FFT's of " + sweepListName;
	
	String FFTWaveName = sweepList[0] + "_FFT_" + FFTNamesWave[0];
	Variable FFTlength = numpnts($FFTWaveName);
	Variable totalwindows = numpnts(FFTNamesWave);
	String tempFFTName;
	Wave tempFFTWave;
	Variable hue,r,g,b;
	Make/O/N=(3) rgbWave;
	for(i=0;i<totalwindows;i+=1)

		Make/O/N=(FFTlength) tempFFTAvg,tempFFTWave;
		for(j=0;j<totalsweeps;j+=1)
			tempFFTName = sweepList[j] + "_FFT_" + FFTNamesWave[i];
			Duplicate/O $tempFFTName tempFFTWave;
			tempFFTAvg += tempFFTWave;
		endfor
		tempFFTAvg /= totalsweeps;

		FFTWaveName = sweepListName + "_FFT_" + FFTNamesWave[i];		
		Duplicate/O tempFFTAvg $FFTWaveName;
		AppendToGraph/W=$windowName $FFTWaveName;
		hue = floor((i / totalwindows) * 65536);
		hslconvert(hue,61680,34672);
		r = rgbWave[0];g = rgbWave[1];b = rgbWave[2];
		ModifyGraph rgb ($FFTWaveName) = (r,g,b);
		
		Variable deltaB=Abs(Bwave[1]-BWave[0]);
		Variable xscale=.5/deltaB;
		SetScale /I x,0,xscale,"",$FFTWaveName;
	endfor
	
	Label left "FFT Amp";Label bottom "1/T";	
	ModifyGraph mode=4,marker=19,msize=2,tick=2,mirror=1,fStyle=1,fSize=16,font="Arial";
	Legend/C/N=text0/F=0;
		
	KillWaves tempFFTAvg,tempFFTWave;
end

// Takes a series of FFTs of sections of the trace targetWave and saves all of the
// values of the resulting FFTs at point FFTpt in one wave to a new wave.  The 
// sections are all of a fixed size and each successive section is shifted one 
// point forward in the targetWave.
//
// targetWave -- the wave from which successive sections are drawn for FFTs
// windowSize -- the size of the sections that are FFTed
// wIndex -- the index ID of the FFT window to be used
// FFTpt -- the point number of the FFT amplitude to be saved to the output wave
Function/S FFTAmpSeries(targetWave,windowSize,wIndex,FFTpt)

	Wave targetWave
	Variable windowSize,wIndex,FFTpt
	
	//Create the output fourier amplitude wave
	Variable outWaveLength=numpnts(targetWave)-windowSize+1
	String outWaveName=NameOfWave(targetWave)+"_FAS"
	Make/O/N=(outWaveLength)/C/D $outWaveName
	Wave/C outWave=$outWaveName
	//Start output wave scaling in the middle of the first Fourier window
	Variable outWaveStart=leftx(targetWave)+.5*(windowSize-1)*deltax(targetWave)
	Variable outWaveDelta=deltax(targetWave)
	SetScale/P x, outWaveStart, outWaveDelta, outWave
	//Append a note to the output wave with information about the wave
	String outWaveNote="DESC:FFT amp series wave;"
	outWaveNote=outWaveNote+"SOURCE:"+NameOfWave(targetWave)+";"
	outWaveNote=outWaveNote+"WINSIZE:"+num2istr(windowSize)+";"
	outWaveNote=outWaveNote+"WINDOW:"+num2istr(wIndex)+";"
	outWaveNote=outWaveNote+"FFTPT:"+num2istr(FFTpt)+";"
	outWaveNote=outWaveNote+"FREQPT:"+num2str(FFTpt/windowSize/deltax(targetwave))+";"
	Note outWave, outWaveNote
	
	//Loop through targetWave taking fourier transforms
	Variable i,phase;String tempFFTname
	for(i=0;i<outWaveLength;i+=1)
		//Take the FFT of the section of targetWave
		tempFFTname=FFTIndexRP(targetWave,wIndex,i,i+windowSize-1)
		Wave/C tempFFT=$tempFFTname
		outWave[i]=tempFFT[FFTpt]
	endfor

	KillWaves tempFFT
	Return outWaveName
end

Function/S FFTIndex(YWave,windowIndex)

	Wave YWave;
	Variable windowIndex;
	Variable coherentgain;
	String windowName;
	
	windowName = FFTWindowName(windowIndex)
	coherentgain = FFTWindowCoherentGainFactor(windowIndex)
	String windowIndexStr = num2str(windowIndex);

	String FFTWaveName=FFTWindow(YWave,windowName,coherentgain,windowIndexStr)

	Return FFTWaveName
end

Function/S FFTIndexRP(YWave,windowIndex,startPt,endPt)

	Wave YWave;
	Variable windowIndex,startPt,endPt
	Variable coherentgain;
	String windowName;
	
	windowName = FFTWindowName(windowIndex)
	coherentgain = FFTWindowCoherentGainFactor(windowIndex)
	String windowIndexStr = num2str(windowIndex);

	String FFTWaveName=FFTWindowRP(YWave,windowName,coherentgain,windowIndexStr,startPt,endPt)

	Return FFTWaveName
end

//Take FFT of targetWave.  Create a second real valued wave with the 
//projections of the complex FFT values onto the unit vector in the 
//complex plane with angle PhasAng to the x axis.  If plotOpt=1, 
//creates a new window named "tempProjWind" and displays the 
//projected FFT wave.  If plotOpt=2, appends the projected wave to a 
//window with this name.
//Returns a list with the name of the FFT wave first and the projected 
//FFT wave second.
Function/S FFTProj(targetWave,PhasAng,plotOpt)

	Wave targetWave
	Variable PhasAng,plotOpt
	
	String FFTWaveName=FFTIndex(targetWave,10)
	Wave/C FFTWave=$FFTWaveName
	
	String ProjFFTWaveName=FFTWaveName+"p"
	Make/D/O/N=(numpnts(FFTWave)) $ProjFFTWaveName
	Wave ProjFFTWave=$ProjFFTWaveName
	
	ProjFFTWave=ProjectCmplx(FFTWave,PhasAng)
	OptionalPlot(ProjFFTWaveName,"tempProjWind",plotOpt)
	
	String FFTNames=FFTWaveName+";"+ProjFFTWaveName
	Return FFTNames
end

Function/S FFTWindow(YWave,windowName,gainfactor,windowAbbrev)

	Wave YWave;
	String windowName,windowAbbrev;
	Variable gainfactor
	String YWaveName, FFTWaveName;
	
	YWaveName=NameofWave(YWave);
	FFTWaveName=YWaveName+"_F"+windowAbbrev;
	Variable paddedLength = 2 * ceil(numpnts(YWave) / 2);
	FFT/PAD=(paddedLength)/WINF=$windowName/DEST=$FFTWaveName YWave;
	
	Wave/C FFToutwave = $FFTWaveName;
	FFToutwave = FFToutwave * 2 / paddedLength / gainfactor
	FFToutwave[0] = FFToutwave[0]/2
	FFToutwave[round(paddedLength / 2)] = FFToutwave[(paddedLength / 2)]/2
	Variable deltaB=DimDelta(YWave,0)
	Variable xscale=(.5/deltaB)//*(paddedLength/(paddedLength-1)  //Can't remember what last factor is for
	SetScale /I x,0,xscale,"",$FFTWaveName
	//FFToutwave=cmplx(cos(2*pi*leftx(YWave)*x),-sin(2*pi*leftx(YWave)*x))*FFToutwave
	
	Return FFTWaveName
end

Function/S FFTWindowRP(YWave,windowName,gainfactor,windowAbbrev,startPt,endPt)

	Wave YWave;
	String windowName,windowAbbrev;
	Variable gainfactor,startPt,endPt
	String YWaveName, FFTWaveName;
	
	YWaveName=NameofWave(YWave);
	FFTWaveName=YWaveName+"_F"+windowAbbrev;
	if (mod(endPt-startPt+1,2)==1)
		endPt=endPt-1
	endif
	Variable paddedLength = endPt-startPt+1
	FFT/RP=[startPt,endPt]/WINF=$windowName/DEST=$FFTWaveName YWave;
	
	Wave/C FFToutwave = $FFTWaveName;
	FFToutwave = FFToutwave * 2 / paddedLength / gainfactor
	FFToutwave[0] = FFToutwave[0]/2
	FFToutwave[(paddedLength / 2 +1)] = FFToutwave[(paddedLength / 2 +1)]/2
	Variable deltaB=DimDelta(YWave,0)
	Variable xscale=(.5/deltaB)*(paddedLength/(paddedLength-1)
	SetScale /I x,0,xscale,"",$FFTWaveName
	FFToutwave=cmplx(cos(2*pi*leftx(YWave)*x),sin(2*pi*leftx(YWave)*x))*FFToutwave
	
	Return FFTWaveName
end


Function/S MFFTIndex(YWave,windowIndex)

	Wave YWave;
	Variable windowIndex;
	Variable coherentgain;
	String windowName;
	
	windowName = FFTWindowName(windowIndex)
	coherentgain = FFTWindowCoherentGainFactor(windowIndex)
	String windowIndexStr = num2str(windowIndex);

	String FFTWaveName=MFFTWindow(YWave,windowName,coherentgain,windowIndexStr)

	Return FFTWaveName
end

Function/S MFFTIndexD(YWave,windowIndex)

	Wave YWave;
	Variable windowIndex;
	Variable coherentgain;
	String windowName;
	
	String YWaveName=NameofWave(YWave);
	String IgorwindowName = "FFT_" + YWaveName;
	DoWindow/K $IgorwindowName;
	Display/N=$IgorwindowName;
	DoWindow/T $IgorwindowName, "FFT of " + YWaveName;
	
	windowName = FFTWindowName(windowIndex)
	coherentgain = FFTWindowCoherentGainFactor(windowIndex)
	String windowIndexStr = num2str(windowIndex);

	String FFTWaveName=MFFTWindow(YWave,windowName,coherentgain,windowIndexStr)

	String winIndexStr;
	sprintf winIndexStr, "%d",windowIndex;
	FFTWaveName = NameofWave(YWave) + "_MF" + winIndexStr;
	AppendToGraph/W=$IgorwindowName $FFTWaveName;
	
	DoGMacro(10)

	Return FFTWaveName
end

Function/S MFFTWindow(YWave,windowName,gainfactor,windowAbbrev)

	Wave YWave;
	String windowName,windowAbbrev;
	Variable gainfactor
	String YWaveName, FFTWaveName;
	
	YWaveName=NameofWave(YWave);
	FFTWaveName=YWaveName+"_MF"+windowAbbrev;
	Variable paddedLength = 2 * ceil(numpnts(YWave) / 2);
	FFT/MAG/PAD=(paddedLength)/WINF=$windowName/DEST=$FFTWaveName YWave;
	
	Wave FFToutwave = $FFTWaveName;
	FFToutwave = FFToutwave * 2 / paddedLength / gainfactor;
	FFToutwave[0] /= 2;
	FFToutwave[(paddedLength / 2 +1)] /= 2;
	Variable deltaB=DimDelta(YWave,0)
	Variable xscale=(.5/deltaB)*(paddedLength/(paddedLength-1)
	SetScale /I x,0,xscale,"",$FFTWaveName;
	
	Return FFTWaveName
end

Function FreqFFTALL(FreqWave,BWave)

	Wave FreqWave,BWave;
	Wave/T FFTNamesWave = root:aux:FFTwindow;
	
	String FreqWaveName;
	
	FreqWaveName=NameofWave(FreqWave);
	String windowName = "FFT_" + FreqWaveName;
	DoWindow/K $windowName;
	Display/N=$windowName;
	DoWindow/T $windowName, "FFT's of " + FreqWaveName;
	
	Variable windowNumber = DimSize(FFTNamesWave,0);
	String FFTWindowName;
	String FFTWaveName, coherentgainStr;
	Variable i,r,g,b,hue,coherentgain;
	for(i=0;i<windowNumber;i+=1)
		FFTWindowName = FFTNamesWave[i][0];
		coherentgainStr = FFTNamesWave[i][1];
		coherentgain = str2num(coherentgainStr);
		FreqFFTIndex(FreqWave,BWave,i);
		FFTWaveName=FreqWaveName+"_FFT_"+FFTWindowName;
		AppendToGraph/W=$windowName $FFTWaveName;
		hue = floor((i / windowNumber) * 65536);
		hslconvert(hue,61680,34672);
		Wave rgbWave;
		r = rgbWave[0];g = rgbWave[1];b = rgbWave[2];
		ModifyGraph rgb ($FFTWaveName) = (r,g,b);
	endfor
	
	Label left "FFT Amp";Label bottom "1/T";	
	ModifyGraph mode=4,marker=19,msize=2,tick=2,mirror=1,fStyle=1,fSize=16,font="Arial";
	Legend/C/N=text0/F=0;
	
	KillWaves rgbWave,imagewave,M_HSL2RGB;

end

Function/S FreqFFTIndex(FreqWave,BWave,windowIndex)

	Wave FreqWave, BWave;
	Variable windowIndex;
	Variable coherentgain;
	String windowName,coherentgainStr;
	
	Wave/T FFTNamesWave = root:aux:FFTwindow;
	windowName = FFTNamesWave[windowIndex][0];
	coherentgainStr = FFTNamesWave[windowIndex][1];
	coherentgain = str2num(coherentgainStr);
	String windowIndexStr = num2str(windowIndex);

	String FFTWaveName=FreqFFTWindow(FreqWave,BWave,windowName,coherentgain,windowIndexStr)

	Return FFTWaveName
end

Function/S FreqFFTIndexD(FreqWave,BWave,windowIndex)

	Wave FreqWave, BWave;
	Variable windowIndex;
	Variable coherentgain;
	String windowName;
	
	String FreqWaveName=NameofWave(FreqWave);
	String IgorwindowName = "FFT_" + FreqWaveName;
	DoWindow/K $IgorwindowName;
	Display/N=$IgorwindowName;
	DoWindow/T $IgorwindowName, "FFT of " + FreqWaveName;
	
	windowName = FFTWindowName(windowIndex)
	coherentgain = FFTWindowCoherentGainFactor(windowIndex)
	String windowIndexStr = num2str(windowIndex);

	String FFTWaveName=FreqFFTWindow(FreqWave,BWave,windowName,coherentgain,windowIndexStr)

	String winIndexStr;
	sprintf winIndexStr, "%d",windowIndex;
	FFTWaveName = NameofWave(FreqWave) + "_F" + winIndexStr;
	AppendToGraph/W=$IgorwindowName $FFTWaveName;
	Label left "FFT Amp (pA)";Label bottom "1/T";
	ModifyGraph prescaleExp (left) = 12;
	ModifyGraph mode=4,marker=19,msize=2,tick=2,mirror=1,fStyle=1,fSize=16,font="Arial";
	Legend/C/N=text0/F=0;

	Return FFTWaveName
end

Function/S FreqFFTWindow(FreqWave,BWave,windowName,gainfactor,windowAbbrev)

	Wave FreqWave, BWave;
	String windowName,windowAbbrev;
	Variable gainfactor
	String FreqWaveName, FFTWaveName;
	
	FreqWaveName=NameofWave(FreqWave);
	FFTWaveName=FreqWaveName+"_F"+windowAbbrev;
	Variable paddedLength = 2 * ceil(numpnts(FreqWave) / 2);
	FFT/MAG/PAD=(paddedLength)/WINF=$windowName/DEST=$FFTWaveName FreqWave;
	
	Wave FFToutwave = $FFTWaveName;
	FFToutwave = FFToutwave * 2 / paddedLength / gainfactor;
	FFToutwave[0] /= 2;
	FFToutwave[(paddedLength / 2 +1)] /= 2;
	Variable deltaB=Abs(Bwave[1]-BWave[0]);
	Variable xscale=(.5/deltaB)*(paddedLength/(paddedLength-1)
	SetScale /I x,0,xscale,"",$FFTWaveName;
	
	Return FFTWaveName
end

Function/S FFTWindowName(windowIndex)

	Variable windowIndex
	
	String FFTwindowNameList="Bartlett;Blackman367;Blackman361;Blackman492;Blackman474;"
	FFTwindowNameList=AddListItem("Cos1;Cos2;Cos3;Cos4;Hamming;Hanning",FFTwindowNameList,";",Inf)
	FFTwindowNameList=AddListItem("KaiserBessel20;KaiserBessel25;KaiserBessel30",FFTwindowNameList,";",Inf)
	FFTwindowNameList=AddListItem("Parzen;Poisson2;Poisson3;Poisson4;Riemann;Square",FFTwindowNameList,";",Inf)
	
	String FFTwindowName=StringFromList(windowIndex,FFTwindowNameList)
	
	Return FFTwindowName
end

Function FFTWindowCoherentGainFactor(windowIndex)

	Variable windowIndex
	
	Make/N=20 FFTTempGainFactorWave
	FFTTempGainFactorWave={.5,.42,.45,.36,.4,.64,.5,.42,.38,.54,.5,.49,.44,.4,.67,.44,.32,.25,.59,1}
	
	Variable gainFactor=FFTTempGainFactorWave[windowIndex]
	
	KillWaves FFTTempGainFactorWave
	
	Return gainFactor
end

Function/S PSDIndex(YWave,windowIndex)

	Wave YWave;
	Variable windowIndex;
	Variable coherentgain;
	String windowName,coherentgainStr;
	
	windowName = FFTWindowName(windowIndex)
	coherentgain = FFTWindowCoherentGainFactor(windowIndex)
	String windowIndexStr = num2str(windowIndex);

	String PSDWaveName=PSDWindow(YWave,windowName,coherentgain,windowIndexStr)

	Return PSDWaveName
end

Function/S PSDIndexD(YWave,windowIndex)

	Wave YWave;
	Variable windowIndex;
	Variable coherentgain;
	String windowName;
	
	String YWaveName=NameofWave(YWave);
	String IgorwindowName = "PSD_" + YWaveName;
	DoWindow/K $IgorwindowName;
	Display/N=$IgorwindowName;
	DoWindow/T $IgorwindowName, "PSD of " + YWaveName;
	
	windowName = FFTWindowName(windowIndex)
	coherentgain = FFTWindowCoherentGainFactor(windowIndex)
	String windowIndexStr = num2str(windowIndex);

	String PSDWaveName=PSDWindow(YWave,windowName,coherentgain,windowIndexStr)

	String winIndexStr;
	sprintf winIndexStr, "%d",windowIndex;
	PSDWaveName = NameofWave(YWave) + "_S" + winIndexStr;
	AppendToGraph/W=$IgorwindowName $PSDWaveName;
	
	Label left "FFT Amp (pA)";Label bottom "1/T";
	ModifyGraph prescaleExp (left) = 12;
	ModifyGraph mode=4,marker=19,msize=2,tick=2,mirror=1,fStyle=1,fSize=16,font="Arial";
	Legend/C/N=text0/F=0;

	Return PSDWaveName
end

Function/S PSDWindow(YWave,windowName,gainfactor,windowAbbrev)

	Wave YWave;
	String windowName,windowAbbrev;
	Variable gainfactor
	String YWaveName, PSDWaveName;
	
	YWaveName=NameofWave(YWave);
	PSDWaveName=YWaveName+"_S"+windowAbbrev;
	Variable paddedLength = 2 * ceil(numpnts(YWave) / 2);
	Variable unpaddedLength = numpnts(YWave)
	FFT/MAGS/PAD=(paddedLength)/WINF=$windowName/DEST=$PSDWaveName YWave;
	
	Wave PSDoutwave = $PSDWaveName;
	PSDoutwave = PSDoutwave * 2 / unpaddedLength / gainfactor;
	PSDoutwave[0] /= 2;
	PSDoutwave[round(paddedLength / 2)] /= 2;
	Variable deltaB=DimDelta(YWave,0)
	PSDoutwave*=deltaB
	Variable xscale=(.5/deltaB)//*(paddedLength/(paddedLength-1) //Not sure why I have this factor
	SetScale /I x,0,xscale,"",$PSDWaveName;
	
	Return PSDWaveName
end

//This script duplicates the wave with waveTag in its title to a subfolder title "CFFT" 
//and then takes the FFT of the wave and saves the real and imaginary parts of 
//the FFT amplitude at the point wavePt in the FFT wave to the waves RealRecWave 
//and ImagRecWave respectively.
Function RecFFTVals(waveTag,wavePt,RealRecWave,ImagRecWave)

	String waveTag
	Variable wavePt
	Wave RealRecWave,ImagRecWave
	
	String targetWaveName=FindWaveWithTag(waveTag)
	Wave targetWave=$targetWaveName
	String origFolder = GetDataFolder(1)
	String subFolderName="CFFT"
	Dup1Wave2Folder(targetWave,subFolderName)
	SetDataFolder $subFolderName
	
	String FFTWaveName=FFTIndex(targetWave,10)
	Wave FFTWave=$FFTWaveName
	RecRealWavePt(FFTWave,wavePt,RealRecWave)
	RecImagWavePt(FFTWave,wavePt,ImagRecWave)
	
	SetDataFolder $origFolder
end