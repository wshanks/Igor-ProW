#pragma rtGlobals=1		// Use modern global access method.

//A collection of functions for averaging and filtering data

Function AddAvgPointEntry(targetWave,Points2Avg)
//Averages Points2Avg points of the wave selected by cursor A on the top graph starting from
//the point selected by pcsr(A).  Adds the averaged value on to the end of targetWave

	Wave targetWave
	Variable Points2Avg
	Wave yWave = CsrWaveRef(A)
	Wave xWave = CsrXWaveRef(A)
	
	Variable direction
	If (waveexists(xWave))
		direction=FindWaveDirection(xWave,pcsr(A))
	else
		direction=1
	endif
	
	Variable i,pointAvg=0
	for(i=0;i<Points2Avg;i+=1)
		pointAvg+=yWave[pcsr(A)+direction*i]
	endfor
	
	pointAvg/=Points2Avg

	AppendPoint2Wave(targetWave,pointAvg)
end

//Create a copy of datawave with every entry deleted for which the corresponding index
//in maskwave is zero.  Return the name of the copied wave.
//Datawave and maskwave should be the same length
Function/S ApplyFilterMask(datawave,maskwave)

	Wave datawave,maskwave;
	
	Variable maskedlength = sum(maskwave);

	String maskedWaveName = NameOfWave(datawave) + "_W";
	Make/D/O/N=(maskedlength) $maskedWaveName
	Wave maskedWave=$maskedWaveName
	
	Variable i;
	Variable j=0;
	for(i=0;i<numpnts(datawave);i+=1)
		if (maskwave[i] == 1)
			maskedwave[j] = datawave[i];
			j +=1;
		endif
	endfor
	
	Return maskedWaveName;
end

//Remove all points from fwave which are more than datawindow away from the corresponding
//entry in trendwave.
//Create a mask wave with the same length as fwave and entries of 1 for points
//not deleted from fwave and 0 for points deleted.  Return mask wave name.
Function/S ApplyWaveFilter(datawave,trendwave,datawindow)

	Wave datawave,trendwave;
	Variable datawindow;
	
	String maskName = "mask" + NameOfWave(datawave)
	Make/D/O/N=(numpnts(datawave)) $maskName
	Wave maskWave=$maskName

	Variable i
	for(i=0;i<numpnts(datawave);i+=1)
		if (abs(datawave[i]-trendwave[i])>datawindow)
			maskwave[i] = 0;
		else
			maskwave[i] = 1;
		endif
	endfor
	
	ApplyFilterMask(datawave,maskwave)

	Return maskName
end

Function AverageSweeps(SweptValueWaveName,w2avgList)
	// Overwrites original waves!

	String w2avgList,SweptValueWaveName;

	Wave SweptValueWave = $SweptValueWaveName;	
	Variable scLen = FindScanLength(SweptValueWave);
	Variable totScans = floor(numpnts(SweptValueWave)/scLen);
	Make/D/O/N=(scLen) tmpAvgdSwpWv;
	Wave tmpAvgdSwpWv;
	
	Variable i,j,k;
	String tmpWaveName;
	for(i=0;i<ItemsInList(w2avgList);i+=1)
		tmpWaveName = StringFromList(i,w2avgList);
		Wave tempWave2Avg = $tmpWaveName;
		tmpAvgdSwpWv = 0;
		for(j=0;j<scLen;j+=1)
			for(k=0;k<totScans;k+=1)
				tmpAvgdSwpWv[j] += tempWave2Avg[j+k*scLen];
			endfor
		endfor
		tmpAvgdSwpWv /= totScans;
		Duplicate/O tmpAvgdSwpWv, tempWave2Avg;
	endfor
	KillWaves tmpAvgdSwpWv;
end

Function/S AvgChunkbyChunk(dataWave,ctrlWave,avgMode)

	Wave dataWave,ctrlWave
	Variable avgMode // 0=normal average,1=tight average
	
	String avgdWaveName = nameofwave(dataWave) + "_avg"
	Make/D/O/N=0 $avgdWaveName
	Wave avgdWave = $avgdWaveName
	
	Variable lastCtrlPnt = numpnts(ctrlWave) - 1
	
	Variable i, chunkLength,newAvgVal
	for (i=0;i<lastCtrlPnt;i+=chunkLength)
		chunkLength = findSecLength(ctrlWave,i)
		Duplicate/O/R=[i,i+chunkLength-1] dataWave,newChunkWave
		if (avgMode == 0)
			newAvgVal = mean(newChunkWave)
		else //if (avgMode == 1)
			newAvgVal = FilteredAverage(newChunkWave)
		endif
		ExtendWavebyOne(avgdWave,newAvgVal)
	endfor
	
	KillWaves newChunkWave
	Return avgdWaveName
end

// Runs TightFreqAvgStdDev on all waves in the current folder that follow the CreateWaveName()
// naming convention for the list given in waveTagList
Function AvgListWaves(fileNum,folderName,waveTagList,nameMode,xWaveIndex)

	Variable fileNum,nameMode,xWaveIndex;
	String folderName,waveTagList;
	
	String xWaveTag = StringFromList(xWaveIndex,waveTagList);
	AvgListWavesTag(fileNum,folderName,waveTagList,nameMode,xWaveTag);
end

Function AvgListWavesbyCtrlWave(targetWaveList,ctrlWaveName)

	String targetWaveList,ctrlWaveName;
	
	Wave tmpCtrlWave = $ctrlWaveName;
	
	String tmpWaveName;
	Variable i;
	for(i=0;i<ItemsInList(targetWaveList);i+=1)
		tmpWaveName = StringFromList(i,targetWaveList);
		Wave tmpWave2Avg = $tmpWaveName;
		AvgWavebyCtrlWave(tmpWave2avg,tmpCtrlWave);
	endfor
end

Function AvgListWavesTag(fileNum,folderName,waveTagList,nameMode,xWaveTag)

	Variable fileNum,nameMode;
	String folderName,waveTagList,xWaveTag;
	
	Variable totWaves = ItemsInList(waveTagList);
	Variable i;
	String tempWaveName;
	
	tempWaveName = CreateWaveName(fileNum,folderName,xWaveTag,nameMode);
	Wave tempWave = $tempWaveName;
	Variable secLength = findSecLength(tempWave,0);
	
	for(i=0;i<totWaves;i+=1)
		xWaveTag = StringFromList(i,waveTagList);
		tempWaveName = CreateWaveName(fileNum,folderName,xWaveTag,nameMode);
		Wave tempWave = $tempWaveName;
		TightFreqAvgStdDev(tempWave,secLength);
	endfor
end

// Average all waves in a folder, using xWave to determine section length for averaging
// Move the averaged waves to a folder called folderName (created if it does not exist)
Function AvgSDev2Folder(xWave,folderName)

	Wave xWave;
	String folderName;
	
	NewDataFolder/O $folderName;
	
	Variable secLength = findSecLength(xWave,0);
	
	Variable totWaves = CountObjects("",1);
	
	Variable i;
	String tempName,newtempName;
	for(i=0;i<totWaves;i+=1)
		tempName = GetIndexedObjName("",1,i);
		Wave tempWave = $tempName;
		
		TightFreqAvgStdDev(tempWave,secLength)
		
		newTempName = tempName + "_tavg";
		Wave newTempWave = $newTempName;
		MvWave(newtempWave,FolderName);
		
		newTempName = tempName + "_tstddev";
		Wave newTempWave = $newTempName;
		MvWave(newtempWave,FolderName);
		//Dup1Wave2Folder(newTempWave,folderName);
		//KillWaves newTempWave;
	endfor

end

Function AvgWavebyCtrlWave(targetWave,ctrlWave)

	Wave targetWave,ctrlWave;
	Variable tolerance = 1e-15 // Problems with == ....
	
	Duplicate/O ctrlWave,tmpCtrlWave;
	SortedNoDupWave(tmpCtrlWave);
	String targetWaveName = nameofwave(targetWave);
	Variable outputPnts = numpnts(tmpCtrlWave);
	String outWaveName = nameofwave(targetWave) + "_avg";
	Make/O/D/N=(outputPnts) $outWaveName;
	Wave outWave = $outWaveName;
	
	Variable i,j,pntCnt,avgVal;
	for(i=0;i<numpnts(tmpCtrlWave);i+=1)
		for(avgVal=0,pntCnt=0,j=0;j<numpnts(targetWave);j+=1)
			if(abs(ctrlWave[j] - tmpCtrlWave[i])<tolerance)
				avgVal+=targetWave[j];
				pntCnt+=1;
			endif
		endfor
		outWave[i]=avgVal/pntCnt;
	endfor
	
	Duplicate/O outWave, $targetWaveName;
	KillWaves outWave,tmpCtrlWave;
end

Function AvgWaves(wavesList)

	Wave/T wavesList;
	
	Variable totalwaves = numpnts(wavesList);
	Variable i,j;
	String currentSweep;
	
	String wavesListName=NameofWave(wavesList);
	String windowName = wavesListName + "_Avgwin";
	DoWindow/K $windowName;
	Display/N=$windowName;
	DoWindow/T $windowName, "Avg of " + wavesListName;
	
	String AvgWaveName = wavesListName + "_Avg";
	String tempName;
	String tempWaveName = wavesList[0];
	Duplicate /O $tempWaveName $AvgWaveName;
	Wave AvgWave = $AvgWaveName;
	AvgWave = 0;
	
	Variable hue,r,g,b;
	Make/O/N=(3) rgbWave;
	for(i=0;i<totalwaves;i+=1)
		tempName = wavesList[i];
		Wave tempWave = $tempName;
		AvgWave += tempWave;

		AppendToGraph/W=$windowName $tempName;
		hue = floor((i / (totalwaves + 2)) * 65536);
		hslconvert(hue,61680,34672);
		r = rgbWave[0];g = rgbWave[1];b = rgbWave[2];
		ModifyGraph rgb ($tempName) = (r,g,b);
	endfor
	AvgWave /= totalwaves
	AppendToGraph/W=$windowName AvgWave;
	hue = floor(((totalwaves + 1) / (totalwaves + 2)) * 65536);
	hslconvert(hue,61680,34672);
	r = rgbWave[0];g = rgbWave[1];b = rgbWave[2];
	ModifyGraph rgb ($AvgWaveName) = (r,g,b);
	
	Label left "FFT Amp";Label bottom "1/T";	
	ModifyGraph mode=4,marker=19,msize=2,tick=2,mirror=1,fStyle=1,fSize=16,font="Arial";
	Legend/C/N=text0/F=0;
	
	KillWaves rgbWave,imagewave,M_HSL2RGB;
end

Function AvgWavesInList(wavesList,avgName)

	String wavesList,avgName;
	
	Variable totalwaves = ItemsInList(wavesList);
	
	String tempWaveName = StringFromList(0,wavesList);
	Wave tempWave = $tempWaveName;
	Duplicate/O tempWave, $avgName;
	Wave avgWave = $avgName;
	
	avgWave = 0;
	Variable i;
	for(i=0;i<totalwaves;i+=1)
		tempWaveName = StringFromList(i,wavesList);
		Wave tempWave = $tempWaveName;
		avgWave += tempWave;
	endfor
	
	avgWave /= totalwaves;
end

Function batchCoerce(fwave,lowsill,highsill,batchsize)

	Wave fwave;
	Variable lowsill,highsill,batchsize;

	Variable numbatches = numpnts(fwave) / batchsize;
	
	String batchName = NameofWave(fwave) + "_b";
	Duplicate/O fwave $batchname;
	
	Wave coercedWave = $batchname;
	
	Make/D/O/N=(batchsize) tempbCoerwave;
	Variable i;
	for(i=0;i<numbatches;i+=1)
		tempbCoerwave = fwave[p+i*batchsize];
		Coerce2Window(tempbCoerwave,lowsill,highsill);
		coercedWave[i*batchsize,i*batchsize+batchsize-1] = tempbCoerwave[p-i*batchsize];
	endfor
	
	KillWaves tempbCoerwave;
end

Function Coerce2Window(fwave,lowsill,highsill)

	Wave fwave;
	Variable lowsill,highsill;
	
	Wave/D checkwave;

	Make/D/O/N=(numpnts(fwave)) checkwave;

	Variable i;
	for(i=0;i<numpnts(fwave);i+=1)
		if ((fwave[i] < lowsill) || (fwave[i] > highsill))
			checkwave[i] = 0;
		else
			checkwave[i] = 1;
		endif
	endfor
	
	Variable goodpoints = sum(checkwave);
	Variable findex = 0;
	if (goodpoints == 0)
		fwave = (lowsill + highsill) / 2;
	else
		Make/D/O/N=(goodpoints) tempfwave;
		for (i=0;i<numpnts(checkwave);i+=1)
			if (checkwave[i] == 1)
				tempfwave[findex] = fwave[i];
				findex+=1;
			endif
		endfor
	endif
	
	Variable badpoints = numpnts(fwave) - goodpoints;
	if ((badpoints > 0) && (goodpoints != 0))
		WaveStats/Q tempfwave;
		for (i=0;i<numpnts(checkwave);i+=1)
			if (checkwave[i] == 0)
				fwave[i] = V_avg;
			endif
		endfor
	endif
	
	KillWaves checkwave,tempfwave;

end

Function FilterChunkWaveGroup(noisyWave,otherWavesList,chunkSize)

	Wave noisyWave;
	String otherWavesList;
	Variable chunkSize;
	
	String tempHistName = "tHW" + nameofwave(noisyWave);
	
	FindBadPointsInChunks(noisyWave,chunkSize,tempHistName);
	
	Wave tempHistoryWave = $tempHistName;
	String tempFilteredName = ApplyFilterMask(noisyWave,tempHistoryWave);
	String noisyName = nameofwave(noisyWave);
	Duplicate/O $tempFilteredName,$noisyName;
	KillWaves $tempFilteredName;
	String tempOtherWaveName;
	Variable i;
	for(i=0;i<ItemsInList(otherWavesList);i+=1)
		tempOtherWaveName = StringFromList(i,otherWavesList);
		Wave tempOtherWave = $tempOtherWaveName;
		tempFilteredName = ApplyFilterMask(tempOtherWave,tempHistoryWave);
		Duplicate/O $tempFilteredName, $tempOtherWaveName;
		KillWaves $tempFilteredName;
	endfor
	
	KillWaves tempHistoryWave;
end	

Function FilteredAverage(opWave)

	Wave opWave;
	Variable midpoint,median,j,opwavelength,numPts,deletedPts;

	WaveStats/Q opWave;
	Sort opWave,opWave;
	midpoint = floor(numpnts(opWave) / 2);
	median = opWave[midpoint];
	deletedPts = 0;
	opwavelength = numpnts(opWave);
	for (j=0;j<opwavelength;j+=1)
		if (opWave[j] > (median + 2 * V_sdev) || opWave[j] < (median - 2 * V_sdev))
			DeletePoints (j-deletedPts), 1, opWave;
			deletedPts += 1;
		endif
	endfor

	do
		WaveStats/Q opWave;
		deletedPts = 0;
		opwavelength = numpnts(opWave);
		for (j=0;j<opwavelength;j+=1)
			if (opWave[j] > (V_avg + 2 * V_sdev) || opWave[j] < (V_avg - 2 * V_sdev))
				DeletePoints (j-deletedPts), 1, opWave;
				deletedPts += 1;
			endif
		endfor
	while (deletedPts != 0)

	Variable avgval = mean(opWave);
	Return avgval;		
end

Function FilteredAvgAndStdDev(opWave)

	Wave opWave;
	Variable midpoint,median,j,opwavelength,numPts,deletedPts;

	WaveStats/Q opWave;
	Sort opWave,opWave;
	midpoint = floor(numpnts(opWave) / 2);
	median = opWave[midpoint];
	deletedPts = 0;
	opwavelength = numpnts(opWave);
	for (j=0;j<opwavelength;j+=1)
		if (opWave[j] > (median + 2 * V_sdev) || opWave[j] < (median - 2 * V_sdev))
			DeletePoints (j-deletedPts), 1, opWave;
			deletedPts += 1;
		endif
	endfor

	do
		WaveStats/Q opWave;
		deletedPts = 0;
		opwavelength = numpnts(opWave);
		for (j=0;j<opwavelength;j+=1)
			if (opWave[j] > (V_avg + 2 * V_sdev) || opWave[j] < (V_avg - 2 * V_sdev))
				DeletePoints (j-deletedPts), 1, opWave;
				deletedPts += 1;
			endif
		endfor
	while (deletedPts != 0)

	Make/D/O/N=1 filteredavgtempwave, filteredstddevtempwave;
	filteredavgtempwave[0] = mean(opWave);
	filteredstddevtempwave[0] = V_sdev;
end

Function FilteredStdDev(opWave)

	Wave opWave;
	Variable midpoint,median,j,opwavelength,numPts,deletedPts;

	WaveStats/Q opWave;
	Sort opWave,opWave;
	midpoint = floor(numpnts(opWave) / 2);
	median = opWave[midpoint];
	deletedPts = 0;
	opwavelength = numpnts(opWave);
	for (j=0;j<opwavelength;j+=1)
		if (opWave[j] > (median + 2 * V_sdev) || opWave[j] < (median - 2 * V_sdev))
			DeletePoints (j-deletedPts), 1, opWave;
			deletedPts += 1;
		endif
	endfor

	do
		WaveStats/Q opWave;
		deletedPts = 0;
		opwavelength = numpnts(opWave);
		for (j=0;j<opwavelength;j+=1)
			if (opWave[j] > (V_avg + 2 * V_sdev) || opWave[j] < (V_avg - 2 * V_sdev))
				DeletePoints (j-deletedPts), 1, opWave;
				deletedPts += 1;
			endif
		endfor
	while (deletedPts != 0)

	Return V_sdev;
end

//See FindBadPointsMLoA
Function FilterPass(dataWave,dataHistory,lowsill,highsill)

	Wave dataWave,dataHistory;
	Variable lowsill, highsill;
	
	Variable i,pointsFiltered=0;
	for(i=0;i<numpnts(dataWave);)
		if ((dataWave[i]<lowsill)||(dataWave[i]>highsill))
			DeletePoints i,1,dataWave;
			UpdateFilterHistory(dataHistory,i);
			pointsFiltered += 1;
		else
			i+=1;
		endif
	endfor
	
	Return pointsFiltered;
end

//See FilterChunkWaveGroup
Function FindBadPointsInChunks(dataWave,chunkSize,dataHistoryName)

	Wave dataWave;
	Variable chunkSize;
	String dataHistoryName;
	
	Make/D/N=0 $dataHistoryName;
	Wave dataHistoryWave = $dataHistoryName;
	Variable numChunks = floor(numpnts(dataWave)/chunkSize);
	Variable i;
	for(i=0;i<numChunks;i+=1)
		Duplicate/O/R=[i*chunkSize,(i+1)*chunkSize-1] dataWave, tempFBPIC;
		Make/D/O/N=(chunkSize) tempHistFBPIC;
		tempHistFBPIC = 1;
		FindBadPointsMLoA(tempFBPIC,tempHistFBPIC);
		Concatenate/NP {tempHistFBPIC}, dataHistoryWave;
	endfor
	
	if (numpnts(dataHistoryWave) != numpnts(dataWave))
		Redimension/N=(numpnts(dataWave)) dataHistoryWave;
	endif
	
	KillWaves tempFBPIC,tempHistFBPIC;
end

//See FindBadPointsInChunks
Function FindBadPointsMLoA(dataWave,dataHistory)
	// Median/Loop on Average Filter
	Wave dataWave,dataHistory;
	
	Variable targetVal = findMedian(dataWave);
	WaveStats/Q dataWave;
	
	FilterPass(dataWave,dataHistory,targetVal - 2 * V_sdev,targetVal + 2 * V_sdev);
	
	Variable deletedPts
	do
		WaveStats/Q dataWave;
		deletedPts = FilterPass(dataWave,dataHistory,v_avg - 2 * V_sdev,v_avg + 2 * V_sdev);
	while (deletedPts != 0)
end

Function IndexedStdDev(ListofWaves,index)

	String ListofWaves;
	Variable Index;
	
	Variable totWaves = ItemsInList(ListofWaves);
	Make/D/O/N=(totWaves) tempIndStdDev;
	
	Variable i;
	String tempListWaveName;
	for(i=0;i<totWaves;i+=1)
		tempListWaveName = StringFromList(i,ListofWaves);
		Wave tempListWave = $tempListWaveName;
		tempIndStdDev[i]= tempListWave[index];
	endfor
	
	WaveStats/Q tempIndStdDev;
	KillWaves tempIndStdDev;
	Return V_sdev;
end

//Creates/modifies the mask wave for targetWave by masking out
//(setting to 0) the points outside of the point interval defined by 
//cursor A and cursor B.  Returns the mask wave name.
//See MaskWavePoints()
Function/S MaskCsrTails()

	Wave yWave = CsrWaveRef(A)
	String origFolder = GetDataFolder(1)
	String yWavePath = GetWavesDataFolder(yWave,1)
	SetDataFolder $yWavePath
	
	String MaskWaveName = MaskTails(yWave,pcsr(A),pcsr(B))
	
	SetDataFolder $origFolder
	Return MaskWaveName
end

//See MaskWavePoints()
//Creates/modifies the mask wave for the wave on the top graph selected 
//by the A cursor so that the point selected by the A cursor is masked out.
Function/S MaskPt()

	Wave yWave = CsrWaveRef(A)
	String origFolder = GetDataFolder(1)
	String yWavePath = GetWavesDataFolder(yWave,1)
	SetDataFolder $yWavePath
	
	String MaskWaveName = MaskWavePoints(yWave,pcsr(A),pcsr(A))
	
	SetDataFolder $origFolder
	Return MaskWaveName
end

//See MaskWavePoints()
//Creates/modifies the mask wave for the wave on the top graph selected 
//by the A cursor so that the points in the interval created by the A and B 
//cursors are masked out (set to 0 in the mask wave).
Function/S MaskPts()

	Wave yWave = CsrWaveRef(A)
	String origFolder = GetDataFolder(1)
	String yWavePath = GetWavesDataFolder(yWave,1)
	SetDataFolder $yWavePath
	
	String MaskWaveName = MaskWavePoints(yWave,pcsr(A),pcsr(B))
	
	SetDataFolder $origFolder
	Return MaskWaveName
end

//Creates/modifies the mask wave for targetWave by masking out
//(setting to 0) the points outside of the point interval defined by 
//pointA and pointB.  Returns the mask wave name.
//See MaskWavePoints()
Function/S MaskTails(targetWave,pointA,pointB)
	
	Wave targetWave
	Variable pointA,pointB
	
	Variable lowPoint,highPoint
	lowPoint=min(pointA,pointB)
	highPoint=max(pointA,pointB)
	
	String MaskWaveName
	if (lowPoint>0)
		MaskWaveName = MaskWavePoints(targetWave,0,lowPoint-1)
	endif
	if (highPoint<numpnts(targetWave)-1)
		MaskWaveName = MaskWavePoints(targetWave,highPoint+1,numpnts(targetWave)-1)
	endif
	Return MaskWaveName
end

//Creates a wave named "m_" + <name of targetWave> with the same 
//scaling and number of points as targetWave with each entry set to 1 
//excpet the points between pointA and pointB inclusive which are set 
//to 0.  If such a wave already exists, it is not overwritten, but the 
//points between pointA and pointB are set to 0.  Returns the name of 
//the created wave.
Function/S MaskWavePoints(targetWave,pointA,pointB)

	Wave targetWave
	Variable pointA,pointB
	
	//Create mask wave if it does not exist
	String maskWaveName="m_"+NameOfWave(targetWave)
	if (waveexists($maskWaveName))
		Wave maskWave=$maskWaveName
	else
		Duplicate/O targetWave, $maskWaveName
		Wave maskWave=$maskWaveName
		maskWave=1
	endif
	
	//Find which of pointA and pointB is lower
	Variable lowPoint,highPoint
	lowPoint=min(pointA,pointB)
	highPoint=max(pointA,pointB)
	
	//Set point range to 0
	Variable i
	for(i=lowPoint;i<highPoint+1;i+=1)
		maskWave[i]=0
	endfor
	
	Return maskWaveName
end

Function PiecewiseAvgStdDev(w,numPtsToAvg)

	Wave w;
	Variable numPtstoAvg;   		//number of data points per field setting
	Variable i;
	String AvgedWaveName,StdDevWaveName;
	Wave tempAvgWave, tempStdDevWave;
	
	MakeDivByN(w,numPtsToAvg);
	String wName = NameofWave(w);
	sprintf AvgedWaveName "%s_%s",wName,"Avg";
	sprintf StdDevWaveName "%s_%s",wName,"StdDev";
	Variable wlength = numpnts(w);

	Variable numPtsAvgedWave = ceil(wlength/numPtsToAvg);
	Make /D/O/N=(numPtsAvgedWave) $AvgedWaveName,$StdDevWaveName;
	Wave tempAvgWave = $AvgedWaveName ;
	Wave tempStdDevWave = $StdDevWaveName;

	for (i=0;i<numPtsAvgedWave;i+=1)
		Duplicate /O/R=[i*numPtsToAvg,numPtsToAvg*i+numPtsToAvg-1] w opWave;
		WaveStats/Q opWave
		tempAvgWave[i]=v_avg
		tempStdDevWave[i]=v_sdev
	endfor
	
	KillWaves opWave
end

//Performs PiecewiseAvgStdDev on all waves in the current data folder
//and all subfolders with waveHead in the name.  num2avg is the number
//of points to average together at a time with PiecewiseAvgStdDev (numPtsToAvg)
Function RecursiveAvg(waveHead,num2avg)

	String waveHead
	Variable num2avg
	
	Variable i
	Variable initialWaveTotal=CountObjects("",1)
	String tempWaveName,listOfMatches=""
	
	//Make a list of all the waves with waveHead in the name
	//before doing any averaging (because averaging might produce more
	//waves with waveHead in the name and cause an infinite loop
	for(i=0;i<initialWaveTotal;i+=1)
		tempWaveName=GetIndexedObjName("",1,i)
		if(strsearch(tempWaveName,waveHead,Inf,1)>-1)
			listOfMatches=AddListItem(tempWaveName,listOfMatches)
		endif
	endfor
	
	//Average all the waves found with waveHead in the title
	for(i=0;i<ItemsInList(listOfMatches);i+=1)
		tempWaveName=StringFromList(i,listOfMatches)
		PiecewiseAvgStdDev($tempWaveName,num2avg)
	endfor
	
	//Recurse RecursiveAvg() on all subfolders of the current data folder
	String tempFolderName
	for(i=0;i<CountObjects("",4);i+=1)
		tempFolderName=GetIndexedObjName("",4,i)
		SetDataFolder $tempFolderName
		RecursiveAvg(waveHead,num2avg)
		SetDataFolder ::
	endfor
	
end

Function StdDevsOfListWaves(ListofWaves,outputWaveName)

	String ListofWaves, outputWaveName;
	
	String tempWaveName = StringFromList(0,ListofWaves);
	Wave tempWave = $tempWaveName;
	Variable LengthofWaves = numpnts(tempWave);
	Make/D/O/N=(LengthofWaves) $outputWaveName;
	Wave outputWave = $outputWaveName;
	
	Variable totWaves = ItemsInList(ListofWaves);
	Variable i;
	for(i=0;i<totWaves;i+=1)
		outputWave[i] = IndexedStdDev(ListofWaves,i);
	endfor
end

Function/S TightFreqAvg(w,numPtsToAvg)

	Wave w;
	Variable numPtstoAvg;   		//number of data points per field setting
	Variable i;
	String AvgedWaveName;
	Wave tempAvgWave;
	
	MakeDivByN(w,numPtsToAvg);
	String wName = NameofWave(w);
	sprintf AvgedWaveName "%s_%s",wName,"tAvg";
	Variable wlength = numpnts(w);

	Variable numPtsAvgedWave = ceil(wlength/numPtsToAvg);
	Make /D/O/N=(numPtsAvgedWave) $AvgedWaveName;
	Duplicate /O $AvgedWaveName tempAvgWave;

	for (i=0;i<numPtsAvgedWave;i+=1)
		Duplicate /O/R=[i*numPtsToAvg,numPtsToAvg*i+numPtsToAvg-1] w opWave;
		tempAvgWave[i]=FilteredAverage(opWave);
	endfor
	
	Duplicate /O tempAvgWave $AvgedWaveName;
	KillWaves tempAvgWave,opWave;
	Return AvgedWaveName
end

Function TightFreqAvgStdDev(w,numPtsToAvg)

	Wave w;
	Variable numPtstoAvg;   		//number of data points per field setting
	Variable i;
	String AvgedWaveName,StdDevWaveName;
	Wave tempAvgWave, tempStdDevWave;
	
	MakeDivByN(w,numPtsToAvg);
	String wName = NameofWave(w);
	sprintf AvgedWaveName "%s_%s",wName,"tAvg";
	sprintf StdDevWaveName "%s_%s",wName,"tStdDev";
	Variable wlength = numpnts(w);

	Variable numPtsAvgedWave = ceil(wlength/numPtsToAvg);
	Make /D/O/N=(numPtsAvgedWave) $AvgedWaveName,$StdDevWaveName;
	Duplicate /O $AvgedWaveName tempAvgWave;
	Duplicate /O $StdDevWaveName tempStdDevWave;

	for (i=0;i<numPtsAvgedWave;i+=1)
		Duplicate /O/R=[i*numPtsToAvg,numPtsToAvg*i+numPtsToAvg-1] w opWave;
		FilteredAvgAndStdDev(opWave);
		Wave filteredavgtempwave,filteredstddevtempwave;
		tempAvgWave[i]=filteredavgtempwave[0];
		tempStdDevWave[i]=filteredstddevtempwave[0];
	endfor
	
	Duplicate /O tempAvgWave $AvgedWaveName;
	Duplicate /O tempStdDevWave $StdDevWaveName;
	KillWaves tempAvgWave,opWave,tempStdDevWave,filteredavgtempwave,filteredstddevtempwave;
end

//See FilterPass()
Function UpdateFilterHistory(historyWave,index)

	Wave historyWave;
	Variable index;
	
	Variable i,j;
	for(i=0,j=0;;j+=1)
		if((i==index)&&(historyWave[j]==1))
			break;
		elseif(historyWave[j] == 1)
			i+=1;
		endif
	endfor
	
	historyWave[j] = 0;
end

//Remove all points from fwave that are below lowsill and above highsill
//Create a mask wave with the same length as fwave and entries of 1 for points
//not deleted from fwave and 0 for points deleted.  Return mask wave name.
Function/S WindowFilter(fwave,lowsill,highsill)

	Wave fwave;
	Variable lowsill,highsill;
	
	String maskName = "mask_" + NameOfWave(fwave);
	Make/O/N=(numpnts(fwave)) $maskName
	Wave maskWave=$maskName

	Variable i;
	for(i=0;i<numpnts(fwave);i+=1)
		if ((fwave[i] < lowsill) || (fwave[i] > highsill))
			maskwave[i] = 0;
		else
			maskwave[i] = 1;
		endif
	endfor
	
	ApplyFilterMask(fwave,maskwave);
	
	Return maskName
end
