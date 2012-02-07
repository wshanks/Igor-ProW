#pragma rtGlobals=1		// Use modern global access method.

#include <Concatenate Waves>

//Functions that chop up parts of waves, insert points in, etc.

Function addPoint2Wave(targetWave,targetIndex,value)
	//Adds an entry with value of value to targetWave at index targetIndex.
	//targetIndex should be either an existing index in targetWave or the index just after the last point
	//Otherwise, an error is printed to command line.
	Wave targetWave
	Variable targetIndex,value
	
	if (targetIndex == numpnts(targetWave))
		AppendPoint2Wave(targetWave,value)
	elseif (targetIndex < numpnts(targetWave))
		InsertPoints targetIndex,1,targetWave
		targetWave[targetIndex]=value
	else
		printf "Error in addPoint2Wave function occurred!\r"
	endif
end

Function AppendPoint2Wave(targetWave,value)
	// Adds a new entry on to the end of targetWave and sets its value to value
	// Target wave needs to be a one dimensional wave
	Wave targetWave
	Variable value
	
	Variable targetPnts = numpnts(targetWave)+1
	Redimension/N=(targetPnts) targetWave
	
	targetWave[targetPnts-1] = value
end

// Chop up all waves in current folder into a series of new folders
// Folders names are foldername+number+postFoldTag
// Waves in the folder have the original wave names plus "_"+number+postWaveTag
// "number" is given by initNum+NumStep*number of the section
// SecLength is the number of points to be chopped out of the original waves for each section
Function ChopUpWaves(foldername,postFoldTag,postWaveTag,initNum,NumStep,SecLength)

	String foldername,postFoldTag,postWaveTag;
	Variable initNum,NumStep,SecLength;

	String origFolder = GetDataFolder(1);
	String tempFolder,curTag;
	Variable i,initSecPnt,finalSecPnt,lastPnt,endFlag=0;
	
	String tempName = GetIndexedObjName("",1,0);
	Wave tempWave = $tempName;
	lastPnt = numpnts(tempWave) - 1;
	
	for(i=0;endFlag<1;i+=1)
		
		curTag = num2str(initNum+i*NumStep)+postFoldTag;
		tempFolder = foldername+curTag;
		DupWaves2Folder(tempFolder);
		SetDataFolder $tempFolder;
		
		initSecPnt = i*SecLength;
		finalSecPnt=(i+1)*SecLength-1;
		if (finalSecPnt >= lastPnt)
			finalSecPnt = lastPnt;
			endFlag = 2;
		endif
		
		TrimAllWavesSec(initSecPnt,finalSecPnt);
		curTag = "_"+num2str(initNum+i*NumStep)+postWaveTag;
		AddTagAllWaves(curTag);
	
		SetDataFolder $origFolder;
	endfor

end

Function ConcatEnumWaves(waveTag,outputName,firstIndex,totalWaves)

	String waveTag,outputName
	Variable firstIndex,totalWaves
	
	String concatStr=""
	Variable i
	for(i=firstIndex;i<firstIndex+totalWaves;i+=1)
		concatStr+=waveTag+num2istr(i)+";"
	endfor
	concatStr=RemoveEnding(concatStr)
	
	ConcatenateWavesInList(outputName,concatStr)
	
	String delWaveName
	for(i=firstIndex;i<firstIndex+totalWaves;i+=1)
		delWaveName=waveTag+num2istr(i)
		Wave delWave=$delWaveName
		KillWaves delWave
	endfor
end

//Stores the last pts2keep points of each chunksize set of points from oldWave in newWave
//oldWave and newWave must already have the correct number of points (no error checking of this).
//See TakeEndChunks()
Function CopyEndChunks(oldWave,newWave,chunksize,pts2keep)

	Wave oldWave,newWave
	Variable chunksize,pts2keep
	
	Variable numchunks = floor(numpnts(oldWave) / chunksize)
	Variable i,j
	for (i=0;i<numchunks;i+=1)
		for (j=0;j<pts2keep;j+=1)
			newWave[i*pts2keep+j] = oldWave[i*chunksize-pts2keep+j]
		endfor
	endfor
end

Function dPtsCtrlWaveNaN(targetWave,ctrlWave)
	//Delete all points in targetWave for which ctrlWave==NaN
	Wave targetWave,ctrlWave
	
	Duplicate/O ctrlWave,tempMaskWave
	tempMaskWave=1
	
	Variable i
	for(i=0;i<numpnts(ctrlWave);i+=1)
		if(numtype(ctrlWave[i])==2)
			tempMaskWave[i]=0
		endif
	endfor
	for(i=numpnts(targetWave);i>=0;i-=1)
		if(tempMaskWave[i]==0)
			DeletePoints i,1,targetWave
		endif
	endfor
	KillWaves tempMaskWave
end

Function DptAll()

	Wave yWave = CsrWaveRef(A)
	String origFolder = GetDataFolder(1)
	String yWavePath = GetWavesDataFolder(yWave,1)
	SetDataFolder $yWavePath

	Variable totWaves = CountObjects("",1);
	String tempName = GetIndexedObjName("",1,0);
	
	Variable point = xcsr(A)
	
	Variable i;
	for(i=0;i<totWaves;i+=1)
		tempName = GetIndexedObjName("",1,i);	
		DeletePoints point,1,$tempName;
	endfor

	SetDataFolder $origFolder
end

Function dPtsCtrlWaveVal(targetWave,ctrlWave,value)
	//Delete all points in targetWave for which ctrlWave==value
	Wave targetWave,ctrlWave
	Variable value
	
	Duplicate/O ctrlWave,tempMaskWave
	tempMaskWave=1
	
	Variable i
	for(i=0;i<numpnts(tempCtrlWave);i+=1)
		if(ctrlWave[i]==value)
			tempMaskWave[i]=0
		endif
	endfor
	for(i=numpnts(targetWave);i>=0;i-=1)
		if(tempMaskWave[i]==0)
			DeletePoints i,1,targetWave
		endif
	endfor
	KillWaves tempMaskWave
end

Function ExciseAllWaves()
	//Deletes all points between cursor A and cursor B from all waves in the 
	//same data folder as the y wave of the trace with cursors A and B on it on the
	//top graph.

	Wave yWave = CsrWaveRef(A)
	String origFolder = GetDataFolder(1)
	String yWavePath = GetWavesDataFolder(yWave,1)
	SetDataFolder $yWavePath
	
	Variable numberofpoints = abs(pcsr(B) - pcsr(A)) + 1;
	
	Variable firstpoint;
	If (pcsr(B) > pcsr(A))
		firstpoint = pcsr(A);
	else
		firstpoint = pcsr(B);
	endif
	
	String tempName = GetIndexedObjName("",1,0);
	Variable totWaves = CountObjects("",1);
	Variable i;
	for(i=0;i<totWaves;i+=1)
		tempName = GetIndexedObjName("",1,i);
		DeletePoints firstpoint,numberofpoints,$tempName;
	endfor
	
	if(firstPoint!=0)
		firstPoint-=1
	endif
	String CsrATrace=CsrTraceName(0,"")
	Cursor A $CsrATrace firstPoint
	String CsrBTrace=CsrTracename(1,"")
	Cursor B $csrBTrace firstPoint

	SetDataFolder $origFolder
end

Function ExtendWavebyOne(targetWave,newVal)

	Wave targetWave;
	Variable newVal;
	
	Variable newPntTot = numpnts(targetWave) + 1;
	Redimension/N=(newPntTot) targetWave;
	targetWave[newPntTot-1] = newVal;
	
end

Function extract1Dfrom4D(indexWave,targetWave,outputWaveName)
	//Extracts a 1D wave from 4D targetWave and copies it to outputWaveName
	//The indices of the 4D targetWave are input in the indexWave.  One of these values should 
	//be nan.  The dimension for that index is the one that will become the row dimension
	//of the output wave
	Wave indexWave,targetWave
	String outputWaveName

	Variable outputIndex=FindNaNIndex(indexWave)
	
	Variable outputLength=DimSize(targetWave,outputIndex)
	Make/O/D/N=(outputLength) $outputWaveName
	Wave outputWave=$outputWaveName
	
	Variable i
	for(i=0;i<numpnts(outputWave);i+=1)
		indexWave[outputIndex]=i
		outputWave[i]=targetWave[indexWave[0]][indexWave[1]][indexWave[2]][indexWave[3]]
	endfor
end

Function InsertEntryMWave(targetWave,insDim,insDimIndex,value)
	//Inserts a new entry in the insDim dimension of targetWave and sets the new entry to value
	//In a multidimensional wave, the new entry will be more than one a point (a row, column, or layer).
	Wave targetWave
	Variable insDim,insDimIndex,value

	Make/O/N=4 insMWaveDims
	//Make the new entry in targetWave
	if (insDimIndex<dimSize(targetWave,insDim))
		InsertPoints/M=(insDim) insDimIndex,1,targetWave
	elseif (dimSize(targetWave,insDim)==insDimIndex)
		insMWaveDims = -1
		insMWaveDims[insDim]=dimSize(targetWave,insDim)+1
		Redimension/N=(insMWaveDims[0],insMWaveDims[1],insMWaveDims[2],insMWaveDims[3]) targetWave
	endif
	
	Make/O/N=4 insMWaveLDims,insMWaveUDims
	//Set the entry in targetWave to value
	Variable i
	for(i=0;i<4;i+=1)
		if((i<(WaveDims(targetWave)-1))&&i!=insDim)
			insMWaveLDims[i]=0
			insMWaveUDims[i]=DimSize(targetWave,i)-1
		elseif (i==insDim)
			insMWaveLDims[i]=insDimIndex
			insMWaveUDims[i]=insDimIndex
		else
			insMWaveLDims[i]=-1
			insMWaveUDims[i]=-1
		endif
	endfor
	SetPointsMWave(targetWave,insMWaveLDims[0],insMWaveUDims[0],insMWaveLDims[1],insMWaveUDims[1],insMWaveLDims[2],insMWaveUDims[2],insMWaveLDims[3],insMWaveUDims[3],value)
	KillWaves insMWaveDims,insMWaveLDims,insMWaveUDims
end

Function MakeDivByN (w,N)
// This function takes in a wave w and removes the appropriate number of 
// points at the end of the wave to make the length of the wave divisible by N
	
	Wave w;
	Variable N;
	
	Variable numPtsToDelete= mod(numpnts(w),N);
	Variable firstPtToDelete=numpnts(w)-numPtsToDelete;
	DeletePoints firstPtToDelete,numPtsToDelete,w;
end

//Finds the complex value of targetWave at the point recPt and takes its projection 
//onto the unit vector in the complex plane with angle phasAng relative to the x axis.  
//Appends this projected value to recWave.
Function RecCmplxProj(targetWave,phasAng,recPt,recWave)

	Wave/C targetWave
	Wave recWave
	Variable phasAng,recPt

	Variable/C Cfft=targetWave[recPt]
	Variable projVal=ProjectCmplx(Cfft,phasAng)
	AppendPoint2Wave(recWave,projVal)
end

//Append the imaginary part of the value of the complex wave targetWave at point
//wavePt to recWave
Function RecImagWavePt(targetWave,wavePt,recWave)

	Wave/C targetWave
	Wave recWave
	Variable wavePt
	
	Variable value=imag(targetWave[wavePt])
	
	AppendPoint2Wave(recWave,value)
end

//Append the magnitude of the value of the complex wave targetWave at point
//wavePt to recWave
Function RecMagWavePt(targetWave,wavePt,recWave)

	Wave/C targetWave
	Wave recWave
	Variable wavePt
	
	Variable value=cabs(targetWave[wavePt])
	
	AppendPoint2Wave(recWave,value)
end

//Append the real part of the value of the complex wave targetWave at point
//wavePt to recWave
Function RecRealWavePt(targetWave,wavePt,recWave)

	Wave/C targetWave
	Wave recWave
	Variable wavePt
	
	Variable value=real(targetWave[wavePt])
	
	AppendPoint2Wave(recWave,value)
end

//Find the wave with waveTag in its name and append its average to recWave
Function RecordTagAvg(waveTag,recWave)

	String waveTag
	Wave recWave
	
	String targetWaveName=FindWaveWithTag(waveTag)
	Wave targetWave=$targetWaveName
	RecordWaveAvg(targetWave,recWave)
end

//Find the average of targetWave and append it to recWave
Function RecordWaveAvg(targetWave,recWave)

	Wave targetWave, recWave
	
	Variable waveAvg=FilteredAverage(targetWave)
	
	AppendPoint2Wave(recWave,waveAvg)
	
end

//Append the value of targetWave at point wavePt to recWave
Function RecWavePt(targetWave,wavePt,recWave)

	Wave targetWave,recWave
	Variable wavePt
	
	Variable value=targetWave[wavePt]
	
	AppendPoint2Wave(recWave,value)
end

Function RemoveBlanks(targetWave)

	Wave targetWave
	
	Variable i
	for(i=0;i<numpnts(targetWave);i+=1)
		if (numtype(targetWave[i])==2)
			DeletePoints i,1,targetWave
			i-=1
		endif
	endfor
end

Function RemoveOffset()

	Wave yWave = CsrWaveRef(A)
	Variable offset = yWave[pCsr(A)]
	yWave -= offset
	
end

Function RmDupEntries(targetWave)
	//For every section of targetWave that contains multiple entries with the same value,
	//all but one of those entries are deleted
	Wave targetWave
	
	Variable i,tempSecLength
	for(i=0;i<numpnts(targetWave);i+=1)
		tempSecLength=findSecLength(targetWave,i)
		if (tempSecLength>1)
			DeletePoints i+1,(tempSecLength-1),targetWave
		endif
	endfor
end

Function SetPointsMWave(targetWave,r1,r2,co1,co2,l1,l2,ch1,ch2,value)
	// Sets all the points between rows r1 and r2, columns co1 and co2,
	// layers l1 and l2 and chunks ch1 and ch2 to value
	// All dimension2 should be greater than dimension1
	// If certain dimensions do not exist, their dimension indices should be set to -1
	Wave targetWave
	Variable r1,r2,co1,co2,l1,l2,ch1,ch2,value
	
	Variable i,j,k,l
	for(i=r1;i<r2+1;i+=1)
		if (co1<0)
			targetWave[i]=value
		else
			for(j=co1;j<co2+1;j+=1)
				if (l1<0)
					targetWave[i][j]=value
				else
					for(k=l1;k<l2+1;k+=1)
						if (ch1<0)
							targetWave[i][j][k]=value
						else
							for(l=ch1;l<ch2+1;l+=1)
								targetWave[i][j][k][l]=value
							endfor
						endif
					endfor
				endif
			endfor
		endif
	endfor
end

Function/S ShiftOverlapWaves(y1,x1,y2,x2)
	// Shift y2 by the average difference between y2 and y1 for the corresponding
	// points in x1 and x2 that are equal
	// Assumes x1 and x2 are single valued
	
	Wave y1,x1,y2,x2
	Variable tolerance = 1e-15 // Problems with ==
	Variable i,j,pntCnt,Shift
	
	for(i=0,pntCnt=0,Shift=0;i<numpnts(x1);i+=1)
		for(j=0;j<numpnts(x2);j+=1)
			if (abs(x1[i] -x2[j]) < tolerance)
				pntCnt+=1
				Shift+=y2[j]-y1[i]
				break
			endif
		endfor
	endfor
	
	if (pntCnt!=0)
		Shift /= pntCnt
	else
		Shift = 0
	endif
	String outWaveName = nameofwave(y2) + "_sh"
	Duplicate/O y2,$outWaveName
	Wave outWave = $outWaveName
	outWave = y2 - Shift
	Return outWaveName
end

//Replace the value of targetWave at pointNum by the average of the values at
//pointNum +/- 1
Function SmoothPoint(targetWave,pointNum)

	Wave targetWave
	Variable pointNum
	
	targetWave[pointNum]=(targetWave[pointNum-1]+targetWave[pointNum+1])/2
end

//Runs SmoothPoint() on the point currently selected by the A cursor
Function SmoothPt()

	Wave targetWave = CsrWaveRef(A)
	
	targetWave[pCsr(A)]=(targetWave[pCsr(A)-1]+targetWave[pCsr(A)+1])/2
end

Function SortedNoDupWave(targetWave)
	// Wave must have at least two points
	Wave targetWave;
	Variable tolerance = 1e-15
	
	Sort targetWave,targetWave;
	
	Variable i;
	for(i=1;i<numpnts(targetWave);)
		if ( abs(targetWave[i] - targetWave[i-1])<tolerance)
			DeletePoints i,1,targetWave;
		else
			i+=1;
		endif
	endfor
end

//wList is a list of wave names with the 0th element being some scanned parameter (with a fixed
//number of points at each value).
//Takes the last pts2keep points from each chunk of points with the same scan wave value from
//every wave in wList and save them as separate waves named as their original waves with an additional
//suffix.
//folderName,fileNum, and fileNameMode are used for forming the correct wave names.
Function TakeEndChunks(wList,folderName,fileNum,fileNameMode,pts2keep,suffix)
	// Assumes first element in wList is the scanned parameter

	String wList,folderName,suffix
	Variable fileNum,fileNameMode,pts2keep
	
	String scanWaveTag = StringFromList(0,wList);
	String scanWaveName = CreateWaveName(fileNum,folderName,scanWaveTag,fileNameMode)
	Wave scanWave = $scanWaveName;
	Variable oldChunkLength = FindSecLength(scanWave,0)
	Variable numChunks = numpnts(scanWave) / oldChunkLength
	Variable newWaveLen = numChunks * pts2keep
	Make/D/N=(newWaveLen) tempECWave
	
	String tempWaveTag,tempWaveName,newWaveName
	Variable i
	for (i=0;i<ItemsInList(wList);i+=1)
		tempWaveTag = StringFromList(i,wList)
		tempWaveName = CreateWaveName(fileNum,folderName,tempWaveTag,fileNameMode)
		Wave tempWave = $tempWaveName
		CopyEndChunks(tempWave,tempECWave,oldChunkLength,pts2keep)
		newWaveName = tempWaveName + suffix
		Duplicate/O tempECWave, $newWaveName
	endfor
	
	KillWaves tempECWave
end

Function TrimAllWaves()
	//Deletes all points not between cursor A and cursor B from all waves in the 
	//same data folder as the y wave of the trace with cursors A and B on it on the
	//top graph.
	Wave yWave = CsrWaveRef(A)
	String origFolder = GetDataFolder(1)
	String yWavePath = GetWavesDataFolder(yWave,1)
	SetDataFolder $yWavePath

	Variable highend = numpnts(yWave) - max(pcsr(A),pcsr(B)) + 1;
	Variable lowend = min(pcsr(A),pcsr(B))
	Variable highpoint = max(pcsr(A),pcsr(B)) + 1;
	Variable Afin,Bfin
	if (pcsr(A)<pcsr(B))
		Afin=0
		Bfin=pcsr(B)-pcsr(A)
	else
		Afin=pcsr(A)-pcsr(B)
		Bfin=0
	endif
	
	Variable totWaves = CountObjects("",1);
	String tempName = GetIndexedObjName("",1,0);
	
	Variable i;
	for(i=0;i<totWaves;i+=1)
		tempName = GetIndexedObjName("",1,i);
		if (highend > 0)
			DeletePoints highpoint,highend,$tempName;
		endif
		if (lowend > 0)
			DeletePoints 0,lowend,$tempName;
		endif
	endfor
	
	String CsrATrace=CsrTraceName(0,"")
	Cursor A $CsrATrace Afin
	String CsrBTrace=CsrTracename(1,"")
	Cursor B $csrBTrace Bfin

	SetDataFolder $origFolder
end

// Delete all points from all waves in the current data folder outside of firstpoint and lastpoint
Function TrimAllWavesSec(firstpoint,lastpoint)

	Variable firstpoint,lastpoint;
	
	Variable totWaves = CountObjects("",1);
	String tempName = GetIndexedObjName("",1,0);
	Variable highend = numpnts($tempName) - lastpoint + 1;
	Variable lowend = firstpoint;
	Variable highpoint = lastpoint + 1;
	
	Variable i;
	for(i=0;i<totWaves;i+=1)
		tempName = GetIndexedObjName("",1,i);
		if (highend > 0)
			DeletePoints highpoint,highend,$tempName;
		endif
		if (lowend > 0)
			DeletePoints 0,lowend,$tempName;
		endif
	endfor

end

//For the wave with the A cursor on it in the active graph, this 
//function deletes all points outside of the A and B cursors.  
//The wave's offset is adjust to account for the points from the
//front of the wave
Function trimGraphWave()
	
	Wave yWave = CsrWaveRef(A)
	String yWavePath = GetWavesDataFolder(yWave,1)

	Variable highend = numpnts(yWave) - max(pcsr(A),pcsr(B)) + 1;
	Variable lowend = min(pcsr(A),pcsr(B))
	Variable highpoint = max(pcsr(A),pcsr(B)) + 1;
	Variable Afin,Bfin
	if (pcsr(A)<pcsr(B))
		Afin=0
		Bfin=pcsr(B)-pcsr(A)
	else
		Afin=pcsr(A)-pcsr(B)
		Bfin=0
	endif
	
	Variable waveDelta=DimDelta(yWave,0)
	Variable waveOffset=DimOffset(yWave,0)
	Variable newOffset=waveOffset+lowend*waveDelta
	
	if (highend > 0)
		DeletePoints highpoint,highend,yWave
	endif
	if (lowend > 0)
		DeletePoints 0,lowend,yWave
	endif
	
	String CsrATrace=CsrTraceName(0,"")
	Cursor A $CsrATrace Afin
	String CsrBTrace=CsrTracename(1,"")
	Cursor B $csrBTrace Bfin
	
	SetScale/P x,newOffset,waveDelta,yWave	
end

//Delete all points from wave outside of firstPoint and lastPoint
//Currently has no error checking for cases when these points lie 
//outside of the wave
Function trimWave(targetWave,firstPoint,lastPoint)

	Wave targetWave
	Variable firstPoint,lastPoint
	
	Variable totPoints=numpnts(targetWave)
	
	Variable highPnts2Cut=totPoints-lastPoint-1 //points to cut off the end
	Variable highPoint=lastPoint+1 //First point to cut
	//Cut high points first so the point numbers aren't changed after
	//cutting the low points
	if (highPnts2Cut > 0)
		DeletePoints highPoint,highPnts2Cut,targetWave
	endif
	//Cut low points
	if (firstPoint > 0)
		DeletePoints 0,firstPoint,targetWave
	endif
end

//Creates a copy of the wave on the top graph with the A cursor on it and then zeroes
//all the points in the copy between the A and B cursors (inclusive).
Function/S zeroGraphWave()

	Wave yWave = CsrWaveRef(A)
	String origFolder = GetDataFolder(1)
	String yWavePath = GetWavesDataFolder(yWave,1)
	SetDataFolder $yWavePath

	String zeroWaveName=NameOfWave(yWave)+"_Z"
	Duplicate/O yWave,$zeroWaveName
	Wave zeroWave=$zeroWaveName
	
	//Correct for pcsr division by 2 for complex waves
	Variable cFactor=1
	if (WaveType(zeroWave) & 0x01)
		cFactor=2
	endif
	
	zeroWavePnts(zeroWave,cFactor*pcsr(A),cFactor*pcsr(B))
	
	String NoteString="Copy of "+GetWavesDataFolder(yWave,2)+" with points between "
	NoteString+=num2str(cFactor*pcsr(A))+" and "+num2str(cFactor*pcsr(B))
	NoteString+="set to zero byzerographwave()"
	Note zeroWave,"---Above this line copied from source wave----"
	Note zeroWave,NoteString

	SetDataFolder $origFolder
	
	Return zeroWaveName
end

//Set all points between firstPoint and lastPoint (inclusive) of targetWave to zero.
Function zeroWavePnts(targetWave,firstPoint,lastPoint)

	Wave targetWave
	Variable firstPoint,lastPoint
	
	Variable step=1
	if (lastPoint<firstPoint)
		step=-1
	endif
	
	Variable i
	for(i=firstPoint;i<lastPoint+1;i+=step)
		targetWave[i]=0
	endfor
end