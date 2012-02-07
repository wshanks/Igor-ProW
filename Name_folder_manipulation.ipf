#pragma rtGlobals=1		// Use modern global access method.

Function AddPrefixAllWaves(waveTag)

	String waveTag;
	
	Variable i,totWaves;
	String oldName,newName;
	
	totWaves = CountObjects("",1);
	for(i=0;i<totWaves;i+=1)
		oldName=GetIndexedObjName("",1,i);
		newName=waveTag+oldName;
		Rename $oldName,$newName;
	endfor
end

Function AddTagAllWaves(waveTag)

	String waveTag;
	
	Variable i,totWaves;
	String oldName,newName;
	
	totWaves = CountObjects("",1);
	for(i=0;i<totWaves;i+=1)
		oldName=GetIndexedObjName("",1,i);
		newName=oldName+waveTag;
		Rename $oldName,$newName;
	endfor
end

Function cleanUpCurveFit()

	KillVariables/Z V_chisq,V_numNaNs,V_numINFs,V_npnts,V_nterms,V_nheld
	KillVariables/Z V_startRow,V_endRow,V_startCol,V_endCol,V_startLayer,V_endLayer
	KillVariables/Z V_startChunk,V_endChunk,S_info
end

Function CleanUpWaveStats()

	KillVariables V_npnts,V_numNaNs,V_numINFs,V_avg,V_Sum,V_sdev;
	KillVariables V_rms,V_adev,V_skew,V_kurt,V_minloc,V_maxloc,V_min;
	KillVariables V_max,V_startRow,V_endRow;
	
end

Function CheckFolderforWave(sourceFolderStr,waveNameStr)

	String sourceFolderStr,waveNamestr; // "" for sourceFolderStr will check current folder
	
	Variable totWaves = CountObjects(sourceFolderStr,1);
	
	Variable i;
	String tempName;
	for(i=0;i<totWaves;i+=1)
		tempName = GetIndexedObjName("",1,i);
		if (stringmatch(tempName,waveNameStr))
			Return 1;
		endif
	endfor
	
	Return 0;
end

Function checkForFolder(folderNameStart)
	String folderNameStart
	
	Variable i
	for(i=0;i<countObjects("",4);i+=1)
		String tempStr=GetIndexedObjName("",4,i)
		if (strsearch(GetIndexedObjName("",4,i),folderNameStart,0,2)==0)
			Return i
		endif
	endfor
	
	Return -1
end

Function/S CoefWaveName(filenum,folder,namemode)
	Variable filenum,namemode
	String folder
	
	Return "coef_"+CreateWaveName(filenum,folder,"",namemode)
end

Function CountWaves(waveTag)

	String waveTag;
	
	Variable i=0;
	String tempName = waveTag + num2str(i);
	if (waveexists($tempName)==0)
		Return 0;
	endif
	for(i=0;waveexists($tempName);i+=1)
		tempName = waveTag + num2str(i);
	endfor
	
	Return i-1;
end

Function/S CreateIgorFolderName(windowsFolder,filenum)

	// Windows folder name is in the form DDMMMYY (eg 06Sep08 for September 6th, 2008)
	String windowsFolder;
	Variable filenum;
	
	String fileNumStr;
	sprintf fileNumStr,  "%.3d",filenum;
	
	String IgorFolderName = windowsFolder[2,4] + windowsFolder[0,1] + " " + fileNumStr;
	Return IgorFolderName;
end

Function/S CreateWaveName(fileNum,folderName,waveNameTag,nameMode)

	Variable fileNum,nameMode;
	String folderName,waveNameTag;
	
	String fileName,fileNumStr;
	String dateStr = folderName[0,1];
	String monthStr = folderName[2,4];
	
	Switch (nameMode)
		Case 0:
			sprintf fileNumStr,  "%.3d",filenum;
			fileName = monthStr+dateStr+"n"+fileNumStr+waveNameTag;
			break;
		default:
			// Insert error handling here
	endswitch

	Return fileName;
end

Function/S CreateWindowName(folder,filenum,idtag)

	String folder,idtag
	Variable filenum
	
	String dateStr = folder[0,1];
	String monthStr = folder[2,4];
	String fileNumStr
	sprintf fileNumStr,  "%.3d",filenum;
	
	String windowName=idtag+"_"+monthStr+"_"+dateStr+"_"+fileNumStr
	
	Return windowName
end

Function Dup1Wave2Folder(w,FolderName)

	Wave w;
	String FolderName;
	
	String origFolder = GetDataFolder(1);
	String wName = NameOfWave(w);
	
	NewDataFolder/O $FolderName
	SetDataFolder $FolderName;
	Duplicate/O w $wName;
	SetDataFolder $origFolder;
	
end

Function DupListWaves2Folder(fileNum,folderName,waveTagList,nameMode,newFolderName,suffixStr)

	Variable fileNum,nameMode;
	String folderName,waveTagList,newFolderName,suffixStr;
	
	Variable totWaves = ItemsInList(waveTagList);
	Variable i;
	String tempWaveName,tempWaveTag;
	
	NewDataFolder/O $newFolderName;
	
	for(i=0;i<totWaves;i+=1)
		tempWaveTag = StringFromList(i,waveTagList);
		tempWaveName = CreateWaveName(fileNum,folderName,tempWaveTag,nameMode);
		tempWaveName += suffixStr;
		Wave tempWave = $tempWaveName;
		Dup1Wave2Folder(tempWave,newFolderName);
	endfor
end

Function DupListWavesRename(fileNum,folderName,waveTagList,nameMode,newFolderName,suffixStr,newNameList)

	Variable fileNum,nameMode;
	String folderName,waveTagList,newFolderName,suffixStr,newNameList;
	
	Variable totWaves = ItemsInList(waveTagList);
	Variable i;
	String tempWaveName,tempWaveTag,newTempName;
	
	NewDataFolder/O $newFolderName;
	
	for(i=0;i<totWaves;i+=1)
		tempWaveTag = StringFromList(i,waveTagList);
		tempWaveName = CreateWaveName(fileNum,folderName,tempWaveTag,nameMode)+suffixStr;
		Wave tempWave = $tempWaveName;
		newTempName = StringFromList(i,newNameList);
		DupWaveRename(tempWave,newFolderName,newTempName);
	endfor
end

Function DupWaves2Folder(newFolderName)

	String newFolderName;
	String origFolder = GetDataFolder(1);

	NewDataFolder /O $newfolderName;
	
	Variable totWaves = CountObjects("",1);
	
	Variable i;
	String tempName,newtempName;
	Wave temp0;
	for(i=0;i<totWaves;i+=1)
		tempName = GetIndexedObjName("",1,i);
		Duplicate/O $tempName temp0;
		SetDataFolder $newfolderName;
		Duplicate /O temp0 $tempName;
		SetDataFolder $origFolder
	endfor
	
	KillWaves temp0;
end

Function DupWaveRename(targetWave,folderName,newWaveName)

	Wave targetWave;
	String folderName,newWaveName;
	
	String origFolder = GetDataFolder(1);
	NewDataFolder/O $folderName;
	SetDataFolder folderName;
	Duplicate/O targetWave,$newWaveName;
	SetDataFolder $origFolder;

end

//Searches the wave names in the current data folder backwards
//for a wave with the string waveTag in its name.  Returns the name
//of the first wave found.
Function/S FindWaveWithTag(waveTag)

	String waveTag
	
	Variable totWaves = CountObjects("",1)
	String tempWaveName,foundWaveName
	
	Variable j
	For(j=0;j<totWaves;j+=1)
		tempWaveName=GetIndexedObjName("",1,j)
		if(strsearch(tempWaveName,waveTag,Inf,1)>0)
			foundWaveName=tempWaveName
			break
		endif
	endfor

	Return foundWaveName
end

Function MoveListWaves2Folder(fileNum,folderName,waveTagList,nameMode,newFolderName,suffixStr)

	Variable fileNum,nameMode;
	String folderName,waveTagList,newFolderName,suffixStr;
	
	Variable totWaves = ItemsInList(waveTagList);
	Variable i;
	String tempWaveName,tempWaveTag;
	
	NewDataFolder/O $newFolderName;
	
	for(i=0;i<totWaves;i+=1)
		tempWaveTag = StringFromList(i,waveTagList);
		tempWaveName = CreateWaveName(fileNum,folderName,tempWaveTag,nameMode)+suffixStr;
		Wave tempWave = $tempWaveName;
		MvWave(tempWave,newFolderName);
	endfor
end

Function MoveListWavesRename(fileNum,folderName,waveTagList,nameMode,newFolderName,suffixStr,newNameList)

	Variable fileNum,nameMode;
	String folderName,waveTagList,newFolderName,suffixStr,newNameList;
	
	Variable totWaves = ItemsInList(waveTagList);
	Variable i;
	String tempWaveName,tempWaveTag,newTempName;
	String origFolder = GetDataFolder(1);
	
	NewDataFolder/O $newFolderName;
	
	for(i=0;i<totWaves;i+=1)
		tempWaveTag = StringFromList(i,waveTagList);
		tempWaveName = CreateWaveName(fileNum,folderName,tempWaveTag,nameMode)+suffixStr;
		Wave tempWave = $tempWaveName;
		newTempName = StringFromList(i,newNameList);
		MvWaveRename(tempWave,newFolderName,newTempName);
	endfor
end

Function MvWave(targetWave,folderName)

	Wave targetWave;
	String folderName;
	
	String origFolder = GetDataFolder(1);
	String targetWaveName = nameofwave(targetWave);
	NewDataFolder/O $folderName;
	SetDataFolder $folderName;
	Duplicate/O targetWave,$targetWaveName;
	SetDataFolder $origFolder;
	
	KillWaves targetWave;

end

Function MvWaveRename(targetWave,folderName,newWaveName)

	Wave targetWave;
	String folderName,newWaveName;
	
	String origFolder = GetDataFolder(1);
	NewDataFolder/O $folderName;
	SetDataFolder $folderName;
	Duplicate/O targetWave,$newWaveName;
	SetDataFolder $origFolder;
	KillWaves targetWave;

end

Function/S ParentFolderName()

	String origFolder = GetDataFolder(1);
	SetDataFolder ::;
	String parentFolder = GetDataFolder(0);
	SetDataFolder $origFolder;
	Return parentFolder;
end

Function/S RenameAnalysisWave(analWaveName)
	
	String analWaveName;
	Wave analWave = $analWaveName;
	String parentName = ParentFolderName();
	String partialParentName = ReplaceString("'",ReplaceString(" ",ParentName,""),"");
	partialParentName = partialParentName[0,7];
	String newName = partialParentName+ nameofwave(analWave);
	RenameWave(analWave,newName);
	Return newName;
end

Function RenameCurrentFolder(newName)
	String newName
	
	String currentFolderPath=GetDataFolder(1)
	
	RenameDataFolder $currentFolderPath,$newName
end

Function RenameWave(oldWave,newName)
	Wave oldWave;
	String newName;
	
	Duplicate/O oldWave,$newName;
	KillWaves oldWave;
end

Function SetNewDataFolder(newFullPath)
	//Sets the current data folder to folder given by newFullPath
	//Any subdirectories from the root directory that do not exist are created
	String newFullPath
	
	SetDataFolder root:
	
	Variable test = ItemsInList(newFullPath,":")
	String currentFolder
	Variable i
	for (i=1;i<ItemsInList(newFullPath,":");i+=1)
		currentFolder=StringFromList(i,newFullPath,":")
		currentFolder=ReplaceString("'",currentFolder,"")
		NewDataFolder/O/S $currentFolder
	endfor
		 
end

Function WaveMatch(wave1,wave2)
	//Returns 1 if wave1 and wave2 are the same wave
	//Returns 0 if they are not
	Wave wave1,wave2
	
	if (stringmatch(GetWavesDataFolder(wave1,2),GetWavesDataFolder(wave2,2)))
		Return 1
	else
		Return 0
	endif
end