#pragma rtGlobals=1		// Use modern global access method.

//There are several functions in this procedure that might need updating:
//one for each of the set of loading functions.

Function GLoadsetup(GLoadmode)

	Variable GLoadmode;

	NVAR GLoadVarExistence;
	Variable varsExist = 0;
	if (NVAR_Exists(GLoadVarExistence) == 1)
		varsExist = 1;
	endif

	Variable/G GLoadVarExistence;
	// Path for April 2008 cooldown
	//String/G gloadfolderpath="C:\\Swelter\\Data\\CooldownApr2008\\SOI7_52\\300mK\\";
	// Path for April 2009 cooldown
	//String/G gloadfolderpath="Y:\\Data\\CooldownApr2009\\SOI7_52\\";
	// Path for May 2009 cooldown
	//String/G gloadfolderpath="Y:\\Data\\Cooldown2May2009\\SOI7_52\\";
	//String/G gloadfolderpath="E:\\My Physics\\Lab Files\\Data Analysis\\Cooldown2May2009\\"
	// Path for September 2008 cooldown
	//String/G gloadfolderpath="C:\\Swelter\\Data\\CooldownSep2008\\SOI7_63\\";
	//Path for November 2008 cooldown on Swelter
	//String/G gloadfolderpath="C:\\Swelter\\Data\\CooldownNov2008\\SOI7_52\\"
	//Path for November 2008 cooldown on Will's laptop
	//String/G gloadfolderpath="Y:\\Data\\CooldownNov2008\\SOI7_52\\"
	//Temp path to load files
	//String/G gloadfolderpath="E:\\My Physics\\Lab Files\\Data Analysis\\2008\\NovemberCooldown\\SOI7_52_CL17\\"
	//Temp path to load files
	//String/G gloadfolderpath="E:\\My Physics\\Lab Files\\Data Analysis\\April 2009\\"
	// Path for April 2010 cooldown 2
	//String/G gloadfolderpath="Y:\\Data\\Cooldown2Apr2010\\CL51\\";
	//String/G gloadfolderpath="E:\\My Physics\\Lab Files\\Data Analysis\\2010\\April Cooldown 2\\"
	//String/G gloadfolderpath="Y:\\Data\\CooldownNov2008\\SOI7_52\\"
	//String/G gloadfolderpath="E:\\My Physics\\Lab Files\\Data Analysis\\Data\\CooldownApr2008\\SOI7_52\\300mK\\"
	String/G gloadfolderpath="E:\\My Physics\\Lab Files\\Data Analysis\Data\\Cooldown2Apr2010\\CL57\\"
	String/G gloadwavelist; // String list of names for waves loaded in by GLoad()
	String/G gloadfilenametag
	Variable/G gloadFirstLine=0
	Variable/G gloadNameMode =0
	
	Variable externAuxExists;
	Switch (GLoadmode)
		Case 0: // Bscan v1
			// Usual wave list for Bscan (as of 9/6/08): Bfield;Frequency;Lock-in Amplitude;Lock-in Phase;34410 Volts;Lock-in ADC2;A;D;B;Time
			externAuxExists = CreateSBscanAuxData();
			SVAR SBscanLoadWaveList,SBscanFileNameTag
			gloadwavelist = SBscanLoadWaveList;
			gloadfilenametag = SBscanFileNameTag
			CleanupSBscanAuxData()
			break;
		Case 1: // Multi-Ampscan v1
			// Usual wave list for Multi-Ampscan (as of 9/6/08): OscOut;Frequency;Lock-in Amplitude;Lock-in Phase;34410 Volts;D;Lock-in 14 Mag;Bfield;Time
			CreateMAmpAuxData();
			SVAR MAmpLoadWavesList,MAmpFileNameTag
			gloadwavelist = MAmpLoadWavesList;
			gloadfilenametag = MAmpFileNameTag
			CleanupMAmpAuxData()
			break;
		Case 2:
			// Thermometer calibration
			// 4.2K rings
			//gloadwavelist = "B;TA;TB;TC;RD;RC;Time";
			gloadwavelist = "B;RD;R1f;Phi1f;TA;TB;TC;Time"
			gloadfilenametag = "bscan";
			break;
		Case 3:
			// Magnetoresistance measurement -- first data format (timescan)
			MRsweepAux()
			SVAR oldMRsweeplist
			SVAR MRsweepfilenametag;
			gloadwavelist = oldMRsweeplist;
			gloadfilenametag = MRsweepfilenametag;
			CleanupMRsweep()
			break;
		Case 4:
			// Magnetoresistance measurement -- second data format (stepped bscan)
			MRsweepAux()
			SVAR MRsweeplist
			SVAR MRsweepfilenametag;
			NVAR MRSweepNamemode
			gloadwavelist = MRsweeplist;
			gloadfilenametag = MRsweepfilenametag;
			gloadNameMode = MRsweepNamemode
			CleanupMRsweep()
			break;
		Case 5:
			//Reserved for quick loading -- adjust list and filename as necessary
			//gloadwavelist = "B;RD;R1f;Phi1f;TA;TB;TC;Time"
			//gloadwavelist = "B;RC;RD;TA;TB;Vab;TC;Phiab;Time"
			//gloadwavelist="B;R1f;Phi1f;freq;R2f;Phi2f;RD;TB;TC;Time"
			//gloadfilenametag = "bscan";
			
			//gloadwavelist="B;freq;R1f;Phi1f;Vref;Vctrl;TC;RD;TB;TA;Time"
			//gloadfilenametag = "RFbscan";
			
			gloadwavelist="gatetime;freq;R1f;Phi1f;Vref;RD;Time"
			gloadfilenametag = "gatescan";
			break;
		Case 6: //ampscan with 1 lock-in
		Case 7: //ampscan with 2 lock-ins
			if (GLoadmode==6)
				LIAmpscanAux(1)
			else
				LIAmpscanAux(2)
			endif
			SVAR LIampNames,LIampFileName
			NVAR LIampFirstLine,LIampNameMode
			gloadNameMode=LIampNameMode
			gloadwavelist=LIampNames
			gloadFirstLine=LIampFirstLine
			gloadfilenametag=LIampFileName
			cleanupLIampscan()
		break
		Case 8: //Timescan bscans early June 2009
			gloadwavelist="B;R1f;Phi1f;freq;R2f;Phi2f;RD;TB;TC;Time"
			gloadfilenametag = "bscan";
			break
		Case 9: //Timescan bscan early June 8 2009
			gloadwavelist="B;R1f;Phi1f;freq;Vctrl;RD;TB;TC;Time"
			gloadfilenametag = "bscan";
			break
		default:
			// Insert error handling here
			break;
	endswitch;

	Return varsExist;
end

Function CleanupGLoad()

	NVAR GLoadVarExistence;
	if (NVAR_Exists(GLoadVarExistence) == 0)
		printf "Error in GLoad functions"
	elseif (GLoadVarExistence>0)
		GLoadVarExistence-=1
	else
		KillStrings/Z gloadfolderpath, gloadWavelist, gloadfilenametag;
		KillVariables/Z GLoadVarExistence,gloadNameMode,gloadFirstLine
	endif
end

Function CreateMAmpAuxData()

	NVAR MAmpVarExistence;
	if (NVAR_Exists(MAmpVarExistence) == 1)
		MAmpVarExistence+=1
	else
		Variable/G MAmpVarExistence=0
	endif

	//String/G MAmpLoadWavesList = "Voscout;B;F;R1f;Phi1f;Vref;RD;R2f;Time";
	//NOTE: some older files have the rings in the following order
	String/G MAmpLoadWavesList = "Voscout;F;R1f;Phi1f;Vref;RD;R2f;B;Time";
	String/G MAmpFileNameTag="ampscan"

	String/G NewMAmpWavesList = "Vpiezo"; // Could add temperature conversion
	String/G x4avgWaveTag = "Voscout";
	Variable/G oscoutDivider = 100;
	
	Variable/G MAmpAnDoPreProcess = 1; 
	// >0 = do preprocessing during loading. <=0 = don't (useful when you want to cut out bad traces from data)
	String/G MAmpAnalysisList = "Vpiezo_tavg;R1f_tavg";
	String/G MAmpAnalNames = "V;R1f";
	String/G MAmpdriveName = StringFromList(0,MAmpAnalNames);
	String/G MAmpresponseName = StringFromList(1,MAmpAnalNames);
	
	String/G MAmpSampList = "F;Vpiezo;Time";
	String/G MAmpSampNames = "F;V;Tim";
	String/G MAmpNoisyName = StringFromList(0,MAmpSampNames);
	String/G MAmpAnalxWave = StringFromList(1,MAmpSampNames);
	String/G MAmpAnaltWave = StringFromList(2,MAmpSampNames);
	String/G MAmpSampOtherList = RemoveFromList(MAmpNoisyName,MAmpSampNames);
	
	Variable/G MAmpfixBesselLinear=1 //Holds linear term in Bessel fit to 0 unless it equals zero

	Return MAmpVarExistence;
end

Function CleanupMAmpAuxData()
	
	Variable deletedVars=0
	NVAR MAmpVarExistence
	if (NVAR_Exists(MAmpVarExistence))
		if(MAmpVarExistence==0)
			KillStrings/Z MAmpLoadWavesList,NewMAmpWavesList,x4avgWaveTag,MAmpAnalysisList,MAmpAnalNames;
			KillStrings/Z MAmpAnalxWave,MAmpNoisyName,MAmpSampOtherlist,MAmpAnaltWave;
			KillStrings/Z MAmpSampList,MAmpSampNames,MAmpdriveName,MAmpresponseName;
			KillStrings/Z MAmpFileNameTag
			KillVariables/Z MAmpVarExistence,oscoutDivider,MAmpAnDoPreProcess,MAmpfixBesselLinear
		else
			MAmpVarExistence-=1
		endif
	else
		printf "Error with CleanupMAmpAuxData()"
	endif
	Return deletedVars
	
end

Function MRsweepAux()

	NVAR MRsweepVarExistence;
	if (NVAR_Exists(MRsweepVarExistence) == 1)
		MRsweepVarExistence+=1
	else
		Variable/G MRsweepVarExistence=0
	endif

	Variable/G MRsampR = 294.9
	Variable/G MRbridgefactor = -94670
	Variable/G MRsweepVarExistence;
	Variable/G MRsweepNamemode=0

	String/G oldMRsweepList = "B;Vab;Phiab;RC;RD;TC;TA;TB;time";
	String/G MRsweepList = "B;Vab;RD;TC;time";
	String/G MRBname = StringFromList(0,MRsweepList)
	String/G MRVabname = StringFromList(1,MRsweepList)
	String/G MRtimename = StringFromList(4,MRsweepList)
	String/G MRSamptimename = "tim"
	String/G MRsweepfilenametag = "bscan";
	String/G MRdRname = "dR";
	String/G MRdRRname = "dRR"

	String/G MRUCForigSampList = MRBname+";"+MRdRRname+";"+MRtimename
	String/G MRUCFnewSampList = MRBname+";"+MRdRRname+";"+MRSamptimename
	String/G MRUCFSampFolder = "Samples"
	
	Return MRsweepVarExistence
end

Function CleanupMRsweep()
	
	Variable deletedVars=0
	NVAR MRsweepVarExistence
	if (NVAR_Exists(MRsweepVarExistence))
		if(MRsweepVarExistence==0)
			KillVariables/Z MRsampR,MRbridgefactor
			KillVariables/Z MRsweepVarExistence,MRsweepNamemode
			KillStrings/Z oldMRsweepList,MRsweepList,MRBname,MRVabname,MRtimename,MRSamptimename
			KillStrings/Z MRsweepfilenametag,MRdRname,MRdRRname,MRUCForigSampList
			KillStrings/Z MRUCFnewSampList,MRUCFSampFolder
		else
			MRsweepVarExistence-=1
		endif
	else
		printf "Error with CleanupMRsweep()"
	endif
	Return deletedVars
end

Function CreateSBscanAuxData()

	NVAR SBscanVarExistence;
	if (NVAR_Exists(SBscanVarExistence) == 1)
		SBscanVarExistence+=1
	else
		Variable/G SBscanVarExistence=0
	endif
	
	String/G SBscanFileNameTag="bscan"

	//String/G SBscanLoadWaveList = "B;F;R1f;Phi1f;Vref;Vctrl;TC;RD;TB;TA;Time";
	String/G SBscanLoadWaveList = "B;F;R1f;Phi1f;Vref;Vctrl;TC;RD;TB;Time";
	//String/G SBscanLoadWaveList = "B;R1f;Phi1f;F;R2f;Phi2f;TB;TC;RD;Time";
	//String/G SBscanLoadWaveList = "B;F;R1f;Phi1f;Vref;Vctrl;TB;TC;RD;TA;Time";
	//String/G SBscanLoadWaveList = "B;F;R1f;Phi1f;Vref;Vctrl;TB;TC;RD;TA;R2f;Time";
	//String/G SBscanLoadWaveList = "B;R1f;Phi1f;Vref;F;Vctrl;TB;TC;RD;Time";
	String/G SBscanAnalList = "B_tavg;F_tavg";
	String/G SBscanAnalNames = "B;F";
	String/G SBscanxWave4avg = StringFromList(0,SBscanLoadWaveList)

	Return SBscanVarExistence;
end

Function CleanupSBscanAuxData()
	
	Variable deletedVars=0
	NVAR SBscanVarExistence
	if (NVAR_Exists(SBscanVarExistence))
		if(SBscanVarExistence==0)
			KillStrings/Z SBscanLoadWaveList,SBscanAnalList,SBscanAnalNames,SBscanxWave4avg;
			KillStrings/Z SBscanFileNameTag
			KillVariables/Z SBscanVarExistence;
			deletedVars=1
		else
			SBscanVarExistence-=1
		endif
	else
		printf "Error with CleanupSBscanAuxData()"
	endif
	Return deletedVars

end

Function LIAmpscanAux(LInum)
	Variable LInum
	
	NVAR LIampVarExistence
	if (NVAR_Exists(LIampVarExistence))
		LIampVarExistence+=1
	else
		Variable/G LIampVarExistence=0
	endif
	
	//Initial guesses for Bessel function fit
	Make/D/O/N=3 LIampCoef={1e-2,30,0}
	
	Variable/G LIampFirstLine=3
	Variable/G LIampNameMode=0
	//Naming conventions
	String/G LIampFileName="ampscan"
	String/G LIampFitPrefix="besselfit_"
	String/G sampleID="CL15"
	String/G LIampFitFunc="BesselPlusLinear"
	String/G LIampFitTextBox="BesselPlusLinearCoef"
	Switch (LInum)
		Case 1:
			Variable/G LIampGloadmode=6
			String/G LIampNames="Drive;Phi1f;Amp1f"
			break
		default:	//Case 2
			Variable/G LIampGloadmode=7
			String/G LIampNames="Drive;Phi1f;Amp1f;Phi2f;Amp2f"
			break
	endswitch
	String/G LIampDriveName=StringFromList(0,LIampNames)
	String/G LIampAmpName=StringFromList(2,LIampNames)
	
	Return LIampVarExistence
end


Function cleanupLIAmpscan()
	Variable deletedVars=0
	NVAR LIampVarExistence
	if (NVAR_Exists(LIampVarExistence))
		if(LIampVarExistence==0)
			KillStrings/Z LIampNames,LIampFileName,LIampDriveName,LIampAmpName
			KillStrings/Z sampleID,LIampFitFunc,LIampFitTextBox,LIampFitPrefix
			KillVariables/Z LIampGloadmode,LIampFirstLine,LIampVarExistence
			KillVariables/Z LIampNameMode
			KillWaves/Z LIampCoef
			deletedVars=1
		else
			LIampVarExistence-=1
		endif
	else
		printf "Error with cleanupLIAmpscan()"
	endif
	Return deletedVars
end

Function GLoad(folder,filenum,GLoadmode)

	String folder;
	Variable filenum,Gloadmode;
	
	String fileNumStr;
	sprintf fileNumStr,  "%.3d",filenum;

	
	String newFolderName = CreateIgorFolderName(folder,filenum)
	NewDataFolder/O $newFolderName;
	SetDataFolder $newFolderName;
	
	GLoadsetup(GLoadmode);
	SVAR gloadfolderpath;
	SVAR gloadwavelist;
	SVAR gloadfilenametag;
	NVAR gloadNameMode,gloadFirstLine
	Wave gloadOptions;
	Variable nameMode=gloadNameMode
	
	String fileName = gloadfolderpath+folder+"\\"+fileNumStr + " " + gloadfilenametag;
	LoadWave/J/D/W/O/L={0,gloadFirstLine,0,0,0}/N=tempname/K=1 fileName;
	// S_waveNames is a semicolon-separated string list of loaded wave names -- could error check loaded file
	
	Variable totWaves = ItemsInList(gloadwavelist);
	Variable i;
	String tempWaveName,newTempWaveName,tempWaveTag;
	for(i=0;i<totWaves;i+=1)
		tempWaveName = "tempname" + num2str(i);
		Wave oldTempWave = $tempWaveName;
		tempWaveTag = StringFromList(i,gloadwavelist);
		newTempWaveName = CreateWaveName(fileNum,folder,tempWaveTag,gloadnameMode)
		RemoveBlanks(oldTempWave)
		Duplicate/O oldTempWave, $newTempWaveName;
		KillWaves oldTempWave;
	endfor
	
	CleanupGLoad()
	Return nameMode
end

//Search folder parent for all files with names matching formatStr and possessing the
//extension "fileExtension." 
//If optionsFlag is equal to 1, GrepString() is used to match formatStr.
//If optionsFlag is equal to 0, StringMatch() is used to match formatStr.
Function/S findAllFiles(parent, formatStr, fileExtension, optionsFlag)

	String parent, formatStr, fileExtension
	Variable optionsFlag
	
	String matchedFilePathList = "";
	
	Variable i
	NewPath/O/Q tempPath, parent
	String subDirList = IndexedDir(tempPath, -1, 1)
	String tempList = ""
	for(i=0;i<ItemsInList(subDirList);i+=1)
		tempLIst = findAllFiles( StringFromList(i, subDirList), formatStr, fileExtension, optionsFlag)
		matchedFilePathList = ListJoin(matchedFilePathList, tempList)
	endfor
	
	NewPath/O/Q tempPath, parent
	String fileList = IndexedFile(tempPath, -1, fileExtension)
	String tempFilePath, tempFileName
	Variable fileMatchBool
	for(i=0;i<ItemsInList(fileList);i+=1)
		tempFilePath = ""
		tempFileName = StringFromList(i, fileList)
		if (GetBinaryDigit(optionsFlag,0))
		//Bit 0 set: use GrepString()
			if (GrepString(tempFileName, formatStr))
				tempFilePath = parent + ":" + tempFileName
			endif
		else
		//Bit 0 not set: use stringmatch()
			if (stringmatch(tempFileName, formatStr))
				tempFilePath = parent + ":" + tempFileName
			endif
		endif
		matchedFilePathList = ListJoin(matchedFilePathList, tempFilePath)
	endfor
	
	KillPath tempPath
	return matchedFilePathList
end