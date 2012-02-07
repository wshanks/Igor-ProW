#pragma rtGlobals=1		// Use modern global access method.

Function AddTextToGraph(windowName,textboxName,textString,noDupBit)
	//Add text to the top graph.  Add to textboxName if that textbox already exists.
	//Otherwise create a textbox with that name.
	String windowName,textboxName,textString
	Variable noDupBit 
	//if noDupBit>=1, textString is not added to textbox if that string is already
	//present in textbox

	if (strlen(windowName) ==  0)
		windowName = WinName(0, 1)	// Name of top graph or table
	endif
	
	//Format string for execute command below
	textString=ReplaceString("\r",textString,"\\r")
	textString="\""+textString+"\""
	
	String currentTextBoxes = AnnotationList(windowName)
	String textboxCmd
	Variable checkForTextBox = WhichListItem(textboxName,currentTextBoxes)
	if (checkForTextBox < 0)
		textboxCmd="TextBox/W="+windowName+"/N="+textboxName+"/A=RB/F=0 "+textString
		printf "%s\r",textBoxCmd
		Execute/Q textboxCmd
	else
		String currentText=StringByKey("TEXT",AnnotationInfo(windowName,textboxName))
		if ((noDupBit<1)||(noDupBit>=1&&strsearch(currentText,textString,0)<0))
			textboxCmd="AppendText/W="+windowName+"/N="+textboxName+" "+textString
			printf "%s\r",textBoxCmd
			Execute/Q textboxCmd
		endif
	endif
end

Function AppendToGraphNoDup(yWave,xWave,windowName,flags)
	//Appends yWave vs xWave to window named windowName ("" = top graph)
	//unless that yWave is already plotted vs that xWave on the graph.
	//If xWave is the same as yWave, yWave is appended vs calculated 
	//(unless it is already on the graph)
	//Any flags desired for the AppendToGraph operation should be entered in the
	//flags string.  For no flags use the "" string.
	
	Wave yWave,xWave
	String windowName,flags
	
	if (strlen(windowName) ==  0)
		windowName = WinName(0, 1)	// Name of top graph or table
	endif
	
	//Determine is the append operation is for a normal trace or a contour trace
	Variable traceBit
	if (wavematch(yWave,xWave))
		traceBit=0 // y vs calculated
	else
		traceBit=1 // y vs x trace
	endif
	
	String ListOfTraces = TraceNameList(windowName,";",1)
	
	Variable i,yMatch,xMatch
	for (i=0;i<ItemsInList(ListOfTraces);i+=1)
		Wave tempTraceYWave=TraceNameToWaveRef(windowName,StringFromList(i,ListOfTraces))
		yMatch=wavematch(yWave,tempTraceYWave)
		
		Wave tempTraceXwave=XWaveRefFromTrace(windowName,StringFromList(i,ListOfTraces))
		if (traceBit)
			xMatch=wavematch(xWave,tempTraceXWave)
		else
			if (waveexists(tempTraceXwave))
				//yWave is plotted vs a wave (for TraceBit==0 we are looking for 
				//a wave plotted vs calculated)
				xMatch=0 
			else
				//yWave is plotted vs calculated
				xMatch=1
			endif
		endif

		if (yMatch&&xMatch)
			Return 0
		endif
	endfor
	
	String appendOperation = "AppendToGraph"
	if (strlen(flags)>0)
		appendOperation+=flags
	endif
	appendOperation+=" "+GetWavesDataFolder(yWave,2)
	if (traceBit==1)
		appendOperation+=" vs "+GetWavesDataFolder(xWave,2)
	endif
	printf "%s\r",appendOperation
	Execute/Q appendOperation
	Return 1
end

Function AppendTraceSections(xwave,ywave,secLength)

	Wave xwave,ywave;
	Variable secLength;
	
	Variable lastpoint = numpnts(ywave) - 1;
	Variable endsec;
	
	Variable i;
	for(i=0;(i-1)*secLength-1<lastpoint;i+=1)
		if ((i+1)*secLength-1>lastpoint)
			endsec = lastpoint;
		else
			endsec = (i+1)*secLength - 1;
		endif
		AppendtoGraph ywave[i*secLength,endsec] vs xwave[i*secLength,endsec];
	endfor
end

// Appends entire ywave vs xwave to graph but does so as separate traces each of size secLength
Function AppendTraceSections2(xwave,ywave,secLength,msize1,mode1,marker1,lsize1)

	// 0 = lines between points, 3 = markers, 4 = markers and lines
	// 19 = marker dots

	Wave xwave,ywave;
	Variable secLength,msize1,mode1,marker1,lsize1;
	
	Variable lastpoint = numpnts(ywave) - 1;
	Variable endsec;
	String temptrace;
	
	Variable i;
	for(i=0;(i-1)*secLength-1<lastpoint;i+=1)
		if ((i+1)*secLength-1>lastpoint)
			endsec = lastpoint;
		else
			endsec = (i+1)*secLength - 1;
		endif
		AppendtoGraph ywave[i*secLength,endsec] vs xwave[i*secLength,endsec];
		temptrace = StringFromList(ItemsInList(TraceNameList("",";",1))-1,TraceNameList("",";",1));
		ModifyGraph msize ($temptrace) = msize1, mode ($temptrace) = mode1, marker ($temptrace) = marker1, lsize ($temptrace) = lsize1;
	endfor
end

Function basicGMac()

	Modifygraph mode=4,msize=2,marker=19
	RainbowifyTraces(.9)
end

Function CheckGraphforTrace(graphName,traceName)

	String graphName,traceName;
	
	String traceList = TraceNameList(graphName,";",1);
	
	String dummyStr="temp";
	Variable i;
	for(i=0;stringmatch(dummyStr,"")==0;i+=1)
		dummyStr = StringFromList(i,traceList,";");
		if (stringmatch(dummyStr,traceName))
			Return 1;
		endif
	endfor
	
	Return 0;
End

Function/S CsrTraceName(csrID,graphNameStr)

	Variable csrID //0==A,1==B
	String graphNameStr //"" for top graph
	
	if (strlen(graphNameStr) ==  0)
		graphNameStr = WinName(0, 1)	// Name of top graph or table
	endif

	if (csrID==0)
		Wave csrTraceWave=CsrWaveRef(A,graphNameStr)
	else
		Wave csrTraceWave=CsrWaveRef(B,graphNameStr)
	endif
	
	String tracesList=TraceNameList(graphNameStr,";",3)
	Variable i
	for(i=0;i<ItemsInList(tracesList);i+=1)
		Wave tempTraceWave=TraceNameToWaveRef(graphNameStr,StringFromList(i,tracesList))
		if(waveMatch(csrTraceWave,tempTraceWave))
			Return StringFromList(i,tracesList)
		endif
	endfor
	
	Return ""
end

Function hslconvert(hue,saturation,lightness)

	Variable hue,lightness,saturation;
	Make/O/N=(1,1,3) imagewave;
	imagewave[0][0][0] = hue;
	imagewave[0][0][1] = saturation;
	imagewave[0][0][2] = lightness;

	ImageTransform hsl2rgb imagewave;
	
	Wave M_HSL2RGB;
	Make/O/N=(3) rgbWave;
	rgbWave[0] = floor(M_HSL2RGB[0][0][0]);
	rgbWave[1]= floor(M_HSL2RGB[0][0][1]);
	rgbWave[2]= floor(M_HSL2RGB[0][0][2]);
	
end

//Does one of three things based on the value of plotOpt:
//plotOpt=0: does nothing
//plotOpt=1: creates a new graph named plotWindowName and displays the wave
//	with name Wave2PlotName
//plotOpt=2: appends wave with Wave2PlotName to graph named plotWindowName. 
//	If the graph does not exist, it is created.
Function OptionalPlot(Wave2PlotName,plotWindowName,plotOpt)

	String Wave2PlotName,plotWindowName
	Variable plotOpt
	
	Switch(plotOpt)
		Case 2:
			if (ItemsInList(WinList(plotWindowName,";","WIN:1"))>0)
				AppendToGraph/W=$plotWindowName $Wave2PlotName
				break
			endif
		Case 1:
			Display/N=$plotWindowName $Wave2PlotName;DoGMacro(10)
		default://Case 0
			break
	Endswitch
end

Function rainbowifyTraces(huefactor)

	// huefactor = scales successive hues so successive traces don't blend together
	// usually huefactor of ~2 works
	// Common graph mode settings: marker = 19,mode = 4 (lines and markers),
	// mode = 3 (markers),msize = 3
	// ModifyGraph mode=4,marker=19,msize=2;
	// ModifyGraph mode=3,marker=19,msize=2;

	Variable huefactor;
	
	String tracelist = TraceNameList("",";",3);
	Variable totTraces = ItemsInList(tracelist);
	Variable i,hue,r,g,b,frac
	String tempTrace;
	for (i=0;i<totTraces;i+=1)
		tempTrace = StringFromList(i,tracelist);
		if(huefactor>0)
			frac=i/totTraces
		else
			frac=1-i/totTraces
		endif
		hue = floor((abs(huefactor)*frac) * 65536);
		hue = mod(hue,65536);
		hslconvert(hue,61680,34672);
		Wave rgbWave;
		r = rgbWave[0];g = rgbWave[1];b = rgbWave[2];
		ModifyGraph rgb ($tempTrace) = (r,g,b);
	endfor
	KillWaves rgbWave,imagewave,M_HSL2RGB;
end

Function rainbowifyTraces2(graphType,huefactor)

	// huefactor = scales successive hues so successive traces don't blend together
	// usually huefactor of ~2 works
	// Common graph mode settings: marker = 19,mode = 4 (lines and markers),
	// mode = 3 (markers),msize = 3
	// ModifyGraph mode=4,marker=19,msize=2;
	// ModifyGraph mode=3,marker=19,msize=2;

	Variable graphType; // 0 for waves plotted vs calculated, 1 for x-y plots
	Variable huefactor;
	
	String tracelist = TraceNameList("",";",graphType);
	Variable totTraces = ItemsInList(tracelist);
	Variable i,hue,r,g,b;
	String tempTrace;
	for (i=0;i<totTraces;i+=1)
		tempTrace = StringFromList(i,tracelist);
		hue = floor((huefactor*i / totTraces) * 65536);
		hue = mod(hue,65536);
		hslconvert(hue,61680,34672);
		Wave rgbWave;
		r = rgbWave[0];g = rgbWave[1];b = rgbWave[2];
		ModifyGraph rgb ($tempTrace) = (r,g,b);
	endfor
	KillWaves rgbWave,imagewave,M_HSL2RGB;
end

Function rainbowifyTraces3(graphType)

	// Common graph mode settings: marker = 19,mode = 4 (lines and markers),
	// mode = 3 (markers),msize = 3
	// ModifyGraph mode=4,marker=19,msize=2;
	// ModifyGraph mode=3,marker=19,msize=2;

	Variable graphType; // 0 for waves plotted vs calculated, 1 for x-y plots
	
	String tracelist = TraceNameList("",";",graphType);
	Variable totTraces = ItemsInList(tracelist);
	Variable i,hue,r,g,b;
	String tempTrace;
	for (i=0;i<totTraces;i+=1)
		tempTrace = StringFromList(i,tracelist);
		hue = floor((.9*i / totTraces) * 65536);
		hslconvert(hue,61680,34672);
		Wave rgbWave;
		r = rgbWave[0];g = rgbWave[1];b = rgbWave[2];
		ModifyGraph rgb ($tempTrace) = (r,g,b);
	endfor
	KillWaves rgbWave,imagewave,M_HSL2RGB;
end

