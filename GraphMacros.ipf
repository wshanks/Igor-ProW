#pragma rtGlobals=1		// Use modern global access method.

Function DoGMacro(macroIndex)

	Variable macroIndex
	
	Switch (macroIndex)
		Case 0:
			Execute "basicXYmarkers()"
			break
		Case 1:
			Execute "CurrentvsBfieldnA()"
			break
		Case 2:
			Execute "CurrentvsBfieldpA()"
			break
		Case 3:
			Execute "FreqNoiseRatio()"
			break
		Case 4:
			Execute "dFreqvsB()"
			break
		Case 5:
			Execute "FreqvsB()"
			break
		Case 6:
			Execute "FreqvsTime()"
			break
		Case 7:
			Execute "FreqvsOscOut()"
			break
		Case 8:
			Execute "AmpScan()"
			break
		Case 9:
			Execute "CurrentvsTemperature()"
			break
		Case 10:
			Execute "FFT_I_pA()"
			break
		Case 11:
			Execute "Poly_Sub()"
			break
		default:
	endswitch
	
end

Proc FreqNoiseRatio() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb=(63719,5624,5624)
	ModifyGraph/Z msize=4
	ModifyGraph/Z fSize=18
	ModifyGraph/Z prescaleExp(bottom)=9
	Label/Z left "\\F'Mathematica1'd\\F'Arial'f\\Bmeasured\\M/\\F'Mathematica1'd\\F'Arial'f\\Bthermal"
	Label/Z bottom "Cantilever Motion Amplitude (\F'Times'nm\\F'Arial')"
	SetAxis/Z left 0,3.0937709
	SetAxis/Z bottom 0,1.5185338e-07
EndMacro

Proc CurrentvsBfieldnA() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb[0]=(63719,5624,5624),rgb[1]=(0,0,65280)
	ModifyGraph/Z msize=2
	ModifyGraph/Z fSize=18
	ModifyGraph/Z prescaleExp(left)=9
	Label/Z left "Current (\\F'Times'nA\\F'Arial')"
	Label/Z bottom "Magnetic field (\\F'Times'T\\F'Arial')"
	ModifyGraph tick=2,mirror=1
EndMacro

Proc CurrentvsBfieldpA() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb[0]=(63719,5624,5624),rgb[1]=(0,0,65280)
	ModifyGraph/Z msize=2
	ModifyGraph/Z fSize=18
	ModifyGraph/Z prescaleExp(left)=12
	Label/Z left "Current (\\F'Times'pA\\F'Arial')"
	Label/Z bottom "Magnetic field (\\F'Times'T\\F'Arial')"
	ModifyGraph tick=2,mirror=1
EndMacro

Proc AutocorrCurrentvsBfieldpA() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode[0]=4
	ModifyGraph/Z marker[0]=8,msize[0]=4
	ModifyGraph/Z mode[1]=0
	ModifyGraph/Z rgb[0]=(65280,21760,0),rgb[1]=(0,0,39168)
	ModifyGraph/Z msize=2
	ModifyGraph/Z fSize=16
	ModifyGraph/Z prescaleExp(left)=24
	ModifyGraph/Z zero(left)=1
	Label/Z left "Current autocorrelation (\\F'Times'pA\S2\M\\F'Arial')"
	Label/Z bottom "Magnetic field lag (\\F'Times'T\\F'Arial')"
	ModifyGraph tick=2,mirror=1
EndMacro

Proc AutocorrCurrentvsBfieldnA() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode[0]=4
	ModifyGraph/Z marker[0]=8,msize[0]=4
	ModifyGraph/Z mode[1]=0
	ModifyGraph/Z rgb[0]=(65280,21760,0),rgb[1]=(0,0,39168)
	ModifyGraph/Z msize=2
	ModifyGraph/Z fSize=16
	ModifyGraph/Z prescaleExp(left)=18
	ModifyGraph/Z lowTrip(left)=0.001
	ModifyGraph/Z zero(left)=1
	Label/Z left "Current autocorrelation (\\F'Times'nA\S2\M\\F'Arial')"
	Label/Z bottom "Magnetic field lag (\\F'Times'T\\F'Arial')"
	ModifyGraph tick=2,mirror=1
EndMacro

Proc PSDCurrentvsBfieldnA() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode[0]=4
	ModifyGraph/Z marker[0]=8,msize[0]=4
	ModifyGraph/Z mode[1]=0
	ModifyGraph/Z rgb[0]=(65280,21760,0),rgb[1]=(0,0,39168)
	ModifyGraph/Z msize=2
	ModifyGraph/Z fSize=16
	ModifyGraph/Z prescaleExp(left)=18
	ModifyGraph/Z lowTrip(left)=0.001
	ModifyGraph/Z zero(left)=1
	Label/Z left "Spectral density (\\F'Times'nA\S2\M T\\F'Arial')"
	Label/Z bottom "\\F'Mathematica1'b\\F'Arial' (\\F'Times'1/T\\F'Arial')"
	ModifyGraph tick=2,mirror=1
EndMacro

Proc basicXYmarkers() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb=(63719,5624,5624)
	ModifyGraph/Z msize=2
	Label/Z left "y"
	Label/Z bottom "x"
EndMacro

Proc dFreqvsB() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb=(63719,5624,5624)
	ModifyGraph/Z msize=2
	ModifyGraph/Z fSize=18
	ModifyGraph/Z prescaleExp(left)=6
	Label/Z left "Frequency Shift (\\F'Mathematica1'm\\F'Times'Hz\\F'Arial')"
	Label/Z bottom "Magnetic field (\\F'Times'T\\F'Arial')"
	ModifyGraph tick=2,mirror=1
EndMacro

Proc FreqvsB() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb=(63719,5624,5624)
	ModifyGraph/Z msize=2
	ModifyGraph/Z fSize=18
	Label/Z left "Frequency (\\F'Times'Hz\\F'Arial')"
	Label/Z bottom "Magnetic field (\\F'Times'T\\F'Arial')"
	ModifyGraph tick=2,mirror=1
EndMacro

Proc FreqvsTime() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z margin(top)=58
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z lSize[1]=2
	ModifyGraph/Z rgb[0]=(63719,5624,5624),rgb[1]=(0,52224,52224)
	ModifyGraph/Z msize=2
	ModifyGraph/Z fSize=16
	ModifyGraph/Z highTrip(bottom)=100000
	Label/Z left "Frequency (\\F'Times'Hz\\F'Arial')"
	Label/Z bottom "Time (\\F'Times's\\F'Arial')"
	ModifyGraph tick=2,mirror=1
EndMacro

Proc FreqvsOscOut() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb[0]=(63719,5624,5624),rgb[1]=(63719,56743,5624),rgb[2]=(19570,63719,5624)
	ModifyGraph/Z rgb[3]=(5624,63719,42803)
	ModifyGraph/Z msize=2
	Label/Z left "Frequency (\\F'Times'Hz\\F'Arial')"
	Label/Z bottom "V\\Bosc out\\M (\\F'Times'mV\\F'Arial')"
EndMacro

Proc AmpScan() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode[0]=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb[0]=(63719,5624,5624),rgb[1]=(0,0,65535)
	ModifyGraph/Z msize=2
	Label/Z left "Cantilever 1st Harmonic Response (\\F'Times'V\\F'Arial')"
	Label/Z bottom "Voscout (\\F'Times'mV\\F'Arial')"
	SetAxis/Z left 0,0.76560003
	SetAxis/Z bottom 0,275
EndMacro

Proc CurrentvsTemperature() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode[3]=3,mode[4]=3
	ModifyGraph/Z marker[3]=19,marker[4]=19
	ModifyGraph/Z lSize[0]=2,lSize[1]=2,lSize[2]=2
	ModifyGraph/Z rgb[0]=(0,52224,52224),rgb[1]=(13056,0,26112),rgb[2]=(0,39168,0),rgb[3]=(0,52224,26368)
	ModifyGraph/Z rgb[4]=(11520,256,65280)
	ModifyGraph/Z msize[3]=3,msize[4]=3
	ModifyGraph/Z log(left)=1
	ModifyGraph/Z fSize=20
	ModifyGraph/Z prescaleExp(left)=9
	Label/Z left "Current (\\F'Times'nA\\F'Arial')"
	Label/Z bottom "Temperature (\\F'Times'K\\F'Arial')"
EndMacro

Proc FFT_I_pA() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z msize=2
	ModifyGraph/Z tick=2
	ModifyGraph/Z mirror=1
	ModifyGraph/Z font="Arial"
	ModifyGraph/Z fSize=16
	ModifyGraph/Z fStyle=1
	ModifyGraph/Z prescaleExp(left)=12
	Label/Z left "FFT Amp (\F'Times'pA\\F'Arial')"
	Label/Z bottom "\F'Times'1/T"
EndMacro

Proc Poly_Sub() : GraphStyle

	ModifyGraph/Z mode[0]=3,mode[1]=0,marker[0]=19,msize[0]=2
	ModifyGraph/Z tick=2,mirror=1,fStyle=1,fSize=16,font="Arial";
	ModifyGraph/Z rgb[0]=(0,0,65535),rgb[1] = (65535,0,0);
EndMacro