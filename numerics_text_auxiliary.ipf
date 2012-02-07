#pragma rtGlobals=1		// Use modern global access method.

Function BitString2Number(bitString)
	// Bit string is a binary string.  Output is the decimal equivalent (as variable, not string)
	String bitString // bitString should be a string consisting of only "1" and "0"

	Variable totBinDigits=strlen(bitString)
	Variable decimalValue=0
	
	String tempChar
	Variable i
	for(i=0;i<totBinDigits;i+=1)
		tempChar=bitString[totBinDigits-1-i]
		decimalValue+=str2num(tempChar)*2^i
	endfor

	Return decimalValue
end

//Removes char2clear from targetString if it is the leading and/or trailing character
//Returns the trimmed version of targetString
Function/S ClearLeadingTrailingChars( targetString, char2clear)

	String targetString, char2clear
	
	if (strlen(char2clear) != 1)
		return targetString
	endif
	
	if ((stringmatch(targetString[strlen(targetString)-1], char2clear) && strlen(targetString)>1))
		targetString = targetString[0,strlen(targetString)-2]
	endif
	
	if (stringmatch(targetString[0], char2clear[0]))
		if(strlen(targetString) > 1)
			targetString = targetString[1,strlen(targetString)-1]
		else
			return ""
		endif
	endif
	
	return targetString
end

Function/S coefTextEntry(coefNameList,coefWave,sigmaWave)

	String coefNameList
	Wave coefWave,sigmaWave
	
	String formattedCoefText=""
	String tempString
	Variable i
	for(i=0;i<itemsinlist(coefNameList);i+=1)
		formattedCoefText+=StringfromList(i,coefNameList)+":\t"
		sprintf tempString, "%.3g",coefWave[i]
		formattedCoefText+=tempString
		formattedCoefText+=" ± "
		sprintf tempString, "%.3g",sigmaWave[i]
		formattedCoefText+=tempString
		formattedCoefText+="\r"
	endfor
	Return formattedCoefText
end

Function/S fitAnnotation(fitFunction,coefWave,sigmaWave)
	String fitFunction
	Wave coefWave,sigmaWave
	
	String coefNameList
	coefNameList=findCoefList(fitFunction)
	
	String annotationText=fitFunction+" coefficients:\r"
	annotationText+=coefTextEntry(coefNameList,coefWave,sigmaWave)
	
	Return annotationText
end

Function GetBinaryDigit(number,digit)
	// Number is the number to be cast into binary
	// digit is the digit to be returned with 0 being the least significant binary digit
	Variable number,digit
	
	String tempBinStr
	sprintf tempBinStr,"%b",number
	Variable totDigits = strlen(tempBinStr)
	String returnDigitStr = tempBinStr[totDigits-1-digit]
	Return str2num(returnDigitStr)
end

//Join Igor-style lists list1 and list2 which use listSepStr as the separating
//string.  listSepStr is an optional parameter.  The default value is ";"
//This function removes any leading or trailing instances of listSepStr from 
//list1 and list2 so that the joined list does not have any empty entries.
Function/S ListJoin(list1, list2,  [listSepStr])
	
	String list1, list2, listSepStr
	
	if ( ParamIsDefault(listSepStr))
		listSepStr = ";"
	endif
	
	list1 = ClearLeadingTrailingChars(list1, listSepStr)
	list2 = ClearLeadingTrailingChars(list2, listSepStr)
	
	if (strlen(list1) > 0)
		return list1 + ";" + list2
	else
		return list2
	endif
end

Function protoFitFunc(w,x)
	Wave w
	Variable x
	
	printf "protoFitFunc called -- procedure error using FuncRef with a FitFunc!\r"
end

Function/S SigmaWaveName(filenum,folder,namemode)
	Variable filenum,namemode
	String folder
	
	Return "sigma_"+CreateWaveName(filenum,folder,"",namemode)
end