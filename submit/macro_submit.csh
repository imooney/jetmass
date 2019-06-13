#!/bin/csh

set ExecPath="~/jetmass/macros"
set shPath="~/jetmass/submit"

set uscore = "_"
set execute= './test' #'root tohist.C'
#set baseprim = ~/jetmass/out/QA/pAu_2015_200_MB_
#if ($1 == 'MB') then
#set type = MB
#endif
#if ($1 == 'HT') then
#set type = HT
#endif 

#set base = "$baseprim$type$uscore"
set basepAu = ~/jetmass/out/QA/pAu_2015_200_HT_
set basepAuMB = ~/jetmass/out/QA/pAu_2015_200_MB_
set basepp = ~/jetmass/out/QA/sum

#PROBLEM! EACH TYPE HAS A DIFFERENT NUMBER OF INPUT FILES. MIGHT NEED TO SCRAP THIS FOR NOW. COME BACK TO IT LATER.

if ( ! -d pAuhists/${type}) then
mkdir -p pAuhists/${type}
endif

if ( ! -d pAulogs/${type}) then
mkdir -p pAulogs/${type}
endif

foreach input ( ${base}* )

set OutBase = `basename $input | sed 's/.root//g'`
set uscore = "_"
set OutBase = "$OutBase$uscore$outFile"

set outLocation = pAuhists/${type}/
set outName = ${OutBase}.root

set Files = ${input}

set LogFile = pAulogs/${outFile}/${OutBase}.log
set ErrFile = pAulogs/${outFile}/${OutBase}.err

set arg = "$outLocation $outName $Files" 

echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${shPath}/qwrap.sh ${ExecPath} $execute $arg
#qsub -V -q erhiq -l mem=4gb -- ${shPath}/qwrap.sh ${ExecPath} root -l ${execute}
