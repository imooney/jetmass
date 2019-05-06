#!/bin/csh
#  A Slightly ugly way to make it as simple as possible to run different analyses,
#  And to make it easier to change variables and keep the book keeping simple
#  Output names and locations are generated by the script, and correspond to the above variables              
# first make sure program is updated and exists 

set ExecPath = `pwd`
#set execute = './bin/analysis'
#set base = /nfs/rhi/STAR/Data/ppHT/picoDst_
set arg = '' 

if ($1 == 'data') then
    make bin/ppdata || exit
    set execute = './bin/ppdata'
    set base = /nfs/rhi/STAR/Data/ppJP2Run12/sum
     # Create the folder name for output
    if ($3 == 'ge') then
	set outFile = data
    endif    
    if ($3 == 'py') then
	set outFile = data
    endif
endif
if ($1 == 'sim') then
    make bin/ppsim || exit
    set execute = './bin/ppsim'
    set base = /nfs/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_pt 
    # Create the folder name for output
    if ($3 == 'ge') then
        set outFile = sim/ge
    endif
    if ($3 == 'py') then
	set outFile = sim/py
    endif
endif
if ($1 == 'matching') then
    make bin/matching || exit
    set execute = './bin/matching'
    set base = /nfs/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_pt
    # Create the folder name for output
    set outFile = matching
endif
if ($1 == 'systematics') then
    make bin/systematics || exit
    set execute = './bin/systematics'
    set base = /nfs/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_pt
    # Create the folder name for output
    set outFile = systematics
endif
if ($1 == 'closure') then
    make bin/closure || exit
    set execute = './bin/closure'
    set base = /nfs/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_pt
    # Create the folder name for output
    set outFile = closure
endif

if ($2 == 'ch') then
    set full_or_ch = "ch"
    echo "Running with only charged jets!"
endif
if ($2 != 'ch') then
    set full_or_ch = "full"
    echo "Running over Ch+Ne jets!"
endif

if ($3 == 'py') then
    set ge_or_py = "py"
    echo "Running Pythia6!"
endif
if ($3 == 'ge') then
    set ge_or_py = "ge"
    echo "Running Pythia6+Geant!"
endif

# Arguments                                                                                                                                                                   
if ( $# < "1") then
        echo 'Error: illegal number of parameters'
        exit
endif

# Create the folder name for output
#set outFile = stock
# Make the directories since they may not exist...                                                                                                                             
if ( ! -d out/${outFile} ) then
mkdir -p out/${outFile}
endif

if ( ! -d log/${outFile} ) then
mkdir -p log/${outFile}
endif

#echo ${base}
# Now Submit jobs for each data file                                                                                                                                           
foreach input ( ${base}* )

# Create the output file base name                                                                                                                                             
set OutBase = `basename $input | sed 's/.root//g'`
set uscore = "_"
set OutBase = "$OutBase$uscore$full_or_ch$uscore$ge_or_py"
    
# Make the output names and path                                                                                                                                               
set outLocation = out/${outFile}/
set outName = ${OutBase}.root

# Input files                                                                                                                                                                  
set Files = ${input}

# Logfiles. Thanks cshell for this "elegant" syntax to split err and out                                                                                                       
set LogFile     = log/${outFile}/${OutBase}.log
set ErrFile     = log/${outFile}/${OutBase}.err

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile

set arg = "$outLocation $outName $full_or_ch $ge_or_py $Files"

echo "now submitting this script: "
echo qsub -V -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

qsub -V -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg   
#make python script later: jet_analysis/submit scripts start with pbs
#add back in a second: -q erhiq
end
