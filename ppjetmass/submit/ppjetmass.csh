#!/bin/csh
#  A Slightly ugly way to make it as simple as possible to run different analyses,                                                                                            
#  And to make it easier to change variables and keep the book keeping simple                                                                                                                                                                             
#  Output names and locations are generated by the script, and correspond to the above variables                                                                              
# first make sure program is updated and exists 

set ExecPath = `pwd`
set execute = './bin/analysis'
set base = /nfs/rhi/STAR/Data/ppHT/picoDst_
set arg = '' 

make clean && make bin/analysis || exit

# Arguments                                                                                                                                                                   

if ( $# != "0") then
        echo 'Error: illegal number of parameters'
        exit
endif

# Create the folder name for output                                                                                                                                            
set outFile = stock
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

set arg = "$outLocation $outName $Files"

qsub -q erhiq -V -l mem=2GB -o $LogFile -e $ErrFile -N test -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
#make python script later: jet_analysis/submit scripts start with pbs
end
