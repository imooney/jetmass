#! /usr/bin/env csh
echo a
set ExecPath = `pwd`# PYTHIA8
echo c
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_0 -o Log/rivet_PY8_0_max30.log -e Log/rivet_PY8_0_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30 --pwd PYTHIA8_5pthat10_1MEvents_200GeV.hepmc -H pythia8_5pthat10_1Mevents_200gev_max30.yoda
echo d
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_1 -o Log/rivet_PY8_1_max30.log -e Log/rivet_PY8_1_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_10PTHAT15 --pwd PYTHIA8_10pthat15_1MEvents_200GeV.hepmc -H pythia8_10pthat15_1Mevents_200gev_max30.yoda
echo e
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_2 -o Log/rivet_PY8_2_max30.log -e Log/rivet_PY8_2_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_15PTHAT20 --pwd PYTHIA8_15pthat20_1MEvents_200GeV.hepmc -H pythia8_15pthat20_1Mevents_200gev_max30.yoda
echo f
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_3 -o Log/rivet_PY8_3_max30.log -e Log/rivet_PY8_3_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_20PTHAT25 --pwd PYTHIA8_20pthat25_1MEvents_200GeV.hepmc -H pythia8_20pthat25_1Mevents_200gev_max30.yoda
echo g
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_4 -o Log/rivet_PY8_4_max30.log -e Log/rivet_PY8_4_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_25PTHAT30 --pwd PYTHIA8_25pthat30_1MEvents_200GeV.hepmc -H pythia8_25pthat30_1Mevents_200gev_max30.yoda
echo h
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_5 -o Log/rivet_PY8_5_max30.log -e Log/rivet_PY8_5_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_30PTHAT35 --pwd PYTHIA8_30pthat35_1MEvents_200GeV.hepmc -H pythia8_30pthat35_1Mevents_200gev_max30.yoda
echo i
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_6 -o Log/rivet_PY8_6_max30.log -e Log/rivet_PY8_6_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_35PTHAT40 --pwd PYTHIA8_35pthat40_1MEvents_200GeV.hepmc -H pythia8_35pthat40_1Mevents_200gev_max30.yoda
echo j
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_7 -o Log/rivet_PY8_7_max30.log -e Log/rivet_PY8_7_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_40PTHAT45 --pwd PYTHIA8_40pthat45_1MEvents_200GeV.hepmc -H pythia8_40pthat45_1Mevents_200gev_max30.yoda
echo k
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_8 -o Log/rivet_PY8_8_max30.log -e Log/rivet_PY8_8_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_45PTHAT50 --pwd PYTHIA8_45pthat50_1MEvents_200GeV.hepmc -H pythia8_45pthat50_1Mevents_200gev_max30.yoda
echo l
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_9 -o Log/rivet_PY8_9_max30.log -e Log/rivet_PY8_9_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_50PTHAT60 --pwd PYTHIA8_50pthat60_1MEvents_200GeV.hepmc -H pythia8_50pthat60_1Mevents_200gev_max30.yoda
echo m
qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_PY8_max30_10 -o Log/rivet_PY8_10_max30.log -e Log/rivet_PY8_10_max30.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_max30_60PTHAT80 --pwd PYTHIA8_60pthat80_1MEvents_200GeV.hepmc -H pythia8_60pthat80_1Mevents_200gev_max30.yoda
echo n
# HERWIG7
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_0 -o Log/rivet_HW7_0.log -e Log/rivet_HW7_0.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP --pwd /nfs/rhi/STAR/Data/herwig/rhic_3_4-S1.hepmc -H herwig7_3pthat4_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_1 -o Log/rivet_HW7_1.log -e Log/rivet_HW7_1.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_4PTHAT5 --pwd /nfs/rhi/STAR/Data/herwig/rhic_4_5-S1.hepmc -H herwig7_4pthat5_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_2 -o Log/rivet_HW7_2.log -e Log/rivet_HW7_2.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_5PTHAT7 --pwd /nfs/rhi/STAR/Data/herwig/rhic_5_7-S1.hepmc -H herwig7_5pthat7_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_3 -o Log/rivet_HW7_3.log -e Log/rivet_HW7_3.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_7PTHAT9 --pwd /nfs/rhi/STAR/Data/herwig/rhic_7_9-S1.hepmc -H herwig7_7pthat9_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_4 -o Log/rivet_HW7_4.log -e Log/rivet_HW7_4.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_9PTHAT11 --pwd /nfs/rhi/STAR/Data/herwig/rhic_9_11-S1.hepmc -H herwig7_9pthat11_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_5 -o Log/rivet_HW7_5.log -e Log/rivet_HW7_5.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_11PTHAT15 --pwd /nfs/rhi/STAR/Data/herwig/rhic_11_15-S1.hepmc -H herwig7_11pthat15_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_6 -o Log/rivet_HW7_6.log -e Log/rivet_HW7_6.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_15PTHAT25 --pwd /nfs/rhi/STAR/Data/herwig/rhic_15_25-S1.hepmc -H herwig7_15pthat25_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_7 -o Log/rivet_HW7_7.log -e Log/rivet_HW7_7.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_25PTHAT35 --pwd /nfs/rhi/STAR/Data/herwig/rhic_25_35-S1.hepmc -H herwig7_25pthat35_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_8 -o Log/rivet_HW7_8.log -e Log/rivet_HW7_8.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_35PTHAT45 --pwd /nfs/rhi/STAR/Data/herwig/rhic_35_45-S1.hepmc -H herwig7_35pthat45_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_9 -o Log/rivet_HW7_9.log -e Log/rivet_HW7_9.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_45PTHAT55 --pwd /nfs/rhi/STAR/Data/herwig/rhic_45_55-S1.hepmc -H herwig7_45pthat55_multiplepthat_500kevents_200gev.yoda
#qsub -V -p 10 -q  erhiq -l mem=4gb -W umask=0022 -N RIVET_HW7_10 -o Log/rivet_HW7_10.log -e Log/rivet_HW7_10.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_55PTHAT65 --pwd /nfs/rhi/STAR/Data/herwig/rhic_55_65-S1.hepmc -H herwig7_55pthat65_multiplepthat_500kevents_200gev.yoda
