#! /usr/bin/env csh #comment to avoid weird line endings
set TYPE = "P8" #options: P8, H7, P6
set ExecPath = `pwd` #files will be written here
set P6Path = /nfs/rhi/STAR/Data/RhicJEWEL #PYTHIA6
set P8Path = /nfs/rhi/STAR/Data/IsaacsPy8 #RhicPythia !  NEED TO CHECK FOR CONSISTENCY!!! #PYTHIA8 //!~~~~~~~~~also need to change back and forth going from had to no had
set H7Path = /nfs/rhi/STAR/Data/RhicHerwig #HERWIG7
set whichdir = woDecay #comment to avoid weird line endings
set Decay = "off"#"noHad" #comment to avoid the weird line endings
set Filename = "PYTHIA8_" #comment to avoid weird line endings
set Fileend = "_1MEvents_200GeV.hepmc" #comment to avoid weird line endings
if ($TYPE == "P6") then #PY6+DECAY
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_PY6_0 -o Log/rivet_PY6_0.log -e Log/rivet_PY6_0.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING --pwd --ignore-beams ${P6Path}/pp200GeV_5pthat80_1Mevt.hepmc -H pp200GeV_5pthat80_p6decayed.yoda #comment to avoid weird line endings
endif #comment to avoid weird line endings
if ($TYPE == "P8") then #PYTHIA8
    if ($Decay == "on") then #comment to avoid weird line endings
	set whichdir = "wDecay" #wDecay #comment to avoid weird line endings
	set Filename = "" #RHIC200_${whichdir}_ #comment to avoid weird line endings
	set Fileend = _1MEvents_200GeV_decays_on.hepmc #-S76456.hepmc #comment to avoid weird line endings
    else if ($Decay == "off") then #comment to avoid weird line endings
	set whichdir = "woDecay" #woDecay #comment to avoid weird line endings
	set Filename = PYTHIA8_ #comment to avoid weird line endings
    	set Fileend = _1MEvents_200GeV_decays_off.hepmc #comment to avoid weird line endings
    else #comment to avoid weird line endings
	set whichdir = "noHad" #noHad #comment to avoid weird line endings
	set Filename = PYTHIA8_ #comment to avoid weird line endings
	set Fileend = _1MEvents_partonLevel_200GeV.hepmc #comment to avoid weird line endings
    endif #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_0 -o Log/rivet_PY8_0_decays_$Decay.log -e Log/rivet_PY8_0_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING --pwd ${P8Path}/${whichdir}/${Filename}5pthat10${Fileend} -H pythia8_5pthat10_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_1 -o Log/rivet_PY8_1_decays_$Decay.log -e Log/rivet_PY8_1_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_10PTHAT15 --pwd ${P8Path}/${whichdir}/${Filename}10pthat15${Fileend} -H pythia8_10pthat15_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_2 -o Log/rivet_PY8_2_decays_$Decay.log -e Log/rivet_PY8_2_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_15PTHAT20 --pwd ${P8Path}/${whichdir}/${Filename}15pthat20${Fileend} -H pythia8_15pthat20_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_3 -o Log/rivet_PY8_3_decays_$Decay.log -e Log/rivet_PY8_3_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_20PTHAT25 --pwd ${P8Path}/${whichdir}/${Filename}20pthat25${Fileend} -H pythia8_20pthat25_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_4 -o Log/rivet_PY8_4_decays_$Decay.log -e Log/rivet_PY8_4_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_25PTHAT30 --pwd ${P8Path}/${whichdir}/${Filename}25pthat30${Fileend} -H pythia8_25pthat30_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_5 -o Log/rivet_PY8_5_decays_$Decay.log -e Log/rivet_PY8_5_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_30PTHAT35 --pwd ${P8Path}/${whichdir}/${Filename}30pthat35${Fileend} -H pythia8_30pthat35_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_6 -o Log/rivet_PY8_6_decays_$Decay.log -e Log/rivet_PY8_6_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_35PTHAT40 --pwd ${P8Path}/${whichdir}/${Filename}35pthat40${Fileend} -H pythia8_35pthat40_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_7 -o Log/rivet_PY8_7_decays_$Decay.log -e Log/rivet_PY8_7_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_40PTHAT45 --pwd ${P8Path}/${whichdir}/${Filename}40pthat45${Fileend} -H pythia8_40pthat45_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_8 -o Log/rivet_PY8_8_decays_$Decay.log -e Log/rivet_PY8_8_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_45PTHAT50 --pwd ${P8Path}/${whichdir}/${Filename}45pthat50${Fileend} -H pythia8_45pthat50_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_9 -o Log/rivet_PY8_9_decays_$Decay.log -e Log/rivet_PY8_9_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_50PTHAT60 --pwd ${P8Path}/${whichdir}/${Filename}50pthat60${Fileend} -H pythia8_50pthat60_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10 -q erhiq   -l mem=4gb -W umask=0022 -N RIVET_PY8_10 -o Log/rivet_PY8_10_decays_$Decay.log -e Log/rivet_PY8_10_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_60PTHAT80 --pwd ${P8Path}/${whichdir}/${Filename}60pthat80${Fileend} -H pythia8_60pthat80_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
endif #comment to avoid weird line endings
if ($TYPE == "H7") then #HERWIG7
    if ($Decay == "on") then #comment to avoid weird line endings
	set whichdir = wDecay #comment to avoid weird line endings
	set Filename = RHIC200_${whichdir}_ #comment to avoid weird line endings
        set Fileend = -S76456.hepmc #comment to avoid weird line endings
    else if ($Decay == "off") then #comment to avoid weird line endings
	set whichdir = woDecay #comment to avoid weird line endings
	set Filename = RHIC200_${whichdir}_ #comment to avoid weird line endings
    	set Fileend = -S76456.hepmc #comment to avoid weird line endings
    else #comment to avoid weird line endings
	set whichdir = noHad #comment to avoid weird line endings
	set Filename = NoHadro_rhic_ #comment to avoid weird line endings
	set Fileend = -S876.hepmc #comment to avoid weird line endings
    endif #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_0 -o Log/rivet_HW7_0_decays_$Decay.log -e Log/rivet_HW7_0_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER --pwd ${H7Path}/${whichdir}/${Filename}5_10${Fileend} -H herwig7_5pthat10_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_1 -o Log/rivet_HW7_1_decays_$Decay.log -e Log/rivet_HW7_1_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_10PTHAT15 --pwd ${H7Path}/${whichdir}/${Filename}10_15${Fileend} -H herwig7_10pthat15_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_2 -o Log/rivet_HW7_2_decays_$Decay.log -e Log/rivet_HW7_2_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_15PTHAT20 --pwd ${H7Path}/${whichdir}/${Filename}15_20${Fileend} -H herwig7_15pthat20_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_3 -o Log/rivet_HW7_3_decays_$Decay.log -e Log/rivet_HW7_3_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_20PTHAT25 --pwd ${H7Path}/${whichdir}/${Filename}20_25${Fileend} -H herwig7_20pthat25_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_4 -o Log/rivet_HW7_4_decays_$Decay.log -e Log/rivet_HW7_4_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_25PTHAT30 --pwd ${H7Path}/${whichdir}/${Filename}25_30${Fileend} -H herwig7_25pthat30_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_5 -o Log/rivet_HW7_5_decays_$Decay.log -e Log/rivet_HW7_5_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_30PTHAT35 --pwd ${H7Path}/${whichdir}/${Filename}30_35${Fileend} -H herwig7_30pthat35_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_6 -o Log/rivet_HW7_6_decays_$Decay.log -e Log/rivet_HW7_6_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_35PTHAT40 --pwd ${H7Path}/${whichdir}/${Filename}35_40${Fileend} -H herwig7_35pthat40_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_7 -o Log/rivet_HW7_7_decays_$Decay.log -e Log/rivet_HW7_7_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_40PTHAT45 --pwd ${H7Path}/${whichdir}/${Filename}40_45${Fileend} -H herwig7_40pthat45_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_8 -o Log/rivet_HW7_8_decays_$Decay.log -e Log/rivet_HW7_8_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_45PTHAT50 --pwd ${H7Path}/${whichdir}/${Filename}45_50${Fileend} -H herwig7_45pthat50_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_9 -o Log/rivet_HW7_9_decays_$Decay.log -e Log/rivet_HW7_9_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_50PTHAT60 --pwd ${H7Path}/${whichdir}/${Filename}50_60${Fileend} -H herwig7_50pthat60_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_10 -o Log/rivet_HW7_10_decays_$Decay.log -e Log/rivet_HW7_10_decays_$Decay.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_HER_60PTHAT80 --pwd ${H7Path}/${whichdir}/${Filename}60_80${Fileend} -H herwig7_60pthat80_1MEvents_200gev_decays_$Decay.yoda #comment to avoid weird line endings
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_1 -o Log/rivet_HW7_1.log -e Log/rivet_HW7_1.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_4PTHAT5 --pwd ${H7Path}/${whichdir}/rhic_4_5-S1.hepmc -H herwig7_4pthat5_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_2 -o Log/rivet_HW7_2.log -e Log/rivet_HW7_2.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_5PTHAT7 --pwd ${H7Path}/${whichdir}/rhic_5_7-S1.hepmc -H herwig7_5pthat7_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_3 -o Log/rivet_HW7_3.log -e Log/rivet_HW7_3.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_7PTHAT9 --pwd ${H7Path}/${whichdir}/rhic_7_9-S1.hepmc -H herwig7_7pthat9_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_4 -o Log/rivet_HW7_4.log -e Log/rivet_HW7_4.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_9PTHAT11 --pwd ${H7Path}/${whichdir}/rhic_9_11-S1.hepmc -H herwig7_9pthat11_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_5 -o Log/rivet_HW7_5.log -e Log/rivet_HW7_5.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_11PTHAT15 --pwd ${H7Path}/${whichdir}/rhic_11_15-S1.hepmc -H herwig7_11pthat15_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_6 -o Log/rivet_HW7_6.log -e Log/rivet_HW7_6.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_15PTHAT25 --pwd ${H7Path}/${whichdir}/rhic_15_25-S1.hepmc -H herwig7_15pthat25_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_7 -o Log/rivet_HW7_7.log -e Log/rivet_HW7_7.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_25PTHAT35 --pwd ${H7Path}/${whichdir}/rhic_25_35-S1.hepmc -H herwig7_25pthat35_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_8 -o Log/rivet_HW7_8.log -e Log/rivet_HW7_8.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_35PTHAT45 --pwd ${H7Path}/${whichdir}/rhic_35_45-S1.hepmc -H herwig7_35pthat45_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_9 -o Log/rivet_HW7_9.log -e Log/rivet_HW7_9.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_45PTHAT55 --pwd ${H7Path}/${whichdir}/rhic_45_55-S1.hepmc -H herwig7_45pthat55_multiplepthat_500kevents_200gev.yoda
    #qsub -V -p 10    -l mem=4gb -W umask=0022 -N RIVET_HW7_10 -o Log/rivet_HW7_10.log -e Log/rivet_HW7_10.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_SPLITTING_DUP_55PTHAT65 --pwd ${H7Path}/${whichdir}/rhic_55_65-S1.hepmc -H herwig7_55pthat65_multiplepthat_500kevents_200gev.yoda
endif #comment to avoid weird line endings
