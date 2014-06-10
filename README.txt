Here is a short (incomplete) documentation overview the package


(1) How to install the L1 Muon Efficiency Tool

export MY_CMSSW_VERSION="CMSSW_6_2_5"
cmsrel $MY_CMSSW_VERSION
cd $MY_CMSSW_VERSION/src

cmsenv

git clone https://github.com/cms-l1-dpg/L1Ntuples.git L1TriggerDPG/L1Ntuples
git clone https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu
git clone https://github.com/battibass/L1MuonTriggerAnalysis.git

scramv1 b -j 9



(2) How to run ntuple production in CRAB

There are already two multicrab.cfg template in:

L1MuonTriggerAnalysis/Running/crab/multicrab_efficiency.cfg
L1MuonTriggerAnalysis/Running/crab/multicrab_rate.cfg

they are setup to run on efficiency or rate samples

Look at the .py cfg it uses and try to customise it for what you 
need to do.
Adapt the multicrab.cfg (including output SEs and directories) 
and run your jobs.

When they are done and if you wrote your jobs in EOS @ CERN,
and you are running from CERN, you can try to use:

L1MuonTriggerAnalysis/Running/crab/getAndMergeRootFiles.py

to automatically retrice and merge the ntuples.
You should adapt the eosBaseDir variable of the script to your 
EOS base directory



(3) How to run the efficiency analysis macro?

cmsenv
cd L1MuonTriggerAnalysis/Efficiency/

EDIT Configuration.C to cope with the correct:
- type of analysis
- GMT quality cut map
- matching windows

root -l

[from root]
.x ../../L1TriggerDPG/L1Ntuples/macros/initL1Analysis.C++
.L GmtEfficiency.C++

goEfficiency("INPUT_FILE.root","OUTPUT_DIR",N_EVENTS<OPTIONAL>)



(4) How to make plot comparisons

You can run an helper macro to superimpose plots from different
studies:

cmsenv
cd L1MuonTriggerAnalysis/Efficiency/

make

./comparePlots.exe FIRST_FILE.root SECOND_FILE.root <ADDITIONAL_FILES.root>


