#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V7_L1L2L3_NoJERsmear_JetByJet_FixCS/
#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V7_L1L2L3_NoJERsmear_JetByJet_FixCS_newSamples
#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V7_L1L2L3_NoJERsmear_JetByJet_FixCS_newSamples_RECO10GeV/
path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoJERsmear_JetByJet_FixCS_newSamples_RECO10GeV_dR02/

#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V7_L1L2L3_NoJERsmear_JetByJet_FixCS_newSamples_RECO10GeV_dR01_dR02/

cd ../src

./main --run C -F -CPPUjets --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat2018_pythia8_2018_AK4CHS.root --outSuffix _MadGraph_noJERsf_dR02  --Generator madgraph

#./main --run C -F -CPRecGenCompare --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat2018_pythia8_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --outSuffix _Flat2018Pythi8MCakaold_MadGraphDATAakanew_noJERsf_dR01

#./main --run C -F -CPRecGenCompare --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS_dR02.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS_dR01.root --outSuffix _dR02akaold_dR01akanew_noJERsf --Generator madgraph

cd ../bash_scripts/