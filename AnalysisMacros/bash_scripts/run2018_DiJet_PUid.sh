#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoJERsmear_DiJet_FixCS_newSamples_RECO10GeV_dR02
#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoJERsmear_DiJet_FixCS_newSamples_RECO10GeV_dR02_addDATApuid/
#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoJERsmear_DiJet_FixCS_newSamples_RECO10GeV_dR02_addDATApuid_addMCpartonPUid/
#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoRes_NoJERsmear_PUidLOOSE3rd_DiJet_FixCS_newSamples_RECO10GeV_dR02_addDATApuid_addMCpartonPUid/
#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoRes_V4ABCJERsmear_PUidTight3rd_DiJet_FixCS_newSamples_RECO10GeV_dR02_addDATApuid_addMCpartonPUid_fixPUidvsPFs/
path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoRes_V4ABCJERsmear_PUidTightALL_DiJet_FixCS_newSamples_RECO10GeV_dR02_addDATApuid_addMCpartonPUid_fixPUidvsPFs/
cd ../src/

./main --run A -F -CPPUjets --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunA_2018_AK4CHS.root --outSuffix _MadGraph_V4JERsf_noL2L3Res_PUidTightALL_PUhists  --Generator madgraph

./main --run A -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunA_2018_AK4CHS.root --outSuffix _V4JERsf_noL2L3Res_PUidTightALL_MadGraph --Generator madgraph

./main --run A -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunA_2018_AK4CHS.root --outSuffix _V4JERsf_noL2L3Res_PUidTightALL_MadGraph --Generator madgraph

#./main --run A -F -CPPUjets --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunA_2018_AK4CHS.root --outSuffix _MadGraph_noJERsf_noL2L3Res_PUhists  --Generator madgraph

# ./main --run A -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --outSuffix _noJERsf_Madgraph_Assym --Generator madgraph

# ./main --run A -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat2018_pythia8_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat2018_pythia8_2018_AK4CHS.root --outSuffix _noJERsf_Pythia8_Assym --Generator madgraph

#./main --run A -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat2018_pythia8_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunA_2018_AK4CHS.root --outSuffix _noJERsf_noL2L3Res_PUidTight3rdJet_Pythia8 



# ./main --run B -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunB_2018_AK4CHS.root --outSuffix _noJERsf_MadGraph --Generator madgraph

# ./main --run B -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat2018_pythia8_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunB_2018_AK4CHS.root --outSuffix _noJERsf_Pythia8 

#./main --run B -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunB_2018_AK4CHS.root --outSuffix _noJERsf_MadGraph --Generator madgraph


# ./main --run C -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunC_2018_AK4CHS.root --outSuffix _noJERsf_MadGraph --Generator madgraph

# ./main --run C -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat2018_pythia8_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunC_2018_AK4CHS.root --outSuffix _noJERsf_Pythia8 

#./main --run C -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunC_2018_AK4CHS.root --outSuffix _noJERsf_MadGraph --Generator madgraph

# ./main --run D -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunD_2018_AK4CHS.root --outSuffix _noJERsf_MadGraph --Generator madgraph

# ./main --run D -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat2018_pythia8_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunD_2018_AK4CHS.root --outSuffix _noJERsf_Pythia8 

#./main --run D -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunD_2018_AK4CHS.root --outSuffix _noJERsf_MadGraph --Generator madgraph

cd ../bash_scripts/