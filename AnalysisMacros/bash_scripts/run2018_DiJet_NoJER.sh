#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V13h_L1L2L3_L2L3Res_NoJERsmear_NoPUid_Lead2jetsPassJetID/
path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V13h_L1L2L3_L2L3Res_NoJERsmear_NoPUid_Lead2jetsPassJetID_addJERSMbool_fixJetSorting/

cd ../src/

#./main --run ABC -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunABC_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

./main --run ABC -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunABC_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

##  ./main --run A -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunA_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

./main --run A -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunA_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

# ./main --run B -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunB_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

./main --run B -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunB_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

# ./main --run C -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunC_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

./main --run C -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunC_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

# ./main --run D -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunD_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

./main --run D -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunD_2018_AK4CHS.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

cd ../bash_scripts/