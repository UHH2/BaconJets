path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoRes_ULtest_2

cd ../src/

./main --run B -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunB_2018_AK4CHS.root --outSuffix _NoJERsf_noL2L3Res_MadGraph_ULtest --Generator madgraph

./main --run D -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunD_2018_AK4CHS.root --outSuffix _NoJERsf_noL2L3Res_MadGraph_ULtest --Generator madgraph


# ./main --run B -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunB_2018_AK4CHS.root --outSuffix _NoJERsf_noL2L3Res_MadGraph_ULtest --Generator madgraph

# ./main --run D -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunD_2018_AK4CHS.root --outSuffix _NoJERsf_noL2L3Res_MadGraph_ULtest --Generator madgraph


cd ../bash_scripts/