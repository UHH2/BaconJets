path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2017/Fall17V32_L1L2L3_NoRes_NoJER_SingleJetTrg_UL17statTest

cd ../src/

#./main --run C -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunC_2017_AK4CHS.root --outSuffix _NoJERsf_noL2L3Res_MadGraph_ULtest --Generator madgraph

./main --run F -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.DATA.DATA_RunF_2017_AK4CHS.root --outSuffix _NoJERsf_noL2L3Res_MadGraph_ULtest --Generator madgraph

cd ../bash_scripts/