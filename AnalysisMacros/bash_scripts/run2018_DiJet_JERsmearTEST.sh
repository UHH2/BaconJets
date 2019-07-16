path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1/2018/Autumn18V8_L1L2L3_NoRes_JERsmearVSnosmear
cd ../src/

./main --run A -F -CPPUjets --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS_NoJERsmear.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS_V4JERsmear.root --outSuffix _MadGraph_noJERsmearMC_V4JERsfDATA_noL2L3Res_PUhists  --Generator madgraph

./main --run A -F -FP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS_NoJERsmear.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS_V4JERsmear.root --outSuffix _MadGraph_noJERsmearMC_V4JERsfDATA_noL2L3Res --Generator madgraph

./main --run A -F -aAP --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS_NoJERsmear.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_HT_2018_AK4CHS_V4JERsmear.root --outSuffix _MadGraph_noJERsmearMC_V4JERsfDATA_noL2L3Res --Generator madgraph


cd ../bash_scripts/