#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1__RelVals_105rel2/2017/WithoutJEC_2/
#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1__RelVals_105rel2/AK4CHS/WithJEC_2017V32_2018V7_2016V11/
#path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1__RelVals_105rel2/AK4CHS_minJet10/WithJEC_2017V32_2018V7_2016V11/
path_Gl=/nfs/dust/cms/user/karavdia/JERC/RunII_102X_v1__RelVals_104/AK4CHS_minJet10/WithJEC_2017V32_2018V7_2016V11/
for thr in 15
do
cd ../src

./main --run F -F -CPRecGenCompare --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_refHCAL_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_newHCAL_AK4CHS.root --outSuffix _2017HCALCheck_oldMC_newDATA --dname RelVals

./main --run F -F -CPRecGenCompare --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_refHCAL_AK4CHS_NoPU.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_newHCAL_AK4CHS_NoPU.root --outSuffix _2017HCALCheck_oldMC_newDATA_NoPU --dname RelVals

./main --run F -F -CPRecGenCompare --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2016_refHCAL_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2016_newHCAL_AK4CHS.root --outSuffix _2016HCALCheck_oldMC_newDATA --dname RelVals


./main --run F -F -CPRecGenCompare --inputMC $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2017_refHCAL_AK4CHS.root --input $path_Gl/uhh2.AnalysisModuleRunner.MC.QCD_Flat_pythia8_2018_newHCAL_AK4CHS.root --outSuffix _2018HCALCheck_old2017MC_new2018DATA --dname RelVals


done
cd ../bash_scripts/