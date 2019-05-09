cd ../src
nevents=10000

root -l -b -q CompareTrees_RelVals.C\(\"N\ matched\ jets\",\"Nmatchedrecojet\",false,20,0,20,$nevents\)
#root -l -b -q CompareTrees_RelVals.C\(\"Njet\",\"Njet\",false,50,0,50,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"Njet\",\"Njet\",false,40,0,40,$nevents\)
#root -l -b -q CompareTrees_RelVals.C\(\"Njet\",\"Njet\",false,20,0,20,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"Ngenjet\",\"Ngenjet\",false,25,0,25,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"nvertices\",\"nvertices\",false,70,0,70,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"nGoodvertices\",\"nGoodvertices\",false,70,0,70,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"jet\_\{1\}\ pt,\ GeV\",\"jet1_pt\",true,100,0,1000,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"jet\_\{1\}\ \#eta\",\"jet1_eta\",true,100,-5.2,5.2,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"jet\_\{1\}\ \#phi\",\"jet1_phi\",true,100,-3.14,3.14,$nevents\)

root -l -b -q CompareTrees_RelVals.C\(\"jet\_\{1\}\ charged\ EM\",\"jet1_chEM\",true,20,0,1,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"jet\_\{1\}\ neutral\ EM\",\"jet1_neutEM\",true,20,0,1,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"jet\_\{1\}\ charged\ HAD\",\"jet1_chHAD\",true,20,0,1,$nevents\)
root -l -b -q CompareTrees_RelVals.C\(\"jet\_\{1\}\ neutral\ HAD\",\"jet1_neutHAD\",true,20,0,1,$nevents\)

cd ../bash_scripts/