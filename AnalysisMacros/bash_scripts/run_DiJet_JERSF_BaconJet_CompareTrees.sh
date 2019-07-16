cd ../src
nevents=-1
root -l -b -q CompareTrees_DiffIDs.C\(\"alpha\",\"alpha\",true,100,0,1,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"rel_r\",\"rel_r\",true,100,0,2,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"rho\",\"rho\",true,60,0,60,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"jet3_eta\",\"jet3_eta\",true,50,0,5,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"Njet\",\"Njet\",false,25,0,25,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"jet3_ptRaw\",\"jet3_ptRaw\",true,50,0,500,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"jet3_pt\",\"jet3_pt\",true,50,0,500,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"jet2_ptRaw\",\"jet2_ptRaw\",true,50,0,500,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"jet2_pt\",\"jet2_pt\",true,50,0,500,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"jet1_ptRaw\",\"jet1_ptRaw\",true,50,0,500,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"jet1_pt\",\"jet1_pt\",true,50,0,500,$nevents\)
root -l -b -q CompareTrees_DiffIDs.C\(\"pt_ave\",\"pt_ave\",true,50,0,500,$nevents\)

cd ../bash_scripts/