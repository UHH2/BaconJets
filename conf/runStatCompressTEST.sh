for Nev in 1k 10k 50k 100k 250k 500k
#for Nev in 500k
do
echo "ZLIB1_${Nev}"
./timeSFrame.sh BaconJets_DiJet_JERsmearing_ComressionTest_ZLIB1_${Nev}.xml
echo ""
#echo "ZLIB1_Robin_${Nev}"
#./timeSFrame.sh BaconJets_DiJet_JERsmearing_ComressionTest_ZLIB1_Robin_${Nev}.xml
#echo ""
echo "LZMA9_${Nev}"
./timeSFrame.sh BaconJets_DiJet_JERsmearing_ComressionTest_LZMA9_${Nev}.xml
echo ""
echo "LZMA9_10kPerFile_${Nev}"
./timeSFrame.sh BaconJets_DiJet_JERsmearing_ComressionTest_LZMA9_10kFile_${Nev}.xml
echo ""
echo "------------------------------------------------------------------------------------"
done