for Nev in 1k 10k 50k 100k 500k
#for Nev in 1k
#for Nev in 500k
do
echo "ZLIB1_${Nev}"
#{ time sleep 1 2> sleep.stderr ; } 2> time.txt
#{ /usr/bin/time -v sframe_main BaconJets_DiJet_JERsmearing_ComressionTest_ZLIB1_${Nev}.xml 1 2>LOG_stderr; } 2> LOG_ZLIB1_${Nev}.txt
/usr/bin/time -v sframe_main BaconJets_DiJet_JERsmearing_ComressionTest_ZLIB1_${Nev}.xml >> LOG_ZLIB1_${Nev}.txt
echo ""
echo "ZLIB1_Robin_${Nev}"
/usr/bin/time -v sframe_main BaconJets_DiJet_JERsmearing_ComressionTest_ZLIB1_Robin_${Nev}.xml >> LOG_ZLIB1_Robin_${Nev}
echo ""
echo "LZMA9_${Nev}"
/usr/bin/time -v sframe_main BaconJets_DiJet_JERsmearing_ComressionTest_LZMA9_${Nev}.xml >> LOG_LZMA9_${Nev}
echo ""
echo "LZMA9_10kPerFile_${Nev}"
/usr/bin/time -v sframe_main BaconJets_DiJet_JERsmearing_ComressionTest_LZMA9_10kFile_${Nev}.xml >> LOG_LZMA9_10kFile_${Nev}
echo ""
echo "------------------------------------------------------------------------------------"
done