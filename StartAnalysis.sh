#! /bin/bash
#
#  Powerweek neutrale Pionen
#  Zur Verfügung gestellt von AG Büsching
#  IKF, Goethe Universität Frankfurt


echo "";

mkdir -p Simulation
mkdir -p Extraction
mkdir -p P_T_Spectra
echo "Starte Reconstruction.C...";
root -q -l -b Reconstruction.C\+\(\)

# echo "Starte Extraction.C...";
# root -q -l -b Extraction.C\+\(\)

# echo "Starte P_TSpectrumExtraction.C...";
# root -q -l -b P_TSpectrumExtraction.C\+\(\)
# Hier wird ein C++ Macro in root aufgerufen
# Das + bedeutet, dass wir den code kompilieren wollen (.C+),
#         die \ werden nur verwendet, da einige symbole in bash
#         "escaped" werden müssen




##################### To Be Modified
#if [ $DoCorrection = 1 ] ; then
#  echo "";
#  echo "Starte Correction.C...";
#  option = "z.B."\,0\,8\,15\,"xD" # irgendeinSetting
#  root -q -l -b  Correction.C\+\(option\)
#fi
#.....
