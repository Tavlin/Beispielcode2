#! /bin/bash
#
#  Powerweek neutrale Pionen
#  Zur Verfügung gestellt von AG Büsching
#  IKF, Goethe Universität Frankfurt


echo "";
echo "Starte Pi0Simulation.C...";
mkdir -p Simulation
root -q -l -b Reconstruction.C\+\(\)
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
