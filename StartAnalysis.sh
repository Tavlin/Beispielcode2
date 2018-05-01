#! /bin/bash
#
#  Powerweek neutrale Pionen
#  Zur Verfügung gestellt von AG Büsching
#  IKF, Goethe Universität Frankfurt


echo "";



# Ordner fuer Bilder machen, Programm ausfuehren, Ordner 1 hoch wegen github
# mkdir -p Reconstructed
# rm Reconstruction_*
# echo "Starte Reconstruction.C...";
# root -q -l -b Reconstruction.C\+\(\)
# cd ..
# rm -r Reconstructed/
# mv  Beispielcode2/Reconstructed/ Reconstructed/
# cd Beispielcode2
#
# mkdir -p Extraction
# rm Extraction_*
# echo "Starte Extraction.C...";
# root -q -l -b Extraction.C\+\(\)
# cd ..
# rm -r Extraction/
# mv  Beispielcode2/Extraction/ Extraction/
# cd Beispielcode2

mkdir -p P_T_Spectra
rm P_TSpectrumExtraction_*
echo "Starte P_TSpectrumExtraction.C...";
root -q -l -b P_TSpectrumExtraction.C\+\(\)
cd ..
rm -r P_T_Spectra/
mv  Beispielcode2/P_T_Spectra/ P_T_Spectra/
cd Beispielcode2
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
