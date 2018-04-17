#include "CommenHeader.h"

void P_TSpectrumExtraction(TString AddName = ""){

  //Open and read file
  TFile* HistoWOBackground_file = new TFile("HistoWOBackground_file.root", "READ");

  if ( HistoWOBackground_file->IsOpen() ) printf("HistoWOBackground_file opened successfully\n");

  if(HistoWOBackground_file->IsZombie()){
    std::cout << "ERROR: HistoWOBackground_file not found" << std::endl;
    return;
  }


  // Erstellen der Canvas
  TCanvas *cP_TSpectrum = new TCanvas("cP_TSpectrum", "",1080,1080);
  SetCanvasStandardSettings(cP_TSpectrum);


  // Erstellen der Latex-Objekte
  TLatex *lP_TSpectrum = new TLatex();
  SetLatexSettings(lP_TSpectrum);


  // Erstelle 1D Histogram fuer Ratio
  TH1D* hP_TSpectrum = new TH1D("hP_TSpectrum","#it{p}_{T} spectrum",27,0.,10.);
  // 27 bins, plus ein overflow bin in welches die high pt messungen eingehen
  // sollen
  SetHistoStandardSettings(hP_TSpectrum);
  TH1D* hSignal[28];

  TH1D* hMinvSpectra[28];
  TF1* fGausFit[28];

  for(int i = 0; i < 28; i++){

    hMinvSpectra[i] = new TH1D(Form("hMinvSpectra[%d]",i),Form("#it{m}_{inv} spectra [%d]",i),150,0.,0.3);
    SetHistoStandardSettings(hMinvSpectra[i]);
    gDirectory->GetObject(Form("hMinvSpectra[%d]",i),hSignal[i]);
    hMinvSpectra[i]->Fit("gaus(0)");
  }




}
