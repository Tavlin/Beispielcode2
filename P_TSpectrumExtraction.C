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

  // Deklarieren der Fit Funktionen
  TF1* fGausFit[28];

  for(int i = 0; i < 28; i++){

    // definieren und auslesen der Histos
    hMinvSpectra[i] = new TH1D(Form("hMinvSpectra[%d]",i),Form("#it{m}_{inv} spectra [%d]",i),150,0.,0.3);
    SetHistoStandardSettings(hMinvSpectra[i]);
    //gDirectory->GetObject(Form("hMinvSpectra[%d]",i),hSignal[i+1]);
    gDirectory->GetObject(Form("hSignal[%d]",i+1),hMinvSpectra[i]);

    // Definieren der Fits und fitten
    fGausFit[i] = new TF1(Form("fGausFit[%d]",i),"gaus", 0, 0.3);
    fGausFit[i]->SetParLimits(1,0.1,0.15);
    fGausFit[i]->SetParLimits(0,0.,10e6);
    hMinvSpectra[i]->Fit(Form("fGausFit[%d]",i),"MR","", 0, 0.3);

    hMinvSpectra[i]->Draw();
    cP_TSpectrum->SaveAs(Form("P_T_Spectra/P_TSpectra(%d).png",i));
  }




}
