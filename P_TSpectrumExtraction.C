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


  TCanvas *cplayceholder = new TCanvas("cplayceholder", "",1080,1080);
  SetCanvasStandardSettings(cplayceholder);


  // Erstellen der Latex-Objekte
  TLatex *lP_TSpectrum = new TLatex();
  SetLatexSettings(lP_TSpectrum);


  // Erstelle 1D Histogram fuer Ratio
  TH1D* hP_TSpectrum = new TH1D("hP_TSpectrum","#it{p}_{T} spectrum",27,0.,8.);
  hP_TSpectrum->SetXTitle("#it{p}_T (GeV/#it{c})");

  // 27 bins, plus ein overflow bin in welches die high pt messungen eingehen
  // sollen
  SetHistoStandardSettings(hP_TSpectrum);
  TH1D* hSignal[28];

  TH1D* hMinvSpectra[28];

  // Deklarieren der Fit Funktionen
  TF1* fGausFit[28];
  TF1* fGausFit_dummy[28];
  Double_t mean[28], sigma[28];
  Double_t integral_value[28];
  Double_t int_error[28];

  cplayceholder->cd();
  for(int i = 0; i < 28; i++){

    // definieren und auslesen der Histos
    hMinvSpectra[i] = new TH1D(Form("hMinvSpectra[%d]",i),Form("#it{m}_{inv} spectra [%d]",i),150,0.,0.3);
    SetHistoStandardSettings(hMinvSpectra[i]);
    //gDirectory->GetObject(Form("hMinvSpectra[%d]",i),hSignal[i+1]);
    gDirectory->GetObject(Form("hSignal[%d]",i+1),hMinvSpectra[i]);

    // Definieren der Fits und fitten
    fGausFit[i] = new TF1(Form("fGausFit[%d]",i),"gaus", 0, 0.15);
    fGausFit_dummy[i] = new TF1(Form("fGausFit_dummy[%d]",i),"gaus", 0, 0.3);
    fGausFit[i]->SetParLimits(1,0.1,0.15);
    fGausFit[i]->SetParLimits(0,0.,10e6);
    hMinvSpectra[i]->Fit(Form("fGausFit[%d]",i),"MR0","", 0, 0.15);
    fGausFit_dummy[i]->SetParameter(0,fGausFit[i]->GetParameter(0));
    fGausFit_dummy[i]->SetParameter(1,fGausFit[i]->GetParameter(1));
    fGausFit_dummy[i]->SetParameter(2,fGausFit[i]->GetParameter(2));

    hMinvSpectra[i]->Draw();
    fGausFit_dummy[i]->Draw("same");
    cplayceholder->SaveAs(Form("P_T_Spectra/P_TSpectra(%d).png",i));

    mean[i] = fGausFit[i]->GetParameter(1);
    sigma[i] = fGausFit[i]->GetParameter(2);

    // 6 Sigma Integration
    integral_value[i] = hMinvSpectra[i]->IntegralAndError(hMinvSpectra[i]->FindBin(mean[i]-6*sigma[i]),hMinvSpectra[i]->FindBin(mean[i]+6*sigma[i]),int_error[i],"");

    // Fill the bins of the p_T spectrum
    hP_TSpectrum->SetBinContent(i,integral_value[i]);
  }
  cP_TSpectrum->cd();
  hP_TSpectrum->SetError(int_error);
  hP_TSpectrum->Draw();
  cP_TSpectrum->SetLeftMargin(0.2);
  gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
  cP_TSpectrum->SetTopMargin(0.025);
  cP_TSpectrum->SetBottomMargin(0.15);
  cP_TSpectrum->SetRightMargin(0.05);
  cP_TSpectrum->SetTickx();
  cP_TSpectrum->SetTicky();
  cP_TSpectrum->SetLogy(0);
  cP_TSpectrum->SetLogx(0);
  cP_TSpectrum->SaveAs(Form("P_T_Spectra/P_TSpectra.png"));



}
