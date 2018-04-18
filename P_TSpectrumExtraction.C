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

  TGaxis::SetMaxDigits(3);

  // Erstellen der Latex-Objekte
  TLatex *laP_TSpectrum = new TLatex();
  SetLatexSettings(laP_TSpectrum);

  // Erstellen der Legende
  TLegend *leP_TSpectrum = new TLegend(0.4,0.75,0.9,0.95);
  SetLegendSettigns(leP_TSpectrum);
  leP_TSpectrum->SetTextSize(0.04);

  // Definieren der Bins fuers pt-Spectrum, vorgegeben durch selektions loop in
  // Extraction.C!

  const Int_t nbins_pt = 29;
  Float_t xbins_pt[nbins_pt+1];

  for (size_t j = 0; j < 26; j++) {
    xbins_pt[j] = j*0.25;
  }
  xbins_pt[26] = 6.5;
  xbins_pt[27] = 7;
  xbins_pt[28] = 8;
  xbins_pt[29] = 10;

  // Erstelle 1D Histogram fuer Ratio
  TH1D* hP_TSpectrum = new TH1D("hP_TSpectrum","#it{p}_{T} spectrum",nbins_pt,xbins_pt);
  hP_TSpectrum->SetXTitle("#it{p}_T (GeV/#it{c})");

  // 27 bins, plus ein overflow bin in welches die high pt messungen eingehen
  // sollen
  SetHistoStandardSettings(hP_TSpectrum);
  TH1D* hSignal[29];

  TH1D* hMinvSpectra[29];

  // Deklarieren der Fit Funktionen
  TF1* fGausFit[29];
  TF1* fGausFit_dummy[29];
  Double_t mean[29], sigma[29];
  Double_t integral_value[29];
  Double_t int_error[29];

  cplayceholder->cd();
  for(int i = 0; i < 29; i++){

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
    hMinvSpectra[i]->Fit(Form("fGausFit[%d]",i),"MR0Q","", 0, 0.15);
    fGausFit_dummy[i]->SetParameter(0,fGausFit[i]->GetParameter(0));
    fGausFit_dummy[i]->SetParameter(1,fGausFit[i]->GetParameter(1));
    fGausFit_dummy[i]->SetParameter(2,fGausFit[i]->GetParameter(2));

    hMinvSpectra[i]->Draw();
    fGausFit_dummy[i]->Draw("same");
    cplayceholder->SaveAs(Form("P_T_Spectra/P_TSpectra(%d).png",i));

    mean[i] = fGausFit[i]->GetParameter(1);
    sigma[i] = fGausFit[i]->GetParameter(2);

    // 6 Sigma Integration (might lower, since out of bounds sometimes!!!!!)
    integral_value[i] = hMinvSpectra[i]->IntegralAndError(hMinvSpectra[i]->FindBin(mean[i]-6*sigma[i]),hMinvSpectra[i]->FindBin(mean[i]+6*sigma[i]),int_error[i],"");

    // Fill the bins of the p_T spectrum
    hP_TSpectrum->SetBinContent(i,integral_value[i]);
    hMinvSpectra[i]->Delete();
  }
  // Draw the actual pt spectrum
  cP_TSpectrum->cd();
  gStyle->SetOptStat(0);
  cP_TSpectrum->SetLeftMargin(0.125);
  cP_TSpectrum->SetTopMargin(0.05);
  cP_TSpectrum->SetBottomMargin(0.15);
  cP_TSpectrum->SetRightMargin(0.05);
  cP_TSpectrum->SetTickx();
  cP_TSpectrum->SetTicky();
  cP_TSpectrum->SetLogy(0);
  cP_TSpectrum->SetLogx(0);

  leP_TSpectrum->AddEntry(hP_TSpectrum,"#splitline{#pi^{0} #it{p}_{T} spectrum}{without background}");
  hP_TSpectrum->SetError(int_error);
  hP_TSpectrum->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hP_TSpectrum->SetMarkerColor(kViolet+7);
  hP_TSpectrum->SetLineColor(kViolet+7);
  hP_TSpectrum->Draw("");
  leP_TSpectrum->Draw("SAME");


  cP_TSpectrum->SaveAs(Form("P_T_Spectra/P_TSpectra.png"));



}
