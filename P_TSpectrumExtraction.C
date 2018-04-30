#include "CommenHeader.h"

void P_TSpectrumExtraction(TString AddName = ""){
  InitStartUp();
  //Open and read file
  TFile* HistoWOBackground_file = new TFile("HistoWOBackground_file.root", "READ");

  if ( HistoWOBackground_file->IsOpen() ) printf("HistoWOBackground_file opened successfully\n");

  if(HistoWOBackground_file->IsZombie()){
    std::cout << "ERROR: HistoWOBackground_file not found" << std::endl;
    return;
  }


  // Erstellen der Canvas
  TCanvas *cP_TSpectrum = new TCanvas("cP_TSpectrum", "",1080,1080);


  TCanvas *cplaceholder[66];

  TGaxis::SetMaxDigits(3);

  // Erstellen der Latex-Objekte
  TLatex *laP_TSpectrum = new TLatex();
  SetLatexSettings(laP_TSpectrum);

  TLatex *ltSignal_pT_projection_clone = new TLatex();
  SetLatexSettings(ltSignal_pT_projection_clone);

  // Erstellen der Legende
  TLegend *leP_TSpectrum = new TLegend(0.5,0.75,0.9,0.95);
  SetLegendSettigns(leP_TSpectrum);
  leP_TSpectrum->SetTextSize(0.04);

  // Definieren der Bins fuers pt-Spectrum, vorgegeben durch selektions loop in
  // Extraction.C!

const Int_t kMaxHit;
  // const Int_t nbins_pt;
  // Float_t xbins_pt[nbins_pt+1];

  TH1D* hP_TSpectrum = new TH1D("hP_TSpectrum","#it{p}_{T} spectrum",nbins_pt,xbins_pt);
  SetHistoStandardSettings(hP_TSpectrum);

  TH1D* hP_TSpectrum_fit = new TH1D("hP_TSpectrum_fit","#it{p}_{T} spectrum",nbins_pt,xbins_pt);
  SetHistoStandardSettings(hP_TSpectrum_fit);

  TH1D* hSignal[66];

  TH1D* hMinvSpectra[66];

  // Deklarieren der Fit Funktionen
  TF1* fGausFit[66];
  TF1* fGausFit_dummy[66];
  Double_t mean[66], sigma[66];
  Double_t integral_value[66];
  Double_t int_error[66];
  Double_t integral_value_fit[66];
  Double_t ymin, ymax;
  for(int i = 0; i < 66; i++){

    // definieren und auslesen der Histos
    hMinvSpectra[i] = new TH1D(Form("hMinvSpectra[%d]",i),Form("#it{m}_{inv} spectra [%d]",i),150,0.,0.3);
    SetHistoStandardSettings(hMinvSpectra[i]);

    cplaceholder[i] = new TCanvas(Form("cplaceholder[%d]",i), "",1080,1080);
    SetCanvasStandardSettings(cplaceholder[i]);

    cplaceholder[i]->cd();

    gDirectory->GetObject(Form("hSignal[%d]",i+1),hMinvSpectra[i]);

    // Definieren der Fits und fitten
    fGausFit[i] = new TF1(Form("fGausFit[%d]",i),"gaus", 0, 0.15);
    fGausFit_dummy[i] = new TF1(Form("fGausFit_dummy[%d]",i),"gaus", 0, 0.3);
    fGausFit[i]->SetParLimits(1,0.1,0.15);
    fGausFit[i]->SetParLimits(0,0.,10e6);
    TFitResultPtr r = hMinvSpectra[i]->Fit(Form("fGausFit[%d]",i),"SIMNRE","", 0, 0.15);
    fGausFit_dummy[i]->SetParameter(0,fGausFit[i]->GetParameter(0));
    fGausFit_dummy[i]->SetParameter(1,fGausFit[i]->GetParameter(1));
    fGausFit_dummy[i]->SetParameter(2,fGausFit[i]->GetParameter(2));

    // hMinvSpectra[i]->Draw();
    // fGausFit_dummy[i]->Draw("same");
    // cplaceholder->SaveAs(Form("P_T_Spectra/P_TSpectra(%d).png",i));

    mean[i] = fGausFit[i]->GetParameter(1);
    sigma[i] = fGausFit[i]->GetParameter(2);

    // 3 Sigma Integration (might lower, since out of bounds sometimes!!!!!)
    integral_value[i] =
    hMinvSpectra[i]->IntegralAndError(hMinvSpectra[i]->FindBin(mean[i]-3*sigma[i]),
    hMinvSpectra[i]->FindBin(mean[i]+3*sigma[i]),int_error[i],"");

    integral_value_fit[i] = fGausFit_dummy[i]->Integral(mean[i]-3*sigma[i],
    mean[i]+3*sigma[i])*150./0.3;

    // Fill the bins of the p_T spectrum
    hP_TSpectrum->SetBinContent(i+1,integral_value[i]/(hP_TSpectrum->GetBinWidth(i+1)));
    hP_TSpectrum->SetBinError(i+1,int_error[i]/(hP_TSpectrum->GetBinWidth(i+1)));
    // cout << "Bin Content = "  << hP_TSpectrum->GetBinContent(i+1) << endl;
    // cout << "Integral value = " << integral_value[i] << endl;
    // cout << "Integral value fit = " << integral_value_fit[i] << endl;
    // cout << "Fehler = " << int_error[i] << endl;
    // cout << "Fehler Histo = " << hP_TSpectrum->GetBinError(i+1) << endl;

    // Float_t ymax = cplaceholder[i]->GetUymax();
    // Float_t ymin = cplaceholder[i]->GetUymin();
    // Float_t xmax = mean[i]+6*sigma[i];
    // Float_t xmin = mean[i]-6*sigma[i];

    const Double_t* GParams, *GMatrixArray;
    GParams = r->GetParams();
    GMatrixArray = r->GetCovarianceMatrix().GetMatrixArray();
    hP_TSpectrum_fit->SetBinContent(i+1,integral_value_fit[i]/(hP_TSpectrum->GetBinWidth(i+1)));
    hP_TSpectrum_fit->SetBinError(i+1,fGausFit_dummy[i]->IntegralError(mean[i]-3*sigma[i],
    mean[i]+3*sigma[i],GParams, GMatrixArray, 1.e-2)*150./0.3);


    hMinvSpectra[i]->Draw("EP");
    fGausFit_dummy[i]->Draw("SAMEP");
    ltSignal_pT_projection_clone->DrawLatexNDC(0.6,0.8,Form("%1.2lf #leq #it{p}_{T} < %1.2lf GeV/#it{c}" ,xbins_pt[i], xbins_pt[i+1]));
    cplaceholder[i]->Update();
    Float_t ymax = cplaceholder[i]->GetUymax();
    Float_t ymin = cplaceholder[i]->GetUymin();
    Float_t xmax = mean[i]+3.*sigma[i];
    Float_t xmin = mean[i]-3.*sigma[i];
    TLine * fitmin = new TLine(xmin,ymin,xmin,ymax);
    TLine * fitmax = new TLine(xmax,ymin,xmax,ymax);

    cout << xmin << " " << ymin << " "  << xmax << " " << ymax << endl;
    fitmin->SetLineColor(kGreen+2);
    fitmax->SetLineColor(kGreen+2);
    fitmin->SetLineWidth(2);
    fitmax->SetLineWidth(2);
    fitmin->DrawLine(xmin,ymin,xmin,ymax);
    fitmax->DrawLine(xmax,ymin,xmax,ymax);

    cplaceholder[i]->SaveAs(Form("P_T_Spectra/P_TSpectra(%d).png",i));
    fitmin->Delete();
    fitmax->Delete();
    hMinvSpectra[i]->Delete();
    cplaceholder[i]->Clear();
    // r->Delete();
  }

  // Normieren auf Anzahl events. normierung auf bin breite bereits beim Fuellen oben!
  hP_TSpectrum->Scale(1./(500.));
  hP_TSpectrum_fit->Scale(1./500.);

  // Draw the actual pt spectrum
  cP_TSpectrum->cd();

  // gStyle->SetOptStat(0);
  cP_TSpectrum->SetLeftMargin(0.15);
  cP_TSpectrum->SetTopMargin(0.05);
  cP_TSpectrum->SetBottomMargin(0.15);
  cP_TSpectrum->SetRightMargin(0.05);
  cP_TSpectrum->SetTickx();
  cP_TSpectrum->SetTicky();
  cP_TSpectrum->SetLogy(1);
  cP_TSpectrum->SetLogx(0);

  leP_TSpectrum->AddEntry(hP_TSpectrum,"#pi^{0} #it{p}_{T} spectrum");
  leP_TSpectrum->AddEntry(hP_TSpectrum_fit,"fit function Integral");
  //hP_TSpectrum->SetError(int_error/(500.*hP_TSpectrum->GetBinWidth(i+1)));
  hP_TSpectrum->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hP_TSpectrum->SetYTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
  hP_TSpectrum->SetMarkerColor(kBlue+1);
  hP_TSpectrum_fit->SetMarkerColor(kRed);
  hP_TSpectrum->SetLineColor(kBlue+1);
  hP_TSpectrum_fit->SetLineColor(kRed);
  hP_TSpectrum->SetMarkerSize(1.2);
  hP_TSpectrum_fit->SetMarkerSize(1.2);
  hP_TSpectrum_fit->SetMarkerStyle(25);
  hP_TSpectrum->SetMarkerStyle(21);

  hP_TSpectrum->Draw("");
  hP_TSpectrum->GetXaxis()->SetRangeUser(1.25, 10.);
  hP_TSpectrum_fit->Draw("SAMEP");
  leP_TSpectrum->Draw("SAME");


  cP_TSpectrum->SaveAs(Form("P_T_Spectra/P_TSpectra.png"));
  hP_TSpectrum->SaveAs(Form("P_TSpectra.root"));

  cP_TSpectrum->Clear();
  laP_TSpectrum->Delete();
  leP_TSpectrum->Delete();



}
