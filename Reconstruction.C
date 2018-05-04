#include "CommenHeader.h"

using namespace std;
void Reconstruction(TString AddName = "") {

  TLatex* poweektex = new TLatex();

  Float_t pi = TMath::Pi();
  std::vector<Double_t> eta_vec;
  std::vector<Double_t> phi_vec;


  // Erstellen der Canvas
  TCanvas *cSignal = new TCanvas("cSignal","",1080,1080);
  SetCanvasStandardSettings(cSignal);

  TCanvas *cSignalmix = new TCanvas("cSignalmix","",1080,1080);
  SetCanvasStandardSettings(cSignalmix);

  TCanvas *cSignal_pT = new TCanvas("cSignal_pT","",1080,1080);
  SetCanvasStandardSettings(cSignal_pT);

  TCanvas *cSignalmix_pT = new TCanvas("cSignalmix_pT","",1080,1080);
  SetCanvasStandardSettings(cSignalmix_pT);

  TCanvas *cRatio = new TCanvas("cRatio", "",1080,1080);
  SetCanvasStandardSettings(cRatio);

  TCanvas *cSignalSubtracted = new TCanvas("cSignalSubtracted", "",1080,1080);
  SetCanvasStandardSettings(cSignalSubtracted);

  TCanvas *cBothSignals = new TCanvas("cBothSignals", "",1080,1080);
  SetCanvasStandardSettings(cBothSignals);

  TCanvas *cWinkelAbdeckung = new TCanvas("cWinkelAbdeckung", "",1080,1080);
  SetCanvasStandardSettings(cWinkelAbdeckung);


  // Erstellen der Latex-Objekte
  TLatex *lSignal = new TLatex();
  SetLatexSettings(lSignal);

  TLatex *lSignalmix = new TLatex();
  SetLatexSettings(lSignalmix);

  TLatex *lSignal_pT = new TLatex();
  SetLatexSettings(lSignal_pT);

  TLatex *lSignalmix_pT = new TLatex();
  SetLatexSettings(lSignalmix_pT);

  TLatex *lRatio = new TLatex();
  SetLatexSettings(lRatio);

  TLatex *lSignalSubtracted = new TLatex();
  SetLatexSettings(lSignalSubtracted);


  // Erstellen der 1D Histos
  TH1F* hSignal = new TH1F("hSignal","invariante Masse",150,0.,0.3);
  SetHistoStandardSettings(hSignal);

  TH1F* hSignalmix = new TH1F("hSignalmix","invariante Masse (mixed events)",150,0.,0.3);
  SetHistoStandardSettings(hSignalmix);

  TH1F* hRatio = new TH1F("hRatio","Ratio",150,0.,0.3);
  SetHistoStandardSettings(hRatio);


  // Erstellen der 2D Histos
  TH2F* hSignal_pT = new TH2F("hSignal_pT","invariante Masse gegen pT",150,0.,0.3,100,0.,10.);
  SetHistoStandardSettings2(hSignal_pT);

  TH2F* hSignalmix_pT = new TH2F("hSignalmix_pT","invariante Masse gegen pT (mixed events)",150,0.,0.3,100,0.,10.);
  SetHistoStandardSettings2(hSignalmix_pT);

  TH2F* hWinkelAbdeckung = new TH2F("hWinkelAbdeckung","Winkelabdeckung same events",100,-1.,1.,180,-1.*pi,1.0001*pi);
  SetHistoStandardSettings2(hWinkelAbdeckung);
  hWinkelAbdeckung->SetXTitle("#it{#eta}");
  hWinkelAbdeckung->SetYTitle("#it{#phi}");
  hWinkelAbdeckung->GetYaxis()->SetBit(TAxis::kLabelsVert);
  hWinkelAbdeckung->GetYaxis()->SetBinLabel(1,"0");
  hWinkelAbdeckung->GetYaxis()->SetBinLabel(hWinkelAbdeckung->GetYaxis()->FindBin(pi), "#pi");
  hWinkelAbdeckung->GetYaxis()->SetBinLabel(hWinkelAbdeckung->GetYaxis()->FindBin(pi/2.),"#frac{-#pi}{2}");
  hWinkelAbdeckung->GetYaxis()->SetBinLabel(hWinkelAbdeckung->GetYaxis()->FindBin(0.),"0");
  hWinkelAbdeckung->GetYaxis()->SetBinLabel(hWinkelAbdeckung->GetYaxis()->FindBin(-pi), "#pi");
  hWinkelAbdeckung->GetYaxis()->SetBinLabel(hWinkelAbdeckung->GetYaxis()->FindBin(-pi/2.),"#frac{-#pi}{2}");
  hWinkelAbdeckung->GetXaxis()->SetNdivisions(5);

  TGaxis::SetMaxDigits(3);

  // read Cluster Tree
  TFile* fDaten = new TFile("pi0_vcal_data.root", "READ");
  if ( fDaten->IsOpen() ) printf("pi0_vcal_data.root opened successfully\n");


  const Int_t iPufferMax = 2;
  //Int_t NMaxEvents = 1000;


  if(fDaten->IsZombie()){
    cout << "ERROR: Data File not found" << endl;
    return;
  }
  DataTree *Daten = new DataTree(fDaten);
  Int_t NMaxEvents = Daten->GetNEvents();

  // Allokiere Arrays in die wir die Cluster schreiben
  Float_t px[iPufferMax][kMaxHit];
  Float_t py[iPufferMax][kMaxHit];
  Float_t pz[iPufferMax][kMaxHit];
  Float_t iCluster[iPufferMax];
  Int_t iNCluster;

  // Setze Anzahl der Cluster in jedem Puffer auf null
  for (Int_t i = 0; i < iPufferMax; i++) {
    iCluster[i] = 0;
  }
  int iPufferAktuell = 0;
  cout<<"Number of events: "<<NMaxEvents<<endl;
  for (Int_t iEvt=0; iEvt < NMaxEvents; iEvt++) {
    if(iEvt%1000 == 0 && iEvt !=0 ) cout << iEvt/1000 << "x10^3 Events have been analyzed" << endl;

    // Anzahl Cluster im Event
    iNCluster = Daten->GetClusterID(iEvt);

    for (Int_t iHit=0; iHit < iNCluster; iHit++) {
      px[iPufferAktuell][iHit] = Daten->GetPX(iEvt, iHit);
      py[iPufferAktuell][iHit] = Daten->GetPY(iEvt, iHit);
      pz[iPufferAktuell][iHit] = Daten->GetPZ(iEvt, iHit);

    }
    iCluster[iPufferAktuell] = iNCluster;

    // starte Analyse des letzten Events
    Float_t pair_pt;
    Float_t minv;
    Float_t px1;
    Float_t py1;
    Float_t pz1;
    Float_t px2;
    Float_t py2;
    Float_t pz2;

    Double_t eta, phi;

    for (int i1 = 0; i1 < iCluster[iPufferAktuell]; i1++) {

      for(int i3 = i1+1; (i3 < iCluster[iPufferAktuell]); i3++){
        // Paare im selben Event
        // ....
        px1 = px[iPufferAktuell][i1];
    	  py1 = py[iPufferAktuell][i1];
    	  pz1 = pz[iPufferAktuell][i1];
    	  px2 = px[iPufferAktuell][i3];
    	  py2 = py[iPufferAktuell][i3];
    	  pz2 = pz[iPufferAktuell][i3];

        // Need to be discussed, phi angeblich 360°!?!
        eta = atanh(pz1/sqrt(px1*px1+py1*py1+pz1*pz1));
        eta_vec.push_back(eta);

        phi = atan2(py1,px1);
        phi_vec.push_back(phi);

        hWinkelAbdeckung->Fill(eta,phi);


        pair_pt = fCalcPT(px1,py1,px2,py2); //pair_pt_same?
        if (pair_pt > 0) {
            minv = fCalcInvMass(px1,py1,pz1,px2,py2,pz2); //minv_same?
            hSignal->Fill(minv);

            hSignal_pT->Fill(minv,pair_pt);
        }
      }
      // Paare aus unterschiedlichen Events (Event Mixing)
      int iPufferAlt;
      if (iPufferAktuell == 0) {
	       iPufferAlt = 1;
      } else {
	       iPufferAlt = 0;
      }
      for (int i2 = 0; i2 < iCluster[iPufferAlt]; i2++) {
  	     px1 = px[iPufferAktuell][i1];
  	     py1 = py[iPufferAktuell][i1];
  	     pz1 = pz[iPufferAktuell][i1];
  	     px2 = px[iPufferAlt][i2];
  	     py2 = py[iPufferAlt][i2];
  	     pz2 = pz[iPufferAlt][i2];

         pair_pt = fCalcPT(px1,py1,px2,py2);
      	 if (pair_pt > 0) {
      	    minv = fCalcInvMass(px1,py1,pz1,px2,py2,pz2);
      	    hSignalmix->Fill(minv);

      	    hSignalmix_pT->Fill(minv,pair_pt);
      	 }

      }
    }

    // Bereite naechstes Event vor

    // Puffer umschalten
    iPufferAktuell++;
    if (iPufferAktuell == iPufferMax) {
      iPufferAktuell = 0;
    }

  }
  fDaten->Close();

  // Wechsle und Zeichne minv same event
  cSignal->cd();
  cSignal->SetTopMargin(0.075);
  // hSignal->Scale(1,"width");
  hSignal->Draw("");
  lSignal->DrawLatex(0.01, 7e3, "#it{m}_{inv} same event");

  // Wechsle und Zeichne minv-pt same event
  cSignal_pT->cd();
  cSignal_pT->SetRightMargin(0.175);
  cSignal_pT->SetBottomMargin(0.125);
  gPad->SetLogz();

  hSignal_pT->GetZaxis()->SetRangeUser(1.e0,1.e4);
  hSignal_pT->Draw("colz");
  lSignal_pT->DrawLatex(0.01, 9., "#it{m}_{inv} #it{p}_{T} same event");

  // Wechsle und Zeichne minv different event
  cSignalmix->cd();
  cSignalmix->SetTopMargin(0.075);
  // hSignalmix->Scale(1,"width");
  hSignalmix->Draw("");
  lSignalmix->DrawLatex(0.01, 14e3, "#it{m}_{inv} mixed event");

  // Wechsle und Zeichne minv-pt different event
  cSignalmix_pT->cd();
  cSignalmix_pT->SetRightMargin(0.175);
  cSignalmix_pT->SetBottomMargin(0.125);
  gPad->SetLogz();
  lSignalmix_pT->DrawLatex(0.01, 9., "#it{m}_{inv} #it{p}_{T} mixed event");

  hSignalmix_pT->GetZaxis()->SetRangeUser(1.e0,1.e4);
  hSignalmix_pT->Draw("colz");


  // Wechsle und zeichne die Ratio
  cRatio->cd();

  // define fit function for ratio
  TF1* ratio_fit = new TF1("ratio_fit","[0]",0.2,0.3);
  ratio_fit->SetLineColor(kRed);



  // wechsle und zeichne Ratio zwischen same und mixed event
  hRatio = (TH1F*)hSignalmix->Clone("hRatio");
  hRatio->Divide(hSignal);
  hRatio->Fit(ratio_fit,"M","",0.2,0.3);
  hRatio->Draw();
  ratio_fit->Draw("l,same");


  cSignalSubtracted->cd();
  cSignalSubtracted->SetTopMargin(0.075);

  // gewichten der mixed events mit der ratio_fit
  TH1F* hSignalmix_clone = (TH1F*)hSignalmix->Clone("hSignalmix_clone");
  TH1F* hSignal_clone = (TH1F*)hSignal->Clone("hSignal_clone");
  hSignalmix_clone->Scale(1/ratio_fit->GetParameter(0));

  hSignal_clone->Add(hSignalmix_clone,-1);


  hSignal_clone->Scale(1,"width");
  hSignal_clone->Draw();




  TFile* HistoFile = new TFile("HistoFile.root", "RECREATE");
  //Lese und speichere in Datei namens HistoFile.root
  if ( HistoFile->IsOpen() ) printf("HistoFile opened successfully\n");
  //gFile = HistoFile;
  hSignal_pT->Write("hSignal_pT");
  hSignalmix_pT->Write("hSignalmix_pT");

  cSignal->SaveAs(Form("Reconstructed/InavrianteMasseSameEvent%s.png", AddName.Data()));
  cSignalmix->SaveAs(Form("Reconstructed/InavrianteMasseDifferentEvent%s.png", AddName.Data()));
  cSignal_pT->SaveAs(Form("Reconstructed/InvarianteMasseTransversalImpulsSameEvent%s.png", AddName.Data()));
  cSignalmix_pT->SaveAs(Form("Reconstructed/InvarianteMasseTransversalImpulsDifferentEvents%s.png", AddName.Data()));
  cRatio->SaveAs(Form("Reconstructed/Ratio%s.png", AddName.Data()));
  cSignalSubtracted->SaveAs(Form("Reconstructed/InavrianteMasseOhneHintergrund%s.png", AddName.Data()));

  // schliesse datei #sauberes Programmieren
  HistoFile->Close();

  // Error plot fuer minv same events
  TCanvas *cErrors = new TCanvas("cErrors","",1080,1080);
  SetCanvasStandardSettings(cErrors);
  Double_t x_Error[150], y_Error[150];
  Double_t x_Error_mms[150], y_Error_mms[150];
  Double_t x_Error_ms_mm[150], y_Error_ms_mm[150];
  Double_t x_Error_ms_mm_root[150], y_Error_ms_mm_root[150];


  Int_t n = 150;
  Double_t fse[150], sse[150]; //first_scaling_error and second_scaling_error
  for (Int_t i = 0; i < n; i++) {
    // same event
    x_Error[i] = hSignal->GetBinContent(i);
    y_Error[i] = hSignal->GetBinError(i);

    // mixed
    x_Error_mms[i] = hSignalmix->GetBinContent(i);
    fse[i] = hSignalmix->GetBinError(i)/ratio_fit->GetParameter(0);
    sse[i] = (hSignalmix->GetBinContent(i)*ratio_fit->GetParError(0)/(ratio_fit->GetParameter(0)*ratio_fit->GetParameter(0)));
    y_Error_mms[i] = sqrt((fse[i]*fse[i])+(sse[i]*sse[i]));

    // same - scaled mixed event haendisch
    x_Error_ms_mm[i] = hSignal_clone->GetBinContent(i);
    y_Error_ms_mm[i] = sqrt(hSignal->GetBinContent(i)+(y_Error_mms[i])*(y_Error_mms[i]))*500;

    // same - mixed event root
    x_Error_ms_mm_root[i] = hSignal_clone->GetBinContent(i);
    y_Error_ms_mm_root[i] = hSignal_clone->GetBinError(i);

  }
  // Error minv same event
  TGraph* grErrors = new TGraph(n,x_Error,y_Error);
  grErrors->SetFillStyle(0);
  grErrors->SetFillColor(0);
  grErrors->SetLineWidth(3);
  grErrors->SetMarkerColor(kRed);
  grErrors->SetMarkerStyle(20);
  grErrors->SetMarkerSize(1.5);

  grErrors->SetTitle(";#frac{d#it{N}_{#gamma #gamma}}{d#it{m}_{inv}} (GeV/#it{c}^{2})^{-1};#it{#sigma} (GeV/#it{c}^{2})^{-1}");

  // Error minv mixed scaled events
  TGraph* grError_mms = new TGraph(n,x_Error_mms,y_Error_mms);
  grError_mms->SetFillStyle(0);
  grError_mms->SetFillColor(0);
  grError_mms->SetLineWidth(3);
  grError_mms->SetLineColor(kRed);
  grError_mms->SetMarkerColor(kRed);
  grError_mms->SetMarkerStyle(20);
  grError_mms->SetMarkerSize(1.5);
  grError_mms->SetTitle(";#frac{d#it{N}_{#gamma #gamma}}{d#it{m}_{inv}} (GeV/#it{c}^{2})^{-1};#it{#sigma} (GeV/#it{c}^{2})^{-1}");

  // Error same - scaled mixed
  TGraph* grErrors_ms_mm = new TGraph(n,x_Error_ms_mm,y_Error_ms_mm);
  grErrors_ms_mm->SetFillStyle(0);
  grErrors_ms_mm->SetFillColor(0);
  grErrors_ms_mm->SetLineColor(kRed);
  grErrors_ms_mm->SetMarkerColor(kRed);
  grErrors_ms_mm->SetMarkerStyle(21);
  grErrors_ms_mm->SetMarkerSize(1.5);
  grErrors_ms_mm->SetTitle(";#frac{d#it{N}_{#gamma #gamma}}{d#it{m}_{inv}} (GeV/#it{c}^{2})^{-1};#it{#sigma} (GeV/#it{c}^{2})^{-1}");

  // Error same - scaled mixed by root functions
  TGraph* grErrors_ms_mm_root = new TGraph(n,x_Error_ms_mm_root,y_Error_ms_mm_root);
  grErrors_ms_mm_root->SetFillStyle(0);
  grErrors_ms_mm_root->SetFillColor(0);
  grErrors_ms_mm_root->SetLineColor(kBlue);
  grErrors_ms_mm_root->SetMarkerColor(kBlue);
  grErrors_ms_mm_root->SetMarkerStyle(25);
  grErrors_ms_mm_root->SetMarkerSize(1.5);
  grErrors_ms_mm_root->SetTitle(";#frac{d#it{N}_{#gamma #gamma}}{d#it{m}_{inv}} (GeV/#it{c}^{2})^{-1}; #it{#sigma} (GeV/#it{c}^{2})^{-1}");

  TLegend* legErrors = new TLegend(0.4,0.3,0.7,0.4);
  SetLegendSettigns(legErrors);
  legErrors->SetTextSize(0.03);
  legErrors->AddEntry(grErrors_ms_mm, "self-calculated uncertainty");
  legErrors->AddEntry(grErrors_ms_mm_root, "root 5.34 calculated uncertainty");

  TF1* sqrtNsame = new TF1("sqrtNsame","sqrt(x)",0.,8.e3);
  sqrtNsame->SetLineColor(kBlue+1);
  sqrtNsame->SetLineWidth(3);
  TLegend* legErrorsSame = new TLegend(0.25,0.75,0.6,0.85);
  SetLegendSettigns(legErrorsSame);
  legErrorsSame->AddEntry(grErrors, "self-calculated uncertainty");
  legErrorsSame->AddEntry(sqrtNsame, "#sqrt{N} fit");


  TF1* sqrtNmixed = new TF1("sqrtNmixed", "sqrt((sqrt(x)/[0])^2+((x*[1])/[2])^2)",0.,16.e3);
  sqrtNmixed->SetLineColor(kBlue+1);
  sqrtNmixed->SetLineWidth(3);
  TLegend* legErrorsMixed = new TLegend(0.25,0.75,0.6,0.85);
  SetLegendSettigns(legErrorsMixed);
  legErrorsMixed->AddEntry(grError_mms, "self-calculated uncertainty");
  legErrorsMixed->AddEntry(sqrtNmixed, "#sqrt{N} fit");


  cErrors->cd();
  cErrors->SetTopMargin(0.075);
  cErrors->SetRightMargin(0.1);
  cErrors->SetBottomMargin(0.18);
  grErrors->GetXaxis()->SetTitleOffset(1.6);
  grErrors->GetYaxis()->SetTitleOffset(1.4);
  grErrors->Draw("AP");
  sqrtNsame->Draw("same");
  legErrorsSame->Draw("same");
  poweektex->SetTextSize(0.04);
  poweektex->DrawLatexNDC(0.45,0.5,poweek_str);
  poweektex->DrawLatexNDC(0.45,0.45,pi0togamma_str);
  poweektex->SetTextSize(0.03);
  poweektex->DrawLatexNDC(0.45,0.4,"same events");
  cErrors->Update();


  cErrors->SaveAs(Form("Reconstructed/ErrorPlot%s.png", AddName.Data()));

  cErrors->Clear();
  cErrors->SetBottomMargin(0.15);
  grError_mms->Draw("AP");
  grError_mms->GetYaxis()->SetTitleOffset(1.4);
  grError_mms->GetXaxis()->SetTitleOffset(1.6);
  grError_mms->Fit(sqrtNmixed,"M");
  sqrtNmixed->Draw("same");
  legErrorsMixed->Draw("same");
  poweektex->SetTextSize(0.04);
  poweektex->DrawLatexNDC(0.45,0.5,poweek_str);
  poweektex->DrawLatexNDC(0.45,0.45,pi0togamma_str);
  poweektex->SetTextSize(0.03);
  poweektex->DrawLatexNDC(0.45,0.4,"mixed events");
  cErrors->SaveAs(Form("Reconstructed/MixedScaledErrorPlot%s.png", AddName.Data()));

  cErrors->Clear();
  cErrors->SetRightMargin(0.07);

  grErrors_ms_mm->Draw("AP");
  grErrors_ms_mm_root->Draw("P");
  grErrors_ms_mm->GetYaxis()->SetTitleOffset(1.4);
  grErrors_ms_mm_root->GetYaxis()->SetTitleOffset(1.4);
  grErrors_ms_mm->GetXaxis()->SetTitleOffset(1.6);
  grErrors_ms_mm_root->GetXaxis()->SetTitleOffset(1.6);
  legErrors->Draw("SAMEP");
  poweektex->SetTextSize(0.04);
  poweektex->DrawLatexNDC(0.4,0.6,poweek_str);
  poweektex->DrawLatexNDC(0.4,0.55,pi0togamma_str);
  poweektex->SetTextSize(0.03);
  poweektex->DrawLatexNDC(0.4,0.5,"same event - scaled mixed event");
  cErrors->Update();

  cErrors->SaveAs(Form("Reconstructed/FurtherErrorPlot%s.png", AddName.Data()));


  grErrors->Delete();
  hSignalmix_pT->Delete();
  hSignal_pT->Delete();
  hRatio->Delete();
  ratio_fit->Delete();
  hSignalmix_clone->Delete();
  hSignal_clone->Delete();
  cErrors->Clear();
  cErrors->Delete();

  cWinkelAbdeckung->cd();
  cWinkelAbdeckung->SetRightMargin(0.175);
  cWinkelAbdeckung->SetBottomMargin(0.125);
  cWinkelAbdeckung->SetTopMargin(0.2);
  gPad->SetLogz();

  hWinkelAbdeckung->GetZaxis()->SetRangeUser(1.e0,1.e4);
  hWinkelAbdeckung->Draw("colz");
  lSignalmix_pT->DrawLatexNDC(0.32, 0.95, "Winkelabdeckung same event");
  cWinkelAbdeckung->Update();

  cWinkelAbdeckung->SaveAs(Form("Reconstructed/Winkelabdeckung%s.png", AddName.Data()));

  cWinkelAbdeckung->Clear();
  hWinkelAbdeckung->Delete();
  cWinkelAbdeckung->Delete();


  std::cout << "eta min = " << *min_element(eta_vec.begin(), eta_vec.end())  << std::endl;
  std::cout << "eta max = " << *max_element(eta_vec.begin(), eta_vec.end())  << std::endl;
  std::cout << "phi min = " << *min_element(phi_vec.begin(), phi_vec.end())  << std::endl;
  std::cout << "phi max = " << *max_element(phi_vec.begin(), phi_vec.end())  << std::endl;









  cout << "finished! :)" << endl;

}
