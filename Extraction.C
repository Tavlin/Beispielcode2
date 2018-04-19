#include "CommenHeader.h"

void Extraction(TString AddName = ""){

  // Erstellen der Canvas
  TCanvas *cRatio = new TCanvas("cRatio", "",1080,1080);
  SetCanvasStandardSettings(cRatio);

  TCanvas *cSignalSubtracted = new TCanvas("cSignalSubtracted", "",1080,1080);
  SetCanvasStandardSettings(cSignalSubtracted);

  TCanvas *cSignal = new TCanvas("cSignal", "",1080,1080);
  SetCanvasStandardSettings(cSignal);


  // Erstellen der Latex-Objekte
  TLatex *lSignal_pT = new TLatex();
  SetLatexSettings(lSignal_pT);

  TLatex *lSignalmix_pT = new TLatex();
  SetLatexSettings(lSignalmix_pT);

  // Erstellen der 2D Histos
  TH2F* hSignal_pT = new TH2F();
  SetHistoStandardSettings2(hSignal_pT);

  TH2F* hSignalmix_pT = new TH2F();
  SetHistoStandardSettings2(hSignalmix_pT);


  // Erstelle 1D Histogram fuer Ratio
  TH1D* hRatio = new TH1D("hRatio","Ratio",150,0.,0.3);
  SetHistoStandardSettings(hRatio);

  TH1D* hSignal[101];
  TH1D* hSignal_pT_projection_clone[101];

  for(int i = 0; i < 100; i++){

    hSignal[i] = new TH1D(Form("hSignal[%d]",i),Form("Blubb[%d]",i),150,0.,0.3);
    SetHistoStandardSettings(hSignal[i]);
  }


  //Open and read file
  TFile* HistoFile = new TFile("HistoFile.root", "READ");
  if ( HistoFile->IsOpen() ) printf("HistoFile opened successfully\n");
  //gFile = HistoFile;

  if(HistoFile->IsZombie()){
    std::cout << "ERROR: Data File not found" << std::endl;
    return;
  }
  // get Histos from file
  gDirectory->GetObject("hSignal_pT",hSignal_pT);
  gDirectory->GetObject("hSignalmix_pT",hSignalmix_pT);


  for(int ip1 = 100; ip1 > 0; ip1--){


    int ip3 = ip1-1;
    int* pip3 = &ip3;

    // make projection
    TH1D* hSignal_pT_projection = hSignal_pT->ProjectionX("hSignal_pT_projection", ip1-1,ip1);
    TH1D* hSignal_pT_mix_projection = hSignalmix_pT->ProjectionX("hSignal_pT_mix_projection", ip1-1,ip1);


    // Wechsle zur Ratio
    cRatio->cd();


    // wechsle und zeichne Ratio zwischen same und mixed event
    TH1D* hRatio = (TH1D*)hSignal_pT_projection->Clone("hRatio");
    hRatio->Divide(hSignal_pT_mix_projection);

    for(int ip2 = ip1-1; hRatio->GetEntries() < 20 && ip2 >= 0; ip2--){

      //cout << "there are " << (int)hRatio->GetEntries() << " Etries in the bin(s)" << endl;

      TH1D* hSignal_pT_projection = hSignal_pT->ProjectionX("hSignal_pT_projection", ip2,ip1);
      TH1D* hSignal_pT_mix_projection = hSignalmix_pT->ProjectionX("hSignal_pT_mix_projection", ip2,ip1);

      // wechsle und zeichne Ratio zwischen same und mixed event
      hRatio = (TH1D*)hSignal_pT_projection->Clone("hRatio");
      hRatio->Divide(hSignal_pT_mix_projection);


      pip3 = &ip2;
      ip3 = ip2;

    }


    // Fit funktionen rechts und links, sowie dann kombiniert
    TF1* ratio_fit = new TF1("ratio_fit","[0]",0,0.3);
    ratio_fit->SetLineColor(kGreen+2);

    SetHistoStandardSettings(hRatio);
    hRatio->SetYTitle("same event/mixed event");

    for (size_t num_bin = 40; num_bin <= 100; num_bin++) {
      hRatio->SetBinContent(num_bin,0);
      hRatio->SetBinError(num_bin,0);
    }

    hRatio->Fit(ratio_fit);
    hRatio->Draw();
    ratio_fit->Draw("L,SAME");


    cSignalSubtracted->cd();
    cSignalSubtracted->SetTopMargin(0.075);

    // Gewichten der mixed events mit der ratio_fit
    TH1D* hSignal_pT_mix_projection_clone = (TH1D*)hSignal_pT_mix_projection->Clone("hSignal_pT_mix_projection_clone");
    hSignal_pT_projection_clone[ip1] = (TH1D*)hSignal_pT_projection->Clone("hSignal_pT_projection_clone");
    hSignal_pT_mix_projection_clone->Multiply(ratio_fit);


    hSignal_pT_projection_clone[ip1]->Add(hSignal_pT_mix_projection_clone,-1);

    SetHistoStandardSettings(hSignal_pT_projection_clone[ip1]);

    // Erstellen der Legenden und Textbox
    TLegend *lSignal_pT_projection_clone = new TLegend(0.65,0.75,0.9,0.95);
    SetLegendSettigns(lSignal_pT_projection_clone);
    lSignal_pT_projection_clone->AddEntry(hSignal_pT_projection_clone[ip1],"#it{m}_{inv} spectrum");

    TLatex *ltSignal_pT_projection_clone = new TLatex();
    SetLatexSettings(ltSignal_pT_projection_clone);

    hSignal_pT_projection_clone[ip1]->SetYTitle("#frac{d#it{N}}{d#it{m}_{inv}} #left(GeV/c^{2}#right)^{-1}");
    hSignal_pT_projection_clone[ip1]->Draw();
    lSignal_pT_projection_clone->Draw("same");
    ltSignal_pT_projection_clone->DrawLatexNDC(0.6,0.8,Form("%1.2lf < #it{p}_{T} < %1.2lf GeV/#it{c}" ,0.1*(double)ip3, 0.1*(double)ip1));

    cSignal->cd();
    // Zeichnen von minv spektrum der einzelnen pt Bereiche mit skaliertem Hintergrund
    hSignal[ip1] = hSignal_pT->ProjectionX(Form("hSignal[%d]",ip1), ip3,ip1);
    TH1D* hSignal_mix = hSignalmix_pT->ProjectionX("hSignal_mix", ip3,ip1);
    hSignal_mix->Multiply(ratio_fit);
    hSignal[ip1]->SetMarkerStyle(34);
    hSignal[ip1]->SetMarkerColor(kBlue);
    hSignal_mix->SetLineColor(kRed);
    hSignal_mix->SetMarkerColor(kRed);
    hSignal_mix->SetMarkerStyle(20);

    Float_t new_y_max = hSignal[ip1]->GetMaximum()*1.5;
    hSignal[ip1]->GetYaxis()->SetRangeUser(0,new_y_max);


    hSignal[ip1]->Draw("EP");
    hSignal_mix->Draw("sameL");


    // Legende fuer minv spektrum der einzelnen pt Bereiche mit skaliertem Hintergrund
    TLegend *lSignal = new TLegend(0.175,0.75,0.9,0.95);
    SetLegendSettigns(lSignal);
    lSignal->AddEntry(hSignal[ip1],"#it{m}_{inv} spectrum same event");
    lSignal->AddEntry(hSignal_mix,"scaled #it{m}_{inv} spectrum mixed event");
    lSignal->Draw("SAME");

    // Latex fuer minv spektrum der einzelnen pt Bereiche mit skaliertem Hintergrund
    TLatex *ltSignal = new TLatex();
    SetLatexSettings(ltSignal);
    ltSignal->DrawLatexNDC(0.175,0.75,Form("%1.2lf < #it{p}_{T} < %1.2lf GeV/#it{c}" ,0.1*(double)ip3, 0.1*(double)ip1));


    cRatio->SaveAs(Form("Extraction/Ratio(%1.2lf-%1.2lf).png", 0.1*(double)ip3, 0.1*(double)ip1));
    cSignalSubtracted->SaveAs(Form("Extraction/InavrianteMasseOhneHintergrund_pt(%1.2lf-%1.2lf).png", 0.1*(double)ip3, 0.1*(double)ip1));
    cSignal->SaveAs(Form("Extraction/InavrianteMasseMitHintergrund(%1.2lf-%1.2lf).png", 0.1*(double)ip3, 0.1*(double)ip1));
    ip1 = ip3+1;

  }

  TFile* HistoWOBackground_file = new TFile("HistoWOBackground_file.root", "UPDATE");
  //Lese und speichere in Datei namens HistoFile.root
  //if ( HistoWOBackground_file->IsOpen() ) printf("HistoWOBackground_file opened successfully\n");

  for (size_t k = 1; k <= 62; k++) {
    hSignal_pT_projection_clone[k]->Write(Form("hSignal[%lu]",k));
  }
  hSignal_pT_projection_clone[65]->Write("hSignal[63]");
  hSignal_pT_projection_clone[70]->Write("hSignal[64]");
  hSignal_pT_projection_clone[79]->Write("hSignal[65]");
  hSignal_pT_projection_clone[100]->Write("hSignal[66]");

  // schliesse datei #sauberes Programmieren
  HistoWOBackground_file->Close();
  cout << "finished! :)" << endl;

















}
