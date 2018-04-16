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

  TH1D* hSignal[40];

  for(int i = 0; i < 40; i++){

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


  for(int ip1 = 40; ip1 > 0; ip1--){


    int ip3 = ip1-1;
    int* pip3 = &ip3;

    // make projection
    TH1D* hSignal_pT_projection = hSignal_pT->ProjectionX("hSignal_pT_projection", ip1-1,ip1);
    TH1D* hSignal_pT_mix_projection = hSignalmix_pT->ProjectionX("hSignal_pT_mix_projection", ip1-1,ip1);


    // Wechsle zur Ratio
    cRatio->cd();
    // Definiere FitFunktion fuer die Ratio
    TF1* ratio_fit = new TF1("ratio_fit","[0]",0,0.3);
    ratio_fit->SetLineColor(kRed);

    // wechsle und zeichne Ratio zwischen same und mixed event
    TH1D* hRatio = (TH1D*)hSignal_pT_projection->Clone("hRatio");
    hRatio->Divide(hSignal_pT_mix_projection);

    for(int ip2 = ip1-1; hRatio->GetEntries() < 20 && ip2 >= 0; ip2--){

      //cout << "there are " << (int)hRatio->GetEntries() << " Etries in the bin(s)" << endl;

      TH1D* hSignal_pT_projection = hSignal_pT->ProjectionX("hSignal_pT_projection", ip2,ip1);
      TH1D* hSignal_pT_mix_projection = hSignalmix_pT->ProjectionX("hSignal_pT_mix_projection", ip2,ip1);

      TF1* ratio_fit = new TF1("ratio_fit","[0]",0,0.3);
      ratio_fit->SetLineColor(kRed);

      // wechsle und zeichne Ratio zwischen same und mixed event
      hRatio = (TH1D*)hSignal_pT_projection->Clone("hRatio");
      hRatio->Divide(hSignal_pT_mix_projection);


      pip3 = &ip2;
      ip3 = ip2;

    }



    hRatio->Fit(ratio_fit,"Q","",0.2,0.3);
    SetHistoStandardSettings(hRatio);
    hRatio->SetYTitle("same event/mixed event");
    hRatio->Draw();
    ratio_fit->Draw("l,same");


    cSignalSubtracted->cd();
    cSignalSubtracted->SetTopMargin(0.075);

    // gewichten der mixed events mit der ratio_fit
    TH1D* hSignal_pT_mix_projection_clone = (TH1D*)hSignal_pT_mix_projection->Clone("hSignal_pT_mix_projection_clone");
    TH1D* hSignal_pT_projection_clone = (TH1D*)hSignal_pT_projection->Clone("hSignal_pT_projection_clone");
    hSignal_pT_mix_projection_clone->Multiply(ratio_fit);


    hSignal_pT_projection_clone->Add(hSignal_pT_mix_projection_clone,-1);

    SetHistoStandardSettings(hSignal_pT_projection_clone);

    // Erstellen der Legenden
    TLegend *lSignal_pT_projection_clone = new TLegend(0.6,0.75,0.9,0.95);
    SetLegendSettigns(lSignal_pT_projection_clone);
    lSignal_pT_projection_clone->AddEntry(hSignal_pT_projection_clone,"#it{m}_{inv} spectrum");

    TLatex *ltSignal_pT_projection_clone = new TLatex();
    SetLatexSettings(ltSignal_pT_projection_clone);

    hSignal_pT_projection_clone->Draw();
    lSignal_pT_projection_clone->Draw("same");
    ltSignal_pT_projection_clone->DrawLatexNDC(0.6,0.8,Form("%1.2lf < #it{p}_{T} < %1.2lf GeV/#it{c}" ,0.25*(double)ip3, 0.25*(double)ip1));

    cSignal->cd();
    // Zeichnen von minv spektrum der einzelnen pt Bereiche mit skaliertem Hintergrund
    hSignal[ip1] = hSignal_pT->ProjectionX(Form("hSignal[%d]",ip1), ip3,ip1);
    TH1D* hSignal_mix = hSignalmix_pT->ProjectionX("hSignal_mix", ip3,ip1);
    hSignal_mix->Multiply(ratio_fit);
    hSignal[ip1]->SetMarkerStyle(34);
    hSignal[ip1]->SetMarkerColor(kBlue);
    hSignal_mix->SetLineColor(kRed);

    Float_t new_y_max = hSignal[ip1]->GetMaximum()*1.5;
    hSignal[ip1]->GetYaxis()->SetRangeUser(0,new_y_max);
    cout << "hSignal[ip1] GetMaximum = " << hSignal[ip1]->GetMaximum()*2 << endl;



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
    ltSignal->DrawLatexNDC(0.175,0.75,Form("%1.2lf < #it{p}_{T} < %1.2lf GeV/#it{c}" ,0.25*(double)ip3, 0.25*(double)ip1));


    cRatio->SaveAs(Form("Extraction/Ratio(%1.2lf-%1.2lf).png", 0.25*(double)ip3, 0.25*(double)ip1));
    cSignalSubtracted->SaveAs(Form("Extraction/InavrianteMasseOhneHintergrund_pt(%1.2lf-%1.2lf).png", 0.25*(double)ip3, 0.25*(double)ip1));
    cSignal->SaveAs(Form("Extraction/InavrianteMasseMitHintergrund(%1.2lf-%1.2lf).png", 0.25*(double)ip3, 0.25*(double)ip1));
    ip1 = ip3+1;

    TFile* HistoWOBackground_file = new TFile("HistoWOBackground_file.root", "UPDATE");
    //Lese und speichere in Datei namens HistoFile.root
    if ( HistoWOBackground_file->IsOpen() ) printf("HistoWOBackground_file opened successfully\n");

    hSignal[ip1]->Write(Form("hSignal[%d]",ip1));

    // schliesse datei #sauberes Programmieren
    HistoWOBackground_file->Close();
    cout << "finished! :)" << endl;

  }


















}
