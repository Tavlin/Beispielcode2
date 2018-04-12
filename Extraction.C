#include "CommenHeader.h"

void Extraction(TString AddName = ""){

  // Erstellen der Canvas  
  TCanvas *cSignal_pT = new TCanvas("cSignal_pT","",1080,1080);
  SetCanvasStandardSettings(cSignal_pT);
  
  TCanvas *cSignalmix_pT = new TCanvas("cSignalmix_pT","",1080,1080);
  SetCanvasStandardSettings(cSignalmix_pT);
  
  TCanvas *cRatio = new TCanvas("cRatio", "",1080,1080);
  SetCanvasStandardSettings(cRatio);
  
  TCanvas *cSignalSubtracted = new TCanvas("cSignalSubtracted", "",1080,1080);
  SetCanvasStandardSettings(cSignalSubtracted);
  
  // Erstellen der Legenden
  TLegend *lSignal_pT_projection_clone = new TLegend(0.5,0.75,0.9,0.95);
  
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
  TH1F* hRatio = new TH1F("hRatio","Ratio",150,0.,0.3);
  SetHistoStandardSettings(hRatio);



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
    TH1F* hRatio = (TH1F*)hSignal_pT_projection->Clone("hRatio");
    hRatio->Divide(hSignal_pT_mix_projection);
    
    for(int ip2 = ip1-1; hRatio->GetEntries() < 20 && ip2 >= 0; ip2--){

      //cout << "there are " << (int)hRatio->GetEntries() << " Etries in the bin(s)" << endl;
      
      TH1D* hSignal_pT_projection = hSignal_pT->ProjectionX("hSignal_pT_projection", ip2,ip1);
      TH1D* hSignal_pT_mix_projection = hSignalmix_pT->ProjectionX("hSignal_pT_mix_projection", ip2,ip1);
      
      TF1* ratio_fit = new TF1("ratio_fit","[0]",0,0.3);
      ratio_fit->SetLineColor(kRed);
      
      // wechsle und zeichne Ratio zwischen same und mixed event
      hRatio = (TH1F*)hSignal_pT_projection->Clone("hRatio");
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
    TH1F* hSignal_pT_mix_projection_clone = (TH1F*)hSignal_pT_mix_projection->Clone("hSignal_pT_mix_projection_clone");
    TH1F* hSignal_pT_projection_clone = (TH1F*)hSignal_pT_projection->Clone("hSignal_pT_projection_clone");
    hSignal_pT_mix_projection_clone->Multiply(ratio_fit);
    
    
    hSignal_pT_projection_clone->Add(hSignal_pT_mix_projection_clone,-1);
    
    SetHistoStandardSettings(hSignal_pT_projection_clone);
    
    string ubin = to_string(ip1);
    string lbin = to_string(ip3);
    
    lSignal_pT_projection_clone->AddEntry(hSignal_pT_projection_clone,)
    
    hSignal_pT_projection_clone->Draw();
    
    
    cRatio->SaveAs(Form("Extraction/Ratio_pt(%1.2lf-%1.2lf).png", 0.25*(double)ip3, 0.25*(double)ip1));
    cSignalSubtracted->SaveAs(Form("Extraction/InavrianteMasseOhneHintergrund_pt(%1.2lf-%1.2lf).png", 0.25*(double)ip3, 0.25*(double)ip1));
    
    ip1 = ip3+1;

  }


















}
