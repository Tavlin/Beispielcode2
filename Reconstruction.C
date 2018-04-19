#include "CommenHeader.h"

void Reconstruction(TString AddName = "") {




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

  TGaxis::SetMaxDigits(3);

  // read Cluster Tree
  TFile* fDaten = new TFile("pi0_vcal_data.root", "READ");


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

    for (int i1 = 0; i1 < iCluster[iPufferAktuell]; i1++) {

      for(int i3 = 0; (i3 != i1) && (i3 < iCluster[iPufferAktuell]); i3++){
        // Paare im selben Event
        // ....
        px1 = px[iPufferAktuell][i1];
    	  py1 = py[iPufferAktuell][i1];
    	  pz1 = pz[iPufferAktuell][i1];
    	  px2 = px[iPufferAktuell][i3];
    	  py2 = py[iPufferAktuell][i3];
    	  pz2 = pz[iPufferAktuell][i3];

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


  // Wechsle und Zeichne minv same event
  cSignal->cd();
  cSignal->SetTopMargin(0.075);
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
  TF1* ratio_fit = new TF1("ratio_fit","[0]",0,0.3);
  ratio_fit->SetLineColor(kRed);



  // wechsle und zeichne Ratio zwischen same und mixed event
  hRatio = (TH1F*)hSignalmix->Clone("hRatio");
  hRatio->Divide(hSignal);
  hRatio->Fit(ratio_fit,"Q");
  hRatio->Draw();
  ratio_fit->Draw("l,same");


  cSignalSubtracted->cd();
  cSignalSubtracted->SetTopMargin(0.075);

  // gewichten der mixed events mit der ratio_fit
  TH1F* hSignalmix_clone = (TH1F*)hSignalmix->Clone("hSignalmix_clone");
  TH1F* hSignal_clone = (TH1F*)hSignal->Clone("hSignal_clone");
  hSignalmix_clone->Divide(ratio_fit);

  hSignal_clone->Add(hSignalmix_clone,-1);

  hSignal_clone->Draw();




  TFile* HistoFile = new TFile("HistoFile.root", "UPDATE");
  //Lese und speichere in Datei namens HistoFile.root
  if ( HistoFile->IsOpen() ) printf("HistoFile opened successfully\n");
  //gFile = HistoFile;
  hSignal_pT->Write("hSignal_pT");
  hSignalmix_pT->Write("hSignalmix_pT");

  cSignal->SaveAs(Form("Simulation/InavrianteMasseSameEvent%s.png", AddName.Data()));
  cSignalmix->SaveAs(Form("Simulation/InavrianteMasseDifferentEvent%s.png", AddName.Data()));
  cSignal_pT->SaveAs(Form("Simulation/InvarianteMasseTransversalImpulsSameEvent%s.png", AddName.Data()));
  cSignalmix_pT->SaveAs(Form("Simulation/InvarianteMasseTransversalImpulsDifferentEvents%s.png", AddName.Data()));
  cRatio->SaveAs(Form("Simulation/Ratio%s.png", AddName.Data()));
  cSignalSubtracted->SaveAs(Form("Simulation/InavrianteMasseOhneHintergrund%s.png", AddName.Data()));

  // schliesse datei #sauberes Programmieren
  HistoFile->Close();
  cout << "finished! :)" << endl;

}
