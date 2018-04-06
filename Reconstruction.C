#include "CommenHeader.h"

void Reconstruction(TString AddName = "") {

  TCanvas *cSignal = new TCanvas("cSignal","",1080,1080);
  SetCanvasStandardSettings(cSignal);
  
  TCanvas *cSignalmix = new TCanvas("cSignalmix","",1080,1080);
  SetCanvasStandardSettings(cSignalmix);
  
  TCanvas *cSignal_pT = new TCanvas("cSignal_pT","",1080,1080);
  SetCanvasStandardSettings(cSignal_pT);
  
  TCanvas *cSignalmix_pT = new TCanvas("cSignalmix_pT","",1080,1080);
  SetCanvasStandardSettings(cSignalmix_pT);
  
  TCanvas *cSignal_subtracted = new TCanvas("cSignal_subtracted","",1080,1080);
  SetCanvasStandardSettings(cSignal_subtracted);
  
  TCanvas *cSignal_pT_subtracted = new TCanvas("cSignal_pT_subtracted","",1080,1080);
  SetCanvasStandardSettings(cSignal_pT_subtracted);
  
  
  TH1F* hSignal = new TH1F("hSignal","invariante Masse",100,0.,0.3);
  SetHistoStandardSettings(hSignal);
  
  TH1F* hSignalmix = new TH1F("hSignalmix","invariante Masse (mixed events)",100,0.,0.3);
  SetHistoStandardSettings(hSignalmix);
  
  TH2F* hSignal_pT = new TH2F("hSignal_pT","invariante Masse gegen pT",100,0.,0.3,40,0.,10.);
  SetHistoStandardSettings2(hSignal_pT);
  
  TH2F* hSignalmix_pT = new TH2F("hSignalmix_pT","invariante Masse gegen pT (mixed events)",100,0.,0.3,40,0.,10.);
  SetHistoStandardSettings2(hSignalmix_pT);

  //TH1F* hSignal_subtracted = new TH1F("hSignal_subtracted","invariante Masse (ohne Hintergrund)",100,0.,0.3);
  //SetHistoStandardSettings(hSignal_subtracted);
  
  //TH2F* hSignal_pT_subtracted = new TH2F("hSignal_pT_subtracted","invariante Masse gegen pT (ohne Hintergrund)",100,0.,0.3,40,0.,10.);
  //SetHistoStandardSettings2(hSignal_pT_subtracted);

  TGaxis::SetMaxDigits(3);
  
  // read Cluster Tree
  TFile* fDaten = new TFile("pi0_vcal_data.root");
  

  const Int_t iPufferMax = 2;
  Int_t NMaxEvents = 1000;


  if(fDaten->IsZombie()){
    cout << "ERROR: Data File not found" << endl;
    return;
  }
  DataTree *Daten = new DataTree(fDaten);
  // NMaxEvents = Daten->GetNEvents()

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

    for (int i1 = 0; i1 < iCluster[iPufferAktuell]; i1++) { //-1 da man wenn man i1+1 benutzt man bereits bei -1 alle Combis durch hat.
    
      for(int i3 = 0; i3 !=i1 && i3 < iCluster[iPufferAktuell]; i3++){
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
      
  // Wechsle und Zeichne minv-pt same event
  cSignal_pT->cd();
  cSignal_pT->SetRightMargin(0.175);
  cSignal_pT->SetBottomMargin(0.125);
  gPad->SetLogz();
      
  hSignal_pT->GetZaxis()->SetRangeUser(1.e0,1.e4);
  hSignal_pT->Draw("colz");
  
  // Wechsle und Zeichne minv different event 
  cSignalmix->cd();
  cSignalmix->SetTopMargin(0.075);
  hSignalmix->Draw("");
          
  // Wechsle und Zeichne minv-pt different event
  cSignalmix_pT->cd();
  cSignalmix_pT->SetRightMargin(0.175);
  cSignalmix_pT->SetBottomMargin(0.125);
  gPad->SetLogz();
         
  hSignalmix_pT->GetZaxis()->SetRangeUser(1.e0,1.e4);
  hSignalmix_pT->Draw("colz");
  
  
  cSignal_subtracted->cd();
  cSignal_subtracted->SetTopMargin(0.075);
  
  TH1F *hSignal_subtracted = (TH1F*) hSignal->Clone();
  
  hSignal_subtracted->Add(hSignalmix,-1);
  hSignal_subtracted->Draw("");
  
  
  cSignal_pT_subtracted->cd();
  cSignal_pT_subtracted->SetRightMargin(0.175);
  cSignal_pT_subtracted->SetBottomMargin(0.125);
  cSignal_pT_subtracted->SetLogz();
  
  //TH2F *hSignal_pT_subtracted = (TH2F*) hSignal_pT->Clone();
  //TH1D *hSignal_pT_projection = hSignal_pT_subtracted->ProjectionX();
  //TH1D *hSignal_pT_projection_mixed = hSignalmix_pT->ProjectionX();
  
  //hSignal_pT_subtracted->Add(hSignal_pT_projection_mixed,-1);
  
  //hSignal_pT_subtracted->Add(hSignalmix_pT,-1);
  //hSignal_pT_subtracted->GetZaxis()->SetRangeUser(1.e-10,1.e1);
  //hSignal_pT_subtracted->Draw("colz");
  
  cSignal->SaveAs(Form("Simulation/InavrianteMasseSameEvent%s.png", AddName.Data()));
  cSignal_pT->SaveAs(Form("Simulation/InvarianteMasseTransversalImpulsSameEvent%s.png", AddName.Data()));
  cSignalmix->SaveAs(Form("Simulation/InvarianteMasseDifferentEvents%s.png", AddName.Data()));
  cSignalmix_pT->SaveAs(Form("Simulation/InvarianteMasseTransversalImpulsDifferentEvents%s.png", AddName.Data()));
  cSignal_subtracted->SaveAs(Form("Simulation/InvarianteMasseSubtracted%s.png", AddName.Data()));
  cSignal_pT_subtracted->SaveAs(Form("Simulation/InvarianteMasseTransversalImpulsSubtracted%s.png", AddName.Data()));
  
}









