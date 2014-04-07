// ******************************
// Control Plotter
// ( Class where to add all control
//   plots to be used during tight
//   selection, hlt matching and 
//   tag selection, before the eff
//   is actually computed
// )
// ******************************

class ControlPlotter
{

public :

  ControlPlotter(TFile* file, std::string name = "") : outFile_(file) , 
						       baseName_(name) { config(); };   
  ~ControlPlotter();

  void fillTight(TriggeredMuon & muon );  
  void fillTrigger(TriggeredMuon & muon );
  void fillGmtDttf(GmtDttfMuon & gmtDttf); 
  
  void config();
  void plot();
  void save() { outFile_->Write(); return ; };
  void printHisto(TCanvas * canvas, std::string tag) {  
    canvas->SaveAs((std::string("plots/")+baseName_+"/"+baseName_+tag+".png").c_str());
    return;
  };
  

protected :

  TFile* outFile_;
  std::string baseName_;

public :

  hTH1Map  hTH1F_;
  
};


ControlPlotter::~ControlPlotter()
{

    hTH1F_.clear();

};

void ControlPlotter::config() 
{ 

  if (baseName_!="")
    {
      outFile_->mkdir(baseName_.c_str());
      outFile_->cd(baseName_.c_str());
      system(std::string("mkdir -p plots/" + baseName_).c_str());
    }

    std::string name  = (baseName_ + "_hTrkHits");
    hTH1F_["hTrkHits"] = new TH1F(name.c_str(),name.c_str(),20,0.5,20.5);

    name  = (baseName_ + "_hChi2");
    hTH1F_["hChi2"] = new TH1F(name.c_str(),name.c_str(),20,0.5,20.5);

    name  = (baseName_ + "_hDTTFDeltaPhi");
    hTH1F_["hDTTFDeltaPhi"] = new TH1F(name.c_str(),name.c_str(),100,0.,TMath::Pi());

    name  = (baseName_ + "_hDTTFDeltaEta");
    hTH1F_["hDTTFDeltaEta"] = new TH1F(name.c_str(),name.c_str(),100,-1.,1.);

    name  = (baseName_ + "_hPrimitivesVsQual");
    hTH1F_["hPrimitivesVsQual"] = new TH2F(name.c_str(),name.c_str(),5,-.5,4.5,8,-0.5,7.5);

    name  = (baseName_ + "_hPrimitiveQualVsChamb");
    hTH1F_["hPrimitiveQualvsChamb"] = new TH2F(name.c_str(),name.c_str(),8,-.5,7.5,4,0.5,4.5);

}

void ControlPlotter::fillTight(TriggeredMuon & muon) 
{
  
  float nTrkHits = muon.mu_->tr_validhits.at(muon.imu_);
  float chi2     = muon.mu_->normchi2.at(muon.imu_);

  hTH1F_["hTrkHits"]->Fill(nTrkHits);
  hTH1F_["hChi2"]->Fill(chi2);

}

void ControlPlotter::fillTrigger(TriggeredMuon & muon) 
{
  
  float dPhi  = muon.deltaPhi();
  float dEta  = muon.deltaEta();

  hTH1F_["hDTTFDeltaPhi"]->Fill(dPhi);
  hTH1F_["hDTTFDeltaEta"]->Fill(dEta);

}

void ControlPlotter::fillGmtDttf(GmtDttfMuon & gmtDttf) 
{

  if (gmtDttf.hasDtTriggerMatch())
    {
      int qual = gmtDttf.dttf_->trQual.at(gmtDttf.idttf_); 

      for (int iSt=1; iSt<=4; ++iSt)
	{
	  if ( gmtDttf.prims_.ScList.at(iSt-1) >= 0)
	    {
	      hTH1F_["hPrimitivesVsQual"]->Fill(iSt,qual);
	    }
	}
      
      std::vector<int>::const_iterator iPrimsIt  = gmtDttf.prims_.PrimitivesId.begin();
      std::vector<int>::const_iterator iPrimsEnd = gmtDttf.prims_.PrimitivesId.end();


      for(;iPrimsIt!=iPrimsEnd; ++iPrimsIt) {
	int primQual = gmtDttf.dttf_->phCode[(*iPrimsIt)];
	int chamb    = gmtDttf.dttf_->phSt[(*iPrimsIt)];
	hTH1F_["hPrimitiveQualvsChamb"]->Fill(primQual,chamb);
      }

    } 
  else
    {
    hTH1F_["hPrimitivesVsQual"]->Fill(0.,0.);
    }
 
}

void ControlPlotter::plot() 
{ 

  setTDRStyle();
  gStyle->SetOptTitle(0);

  TCanvas *cTrkHits = new TCanvas((baseName_+"cTrkHits").c_str(),
				  (baseName_+"cTrkHits").c_str(),500,500);
  cTrkHits->cd();
  cTrkHits->SetGrid(); 
  cTrkHits->SetLogy();  

  hTH1F_["hTrkHits"]->Draw();
  hTH1F_["hTrkHits"]->GetXaxis()->SetTitleSize(0.04);
  hTH1F_["hTrkHits"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cTrkHits,"nTrkHits");

  TCanvas *cChi2 = new TCanvas((baseName_+"cChi2").c_str(),
				  (baseName_+"cChi2").c_str(),500,500);
  cChi2->cd();
  cChi2->SetGrid();  
  cChi2->SetLogy();  

  hTH1F_["hChi2"]->Draw();
  hTH1F_["hChi2"]->GetXaxis()->SetTitleSize(0.04);
  hTH1F_["hChi2"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cChi2,"nChi2");

  TCanvas *cDTTFDeltaPhi = new TCanvas((baseName_+"cDTTFDeltaPhi").c_str(),
				  (baseName_+"cDTTFDeltaPhi").c_str(),500,500);
  cDTTFDeltaPhi->cd();
  cDTTFDeltaPhi->SetGrid();  

  hTH1F_["hDTTFDeltaPhi"]->Draw();
  hTH1F_["hDTTFDeltaPhi"]->GetXaxis()->SetTitleSize(0.04);
  hTH1F_["hDTTFDeltaPhi"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cDTTFDeltaPhi,"nDTTFDeltaPhi");

  TCanvas *cDTTFDeltaEta = new TCanvas((baseName_+"cDTTFDeltaEta").c_str(),
				  (baseName_+"cDTTFDeltaEta").c_str(),500,500);
  cDTTFDeltaEta->cd();
  cDTTFDeltaEta->SetGrid();  

  hTH1F_["hDTTFDeltaEta"]->Draw();
  hTH1F_["hDTTFDeltaEta"]->GetXaxis()->SetTitleSize(0.04);
  hTH1F_["hDTTFDeltaEta"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cDTTFDeltaEta,"nDTTFDeltaEta");

  TCanvas *cPrimitivesVsQual= new TCanvas((baseName_+"cPrimitivesVsQual").c_str(),
					  (baseName_+"cPrimitivesVsQual").c_str(),500,500);
  cPrimitivesVsQual->cd();
  cPrimitivesVsQual->SetGrid();  

  hTH1F_["hPrimitivesVsQual"]->Draw("colztext");
  hTH1F_["hPrimitivesVsQual"]->GetXaxis()->SetTitleSize(0.04);
  hTH1F_["hPrimitivesVsQual"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cPrimitivesVsQual,"hPrimitivesVsQual");

  TCanvas *cPrimitiveQualvsChamb= new TCanvas((baseName_+"cPrimitiveQualvsChamb").c_str(),
					      (baseName_+"cPrimitiveQualvsChamb").c_str(),500,500);
  cPrimitiveQualvsChamb->cd();
  cPrimitiveQualvsChamb->SetGrid();  

  hTH1F_["hPrimitiveQualvsChamb"]->Draw("colztext");
  hTH1F_["hPrimitiveQualvsChamb"]->GetXaxis()->SetTitleSize(0.04);
  hTH1F_["hPrimitiveQualvsChamb"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cPrimitiveQualvsChamb,"hPrimitiveQualvsChamb");


  save();

  return;

}
