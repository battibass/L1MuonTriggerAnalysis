// ******************************
// Control Plotter
// ( Class where to add all control plots to be used during tight
// selection, hlt matching and tag selection, before the eff is 
// actually computed )
// ******************************

class ControlPlotter
{

public :

  ControlPlotter(TFile* file, std::string name = "") : outFile_(file) , 
						       baseName_(name) { config(); };   
  ~ControlPlotter();

  void fillTight(TriggeredMuon & muon );  
  void fillTrigger(TriggeredMuon & muon );
  
  void config();
  void plotAndSave();

  void save() { outFile_->Write(); return ; };
  void printHisto(TCanvas * canvas, std::string tag) 
  {  
    canvas->SaveAs((std::string("plots/")+baseName_+"/"+baseName_+tag+".png").c_str());
    return;
  };
  

protected :

  TFile* outFile_;
  std::string baseName_;

public :

  hTH1Map  hTH1_;
  
};


ControlPlotter::~ControlPlotter()
{
    hTH1_.clear();
};

void ControlPlotter::config() 
{ 
  if (baseName_!="")
    {
      outFile_->mkdir(baseName_.c_str());
      outFile_->cd(baseName_.c_str());
      system(string("mkdir -p plots/" + baseName_).c_str());
    }

  string name  = (baseName_ + "_hTrkHits");
  hTH1_["hTrkHits"] = new TH1F(name.c_str(),name.c_str(),20,0.5,20.5);

  name  = (baseName_ + "_hChi2");
  hTH1_["hChi2"] = new TH1F(name.c_str(),name.c_str(),20,0.5,20.5);

  name  = (baseName_ + "_hGmtDeltaPhi");
  hTH1_["hGmtDeltaPhi"] = new TH1F(name.c_str(),name.c_str(),100,0.,TMath::Pi());

  name  = (baseName_ + "_hGmtDeltaEta");
  hTH1_["hGmtDeltaEta"] = new TH1F(name.c_str(),name.c_str(),100,-1.,1.);

  name  = (baseName_ + "_hGmtDeltaPhiVsEta");
  hTH1_["hGmtDeltaPhiVsEta"] = new TH2F(name.c_str(),name.c_str(),24,-2.4,2.4,100,0.,TMath::Pi());

  name  = (baseName_ + "_hGmtDeltaEtaVsEta");
  hTH1_["hGmtDeltaEtaVsEta"] = new TH2F(name.c_str(),name.c_str(),24,-2.4,2.4,100,-1.,1.);
  
}

void ControlPlotter::fillTight(TriggeredMuon & muon) 
{
  
  float nTrkHits = muon.mu_->tr_validhits.at(muon.imu_);
  float chi2     = muon.mu_->normchi2.at(muon.imu_);

  hTH1_["hTrkHits"]->Fill(nTrkHits);
  hTH1_["hChi2"]->Fill(chi2);

}

void ControlPlotter::fillTrigger(TriggeredMuon & muon) 
{
  
  float dPhi  = muon.deltaPhi();
  float dEta  = muon.deltaEta();

  float eta   = muon.my_gmt->eta(muon.my_igmt);

  hTH1_["hGmtDeltaPhi"]->Fill(dPhi);
  hTH1_["hGmtDeltaEta"]->Fill(dEta);

  hTH1_["hGmtDeltaPhi"]->Fill(eta,dPhi);
  hTH1_["hGmtDeltaEta"]->Fill(eta,dEta);

}

void ControlPlotter::plotAndSave() 
{ 

  setTDRStyle();
  gStyle->SetOptTitle(0);

  TCanvas *cTrkHits = new TCanvas((baseName_+"cTrkHits").c_str(),
				  (baseName_+"cTrkHits").c_str(),500,500);
  cTrkHits->cd();
  cTrkHits->SetGrid(); 
  cTrkHits->SetLogy();  

  hTH1_["hTrkHits"]->Draw();
  hTH1_["hTrkHits"]->GetXaxis()->SetTitleSize(0.04);
  hTH1_["hTrkHits"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cTrkHits,"nTrkHits");

  TCanvas *cChi2 = new TCanvas((baseName_+"cChi2").c_str(),
				  (baseName_+"cChi2").c_str(),500,500);
  cChi2->cd();
  cChi2->SetGrid();  
  cChi2->SetLogy();  

  hTH1_["hChi2"]->Draw();
  hTH1_["hChi2"]->GetXaxis()->SetTitleSize(0.04);
  hTH1_["hChi2"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cChi2,"nChi2");

  TCanvas *cGmtDeltaPhi = new TCanvas((baseName_+"cGmtDeltaPhi").c_str(),
				  (baseName_+"cGmtDeltaPhi").c_str(),500,500);
  cGmtDeltaPhi->cd();
  cGmtDeltaPhi->SetGrid();  

  hTH1_["hGmtDeltaPhi"]->Draw();
  hTH1_["hGmtDeltaPhi"]->GetXaxis()->SetTitleSize(0.04);
  hTH1_["hGmtDeltaPhi"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cGmtDeltaPhi,"nGmtDeltaPhi");

  TCanvas *cGmtDeltaEta = new TCanvas((baseName_+"cGmtDeltaEta").c_str(),
				  (baseName_+"cGmtDeltaEta").c_str(),500,500);
  cGmtDeltaEta->cd();
  cGmtDeltaEta->SetGrid();  

  hTH1_["hGmtDeltaEta"]->Draw();
  hTH1_["hGmtDeltaEta"]->GetXaxis()->SetTitleSize(0.04);
  hTH1_["hGmtDeltaEta"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cGmtDeltaEta,"nGmtDeltaEta");

  TCanvas *cGmtDeltaPhiVsEta = new TCanvas((baseName_+"cGmtDeltaPhiVsEta").c_str(),
				  (baseName_+"cGmtDeltaPhiVsEta").c_str(),500,500);
  cGmtDeltaPhiVsEta->cd();
  cGmtDeltaPhiVsEta->SetGrid();  

  hTH1_["hGmtDeltaPhiVsEta"]->Draw("colz");
  hTH1_["hGmtDeltaPhiVsEta"]->GetXaxis()->SetTitleSize(0.04);
  hTH1_["hGmtDeltaPhiVsEta"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cGmtDeltaPhiVsEta,"nGmtDeltaPhiVsEta");

  TCanvas *cGmtDeltaEtaVsEta = new TCanvas((baseName_+"cGmtDeltaEtaVsEta").c_str(),
				  (baseName_+"cGmtDeltaEtaVsEta").c_str(),500,500);
  cGmtDeltaEtaVsEta->cd();
  cGmtDeltaEtaVsEta->SetGrid();  

  hTH1_["hGmtDeltaEtaVsEta"]->Draw("colz");
  hTH1_["hGmtDeltaEtaVsEta"]->GetXaxis()->SetTitleSize(0.04);
  hTH1_["hGmtDeltaEtaVsEta"]->GetYaxis()->SetTitleSize(0.04);
  
  printHisto(cGmtDeltaEta,"nGmtDeltaEta");

  save();

  return;

}
