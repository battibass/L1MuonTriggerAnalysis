// ******************************
// Efficiency Plotter
// ( Provide plots for efficiency and compute efficiency using
// TEfficiency from ROOT itself )
// ******************************

class EfficiencyPlotter
{

public :

  EfficiencyPlotter(TFile* file, 
		    std::string name = "", float minPt = 0) : outFile_(file) , 
							      baseName_(name),
							      minPt_(minPt) { config(); };   
  ~EfficiencyPlotter();

  void fill(triggeredMuonsIt & muon, 
	    L1Analysis::L1AnalysisEventDataFormat * event,
	    L1Analysis::L1AnalysisRecoVertexDataFormat * vtx);
  
  void config();
  void plotAndSave();
  void save() { outFile_->Write(); return ; };
  void printHisto(TCanvas * canvas, std::string tag) {  
    canvas->SaveAs((std::string("plots/")+baseName_+"/"+baseName_+tag+".png").c_str());
    return;
  };
  float minMuPt() { return minPt_+10; };
  

protected :

  TFile* outFile_;
  std::string baseName_;
  float minPt_;

public :

  histoMap histos_;
  hTH1Map  hTH1F_;
  
};


EfficiencyPlotter::~EfficiencyPlotter()
{

    histos_.clear();

};

void EfficiencyPlotter::config() 
{ 

  if (baseName_!="")
    {
      outFile_->mkdir(baseName_.c_str());
      outFile_->cd(baseName_.c_str());
      system(string("mkdir -p plots/" + baseName_).c_str());
    }

    string name  = (baseName_ + "_hGmtPtVsQual");
    string title = name + ";tight muon p_{t} [GeV/c];efficiency";
    hTH1F_["hGmtPtVsQual"] = new TH2F(name.c_str(),name.c_str(),60,0.5,60.5,7,0.5,7.5);

    const char *etaTags[7] = { "", "Fine", "Coarse", "1", "2", "3", "4" };

    name  = baseName_ + "_hEffVsPt";
    histos_["hEffVsPt"] = new TEfficiency(name.c_str(),title.c_str(),60,0.5,60.5);

    name  = baseName_ + "_hEffVsEta";
    title = name + ";tight muon #eta;efficiency";
    histos_["hEffVsEta"] = new TEfficiency(name.c_str(),title.c_str(),56,-1.05,1.05);

	name  = (baseName_ + "_hGmtEtaVsQual" + tag.str());
	hTH1F_["hGmtEtaVsQual" + tag.str()] = new TH2F(name.c_str(),name.c_str(),56,-1.05,1.05,7,0.5,7.5);

	name  = (baseName_ + "_hGmtEtaVsChamb" + tag.str());
	hTH1F_["hGmtEtaVsChamb" + tag.str()] = new TH2F(name.c_str(),name.c_str(),56,-1.05,1.05,11,0.5,11.5);

	name  = (baseName_ + "_hGmtBinEta" + tag.str());
	hTH1F_["hGmtBinEta" + tag.str()] = new TH1F(name.c_str(),name.c_str(),64,-1.2,1.2); // check Gmt binning in GMT DQM

	name  = (baseName_ + "_hGmtBinEtaVsQual" + tag.str());
	hTH1F_["hGmtBinEtaVsQual" + tag.str()] = new TH2F(name.c_str(),name.c_str(),64,-1.2,1.2,7,0.5,7.5); // check Gmt binning in GMT DQM

	name  = (baseName_ + "_hGmtBinEtaVsChamb" + tag.str());
	hTH1F_["hGmtBinEtaVsChamb" + tag.str()] = new TH2F(name.c_str(),name.c_str(),64,-1.2,1.2,11,0.5,11.5); // check Gmt binning in GMT DQM

      }

    name  = (baseName_ + "_hEffVsPhi");
    title = name + ";tight muon #phi [rad];efficiency";
    histos_["hEffVsPhi"] = new TEfficiency(name.c_str(),title.c_str(),48,-TMath::Pi(),TMath::Pi());

    name  = (baseName_ + "_hGmtPhiVsQual");
    hTH1F_["hGmtPhiVsQual"] = new TH2F(name.c_str(),name.c_str(),48,-TMath::Pi(),TMath::Pi(),7,0.5,7.5);

    name  = (baseName_ + "_hGmtPhiVsChamb");
    hTH1F_["hGmtPhiVsChamb"] = new TH2F(name.c_str(),name.c_str(),48,-TMath::Pi(),TMath::Pi(),11,0.5,11.5);

    name  = (baseName_ + "_hEffPhiVsEta");
    title = name + ";tight muon #phi [rad];tight muon #eta;efficiency";
    histos_["hEffPhiVsEta"] = new TEfficiency(name.c_str(),title.c_str(),48,-TMath::Pi(),TMath::Pi(),56,-1.05,1.05);

    name  = (baseName_ + "_hEffVsLumi");
    title = name + ";inst luminosity [10^{30}cm^{-2}s^{-1}];efficiency";
    histos_["hEffVsLumi"] = new TEfficiency(name.c_str(),title.c_str(),50,2000.,7500.);

    name  = (baseName_ + "_hGmtLumiVsQual");
    hTH1F_["hGmtLumiVsQual"] = new TH2F(name.c_str(),name.c_str(),50,2000.,7500.,7,0.5,7.5);

    name  = (baseName_ + "_hEffVsVtx");
    title = name + ";N. reco vtx ;efficiency";
    histos_["hEffVsVtx"] = new TEfficiency(name.c_str(),title.c_str(),25,0.,50.);

    name  = (baseName_ + "_hGmtVtxVsQual");
    hTH1F_["hGmtVtxVsQual"] = new TH2F(name.c_str(),name.c_str(),25,0.,50.,7,0.5,7.5);


}

void EfficiencyPlotter::fill(triggeredMuonsIt & muon, 
			     L1Analysis::L1AnalysisEventDataFormat * event,
			     L1Analysis::L1AnalysisRecoVertexDataFormat * vtx) 
{

  // CB find a better palce for this
  qualToFine chambToIndex; 

  chambToIndex[qualPair(1,4)] = 1;
  chambToIndex[qualPair(2,4)] = 2;
  chambToIndex[qualPair(2,3)] = 3;
  chambToIndex[qualPair(3,4)] = 4;
  chambToIndex[qualPair(3,3)] = 5;
  chambToIndex[qualPair(3,2)] = 6;
  chambToIndex[qualPair(4,4)] = 7;
  chambToIndex[qualPair(5,4)] = 8;
  chambToIndex[qualPair(6,4)] = 9;
  chambToIndex[qualPair(6,3)] = 10;
  chambToIndex[qualPair(7,4)] = 11;
  
  float pt     = muon->mu_->pt.at(muon->imu_);
  float eta    = muon->mu_->eta.at(muon->imu_);
  float phi    = muon->mu_->phi.at(muon->imu_);
  int   phiBin = (int)((phi+TMath::Pi())/(TMath::Pi()/24.))%4+1;

  int  nVtx  = vtx ? vtx->nVtx : 0;
  float lumi = event->instLumi;

  stringstream phiTagBin;
  phiTagBin << phiBin;
  std::string phiTag = phiTagBin.str();

  bool hasTrigger = muon->hasDtTriggerMatch() && 
  ((muon->gmt()->Ptdt.at(muon->igmt())) +0.01 > minPt_);

  //  float qualDt = muon->DTTriggerQual();

  float dttfEta     = 999.;
  bool  dttfFineEta = false;
  int dttfQual      = 999;
  int dttfChambId   = 999;
  
  if (muon->hasDtTriggerMatch()) 
    {
      dttfEta  = muon->gmt()->Etadt.at(muon->igmt());
      dttfQual = muon->gmt()->Qualdt.at(muon->igmt());
      dttfFineEta = muon->gmt()->FineEtadt.at(muon->igmt()) ? false : true ;
      
      int dtOuterCh = 0;
      for (int iSt=1; iSt<=4; ++iSt)
	if ( muon->gmtDttf_->prims_.ScList.at(iSt-1) >= 0)
	  dtOuterCh = iSt; 
      
      dttfChambId = chambToIndex.find(qualPair(dttfQual,dtOuterCh))->second;
    }
	
  if (fabs(eta) < MAX_MU_ETA)
    {
      histos_["hEffVsPt"]->Fill(hasTrigger,pt);
      histos_["hEffVsPtFine"]->Fill(hasTrigger&&dttfFineEta,pt);
      histos_["hEffVsPtCoarse"]->Fill(hasTrigger&&(!dttfFineEta),pt);

      if (hasTrigger) 
	{
	  hTH1F_["hGmtPtVsQual"]->Fill(pt,dttfQual);
	  hTH1F_["hGmtPtVsChamb"]->Fill(pt,dttfChambId);
	}
    }  

  if (pt > minMuPt())
    {
      histos_["hEffVsEta"]->Fill(hasTrigger,eta);

      histos_["hEffVsEta"+phiTag]->Fill(hasTrigger,eta);
 
      histos_["hEffVsEtaFine"]->Fill(hasTrigger&&dttfFineEta,eta);
      histos_["hEffVsEtaCoarse"]->Fill(hasTrigger&&(!dttfFineEta),eta);

      histos_["hEffPhiVsEta"]->Fill(hasTrigger,phi,eta);

      if (hasTrigger) 
	{
	  hTH1F_["hGmtBinEta"]->Fill(dttfEta);
	  hTH1F_["hGmtEtaVsQual"]->Fill(eta,dttfQual);
	  hTH1F_["hGmtEtaVsChamb"]->Fill(eta,dttfChambId);
	  hTH1F_["hGmtBinEtaVsQual"]->Fill(dttfEta,dttfQual);
	  hTH1F_["hGmtBinEtaVsChamb"]->Fill(dttfEta,dttfChambId);

	  hTH1F_["hGmtBinEta"+phiTag]->Fill(dttfEta);
          hTH1F_["hGmtEtaVsQual"+phiTag]->Fill(eta,dttfQual);
          hTH1F_["hGmtEtaVsChamb"+phiTag]->Fill(eta,dttfChambId);
          hTH1F_["hGmtBinEtaVsQual"+phiTag]->Fill(dttfEta,dttfQual);
	  hTH1F_["hGmtBinEtaVsChamb"+phiTag]->Fill(dttfEta,dttfChambId);

	  if( dttfFineEta ) 
	    {
	      hTH1F_["hGmtBinEtaFine"]->Fill(dttfEta);
	      hTH1F_["hGmtEtaVsQualFine"]->Fill(eta,dttfQual);
	      hTH1F_["hGmtBinEtaVsQualFine"]->Fill(dttfEta,dttfQual);
	      hTH1F_["hGmtEtaVsChambFine"]->Fill(eta,dttfChambId);
	      hTH1F_["hGmtBinEtaVsChambFine"]->Fill(dttfEta,dttfChambId);
	    }
	  else 
	    {
	      hTH1F_["hGmtBinEtaCoarse"]->Fill(dttfEta);
	      hTH1F_["hGmtEtaVsQualCoarse"]->Fill(eta,dttfQual);
	      hTH1F_["hGmtBinEtaVsQualCoarse"]->Fill(dttfEta,dttfQual);
	      hTH1F_["hGmtEtaVsChambCoarse"]->Fill(eta,dttfChambId);
	      hTH1F_["hGmtBinEtaVsChambCoarse"]->Fill(dttfEta,dttfChambId);
	    }
	}
      
    }

  if (pt > minMuPt() && fabs(eta) < MAX_MU_ETA)
    {
      histos_["hEffVsPhi"]->Fill(hasTrigger,phi);
      if (hasTrigger) 
	{
	  hTH1F_["hGmtPhiVsQual"]->Fill(phi,dttfQual);
	  hTH1F_["hGmtPhiVsChamb"]->Fill(phi,dttfChambId);
	}

      histos_["hEffVsLumi"]->Fill(hasTrigger,lumi);
      if (hasTrigger) 
	{
	  hTH1F_["hGmtLumiVsQual"]->Fill(lumi,dttfQual);
	}

      histos_["hEffVsVtx"]->Fill(hasTrigger,nVtx);
      if (hasTrigger) 
	{
	  hTH1F_["hGmtVtxVsQual"]->Fill(nVtx,dttfQual);
	}
    }  


  return;

}

void EfficiencyPlotter::plotAndSave() 
{ 

  setTDRStyle();
  gStyle->SetOptTitle(0);

  TCanvas *cEffvsPt = new TCanvas((baseName_+"cEffvsPt").c_str(),
				  (baseName_+"cEffvsPt").c_str(),500,500);
  cEffvsPt->cd();
  cEffvsPt->SetGrid();  

  histos_["hEffVsPt"]->Draw();

    TH1 const *hEffVsPtTotal = histos_["hEffVsPt"]->GetTotalHistogram(); 
  
  THStack * hEffvsPtStack = new THStack("hEffvsPtStack", "hEffvsPtStack");
  for(int ybin=1; ybin<=hTH1F_["hGmtPtVsQual"]->GetNbinsY(); ++ybin)
    {
      stringstream px ;
      px<<ybin<<"QualityPt";
      TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtPtVsQual"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

      projection->SetFillColor(colorMap[ybin-1]);
      projection->SetLineColor(colorMap[ybin-1]);
      projection->Divide(hEffVsPtTotal);
      hEffvsPtStack->Add(projection);
    }

  hEffvsPtStack->Draw("same");
  histos_["hEffVsPt"]->Draw("same");
  
  cEffvsPt->Update();
  histos_["hEffVsPt"]->GetPaintedGraph()->GetXaxis()->SetRangeUser(0.,60.);
  histos_["hEffVsPt"]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.,1.1);
  histos_["hEffVsPt"]->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.04);
  histos_["hEffVsPt"]->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.04);

  printHisto(cEffvsPt,"EffvsPt");


  TCanvas *cEffvsPtChamb = new TCanvas((baseName_+"cEffvsPtChamb").c_str(),
				       (baseName_+"cEffvsPtChamb").c_str(),500,500);
  cEffvsPtChamb->cd();
  cEffvsPtChamb->SetGrid();  

  histos_["hEffVsPt"]->Draw();

  THStack * hEffvsPtChambStack = new THStack("hEffvsPtChambStack", "hEffvsPtChambStack");
  for(int ybin=1; ybin<=hTH1F_["hGmtPtVsChamb"]->GetNbinsY(); ++ybin)
    {
      stringstream px ;
      px<<ybin<<"Chamb";
      TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtPtVsChamb"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

      projection->SetFillColor(colorMapFine[ybin-1]);
      projection->SetLineColor(colorMapFine[ybin-1]);
      projection->Divide(hEffVsPtTotal);
      hEffvsPtChambStack->Add(projection);
    }

  hEffvsPtChambStack->Draw("same");
  histos_["hEffVsPt"]->Draw("same");
  
  cEffvsPtChamb->Update();
  histos_["hEffVsPt"]->GetPaintedGraph()->GetXaxis()->SetRangeUser(0.,60.);
  histos_["hEffVsPt"]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.,1.1);
  histos_["hEffVsPt"]->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.04);
  histos_["hEffVsPt"]->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.04);

  printHisto(cEffvsPtChamb,"EffvsPtChamb");


  const char *etaTags[7] = { "", "Fine", "Coarse", "1", "2", "3", "4" };

  for (int iTag=0; iTag<7; ++ iTag) 
    {
      
      stringstream tag; tag << etaTags[iTag];
      
      TCanvas *cEffvsEta = new TCanvas((baseName_+"cEffvsEta"+tag.str()).c_str(),
				       (baseName_+"cEffvsEta"+tag.str()).c_str(),500,500);
      cEffvsEta->cd();
      cEffvsEta->SetGrid();  

      histos_["hEffVsEta"+tag.str()]->Draw();

      TH1 const *hEffVsEtaTotal = histos_["hEffVsEta"+tag.str()]->GetTotalHistogram(); 

      THStack * hEffvsEtaStack = new THStack(("hEffvsEtaStack"+tag.str()).c_str(), 
					     ("hEffvsEtaStack"+tag.str()).c_str());
      for(int ybin=1; ybin<=hTH1F_["hGmtEtaVsQual"+tag.str()]->GetNbinsY(); ++ybin)
      {
        stringstream px ;
        px<<ybin<<"QualityEta"<<etaTags[iTag];
        TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtEtaVsQual"+tag.str()])->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMap[ybin-1]);
        projection->SetLineColor(colorMap[ybin-1]);
	projection->Divide(hEffVsEtaTotal);
	hEffvsEtaStack->Add(projection);
      }

      hEffvsEtaStack->Draw("same");
      histos_["hEffVsEta"+tag.str()]->Draw("same");

      cEffvsEta->Update();
      histos_["hEffVsEta"+tag.str()]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.1);
      histos_["hEffVsEta"+tag.str()]->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.04);
      histos_["hEffVsEta"+tag.str()]->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.04);

      printHisto(cEffvsEta,"EffvsEta"+tag.str());


      TCanvas *cEffvsEtaGmt = new TCanvas((baseName_+"cEffvsEtaGmt"+tag.str()).c_str(),
					   (baseName_+"cEffvsEtaGmt"+tag.str()).c_str(),500,500);
      cEffvsEtaGmt->cd();
      cEffvsEtaGmt->SetGrid();  

      TH1 *hGmtBinEta = hTH1F_["hGmtBinEta"+tag.str()]; 

      THStack * hEffvsGmtEtaStack = new THStack(("hEffvsGmtEtaStack"+tag.str()).c_str(), 
						 ("hEffGmtvsEtaStack"+tag.str()).c_str());
      for(int ybin=1; ybin<=hTH1F_["hGmtBinEtaVsQual"+tag.str()]->GetNbinsY(); ++ybin)
	{
	  stringstream px ;
	  px<<ybin<<"QualityGmtEta"<<etaTags[iTag];
	  TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtBinEtaVsQual"+tag.str()]) ->ProjectionX(px.str().c_str(),ybin,ybin);
	  
	  projection->SetFillColor(colorMap[ybin-1]);
	  projection->SetLineColor(colorMap[ybin-1]);
	  projection->Divide(hGmtBinEta);
	  hEffvsGmtEtaStack->Add(projection);
	}
      
      hEffvsGmtEtaStack->Draw();
      hEffvsGmtEtaStack->GetXaxis()->SetTitle("Gmt #eta");
      hEffvsGmtEtaStack->GetYaxis()->SetTitle("quality contrib. to efficiency");
      hEffvsGmtEtaStack->GetXaxis()->SetTitleSize(0.04);
      hEffvsGmtEtaStack->GetYaxis()->SetTitleSize(0.04);
      
      
      printHisto(cEffvsEtaGmt,"EffvsGmtEta"+tag.str());

      TCanvas *cEffvsEtaChamb = new TCanvas((baseName_+"cEffvsEtaChamb"+tag.str()).c_str(),
					    (baseName_+"cEffvsEtaChamb"+tag.str()).c_str(),500,500);
      cEffvsEtaChamb->cd();
      cEffvsEtaChamb->SetGrid();  

      histos_["hEffVsEta"+tag.str()]->Draw();

      THStack * hEffvsEtaChambStack = new THStack(("hEffvsEtaChambStack"+tag.str()).c_str(), 
						  ("hEffvsEtaChambStack"+tag.str()).c_str());
      for(int ybin=1; ybin<=hTH1F_["hGmtEtaVsChamb"+tag.str()]->GetNbinsY(); ++ybin)
	{
	  stringstream px ;
	  px<<ybin<<"ChambEta"<<etaTags[iTag];
	  TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtEtaVsChamb"+tag.str()])->ProjectionX(px.str().c_str(),ybin,ybin);

	  projection->SetFillColor(colorMapFine[ybin-1]);
	  projection->SetLineColor(colorMapFine[ybin-1]);
	  projection->Divide(hEffVsEtaTotal);
	  hEffvsEtaChambStack->Add(projection);
	}
      
      hEffvsEtaChambStack->Draw("same");
      histos_["hEffVsEta"+tag.str()]->Draw("same");
      
      cEffvsEtaChamb->Update();
      histos_["hEffVsEta"+tag.str()]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.1);
      histos_["hEffVsEta"+tag.str()]->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.04);
      histos_["hEffVsEta"+tag.str()]->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.04);
      
      printHisto(cEffvsEtaChamb,"EffvsEtaChamb"+tag.str());
      
      
      TCanvas *cEffvsEtaChambGmt = new TCanvas((baseName_+"cEffvsEtaChambGmt"+tag.str()).c_str(),
						(baseName_+"cEffvsEtaChambGmt"+tag.str()).c_str(),500,500);
      cEffvsEtaChambGmt->cd();
      cEffvsEtaChambGmt->SetGrid();  
      
      THStack * hEffvsGmtEtaChambStack = new THStack(("hEffvsGmtEtaChambStack"+tag.str()).c_str(), 
						      ("hEffGmtvsEtaChambStack"+tag.str()).c_str());
      for(int ybin=1; ybin<=hTH1F_["hGmtBinEtaVsChamb"+tag.str()]->GetNbinsY(); ++ybin)
	{
	  stringstream px ;
	  px<<ybin<<"ChambGmtEta"<<etaTags[iTag];
	  TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtBinEtaVsChamb"+tag.str()]) ->ProjectionX(px.str().c_str(),ybin,ybin);
	  
	  projection->SetFillColor(colorMapFine[ybin-1]);
	  projection->SetLineColor(colorMapFine[ybin-1]);
	  projection->Divide(hGmtBinEta);
	  hEffvsGmtEtaChambStack->Add(projection);
	}
      
      hEffvsGmtEtaChambStack->Draw();
      hEffvsGmtEtaChambStack->GetXaxis()->SetTitle("Gmt #eta");
      hEffvsGmtEtaChambStack->GetYaxis()->SetTitle("quality contrib. to efficiency");
      hEffvsGmtEtaChambStack->GetXaxis()->SetTitleSize(0.04);
      hEffvsGmtEtaChambStack->GetYaxis()->SetTitleSize(0.04);
      
      
      printHisto(cEffvsEtaChambGmt,"EffvsGmtEtaChamb"+tag.str());
      
    }

    TCanvas *cEffvsPhi = new TCanvas((baseName_+"cEffvsPhi").c_str(),
				     (baseName_+"cEffvsPhi").c_str(),500,500);
    cEffvsPhi->cd();
    cEffvsPhi->SetGrid();  

    histos_["hEffVsPhi"]->Draw();

    TH1 const *hEffVsPhiTotal = histos_["hEffVsPhi"]->GetTotalHistogram(); 

    THStack * hEffvsPhiStack = new THStack("hEffvsPhiStack", "hEffvsPhiStack");
    for(int ybin=1; ybin<=hTH1F_["hGmtPhiVsQual"]->GetNbinsY(); ++ybin)
      {
        stringstream px ;
        px<<ybin<<"QualityPhi";
        TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtPhiVsQual"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMap[ybin-1]);
        projection->SetLineColor(colorMap[ybin-1]);
	projection->Divide(hEffVsPhiTotal);
	hEffvsPhiStack->Add(projection);
      }

    hEffvsPhiStack->Draw("same");
    histos_["hEffVsPhi"]->Draw("same");

    cEffvsPhi->Update();
    histos_["hEffVsPhi"]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.,1.1);
    histos_["hEffVsPhi"]->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.04);
    histos_["hEffVsPhi"]->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.04);

    printHisto(cEffvsPhi,"EffvsPhi");


    TCanvas *cEffvsPhiChamb = new TCanvas((baseName_+"cEffvsPhiChamb").c_str(),
					  (baseName_+"cEffvsPhiChamb").c_str(),500,500);
    cEffvsPhiChamb->cd();
    cEffvsPhiChamb->SetGrid();  

    histos_["hEffVsPhi"]->Draw();

    THStack * hEffvsPhiChambStack = new THStack("hEffvsPhiChambStack", "hEffvsPhiChambStack");
    for(int ybin=1; ybin<=hTH1F_["hGmtPhiVsChamb"]->GetNbinsY(); ++ybin)
      {
        stringstream px ;
        px<<ybin<<"ChambPhi";
        TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtPhiVsChamb"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMapFine[ybin-1]);
        projection->SetLineColor(colorMapFine[ybin-1]);
	projection->Divide(hEffVsPhiTotal);
	hEffvsPhiChambStack->Add(projection);
      }

    hEffvsPhiChambStack->Draw("same");
    histos_["hEffVsPhi"]->Draw("same");

    cEffvsPhiChamb->Update();
    histos_["hEffVsPhi"]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.,1.1);
    histos_["hEffVsPhi"]->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.04);
    histos_["hEffVsPhi"]->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.04);

    printHisto(cEffvsPhiChamb,"EffvsPhiChamb");


    TCanvas *cEffPhivsEta = new TCanvas((baseName_+"cEffPhivsEta").c_str(),
					(baseName_+"cEffPhivsEta").c_str(),500,500);
    cEffPhivsEta->cd();
    cEffPhivsEta->SetGrid();  

    histos_["hEffPhiVsEta"]->Draw("colz");
    //    std::cout << histos_["hEffPhiVsEta"]->GetPaintedGraph() << std::endl;
//     histos_["hEffPhiVsEta"]->GetPaintedHistogram()->SetMaximum(1.);
//     histos_["hEffPhiVsEta"]->GetPaintedHistogram()->GetXaxis()->SetTitleSize(0.04);
//     histos_["hEffPhiVsEta"]->GetPaintedHistogram()->GetYaxis()->SetTitleSize(0.04);

    printHisto(cEffPhivsEta,"EffPhivsEta");


    TCanvas *cEffvsVtx = new TCanvas((baseName_+"cEffvsVtx").c_str(),
				     (baseName_+"cEffvsVtx").c_str(),500,500);
    cEffvsVtx->cd();
    cEffvsVtx->SetGrid();  
    
    histos_["hEffVsVtx"]->Draw();

    TH1 const *hEffVsVtxTotal = histos_["hEffVsVtx"]->GetTotalHistogram(); 
  
    THStack * hEffvsVtxStack = new THStack("hEffvsVtxStack", "hEffvsVtxStack");
    for(int ybin=1; ybin<=hTH1F_["hGmtVtxVsQual"]->GetNbinsY(); ++ybin)
      {
	stringstream px ;
	px<<ybin<<"QualityVtx";
	TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtVtxVsQual"]) ->ProjectionX(px.str().c_str(),ybin,ybin);
	
	projection->SetFillColor(colorMap[ybin-1]);
	projection->SetLineColor(colorMap[ybin-1]);
	projection->Divide(hEffVsVtxTotal);
	hEffvsVtxStack->Add(projection);
      }
    
    hEffvsVtxStack->Draw("same");
    histos_["hEffVsVtx"]->Draw("same");
    
    cEffvsVtx->Update();
    histos_["hEffVsVtx"]->GetPaintedGraph()->GetXaxis()->SetRangeUser(5.,30.);
    histos_["hEffVsVtx"]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.,1.1);
    histos_["hEffVsVtx"]->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.04);
    histos_["hEffVsVtx"]->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.04);
    
    printHisto(cEffvsVtx,"EffvsVtx");


    TCanvas *cEffvsLumi = new TCanvas((baseName_+"cEffvsLumi").c_str(),
				      (baseName_+"cEffvsLumi").c_str(),500,500);
    cEffvsLumi->cd();
    cEffvsLumi->SetGrid();  
    
     histos_["hEffVsLumi"]->Draw();

    TH1 const *hEffVsLumiTotal = histos_["hEffVsLumi"]->GetTotalHistogram(); 
  
    THStack * hEffvsLumiStack = new THStack("hEffvsLumiStack", "hEffvsLumiStack");
    for(int ybin=1; ybin<=hTH1F_["hGmtLumiVsQual"]->GetNbinsY(); ++ybin)
      {
	stringstream px ;
	px<<ybin<<"QualityLumi";
	TH1D* projection = static_cast<TH2*>(hTH1F_["hGmtLumiVsQual"]) ->ProjectionX(px.str().c_str(),ybin,ybin);
	
	projection->SetFillColor(colorMap[ybin-1]);
	projection->SetLineColor(colorMap[ybin-1]);
	projection->Divide(hEffVsLumiTotal);
	hEffvsLumiStack->Add(projection);
      }
    
    hEffvsLumiStack->Draw("same");
    histos_["hEffVsLumi"]->Draw("same");
    
    cEffvsLumi->Update();
    histos_["hEffVsLumi"]->GetPaintedGraph()->GetXaxis()->SetRangeUser(2000.,6000.);
    histos_["hEffVsLumi"]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.,1.1);
    histos_["hEffVsLumi"]->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.04);
    histos_["hEffVsLumi"]->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.04);
    
    printHisto(cEffvsLumi,"EffvsLumi");
    
    save();

    return;

}

//////////////////////////////////////////////////////////
// Helper classfunction to plot togheter efficiency plots
//////////////////////////////////////////////////////////

void plotAndSaveAll(std::vector<EfficiencyPlotter*> & plotters, std::string plotName)
{
  
  TCanvas *canvas = new TCanvas(("c"+plotName).c_str(),
				("c"+plotName).c_str(),500,500);
  canvas->cd();
  canvas->SetGrid();  
  
  std::vector<EfficiencyPlotter*>::const_iterator plotterIt  = plotters.begin();
  std::vector<EfficiencyPlotter*>::const_iterator plotterEnd = plotters.end();

  int iPlot = 0;
  TF1 *fitFunc = 0;
  for(;plotterIt!=plotterEnd;++plotterIt,++iPlot)
    {
      TEfficiency * plot = (*plotterIt)->histos_["h"+plotName];
      plot->SetLineColor(iPlot+1);
      plot->Draw(iPlot == 0 ? "" : "same");

      canvas->Update();

      if (plotName.find("Pt") == string::npos) 
	plot->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.8,1.1);

      if (iPlot == 0 && plotName.find("Lumi") != string::npos)
	{ 
	  plot->GetPaintedGraph()->Fit("pol1");
	  fitFunc = plot->GetPaintedGraph()->GetFunction("pol1");
	}
    }

  if (fitFunc) 
    {
      fitFunc->Draw("same");
      canvas->Update();
    }
    
  system(string("mkdir -p plots/All/").c_str());
  canvas->SaveAs(("plots/All/All"+plotName+".png").c_str());
  
  return;

}

