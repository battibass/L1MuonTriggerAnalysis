#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "THStack.h"
#include "TColor.h"

#include "TStyle.h"

#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <map>

using namespace std;

// ********************************
// (1)
// Basic plotter class 
// (virtual not possible to make
// an  instance of this class)
// ********************************

class Plotter
{

    public :
        Plotter() : outFile_(0) {};   
        ~Plotter() {};

        virtual void fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
			  int igmt, int idt, int irpcb, int icsc, 
			  int irpcf, TriggeredMuon tMu) = 0;
        virtual void scale(float scaleFactor) = 0;
        virtual void config(std::string baseName = "") = 0;

        virtual void plot() = 0;

        void printHisto(TCanvas * canvas, std::string tag) {  
	  canvas->SaveAs((std::string("plots/")+baseName_+"/"+baseName_+tag+".png").c_str());
	  return;
	};

    protected :
        TFile* outFile_;
        hTH1Map histos_;
        bool isCrossSection_;
        std::string baseName_;

};


// ********************************
// (2)
// DT Upgrade Plotter class 
// (pt dependent plots and rates)
// ********************************

class DTRatePtPlotter : public Plotter {

    public:
        DTRatePtPlotter(TFile* file, bool isCrossSection) ;
        ~DTRatePtPlotter(); 

        void fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
		  int igmt, int idt, int irpcb, int icsc, 
		  int irpcf, TriggeredMuon tMu) ;
        void config(std::string baseName = "");
        void scale(float scaleFactor) ;
        void plot();

};


DTRatePtPlotter::DTRatePtPlotter(TFile* file, bool isCrossSection)
{ 

    outFile_ = file;
    isCrossSection_ = isCrossSection;

};

DTRatePtPlotter::~DTRatePtPlotter()
{

    histos_.clear();

};

void DTRatePtPlotter::config(std::string baseName) 
{ 

    baseName_ = baseName;

    if (baseName_!="")
    {
        outFile_->mkdir(baseName_.c_str());
        outFile_->cd(baseName_.c_str());
	system(string("mkdir -p plots/" + baseName_).c_str());
    }

    string name = (baseName_ + "_hRatevsPt");  
    histos_["hRatevsPt"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5);

    name = (baseName_ + "_hRatevsPtGlobal");  
    histos_["hRatevsPtGlobal"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5);

    name = (baseName_ + "_hRatevsPtTight");  
    histos_["hRatevsPtTight"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5);

    name = (baseName_ + "_hRatevsPtTightPt");  
    histos_["hRatevsPtTightPt"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5);

    name = (baseName_ + "_hRatePtvsEta");  
    histos_["hRatePtvsEta"] = new TH2F(name.c_str(),name.c_str(),64,-1.2,1.2,81,-0.5,80.5);

    name = (baseName_ + std::string("_hRateQualityvsPt"));  
    histos_["hRateQualityvsPt"] = new TH2F(name.c_str(),name.c_str(),81,-0.5,80.5,7,0.5,7.5);

    name = (baseName_ + std::string("_hRateChambvsPt"));  
    histos_["hRateChambvsPt"] = new TH2F(name.c_str(),name.c_str(),81,-0.5,80.5,11,0.5,11.5);

    name = (baseName_ + std::string("_hGmtPtvsMuPt"));  
    histos_["hGmtPtvsMuPt"] = new TH2F(name.c_str(),name.c_str(),70,0.5,140.5,70,0.5,140.5);

}

void DTRatePtPlotter::fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
			   int igmt, int idt, int irpcb, int icsc, 
			   int irpcf, TriggeredMuon tMu) 
{

    if (igmt != -1 || irpcb != -1 || icsc != -1 || irpcf != -1)
        std::cout << "Hey DtRatePtPlotter does not use igmt.\nCheck your code!" << std::endl;

    float pt   = tMu.gmt()->Ptdt[tMu.igmt()];
    float eta  = tMu.gmt()->Etadt[tMu.igmt()];
    float qual = tMu.gmt()->Qualdt[tMu.igmt()];

    bool  hasGlbMuon    = tMu.hasMuon();
    bool  hasTightMuon  = tMu.hasTightMuon();
    float muPt          = hasTightMuon ? tMu.mu_->pt.at(tMu.imu_) : -1 ;

    int dttfChambId = 0;

    int dtOuterCh = 0;
    for (int iSt=1; iSt<=4; ++iSt)
      if ( tMu.gmtDttf_->prims_.ScList.at(iSt-1) >= 0)
	dtOuterCh = iSt; 
      
    dttfChambId = chambToIndex(qual,dtOuterCh);

    if (hasGlbMuon) 
      {
	histos_["hGmtPtvsMuPt"]->Fill(pt,muPt);
      }

    for (int ipt=0; ipt<=pt;++ipt) 
    {
        histos_["hRatevsPt"]->Fill(ipt);
        if (hasGlbMuon) histos_["hRatevsPtGlobal"]->Fill(ipt);
        if (hasTightMuon) histos_["hRatevsPtTight"]->Fill(ipt);
        histos_["hRatePtvsEta"]->Fill(eta,ipt);
        histos_["hRateQualityvsPt"]->Fill(ipt,qual);
        histos_["hRateChambvsPt"]->Fill(ipt,dttfChambId);
    }

    
    if (hasTightMuon)
      {
	for (int ipt=0; ipt<=muPt;++ipt) 
	  {
	    histos_["hRatevsPtTightPt"]->Fill(ipt);
	  }
      }

    return; 

}

void DTRatePtPlotter::plot() 
{ 

    setTDRStyle();
    gStyle->SetOptTitle(0);

    TCanvas *cRatevsPt = new TCanvas((baseName_+"cRatevsPt").c_str(),
            (baseName_+"cRatevsPt").c_str(),500,500);
    cRatevsPt->cd();
    cRatevsPt->SetLogy();
    cRatevsPt->SetGrid();  

    histos_["hRatevsPt"]->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
    histos_["hRatevsPt"]->GetXaxis()->SetTitleOffset(0.9);
    histos_["hRatevsPt"]->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    histos_["hRatevsPt"]->GetYaxis()->SetRangeUser(100.,50000.);
    histos_["hRatevsPt"]->SetLineColor(kRed);
    histos_["hRatevsPt"]->SetLineWidth(2);
    histos_["hRatevsPt"]->Draw();

    histos_["hRatevsPtGlobal"]->SetLineColor(kGreen);
    histos_["hRatevsPtGlobal"]->SetLineWidth(2);
    histos_["hRatevsPtGlobal"]->Draw("same");

    histos_["hRatevsPtTight"]->SetLineColor(kBlue);
    histos_["hRatevsPtTight"]->SetLineWidth(2);
    histos_["hRatevsPtTight"]->Draw("same");

    histos_["hRatevsPtTightPt"]->SetLineColor(kCyan);
    histos_["hRatevsPtTightPt"]->SetLineWidth(2);
    histos_["hRatevsPtTightPt"]->Draw("same");

    printHisto(cRatevsPt,"RatevsPt");


    TCanvas *cRatevsPtRatio = new TCanvas((baseName_+"cRatevsPtRatio").c_str(),
					  (baseName_+"cRatevsPtRatio").c_str(),500,500);
    cRatevsPtRatio->cd();
    cRatevsPtRatio->SetGrid();  


    TH1F * hRatevsPtGlobalRatio = static_cast<TH1F*>(histos_["hRatevsPtGlobal"]->Clone("hRatevsPtGlobalRatio"));
    hRatevsPtGlobalRatio->Divide(histos_["hRatevsPt"]);

    TH1F * hRatevsPtTightRatio = static_cast<TH1F*>(histos_["hRatevsPtTight"]->Clone("hRatevsPtTightRatio"));
    hRatevsPtTightRatio->Divide(histos_["hRatevsPt"]);

    TH1F * hRatevsPtTightPtRatio = static_cast<TH1F*>(histos_["hRatevsPtTightPt"]->Clone("hRatevsPtTightPtRatio"));
    hRatevsPtTightPtRatio->Divide(histos_["hRatevsPt"]);

    hRatevsPtGlobalRatio->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
    hRatevsPtGlobalRatio->GetXaxis()->SetTitleOffset(0.9);
    hRatevsPtGlobalRatio->GetYaxis()->SetTitle("Rate Ratio");
    hRatevsPtGlobalRatio->GetYaxis()->SetRangeUser(0.,1.);

    hRatevsPtGlobalRatio->SetLineColor(kGreen);
    hRatevsPtGlobalRatio->SetLineWidth(2);
    hRatevsPtGlobalRatio->Draw();

    hRatevsPtTightRatio->SetLineColor(kBlue);
    hRatevsPtTightRatio->SetLineWidth(2);
    hRatevsPtTightRatio->Draw("same");

    hRatevsPtTightPtRatio->SetLineColor(kCyan);
    hRatevsPtTightPtRatio->SetLineWidth(2);
    hRatevsPtTightPtRatio->Draw("same");

    printHisto(cRatevsPtRatio,"RatevsPtRatio");


    TCanvas *cGmtPtvsMuPt = new TCanvas((baseName_+"cGmtPtvsMuPt").c_str(),
            (baseName_+"cGmtPtvsMuPt").c_str(),500,500);
    cGmtPtvsMuPt->cd();
    cGmtPtvsMuPt->SetLogz();
    cGmtPtvsMuPt->SetGrid();  

    histos_["hGmtPtvsMuPt"]->GetXaxis()->SetTitle("GMT p_{T}");
    histos_["hGmtPtvsMuPt"]->GetYaxis()->SetTitle("Glb Mu p_{T}");
    histos_["hGmtPtvsMuPt"]->Draw("colz");

    printHisto(cRatevsPt,"RatevsPt");

    int ptBins[5] = {12, 16, 20, 25, 30}; // CB adapt pt thresholds if needed

    for (int bin=0; bin<5; ++bin) 
    {
      std::cout << "pt : " << ptBins[bin] << " rate : " 
                <<   histos_["hRatevsPt"]->GetBinContent(histos_["hRatevsPt"]->FindBin(ptBins[bin]))
		<< std::endl;
    }

    TCanvas *cRatePtvsEta = new TCanvas((baseName_+"cRatePtvsEta").c_str(),
            (baseName_+"cRatePtvsEta").c_str(),500,500);
    cRatePtvsEta->cd();
    cRatePtvsEta->SetGrid();

    histos_["hRatePtvsEta"]->GetXaxis()->SetTitle("DT Trigger #eta");
    histos_["hRatePtvsEta"]->GetXaxis()->SetTitleOffset(0.9);
    histos_["hRatePtvsEta"]->GetYaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");

    for(int ybin=1; ybin<=histos_["hRatePtvsEta"]->GetNbinsY(); ++ybin)
    {
        float integral = static_cast<TH2*>(histos_["hRatePtvsEta"])->Integral(0,
                histos_["hRatePtvsEta"]->GetNbinsX()+1,ybin,ybin);
        for(int xbin=1; xbin<=histos_["hRatePtvsEta"]->GetNbinsX(); ++xbin)
        {
            float value = histos_["hRatePtvsEta"]->GetBinContent(xbin,ybin) / integral;
            histos_["hRatePtvsEta"]->SetBinContent(xbin,ybin,value);
        }
    }

    histos_["hRatePtvsEta"]->SetMaximum(0.25);
    histos_["hRatePtvsEta"]->Draw("colz");

    printHisto(cRatePtvsEta,"RatePtvsEta");


    TCanvas *cRateQualityvsPt = new TCanvas((baseName_+"cRateQualityvsPt").c_str(),
            (baseName_+"cRateQualityvsPt").c_str(),500,500);
    cRateQualityvsPt->cd();
    cRateQualityvsPt->SetLogy();
    cRateQualityvsPt->SetGrid();


    THStack * hRatevsPtQualityStack = new THStack("hRatevsPtQualityStack", "hRatevsPtQualityStack");
    THStack * hRatevsPtQualityContriStack = new THStack("hRatevsPtQualityContriStack", "hRatevsPtQualityContriStack");

    for(int ybin=1; ybin<=histos_["hRateQualityvsPt"]->GetNbinsY(); ++ybin)
    {
        stringstream px ;
        px<<ybin<<"Quality";
        TH1D* projection = static_cast<TH2*>(histos_["hRateQualityvsPt"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMap[ybin-1]);
        projection->SetLineColor(colorMap[ybin-1]);
        projection->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
        projection->GetYaxis()->SetTitle(isCrossSection_ ?
                "DT Trigger x-section [#mub]" :
                "rate [arb. units]");

        hRatevsPtQualityStack->Add(projection);

        px<<ybin<<"Contri";
        TH1D* contriProj = new TH1D(px.str().c_str(),px.str().c_str(),81,-0.5,80.5); 
        for(int xbin=1; xbin<=81; ++xbin)
        {
//            contriProj->SetBinContent(xbin, 1.*projection->GetBinContent(xbin)/histos_["hRatevsPtNoHorns"]->GetBinContent(xbin));
            contriProj->SetBinContent(xbin, 1.*projection->GetBinContent(xbin)/histos_["hRatevsPt"]->GetBinContent(xbin));
        }
        contriProj->SetFillColor(colorMap[ybin-1]);
        contriProj->SetLineColor(colorMap[ybin-1]);
        contriProj->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
        contriProj->GetYaxis()->SetTitle("Contribution to total");

        hRatevsPtQualityContriStack->Add(contriProj);

    }

    hRatevsPtQualityStack->Draw();
    hRatevsPtQualityStack->GetHistogram()->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
    hRatevsPtQualityStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsPtQualityStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsPtQualityStack->Draw("Asame");

    cRateQualityvsPt->Update();
    cRateQualityvsPt->Modified();

    printHisto(cRateQualityvsPt,"RateQualityvsPt");

    TCanvas *cRateQualityvsPtContri = new TCanvas((baseName_+"cRateQualityvsPtContri").c_str(),
						  (baseName_+"cRateQualityvsPtContri").c_str(),
						  500,500);
    cRateQualityvsPtContri->cd();
    cRateQualityvsPtContri->SetGrid();

    hRatevsPtQualityContriStack->Draw();
    hRatevsPtQualityContriStack->GetHistogram()->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
    hRatevsPtQualityContriStack->GetHistogram()->GetYaxis()->SetTitle("Contribution to total");
    hRatevsPtQualityContriStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsPtQualityContriStack->GetHistogram()->GetXaxis()->SetRangeUser(0,40);
    hRatevsPtQualityContriStack->Draw("Asame");

    cRateQualityvsPtContri->Update();
    cRateQualityvsPtContri->Modified();

    printHisto(cRateQualityvsPtContri,"RateQualityvsPtContri");


    TCanvas *cRateChambvsPt = new TCanvas((baseName_+"cRateChambvsPt").c_str(),
            (baseName_+"cRateChambvsPt").c_str(),500,500);
    cRateChambvsPt->cd();
    cRateChambvsPt->SetLogy();
    cRateChambvsPt->SetGrid();


    THStack * hRatevsPtChambStack = new THStack("hRatevsPtChambStack", "hRatevsPtChambStack");
    THStack * hRatevsPtChambContriStack = new THStack("hRatevsPtChambContriStack", "hRatevsPtChambContriStack");

    for(int ybin=1; ybin<=histos_["hRateChambvsPt"]->GetNbinsY(); ++ybin)
    {
        stringstream px ;
        px<<ybin<<"Chamb";
        TH1D* projection = static_cast<TH2*>(histos_["hRateChambvsPt"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMapFine[ybin-1]);
        projection->SetLineColor(colorMapFine[ybin-1]);
        projection->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
        projection->GetYaxis()->SetTitle(isCrossSection_ ?
                "DT Trigger x-section [#mub]" :
                "rate [arb. units]");

        hRatevsPtChambStack->Add(projection);

        px<<ybin<<"Contri";
        TH1D* contriProj = new TH1D(px.str().c_str(),px.str().c_str(),81,-0.5,80.5); 
        for(int xbin=1; xbin<=81; ++xbin)
        {
//            contriProj->SetBinContent(xbin, 1.*projection->GetBinContent(xbin)/histos_["hRatevsPtNoHorns"]->GetBinContent(xbin));
            contriProj->SetBinContent(xbin, 1.*projection->GetBinContent(xbin)/histos_["hRatevsPt"]->GetBinContent(xbin));
        }
        contriProj->SetFillColor(colorMapFine[ybin-1]);
        contriProj->SetLineColor(colorMapFine[ybin-1]);
        contriProj->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
        contriProj->GetYaxis()->SetTitle("Contribution to total");

        hRatevsPtChambContriStack->Add(contriProj);

    }

    hRatevsPtChambStack->Draw();
    hRatevsPtChambStack->GetHistogram()->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
    hRatevsPtChambStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsPtChambStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsPtChambStack->Draw("Asame");

    cRateChambvsPt->Update();
    cRateChambvsPt->Modified();

    printHisto(cRateChambvsPt,"RateChambvsPt");

    TCanvas *cRateChambvsPtContri = new TCanvas((baseName_+"cRateChambvsPtContri").c_str(),
						  (baseName_+"cRateChambvsPtContri").c_str(),
						  500,500);
    cRateChambvsPtContri->cd();
    cRateChambvsPtContri->SetGrid();

    hRatevsPtChambContriStack->Draw();
    hRatevsPtChambContriStack->GetHistogram()->GetXaxis()->SetTitle("DT Trigger p_{t} [GeV/c]");
    hRatevsPtChambContriStack->GetHistogram()->GetYaxis()->SetTitle("Contribution to total");
    hRatevsPtChambContriStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsPtChambContriStack->GetHistogram()->GetXaxis()->SetRangeUser(0,40);
    hRatevsPtChambContriStack->Draw("Asame");

    cRateChambvsPtContri->Update();
    cRateChambvsPtContri->Modified();

    printHisto(cRateChambvsPtContri,"RateChambvsPtContri");

    return;

}

void DTRatePtPlotter::scale(float scaleFactor) 
{ 

    hTH1MapIt histoIt  = histos_.begin();
    hTH1MapIt histoEnd = histos_.end();

    for(;histoIt!=histoEnd;++histoIt)
        histoIt->second->Scale(scaleFactor);

    return;

}


// ****************************************
// (3)
// DT Upgrade Plotter class
// (eta and phi pt_cut dependent plots)
// ****************************************

class DTRateEtaPhiPlotter : public Plotter {

    public:
        DTRateEtaPhiPlotter(TFile* file, bool isCrossSection) ;
        ~DTRateEtaPhiPlotter(); 

        void chamberCombLabels(TH1* h);
        int phiToBin(float phi); 

        void fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
		  int igmt, int idt, int irpcb, int icsc, 
		  int irpcf, TriggeredMuon tMu) ;
        void config(std::string baseName = "");
        void scale(float scaleFactor) ;
        void plot();

    private:
  map<float,int> phiBins;

};

DTRateEtaPhiPlotter::DTRateEtaPhiPlotter(TFile* file, bool isCrossSection)
{ 

    outFile_ = file;
    isCrossSection_ = isCrossSection;

};

DTRateEtaPhiPlotter::~DTRateEtaPhiPlotter()
{

    histos_.clear();

};


void DTRateEtaPhiPlotter::chamberCombLabels(TH1* h)
{

  h->GetXaxis()->SetBinLabel(1,"MB3/4");
  h->GetXaxis()->SetBinLabel(2,"MB2/4");
  h->GetXaxis()->SetBinLabel(3,"MB2/3");
  h->GetXaxis()->SetBinLabel(4,"MB1/4");
  h->GetXaxis()->SetBinLabel(5,"MB1/3");
  h->GetXaxis()->SetBinLabel(6,"MB1/2");
  h->GetXaxis()->SetBinLabel(7,"MB2/3/4");
  h->GetXaxis()->SetBinLabel(8,"MB1/3/4");
  h->GetXaxis()->SetBinLabel(9,"MB1/2/4");
  h->GetXaxis()->SetBinLabel(10,"MB1/2/3");
  h->GetXaxis()->SetBinLabel(11,"MB1/2/3/4");

}

int DTRateEtaPhiPlotter::phiToBin(float phi) 
{

  map<float,int>:: const_iterator phiBinIt  = phiBins.begin();
  map<float,int>:: const_iterator phiBinEnd = phiBins.end();
 
  int iBin = 0;
  for(;phiBinIt!=phiBinEnd;++phiBinIt)
    {
      if (fabs(phiBinIt->first - phi)<0.01)
	{
	  iBin = phiBinIt->second;
	  break;
	}
    }

  return (iBin+3)/6;

}


void DTRateEtaPhiPlotter::config(std::string baseName) 
{ 

	phiBins[-0.0001]=0;
	phiBins[0.0435332]=1;
	phiBins[0.0871665]=2;
	phiBins[0.1308]=3;
	phiBins[0.174433]=4;
	phiBins[0.218066]=5;
	phiBins[0.261699]=6;
	phiBins[0.305333]=7;
	phiBins[0.348966]=8;
	phiBins[0.392599]=9;
	phiBins[0.436232]=10;
	phiBins[0.479866]=11;
	phiBins[0.523499]=12;
	phiBins[0.567132]=13;
	phiBins[0.610765]=14;
	phiBins[0.654398]=15;
	phiBins[0.698032]=16;
	phiBins[0.741665]=17;
	phiBins[0.785298]=18;
	phiBins[0.828931]=19;
	phiBins[0.872565]=20;
	phiBins[0.916198]=21;
	phiBins[0.959831]=22;
	phiBins[1.00346]=23;
	phiBins[1.0471]=24;
	phiBins[1.09073]=25;
	phiBins[1.13436]=26;
	phiBins[1.178]=27;
	phiBins[1.22163]=28;
	phiBins[1.26526]=29;
	phiBins[1.3089]=30;
	phiBins[1.35253]=31;
	phiBins[1.39616]=32;
	phiBins[1.4398]=33;
	phiBins[1.48343]=34;
	phiBins[1.52706]=35;
	phiBins[1.5707]=36;
	phiBins[1.61433]=37;
	phiBins[1.65796]=38;
	phiBins[1.7016]=39;
	phiBins[1.74523]=40;
	phiBins[1.78886]=41;
	phiBins[1.8325]=42;
	phiBins[1.87613]=43;
	phiBins[1.91976]=44;
	phiBins[1.9634]=45;
	phiBins[2.00703]=46;
	phiBins[2.05066]=47;
	phiBins[2.0943]=48;
	phiBins[2.13793]=49;
	phiBins[2.18156]=50;
	phiBins[2.22519]=51;
	phiBins[2.26883]=52;
	phiBins[2.31246]=53;
	phiBins[2.35609]=54;
	phiBins[2.39973]=55;
	phiBins[2.44336]=56;
	phiBins[2.48699]=57;
	phiBins[2.53063]=58;
	phiBins[2.57426]=59;
	phiBins[2.61789]=60;
	phiBins[2.66153]=61;
	phiBins[2.70516]=62;
	phiBins[2.74879]=63;
	phiBins[2.79243]=64;
	phiBins[2.83606]=65;
	phiBins[2.87969]=66;
	phiBins[2.92333]=67;
	phiBins[2.96696]=68;
	phiBins[3.01059]=69;
	phiBins[3.05423]=70;
	phiBins[3.09786]=71;
	phiBins[3.14149]=72;
	phiBins[3.18513]=73;
	phiBins[3.22876]=74;
	phiBins[3.27239]=75;
	phiBins[3.31603]=76;
	phiBins[3.35966]=77;
	phiBins[3.40329]=78;
	phiBins[3.44693]=79;
	phiBins[3.49056]=80;
	phiBins[3.53419]=81;
	phiBins[3.57783]=82;
	phiBins[3.62146]=83;
	phiBins[3.66509]=84;
	phiBins[3.70872]=85;
	phiBins[3.75236]=86;
	phiBins[3.79599]=87;
	phiBins[3.83962]=88;
	phiBins[3.88326]=89;
	phiBins[3.92689]=90;
	phiBins[3.97052]=91;
	phiBins[4.01416]=92;
	phiBins[4.05779]=93;
	phiBins[4.10142]=94;
	phiBins[4.14506]=95;
	phiBins[4.18869]=96;
	phiBins[4.23232]=97;
	phiBins[4.27596]=98;
	phiBins[4.31959]=99;
	phiBins[4.36322]=100;
	phiBins[4.40686]=101;
	phiBins[4.45049]=102;
	phiBins[4.49412]=103;
	phiBins[4.53776]=104;
	phiBins[4.58139]=105;
	phiBins[4.62502]=106;
	phiBins[4.66866]=107;
	phiBins[4.71229]=108;
	phiBins[4.75592]=109;
	phiBins[4.79956]=110;
	phiBins[4.84319]=111;
	phiBins[4.88682]=112;
	phiBins[4.93046]=113;
	phiBins[4.97409]=114;
	phiBins[5.01772]=115;
	phiBins[5.06135]=116;
	phiBins[5.10499]=117;
	phiBins[5.14862]=118;
	phiBins[5.19225]=119;
	phiBins[5.23589]=120;
	phiBins[5.27952]=121;
	phiBins[5.32315]=122;
	phiBins[5.36679]=123;
	phiBins[5.41042]=124;
	phiBins[5.45405]=125;
	phiBins[5.49769]=126;
	phiBins[5.54132]=127;
	phiBins[5.58495]=128;
	phiBins[5.62859]=129;
	phiBins[5.67222]=130;
	phiBins[5.71585]=131;
	phiBins[5.75949]=132;
	phiBins[5.80312]=133;
	phiBins[5.84675]=134;
	phiBins[5.89039]=135;
	phiBins[5.93402]=136;
	phiBins[5.97765]=137;
	phiBins[6.02129]=138;
	phiBins[6.06492]=139;
	phiBins[6.10855]=140;
	phiBins[6.15219]=-3;
	phiBins[6.19582]=-2;
	phiBins[6.23945]=-1;

    baseName_ = baseName;

    if (baseName!="")
    {
        outFile_->mkdir(baseName.c_str());
        outFile_->cd(baseName.c_str());
	system(string("mkdir -p plots/" + baseName_).c_str());
    }

    string name = (baseName_ + std::string("_hRatevsPhi"));  
    histos_["hRatevsPhi"] = new TH1F(name.c_str(),name.c_str(),24,-0.5,23.5);

    name = (baseName + std::string("_hRateChambvsPhi"));  
    histos_["hRateChambvsPhi"] = new TH2F(name.c_str(),name.c_str(),24,-0.5,23.5,11,0.5,11.5);

    name = (baseName_ + std::string("_hRatevsEta"));  
    histos_["hRatevsEta"] = new TH1F(name.c_str(),name.c_str(),32,-1.2,1.2);

    name = (baseName_ + std::string("_hRatevsQuality"));  
    histos_["hRatevsQuality"] = new TH1F(name.c_str(),name.c_str(),7,0.5,7.5);

    name = (baseName + std::string("_hRateQualityvsEta"));  
    histos_["hRateQualityvsEta"] = new TH2F(name.c_str(),name.c_str(),32,-1.2,1.2,7,0.5,7.5);

    name = (baseName_ + std::string("_hRatevsChamb"));  
    histos_["hRatevsChamb"] = new TH1F(name.c_str(),name.c_str(),11,0.5,11.5);

    name = (baseName_ + std::string("_hRatevsChambGlobal"));  
    histos_["hRatevsChambGlobal"] = new TH1F(name.c_str(),name.c_str(),11,0.5,11.5);

    name = (baseName_ + std::string("_hRatevsChambTight"));  
    histos_["hRatevsChambTight"] = new TH1F(name.c_str(),name.c_str(),11,0.5,11.5);

    name = (baseName_ + std::string("_hRatevsChambTightPt"));  
    histos_["hRatevsChambTightPt"] = new TH1F(name.c_str(),name.c_str(),11,0.5,11.5);

    name = (baseName + std::string("_hRatePrimQualvsChamb"));  
    histos_["hRatePrimQualvsChamb"] = new TH2F(name.c_str(),name.c_str(),11,0.5,11.5,5,-0.5,4.5);

    name = (baseName + std::string("_hRatePrimQualvsChambGlobal"));  
    histos_["hRatePrimQualvsChambGlobal"] = new TH2F(name.c_str(),name.c_str(),11,0.5,11.5,5,-0.5,4.5);

    name = (baseName + std::string("_hRatePrimQualvsChambTight"));  
    histos_["hRatePrimQualvsChambTight"] = new TH2F(name.c_str(),name.c_str(),11,0.5,11.5,5,-0.5,4.5);

    name = (baseName + std::string("_hRatePrimQualvsChambTightPt"));  
    histos_["hRatePrimQualvsChambTightPt"] = new TH2F(name.c_str(),name.c_str(),11,0.5,11.5,5,-0.5,4.5);

    name = (baseName + std::string("_hRateChambvsEta"));  
    histos_["hRateChambvsEta"] = new TH2F(name.c_str(),name.c_str(),32,-1.2,1.2,11,0.5,11.5);

    name = (baseName_ + std::string("_hRatevsMuonKind"));  
    histos_["hRatevsMuonKind"] = new TH1F(name.c_str(),name.c_str(),4,0.5,4.5);


}

void DTRateEtaPhiPlotter::fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
			       int igmt, int idt, int irpcb, int icsc, 
			       int irpcf, TriggeredMuon tMu) 
{

  if (igmt != -1 || irpcb != -1 || icsc != -1 || irpcf != -1)
        std::cout << "Hey DtRateEtaPhiPlotter does not use igmt, irpcb, icsc or irpcf.\nCheck your code!" << std::endl;

    float eta  = tMu.gmt()->Etadt[tMu.igmt()];
    float qual = tMu.gmt()->Qualdt[tMu.igmt()];
    float phi  = tMu.gmt()->Phi[tMu.igmt()] - 0.0001;

    int phiBin = phiToBin(phi);

    int dttfChambId = 0;

    int dtOuterCh = 0;
    int nCh = 0;
    for (int iSt=1; iSt<=4; ++iSt)
      {
	nCh++;
	if ( tMu.gmtDttf_->prims_.ScList.at(iSt-1) >= 0)
	dtOuterCh = iSt; 
      }
      
    dttfChambId = chambToIndex(qual,dtOuterCh);

    bool hasGlbMuon   = tMu.hasMuon();
    bool hasTightMuon = tMu.hasTightMuon();
    bool isGoodPt     = tMu.hasTightMuon() && tMu.mu_->pt.at(tMu.imu_) > 16.;

    int nCorrPrims = 0;
    
    vector<int>::const_iterator iPrimsIt  = tMu.gmtDttf_->prims_.PrimitivesId.begin();
    vector<int>::const_iterator iPrimsEnd = tMu.gmtDttf_->prims_.PrimitivesId.end();


    for(;iPrimsIt!=iPrimsEnd; ++iPrimsIt) {
      int primQual = tMu.gmtDttf_->dttf_->phCode[(*iPrimsIt)];
      (primQual >= 4) && nCorrPrims++;
    }

//     if (tMu.gmtDttf_->prims_.PrimitivesId.size() == 1)
//       {
// 	cout << "qual : " << qual
// 	     << " nCh : "  << nCh
// 	     << " wh : " << tMu.gmtDttf_->dttf_->trWh.at(tMu.gmtDttf_->idttf_)
// 	     << " size : " << tMu.gmtDttf_->prims_.PrimitivesId.size()
// 	     << " outer : "<< dtOuterCh
// 	     << " chId : " << dttfChambId
// 	     << " corr : "<< nCorrPrims 
// 	     << endl;
//       }

//     if (tMu.gmtDttf_->prims_.PrimitivesId.size() == 3 && qual == 3)
//       {
// 	cout << "qual : " << qual
// 	     << " nCh : "  << nCh
// 	     << " wh : " << tMu.gmtDttf_->dttf_->trWh.at(tMu.gmtDttf_->idttf_)
// 	     << " size : " << tMu.gmtDttf_->prims_.PrimitivesId.size()
// 	     << " outer : "<< dtOuterCh
// 	     << " chId : " << dttfChambId
// 	     << " corr : "<< nCorrPrims 
// 	     << endl;
//       }


    histos_["hRatevsEta"]->Fill(eta);
    histos_["hRatevsQuality"]->Fill(qual);
    histos_["hRatevsChamb"]->Fill(dttfChambId);
    if (hasGlbMuon) histos_["hRatevsChambGlobal"]->Fill(dttfChambId);
    if (hasTightMuon) histos_["hRatevsChambTight"]->Fill(dttfChambId);
    if (isGoodPt) histos_["hRatevsChambTightPt"]->Fill(dttfChambId);
    histos_["hRatePrimQualvsChamb"]->Fill(dttfChambId,nCorrPrims);
    if (hasGlbMuon) histos_["hRatePrimQualvsChambGlobal"]->Fill(dttfChambId,nCorrPrims);
    if (hasTightMuon) histos_["hRatePrimQualvsChambTight"]->Fill(dttfChambId,nCorrPrims);
    if (isGoodPt) histos_["hRatePrimQualvsChambTightPt"]->Fill(dttfChambId,nCorrPrims);
    histos_["hRatevsMuonKind"]->Fill(isGoodPt ? 4 : hasTightMuon ? 3 : hasGlbMuon ? 2 : 1); 
    histos_["hRateQualityvsEta"]->Fill(eta,qual);
    histos_["hRateChambvsEta"]->Fill(eta,dttfChambId);
    histos_["hRatevsPhi"]->Fill(phiBin);
    histos_["hRateChambvsPhi"]->Fill(phiBin,dttfChambId);

    return; 

}

void DTRateEtaPhiPlotter::plot() 
{ 

  // map<float,float>:: const_iterator phiBinIt  = phiBins.begin();
  // map<float,float>:: const_iterator phiBinEnd = phiBins.end();
 
  // for (int iBin=0; phiBinIt != phiBinEnd; ++phiBinIt, ++iBin)
  //   {
  //     cout << iBin << "\tphiBins[" << phiBinIt->first << "]=" << iBin <<";" << endl;
  //   }

    setTDRStyle();
    gStyle->SetOptTitle(0);

    TCanvas *cRatevsPhi = new TCanvas((baseName_+"cRatevsPhi").c_str(),
            (baseName_+"cRatevsPhi").c_str(),500,500);
    cRatevsPhi->cd();

    histos_["hRatevsPhi"]->GetXaxis()->SetTitle("DT Trigger #phi");
    histos_["hRatevsPhi"]->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
//    histos_["hRatevsPhi"]->GetYaxis()->SetTitleOffset(0.9);
    histos_["hRatevsPhi"]->Draw();

    printHisto(cRatevsPhi,"RatevsPhi");

    TCanvas *cRatevsQuality = new TCanvas((baseName_+"cRatevsQuality").c_str(),
            (baseName_+"cRatevsQuality").c_str(),500,500);
    cRatevsQuality->cd();
    cRatevsQuality->SetGrid();

    histos_["hRatevsQuality"]->GetXaxis()->SetTitle("DT Trigger quality");
    histos_["hRatevsQuality"]->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
//    histos_["hRatevsQuality"]->GetYaxis()->SetTitleOffset(0.9);
    histos_["hRatevsQuality"]->SetFillColor(4);
    histos_["hRatevsQuality"]->SetLineColor(4);
    histos_["hRatevsQuality"]->Draw();

    printHisto(cRatevsQuality,"RatevsQuality");

    TCanvas *cRatevsChamb = new TCanvas((baseName_+"cRatevsChamb").c_str(),
            (baseName_+"cRatevsChamb").c_str(),500,500);
    cRatevsChamb->cd();
    cRatevsChamb->SetGrid();

    float entries = histos_["hRatevsChamb"]->Integral();

    histos_["hRatevsChamb"]->GetXaxis()->SetTitle("Chamber Combination");
    histos_["hRatevsChamb"]->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
//    histos_["hRatevsChamb"]->GetYaxis()->SetTitleOffset(0.9);
    histos_["hRatevsChamb"]->SetFillColor(kRed);
    histos_["hRatevsChamb"]->SetLineColor(kRed);
    chamberCombLabels(histos_["hRatevsChamb"]);
    histos_["hRatevsChamb"]->Scale(1./entries);
    histos_["hRatevsChamb"]->Draw();

    histos_["hRatevsChambGlobal"]->SetFillColor(kGreen);
    histos_["hRatevsChambGlobal"]->SetLineColor(kGreen);
    histos_["hRatevsChambGlobal"]->Scale(1./entries);
    histos_["hRatevsChambGlobal"]->Draw("same");

    histos_["hRatevsChambTight"]->SetFillColor(kBlue);
    histos_["hRatevsChambTight"]->SetLineColor(kBlue);
    histos_["hRatevsChambTight"]->Scale(1./entries);
    histos_["hRatevsChambTight"]->Draw("same");

    histos_["hRatevsChambTightPt"]->SetFillColor(kCyan);
    histos_["hRatevsChambTightPt"]->SetLineColor(kCyan);
    histos_["hRatevsChambTightPt"]->Scale(1./entries);
    histos_["hRatevsChambTightPt"]->Draw("same");

    printHisto(cRatevsChamb,"RatevsChamb");

    TCanvas *cRatevsMuonKind = new TCanvas((baseName_+"cRatevsMuonKind").c_str(),
            (baseName_+"cRatevsMuonKind").c_str(),500,500);
    cRatevsMuonKind->cd();
    cRatevsMuonKind->SetGrid();

    histos_["hRatevsMuonKind"]->GetXaxis()->SetTitle("Muon Kind");
    histos_["hRatevsMuonKind"]->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
//    histos_["hRatevsMuonKind"]->GetYaxis()->SetTitleOffset(0.9);
    histos_["hRatevsMuonKind"]->SetFillColor(4);
    histos_["hRatevsMuonKind"]->SetLineColor(4);
    histos_["hRatevsMuonKind"]->Draw();

    printHisto(cRatevsMuonKind,"RatevsMuonKind");


    TCanvas *cRateQualityvsEta = new TCanvas((baseName_+"cRateQualityvsEta").c_str(),
            (baseName_+"cRateQualityvsEta").c_str(),500,500);
    cRateQualityvsEta->cd();
    cRateQualityvsEta->SetGridy();


    THStack * hRatevsEtaQualityStack = new THStack((baseName_+"hRatevsEtaQualityStack").c_str(),
            (baseName_+"hRatevsEtaQualityStack").c_str());
    for(int ybin=1; ybin<=histos_["hRateQualityvsEta"]->GetNbinsY(); ++ybin)
    {
        stringstream px ;
        px<<ybin<<baseName_;
        TH1D* projection = static_cast<TH2*>(histos_["hRateQualityvsEta"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMap[ybin-1]);
        projection->SetLineColor(colorMap[ybin-1]);
        projection->GetXaxis()->SetTitle("DT Trigger #eta");
        projection->GetYaxis()->SetTitle(isCrossSection_ ? 
                "x-section [#mub]" : "rate [arb. units]");
        hRatevsEtaQualityStack->Add(projection);
    }

    hRatevsEtaQualityStack->Draw();
    hRatevsEtaQualityStack->GetHistogram()->GetXaxis()->SetTitle("DT Trigger #eta");
    hRatevsEtaQualityStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsEtaQualityStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsEtaQualityStack->Draw("Asame");

    cRateQualityvsEta->Update();
    cRateQualityvsEta->Modified();

    printHisto(cRateQualityvsEta,"RateQualityvsEta");


    TCanvas *cRateChambvsEta = new TCanvas((baseName_+"cRateChambvsEta").c_str(),
            (baseName_+"cRateChambvsEta").c_str(),500,500);
    cRateChambvsEta->cd();
    cRateChambvsEta->SetGridy();


    THStack * hRatevsEtaChambStack = new THStack((baseName_+"hRatevsEtaChambStack").c_str(),
            (baseName_+"hRatevsEtaChambStack").c_str());
    for(int ybin=1; ybin<=histos_["hRateChambvsEta"]->GetNbinsY(); ++ybin)
    {
        stringstream px ;
        px<<ybin<<baseName_;
        TH1D* projection = static_cast<TH2*>(histos_["hRateChambvsEta"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMapFine[ybin-1]);
        projection->SetLineColor(colorMapFine[ybin-1]);
        projection->GetXaxis()->SetTitle("DT Trigger #eta");
        projection->GetYaxis()->SetTitle(isCrossSection_ ? 
                "x-section [#mub]" : "rate [arb. units]");
        hRatevsEtaChambStack->Add(projection);
    }

    hRatevsEtaChambStack->Draw();
    hRatevsEtaChambStack->GetHistogram()->GetXaxis()->SetTitle("DT Trigger #eta");
    hRatevsEtaChambStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsEtaChambStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsEtaChambStack->Draw("Asame");

    cRateChambvsEta->Update();
    cRateChambvsEta->Modified();

    printHisto(cRateChambvsEta,"RateChambvsEta");


    TCanvas *cRateChambvsPhi = new TCanvas((baseName_+"cRateChambvsPhi").c_str(),
            (baseName_+"cRateChambvsPhi").c_str(),500,500);
    cRateChambvsPhi->cd();
    cRateChambvsPhi->SetGridy();


    THStack * hRatevsPhiChambStack = new THStack((baseName_+"hRatevsPhiChambStack").c_str(),
            (baseName_+"hRatevsPhiChambStack").c_str());
    for(int ybin=1; ybin<=histos_["hRateChambvsPhi"]->GetNbinsY(); ++ybin)
    {
        stringstream px ;
        px<<ybin<<baseName_<<"Phi";
        TH1D* projection = static_cast<TH2*>(histos_["hRateChambvsPhi"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMapFine[ybin-1]);
        projection->SetLineColor(colorMapFine[ybin-1]);
        projection->GetXaxis()->SetTitle("DT Trigger #phi");
        projection->GetYaxis()->SetTitle(isCrossSection_ ? 
                "x-section [#mub]" : "rate [arb. units]");
        hRatevsPhiChambStack->Add(projection);
    }

    hRatevsPhiChambStack->Draw();
    hRatevsPhiChambStack->GetHistogram()->GetXaxis()->SetTitle("DT Trigger #phi");
    hRatevsPhiChambStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsPhiChambStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsPhiChambStack->Draw("Asame");

    cRateChambvsPhi->Update();
    cRateChambvsPhi->Modified();

    printHisto(cRateChambvsPhi,"RateChambvsPhi");

    TCanvas *cRatePrimQualvsChamb = new TCanvas((baseName_+"cRatePrimQualvsChamb").c_str(),
						(baseName_+"cRatePrimQualvsChamb").c_str(),500,500);
    cRatePrimQualvsChamb->cd();
    cRatePrimQualvsChamb->SetGridy();

    
    THStack * hRatevsChambPrimQualStack = new THStack((baseName_+"hRatevsChambPrimQualStack").c_str(),
						      (baseName_+"hRatevsChambPrimQualStack").c_str());

    histos_["hRatePrimQualvsChamb"]->Scale(1./histos_["hRatePrimQualvsChamb"]->Integral());    

    for(int ybin=1; ybin<=histos_["hRatePrimQualvsChamb"]->GetNbinsY(); ++ybin)
    {
      stringstream px ;
      px<<ybin<<baseName_ + "ChambPrimQual";
      TH1D* projection = static_cast<TH2*>(histos_["hRatePrimQualvsChamb"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMap[ybin-1]);
        projection->SetLineColor(colorMap[ybin-1]);
	chamberCombLabels(projection);
        projection->GetXaxis()->SetTitle("Chamber Combination");
        projection->GetYaxis()->SetTitle(isCrossSection_ ? 
                "x-section [#mub]" : "rate [arb. units]");
        hRatevsChambPrimQualStack->Add(projection);
    }

    hRatevsChambPrimQualStack->Draw();
    chamberCombLabels(hRatevsChambPrimQualStack->GetHistogram());
    hRatevsChambPrimQualStack->GetHistogram()->GetXaxis()->SetTitle("Chamber Combination");
    
    hRatevsChambPrimQualStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsChambPrimQualStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsChambPrimQualStack->Draw("Asame");

    cRatePrimQualvsChamb->Update();
    cRatePrimQualvsChamb->Modified();

    printHisto(cRatePrimQualvsChamb,"RatePrimQualvsChamb");


    TCanvas *cRatePrimQualvsChambGlobal = new TCanvas((baseName_+"cRatePrimQualvsChambGlobal").c_str(),
						(baseName_+"cRatePrimQualvsChambGlobal").c_str(),500,500);
    cRatePrimQualvsChambGlobal->cd();
    cRatePrimQualvsChambGlobal->SetGridy();

    
    THStack * hRatevsChambGlobalPrimQualStack = new THStack((baseName_+"hRatevsChambGlobalPrimQualStack").c_str(),
						      (baseName_+"hRatevsChambGlobalPrimQualStack").c_str());

    histos_["hRatePrimQualvsChambGlobal"]->Scale(1./histos_["hRatePrimQualvsChambGlobal"]->Integral());    

    for(int ybin=1; ybin<=histos_["hRatePrimQualvsChambGlobal"]->GetNbinsY(); ++ybin)
    {
      stringstream px ;
      px<<ybin<<baseName_ + "ChambGlobalPrimQual";
      TH1D* projection = static_cast<TH2*>(histos_["hRatePrimQualvsChambGlobal"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMap[ybin-1]);
        projection->SetLineColor(colorMap[ybin-1]);
        projection->GetXaxis()->SetTitle("Chamber Combination");
	chamberCombLabels(projection);
        projection->GetYaxis()->SetTitle(isCrossSection_ ? 
                "x-section [#mub]" : "rate [arb. units]");
        hRatevsChambGlobalPrimQualStack->Add(projection);
    }

    hRatevsChambGlobalPrimQualStack->Draw();
    chamberCombLabels(hRatevsChambGlobalPrimQualStack->GetHistogram());
    hRatevsChambGlobalPrimQualStack->GetHistogram()->GetXaxis()->SetTitle("Chamber Combination");
    hRatevsChambGlobalPrimQualStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsChambGlobalPrimQualStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsChambGlobalPrimQualStack->Draw("Asame");

    cRatePrimQualvsChambGlobal->Update();
    cRatePrimQualvsChambGlobal->Modified();

    printHisto(cRatePrimQualvsChambGlobal,"RatePrimQualvsChambGlobal");


    TCanvas *cRatePrimQualvsChambTight = new TCanvas((baseName_+"cRatePrimQualvsChambTight").c_str(),
						(baseName_+"cRatePrimQualvsChambTight").c_str(),500,500);
    cRatePrimQualvsChambTight->cd();
    cRatePrimQualvsChambTight->SetGridy();

    
    THStack * hRatevsChambTightPrimQualStack = new THStack((baseName_+"hRatevsChambTightPrimQualStack").c_str(),
						      (baseName_+"hRatevsChambTightPrimQualStack").c_str());

    histos_["hRatePrimQualvsChambTight"]->Scale(1./histos_["hRatePrimQualvsChambTight"]->Integral());    

    for(int ybin=1; ybin<=histos_["hRatePrimQualvsChambTight"]->GetNbinsY(); ++ybin)
    {
      stringstream px ;
      px<<ybin<<baseName_ + "ChambTightPrimQual";
      TH1D* projection = static_cast<TH2*>(histos_["hRatePrimQualvsChambTight"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMap[ybin-1]);
        projection->SetLineColor(colorMap[ybin-1]);
	chamberCombLabels(projection);
	projection->GetXaxis()->SetTitle("Chamber Combination");
        projection->GetYaxis()->SetTitle(isCrossSection_ ? 
                "x-section [#mub]" : "rate [arb. units]");
        hRatevsChambTightPrimQualStack->Add(projection);
    }

    hRatevsChambTightPrimQualStack->Draw();
    chamberCombLabels(hRatevsChambTightPrimQualStack->GetHistogram());
    hRatevsChambTightPrimQualStack->GetHistogram()->GetXaxis()->SetTitle("Chamber Combination");
    hRatevsChambTightPrimQualStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsChambTightPrimQualStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsChambTightPrimQualStack->Draw("Asame");

    cRatePrimQualvsChambTight->Update();
    cRatePrimQualvsChambTight->Modified();

    printHisto(cRatePrimQualvsChambTight,"RatePrimQualvsChambTight");


    TCanvas *cRatePrimQualvsChambTightPt = new TCanvas((baseName_+"cRatePrimQualvsChambTightPt").c_str(),
						(baseName_+"cRatePrimQualvsChambTightPt").c_str(),500,500);
    cRatePrimQualvsChambTightPt->cd();
    cRatePrimQualvsChambTightPt->SetGridy();

    
    THStack * hRatevsChambTightPtPrimQualStack = new THStack((baseName_+"hRatevsChambTightPtPrimQualStack").c_str(),
						      (baseName_+"hRatevsChambTightPtPrimQualStack").c_str());

    histos_["hRatePrimQualvsChambTightPt"]->Scale(1./histos_["hRatePrimQualvsChambTightPt"]->Integral());    


    for(int ybin=1; ybin<=histos_["hRatePrimQualvsChambTightPt"]->GetNbinsY(); ++ybin)
    {
      stringstream px ;
      px<<ybin<<baseName_ + "ChambTightPtPrimQual";
      TH1D* projection = static_cast<TH2*>(histos_["hRatePrimQualvsChambTightPt"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

        projection->SetFillColor(colorMap[ybin-1]);
        projection->SetLineColor(colorMap[ybin-1]);
	chamberCombLabels(projection);
	projection->GetXaxis()->SetTitle("Chambter Combination");
        projection->GetYaxis()->SetTitle(isCrossSection_ ? 
                "x-section [#mub]" : "rate [arb. units]");
        hRatevsChambTightPtPrimQualStack->Add(projection);
    }

    hRatevsChambTightPtPrimQualStack->Draw();
    chamberCombLabels(hRatevsChambTightPtPrimQualStack->GetHistogram());
    hRatevsChambTightPtPrimQualStack->GetHistogram()->GetXaxis()->SetTitle("Chamber Combination");
    hRatevsChambTightPtPrimQualStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    hRatevsChambTightPtPrimQualStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
    hRatevsChambTightPtPrimQualStack->Draw("Asame");

    cRatePrimQualvsChambTightPt->Update();
    cRatePrimQualvsChambTightPt->Modified();

    printHisto(cRatePrimQualvsChambTightPt,"RatePrimQualvsChambTightPt");
    
    return;

}

void DTRateEtaPhiPlotter::scale(float scaleFactor) 
{ 

    hTH1MapIt histoIt  = histos_.begin();
    hTH1MapIt histoEnd = histos_.end();

    for(;histoIt!=histoEnd;++histoIt)
        histoIt->second->Scale(scaleFactor);

    return;

}


//*************************************
// (4) 
// Plotter for GMT rate plots
// (plot all that is not only related
// to DTTF)
//*************************************

class GMTRatePlotter : public Plotter {

    public:
        GMTRatePlotter(TFile *file, bool isCrossSection) ;
        ~GMTRatePlotter(); 

        void fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
		  int igmt, int idt, int irpcb, int icsc,
		  int irpcf, TriggeredMuon tMu) ; // CB add RPCb RPCf CSC
        void config(std::string baseName = "");
        void scale(float scaleFactor) ;
        void plot();

};

GMTRatePlotter::GMTRatePlotter(TFile *file, bool isCrossSection)
{ 

    outFile_ = file;
    isCrossSection_ = isCrossSection;

};

GMTRatePlotter::~GMTRatePlotter()
{

    histos_.clear();

};

void GMTRatePlotter::config(std::string baseName) 
{ 

    baseName_ = baseName;

    if (baseName!="")
    {
        outFile_->mkdir(baseName.c_str());
        outFile_->cd(baseName.c_str());
	system(string("mkdir -p plots/" + baseName_).c_str());
    }

    string name = (baseName + "_hGMTRates"); //CB rates for GMT, GMT forward, GMT barrel, DT  
    histos_["hGMTRates"] = new TH1F(name.c_str(),name.c_str(),7,0.5,7.5); // add RPCb, RPCf, CSC

    name = (baseName + "_hGMTRatevsPt" ); 
    histos_["hGMTRatevsPt"]  = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5); // add GMT forward, GMT barrel, DT, RPCb, RPCf, CSC
    name = (baseName + "_hGMTfRatevsPt"); 
    histos_["hGMTfRatevsPt"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5); 
    name = (baseName + "_hGMTbRatevsPt"); 
    histos_["hGMTbRatevsPt"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5); 
    name = (baseName + "_hDTRatevsPt"  ); 
    histos_["hDTRatevsPt"]   = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5); 
    name = (baseName + "_hRPCfRatevsPt"); 
    histos_["hRPCfRatevsPt"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5); 
    name = (baseName + "_hRPCbRatevsPt"); 
    histos_["hRPCbRatevsPt"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5); 
    name = (baseName + "_hCSCRatevsPt" ); 
    histos_["hCSCRatevsPt"]  = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5); 

}

void GMTRatePlotter::fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
			  int igmt, int idt, int irpcb, int icsc,
			  int irpcf, TriggeredMuon tMu) 
{
  
  if (idt != -1) {
    
    histos_["hGMTRates"]->Fill(4);  // 4 is DT
    float dtPt = gmt->Ptdt[idt];
    
    for (int ipt=0; ipt<=dtPt;++ipt) 
      {
	histos_["hDTRatevsPt"]  ->Fill(ipt); 
      }

  }
  
  if (irpcb != -1) {
    
    histos_["hGMTRates"]->Fill(5);  // 5 is RPCB
    float rpcbPt = gmt->Ptrpcb[irpcb];
    
    for (int ipt=0; ipt<=rpcbPt;++ipt) 
      {
	histos_["hRPCbRatevsPt"]  ->Fill(ipt); 
      }

  }
  
  if (icsc != -1) {
    
    histos_["hGMTRates"]->Fill(7);  // 7 is CSC
    float cscPt = gmt->Ptcsc[icsc];
    
    for (int ipt=0; ipt<=cscPt;++ipt) 
      {
	histos_["hCSCRatevsPt"]  ->Fill(ipt); 
      }

  }
  
  if (irpcf != -1) {
    
    histos_["hGMTRates"]->Fill(6);  // 6 is RPCF
    float rpcfPt = gmt->Ptrpcf[irpcf];
    
    for (int ipt=0; ipt<=rpcfPt;++ipt) 
      {
	histos_["hRPCfRatevsPt"]  ->Fill(ipt); 
      }

  }
  
  // Then filling GMT bins for events that fired GMT
  if (igmt != -1) 
    {

      bool isDT   = gmt->IdxDTBX[igmt] > -1 ? true : false;
      bool isRPCb = gmt->IdxRPCb[igmt] > -1 ? true : false;
      bool isBarrel = isDT || isRPCb;

      bool isCSC  = gmt->IdxCSC[igmt]  > -1 ? true : false;
      bool isRPCf = gmt->IdxRPCf[igmt] > -1 ? true : false;
      bool isEndcap = isCSC || isRPCf;

      float gmtPt = gmt->Pt[igmt];
      
      histos_["hGMTRates"]->Fill(1);                 // 1 is GMT
      if (isBarrel)  histos_["hGMTRates"]->Fill(2);  // 2 is GMT barrel
      if (isEndcap)  histos_["hGMTRates"]->Fill(3);  // 3 is GMT encap

      for (int ipt=0; ipt<=gmtPt;++ipt) 
        {
	  histos_["hGMTRatevsPt"] ->Fill(ipt);
	  if (isEndcap)  histos_["hGMTfRatevsPt"]->Fill(ipt); 
	  if (isBarrel)  histos_["hGMTbRatevsPt"]->Fill(ipt); 
        }
      
    }

}

void GMTRatePlotter::plot() 
{ 

    setTDRStyle();
    gStyle->SetOptTitle(0);

    TCanvas *cRates = new TCanvas((baseName_+"cGMTRates").c_str(),
            (baseName_+"cGMTRates").c_str(),500,500);
    cRates->cd();
    cRates->SetGrid();  

    histos_["hGMTRates"]->GetXaxis()->SetBinLabel(1,"GMT");
    histos_["hGMTRates"]->GetXaxis()->SetBinLabel(2,"GMT barrel");
    histos_["hGMTRates"]->GetXaxis()->SetBinLabel(3,"GMT forward");
    histos_["hGMTRates"]->GetXaxis()->SetBinLabel(4,"DT");
    histos_["hGMTRates"]->GetXaxis()->SetBinLabel(5,"RPC barrel");
    histos_["hGMTRates"]->GetXaxis()->SetBinLabel(6,"RPC forward");
    histos_["hGMTRates"]->GetXaxis()->SetBinLabel(7,"CSC");
    histos_["hGMTRates"]->GetYaxis()->SetTitle(isCrossSection_ ? 
            "x-section [#mub]" : "rate [arb. units]");
    histos_["hGMTRates"]->GetYaxis()->SetTitleOffset(0.9);
    histos_["hGMTRates"]->GetYaxis()->SetRange(0.,-1);

    histos_["hGMTRates"]->SetFillColor(4);
    
    histos_["hGMTRates"]->Draw();

    printHisto(cRates,"nGMTRates");


    TCanvas *cGMTRatevsPt = new TCanvas((baseName_+"cGMTRatevsPt").c_str(), (baseName_+"cGMTRatevsPt").c_str(),500,500);
    cGMTRatevsPt->SetLogy();
    cGMTRatevsPt->SetGrid();
    cGMTRatevsPt->cd();   
    histos_["hGMTRatevsPt"] ->SetLineWidth(2);
    histos_["hGMTRatevsPt"] ->GetXaxis()->SetTitle("pt [GeV]");
    histos_["hGMTRatevsPt"] ->GetYaxis()->SetTitle(isCrossSection_ ?  "x-section [#mub]" : "rate [arb. units]");
    histos_["hGMTRatevsPt"] ->Draw(); //cGMTRatevsPt ->SaveAs((baseName_+"nGMTRatevsPt.png" ).c_str());

    histos_["hGMTfRatevsPt"]->SetLineColor(kRed);
    histos_["hGMTfRatevsPt"]->SetLineWidth(2);
    histos_["hGMTfRatevsPt"]->GetXaxis()->SetTitle("pt [GeV]");
    histos_["hGMTfRatevsPt"]->GetYaxis()->SetTitle(isCrossSection_ ?  "x-section [#mub]" : "rate [arb. units]");
    histos_["hGMTfRatevsPt"]->Draw("same"); //cGMTfRatevsPt->SaveAs((baseName_+"nGMTfRatevsPt.png").c_str());
    
    histos_["hGMTbRatevsPt"]->SetLineColor(kOrange);
    histos_["hGMTbRatevsPt"]->SetLineWidth(2);
    histos_["hGMTbRatevsPt"]->GetXaxis()->SetTitle("pt [GeV]");
    histos_["hGMTbRatevsPt"]->GetYaxis()->SetTitle(isCrossSection_ ?  "x-section [#mub]" : "rate [arb. units]");
    histos_["hGMTbRatevsPt"]->Draw("same"); //cGMTbRatevsPt->SaveAs((baseName_+"nGMTbRatevsPt.png").c_str());
    
    histos_["hDTRatevsPt"]  ->SetLineColor(kBlue);
    histos_["hDTRatevsPt"]  ->SetLineWidth(2);
    histos_["hDTRatevsPt"]  ->GetXaxis()->SetTitle("pt [GeV]");
    histos_["hDTRatevsPt"]  ->GetYaxis()->SetTitle(isCrossSection_ ?  "x-section [#mub]" : "rate [arb. units]");
    histos_["hDTRatevsPt"]  ->Draw("same"); //cDTRatevsPt  ->SaveAs((baseName_+"nDTRatevsPt.png"  ).c_str());
    
//     histos_["hRPCfRatevsPt"]->SetLineColor(kYellow);
//     histos_["hRPCfRatevsPt"]->GetXaxis()->SetTitle("pt [GeV]");
//     histos_["hRPCfRatevsPt"]->GetYaxis()->SetTitle(isCrossSection_ ?  "x-section [#mub]" : "rate [arb. units]");
//     histos_["hRPCfRatevsPt"]->Draw("same"); //cRPCfRatevsPt->SaveAs((baseName_+"nRPCfRatevsPt.png").c_str());
    
//     histos_["hRPCbRatevsPt"]->SetLineColor(kPink);
//     histos_["hRPCbRatevsPt"]->GetXaxis()->SetTitle("pt [GeV]");
//     histos_["hRPCbRatevsPt"]->GetYaxis()->SetTitle(isCrossSection_ ?  "x-section [#mub]" : "rate [arb. units]");
//     histos_["hRPCbRatevsPt"]->Draw("same"); //cRPCbRatevsPt->SaveAs((baseName_+"nRPCbRatevsPt.png").c_str());
    
//     histos_["hCSCRatevsPt"] ->SetLineColor(kGreen);
//     histos_["hCSCRatevsPt"] ->GetXaxis()->SetTitle("pt [GeV]");
//     histos_["hCSCRatevsPt"] ->GetYaxis()->SetTitle(isCrossSection_ ?  "x-section [#mub]" : "rate [arb. units]");
//     histos_["hCSCRatevsPt"] ->Draw("same"); //cCSCRatevsPt ->SaveAs((baseName_+"nCSCRatevsPt.png" ).c_str());

    printHisto(cGMTRatevsPt,"nGMTRatevsPt" );    

}


void GMTRatePlotter::scale(float scaleFactor) 
{ 

    hTH1MapIt histoIt  = histos_.begin();
    hTH1MapIt histoEnd = histos_.end();

    for(;histoIt!=histoEnd;++histoIt)
        histoIt->second->Scale(scaleFactor);

    return;

}

// ********************************
// (5)
// GMT Barrel Upgrade Plotter class 
// (pt dependent plots and rates in
// barrel)
// ********************************

class GMTBarrelRatePtPlotter : public Plotter {

public:
  GMTBarrelRatePtPlotter(TFile* file, bool isCrossSection) ;
  ~GMTBarrelRatePtPlotter(); 
  
  void fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
	    int igmt, int idt, int irpcb, int icsc, 
	    int irpcf, TriggeredMuon tMu) ;
  void config(std::string baseName = "");
  void scale(float scaleFactor) ;
  void plot();
  
};


GMTBarrelRatePtPlotter::GMTBarrelRatePtPlotter(TFile* file, bool isCrossSection)
{ 
  
  outFile_ = file;
  isCrossSection_ = isCrossSection;

};

GMTBarrelRatePtPlotter::~GMTBarrelRatePtPlotter()
{

  histos_.clear();
  
};

void GMTBarrelRatePtPlotter::config(std::string baseName) 
{ 
  
  baseName_ = baseName;
  
  if (baseName_!="")
    {
      outFile_->mkdir(baseName_.c_str());
      outFile_->cd(baseName_.c_str());
      system(string("mkdir -p plots/" + baseName_).c_str());
    }

  string name = (baseName_ + "_hRatevsPt");  
  histos_["hRatevsPt"] = new TH1F(name.c_str(),name.c_str(),81,-0.5,80.5);

  name = (baseName_ + std::string("_hRateQualityvsPt"));  
  histos_["hRateQualityvsPt"] = new TH2F(name.c_str(),name.c_str(),81,-0.5,80.5,7,0.5,7.5);

}

void GMTBarrelRatePtPlotter::fill(L1Analysis::L1AnalysisGMTDataFormat * gmt, 
				  int igmt, int idt, int irpcb, int icsc, 
				  int irpcf, TriggeredMuon tMu) 
{

  float pt   = gmt->Pt[igmt];
  //float eta  = gmt->Eta[igmt];
  float qual = gmt->Qual[igmt];

  for (int ipt=0; ipt<=pt;++ipt) 
    {
      histos_["hRatevsPt"]->Fill(ipt);
      histos_["hRateQualityvsPt"]->Fill(ipt,qual);
    }
  
  return; 
  
}

void GMTBarrelRatePtPlotter::plot() 
{ 
  
  setTDRStyle();
  gStyle->SetOptTitle(0);
  
  TCanvas *cRatevsPt = new TCanvas((baseName_+"cRatevsPt").c_str(),
				   (baseName_+"cRatevsPt").c_str(),500,500);
  cRatevsPt->cd();
  cRatevsPt->SetLogy();
  cRatevsPt->SetGrid();  

  histos_["hRatevsPt"]->GetXaxis()->SetTitle("GMTBarrel Trigger p_{t} [GeV/c]");
  histos_["hRatevsPt"]->GetXaxis()->SetTitleOffset(0.9);
  histos_["hRatevsPt"]->GetYaxis()->SetTitle(isCrossSection_ ? 
					     "x-section [#mub]" : "rate [arb. units]");
  histos_["hRatevsPt"]->SetLineColor(kBlue);
  histos_["hRatevsPt"]->SetLineWidth(2);
  histos_["hRatevsPt"]->Draw();
  
  printHisto(cRatevsPt,"RatevsPt");

  int ptBins[5] = {12, 16, 20, 25, 30}; // CB adapt pt thresholds if needed
  
  for (int bin=0; bin<5; ++bin) 
    {
      std::cout << "pt : " << ptBins[bin]
		<< " rate (full GMTBarrel Trigger eta) : " 
                << histos_["hRatevsPt"]->GetBinContent(histos_["hRatevsPt"]->FindBin(ptBins[bin]))
		<< std::endl;
    }

  TCanvas *cRateQualityvsPt = new TCanvas((baseName_+"cRateQualityvsPt").c_str(),
					  (baseName_+"cRateQualityvsPt").c_str(),500,500);
  cRateQualityvsPt->cd();
  cRateQualityvsPt->SetLogy();
  cRateQualityvsPt->SetGrid();


  THStack * hRatevsPtQualityStack = new THStack("hRatevsPtQualityStack", "hRatevsPtQualityStack");
  THStack * hRatevsPtQualityContriStack = new THStack("hRatevsPtQualityContriStack", "hRatevsPtQualityContriStack");
  for(int ybin=1; ybin<=histos_["hRateQualityvsPt"]->GetNbinsY(); ++ybin)
    {
      stringstream px ;
      px<<ybin<<"QualityGmtBarrel";
      TH1D* projection = static_cast<TH2*>(histos_["hRateQualityvsPt"]) ->ProjectionX(px.str().c_str(),ybin,ybin);

      projection->SetFillColor(colorMap[ybin-1]);
      projection->SetLineColor(colorMap[ybin-1]);
      projection->GetXaxis()->SetTitle("GMTBarrel Trigger p_{t} [GeV/c]");
      projection->GetYaxis()->SetTitle(isCrossSection_ ?
				       "GMTBarrel Trigger x-section [#mub]" :
				       "rate [arb. units]");

      hRatevsPtQualityStack->Add(projection);

      px<<"ContriGmtBarrel";
      TH1D* contriProj = new TH1D(px.str().c_str(),px.str().c_str(),81,-0.5,80.5); 
      for(int xbin=1; xbin<=81; ++xbin)
        {
	  contriProj->SetBinContent(xbin, 1.*projection->GetBinContent(xbin)/histos_["hRatevsPt"]->GetBinContent(xbin));
        }
      contriProj->SetFillColor(colorMap[ybin-1]);
      contriProj->SetLineColor(colorMap[ybin-1]);
      contriProj->GetXaxis()->SetTitle("GMTBarrel Trigger p_{t} [GeV/c]");
      contriProj->GetYaxis()->SetTitle("Contribution to total");
      
      hRatevsPtQualityContriStack->Add(contriProj);
      
    }

  hRatevsPtQualityStack->Draw();
  hRatevsPtQualityStack->GetHistogram()->GetXaxis()->SetTitle("GMTBarrel Trigger p_{t} [GeV/c]");
  hRatevsPtQualityStack->GetHistogram()->GetYaxis()->SetTitle(isCrossSection_ ? 
							      "x-section [#mub]" : "rate [arb. units]");
  hRatevsPtQualityStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
  hRatevsPtQualityStack->Draw("Asame");
  
  cRateQualityvsPt->Update();
  cRateQualityvsPt->Modified();
  
  printHisto(cRateQualityvsPt,"RateQualityvsPt");
  
  TCanvas *cRateQualityvsPtContri = new TCanvas((baseName_+"cRateQualityvsPtContri").c_str(),
						(baseName_+"cRateQualityvsPtContri").c_str(),500,500);
  cRateQualityvsPtContri->cd();
  cRateQualityvsPtContri->SetGrid();
  
  hRatevsPtQualityContriStack->Draw();
  hRatevsPtQualityContriStack->GetHistogram()->GetXaxis()->SetTitle("GMTBarrel Trigger p_{t} [GeV/c]");
  hRatevsPtQualityContriStack->GetHistogram()->GetYaxis()->SetTitle("Contribution to total");
  hRatevsPtQualityContriStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
  hRatevsPtQualityContriStack->GetHistogram()->GetXaxis()->SetRangeUser(0,40);
  hRatevsPtQualityContriStack->Draw("Asame");
  
  cRateQualityvsPtContri->Update();
  cRateQualityvsPtContri->Modified();
  
  printHisto(cRateQualityvsPtContri,"RateQualityvsPtContri");
  
  return;

}


void GMTBarrelRatePtPlotter::scale(float scaleFactor) 
{ 

    hTH1MapIt histoIt  = histos_.begin();
    hTH1MapIt histoEnd = histos_.end();

    for(;histoIt!=histoEnd;++histoIt)
        histoIt->second->Scale(scaleFactor);

    return;

}






