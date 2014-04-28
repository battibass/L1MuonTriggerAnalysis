#include "../../L1TriggerDPG/L1Ntuples/macros/L1Ntuple.h"
#include "Configuration.C"
#include "CommonUtils.C"
#include "../CommonUtils/TriggeredMuon.C"
#include "ControlPlotter.C"
#include "TriggeredMuons.C"
#include "EfficiencyPlotter.C"

// --------------------------------------------------------------------
//                       GmtEfficiency macro definition
// --------------------------------------------------------------------

class GmtEfficiency : public L1Ntuple
{
public :
  
  //constructor    
  GmtEfficiency(std::string filename) : L1Ntuple(filename) {}
  GmtEfficiency()  {}
  ~GmtEfficiency() {}

  void runEfficiency(int nEvents,  std::string outFileName);
  

};
        

// --------------------------------------------------------------------
//                             run function 
// --------------------------------------------------------------------


void GmtEfficiency::runEfficiency(int nEvents, std::string outFileName) {

  system("mkdir -p results");
  TFile *outFile = new TFile(std::string("results/" + outFileName + ".root").c_str() ,"recreate");

  std::vector<EfficiencyPlotter *> plotters;
  plotters.push_back(new EfficiencyPlotter(outFile,outFileName,"Pt16",16));
  plotters.push_back(new EfficiencyPlotter(outFile,outFileName,"Pt20",20));
  plotters.push_back(new EfficiencyPlotter(outFile,outFileName,"Pt25",25));
  plotters.push_back(new EfficiencyPlotter(outFile,outFileName,"Pt30",30));

  ControlPlotter *controlPlots = new ControlPlotter(outFile,outFileName,"ControlPlots");

  int nevents = nEvents == 0 ? GetEntries() : nEvents;
        
  std::cout << "File size is " << nevents << " events." << std::endl;
  for (Long64_t event=0; event<nevents; ++event)
    { 
      Long64_t eventEntry = LoadTree(event); 
      if (eventEntry < 0) break;
      GetEntry(event);

      if (event%50000 == 0) 
	std::cout << "Processed " << event << " events." << std::endl;

      TriggeredMuons trigMuons(recoMuon_,gmt_);
      
      trigMuons.findTightMuons(controlPlots)     && 
      trigMuons.findProbes()                     &&
      trigMuons.runTriggerMatching(controlPlots);
  
      triggeredMuonsIt trigMuonsIt  = trigMuons.my_triggered_muons.begin();
      triggeredMuonsIt trigMuonsEnd = trigMuons.my_triggered_muons.end();

      for(;trigMuonsIt!=trigMuonsEnd;++trigMuonsIt)
	{
	  std::vector<EfficiencyPlotter*>::const_iterator plotterIt  = plotters.begin();
	  std::vector<EfficiencyPlotter*>::const_iterator plotterEnd = plotters.end();
	  
	  for(;plotterIt!=plotterEnd;++plotterIt)
	    {
	      (*plotterIt)->fill(trigMuonsIt,(doreco ? recoVertex_ : 0));
	    }

	}

    } // end event loop

  std::vector<EfficiencyPlotter*>::const_iterator plotterIt  = plotters.begin();
  std::vector<EfficiencyPlotter*>::const_iterator plotterEnd = plotters.end();
  
  for(;plotterIt!=plotterEnd;++plotterIt)
    {
      (*plotterIt)->plotAndSave();
    }

  plotAndSaveAll(plotters,(outFileName + "/All").c_str(),"EffVsPt");
  plotAndSaveAll(plotters,(outFileName + "/All").c_str(),"EffVsPtBarrel");
  plotAndSaveAll(plotters,(outFileName + "/All").c_str(),"EffVsPtOverlap");
  plotAndSaveAll(plotters,(outFileName + "/All").c_str(),"EffVsPtEndcap");
  plotAndSaveAll(plotters,(outFileName + "/All").c_str(),"EffVsEta");
  plotAndSaveAll(plotters,(outFileName + "/All").c_str(),"EffVsPhi");

  controlPlots->plotAndSave();

  outFile->Write();

}

// --------------------------------------------------------------------

void goEfficiency(std::string ntupleFilePath, std::string outFileName, int nEvents = 0) 
{
  
  GmtEfficiency analyzerPhysics(ntupleFilePath.c_str()); 
  analyzerPhysics.runEfficiency(nEvents, outFileName);
  
}

