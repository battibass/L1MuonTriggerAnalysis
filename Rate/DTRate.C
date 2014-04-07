#include "../../UserCode/L1TriggerDPG/macros/L1Ntuple.h" // CB fix path
#include "CommonUtils.C" // CB fix path
#include "../CommonUtils/TriggeredMuon.C" // CB fix path
#include "DTControlPlotter.C"
#include "TriggeredMuons.C" // CB fix path
#include "DTRateAlgo.C"

// --------------------------------------------
// Configuration parameters to perform scaling 
// --------------------------------------------

// Number of processed lumisections
#define N_LUMIS 70

// Lumi from lumiCalc2.py corrected for PU (1^e33)
#define LUMI_REF .0185 

// the seed bit prescale
#define SEED_PRESCALE 23.

//the prescale value of the stream
#define STREAM_PRESCALE 1


// ----------------------------------------------
// DTRate steering (configuration and run) class 
// ----------------------------------------------

class DTRate : public L1Ntuple
{
    public :

        //constructor    
        DTRate(std::string filename) : L1Ntuple(filename) {}
        DTRate()  {}
        ~DTRate() {}

        void runRates(int nEvents = 0);

    private :

        void FillBits();
        float ScaleFactor(float targetLumi); // targetLumi is the reference luminosity 
                                             // (in units of 1e33) to which you want to
                                             // scale the rates

        bool PhysicsBits[128];

};

void DTRate::FillBits() {

    for (int ibit=0; ibit < 128; ibit++) 
    {
        PhysicsBits[ibit] = 0;
        if (ibit<64) 
        {
            PhysicsBits[ibit] = (gt_->tw1[2]>>ibit)&1;
        }
        else 
        {
            PhysicsBits[ibit] = (gt_->tw2[2]>>(ibit-64))&1;
        }
    }

}       

// ------------------------------------------------------------------

float DTRate::ScaleFactor(float targetLumi) {

    float scal = 1./(23.3) ;       // 1 LS
    scal = scal/(float)N_LUMIS ;  //  the ntuple has NLUMIS LSs
    scal = scal * STREAM_PRESCALE;      // nanoDST is p'ed by 10
    scal = scal * targetLumi / LUMI_REF;
    scal = scal / 1000. ;  // rate in kHz / mub 
    scal = scal * SEED_PRESCALE; // scaling for the prescale of the seed bit 

    return scal;

}

// --------------------------------------------------------------------
//                             run function 
// --------------------------------------------------------------------


void DTRate::runRates(int nEvents) {

    system("mkdir -p results");
    TFile *outFile = new TFile("results/ratePlots.root","recreate");

    Plotter* plotterDtPt  = new DTRatePtPlotter(outFile,false);
    plotterDtPt->config("DT_Pt");

    Plotter* plotterDtEtaPhi16  = new DTRateEtaPhiPlotter(outFile,false);
    plotterDtEtaPhi16->config("DT_EtaPhi_Pt16");

    Plotter* plotterDtEtaPhi20  = new DTRateEtaPhiPlotter(outFile,false);
    plotterDtEtaPhi20->config("DT_EtaPhi_Pt20");

    Plotter* plotterDtEtaPhi25  = new DTRateEtaPhiPlotter(outFile,false);
    plotterDtEtaPhi25->config("DT_EtaPhi_Pt25");

    Plotter* plotterGMTPt = new GMTRatePlotter(outFile,false);
    plotterGMTPt->config("GMT_Pt");

    Plotter* plotterGMTPt16 = new GMTRatePlotter(outFile,false);
    plotterGMTPt16->config("GMT_Pt16");

    Plotter* plotterGMTPt20 = new GMTRatePlotter(outFile,false);
    plotterGMTPt20->config("GMT_Pt20");

    Plotter* plotterGMTPt25 = new GMTRatePlotter(outFile,false);
    plotterGMTPt25->config("GMT_Pt25");

    Plotter* plotterGMTBarrelPt = new GMTBarrelRatePtPlotter(outFile,false);
    plotterGMTBarrelPt->config("GMTBarrel_Pt");

    float scaleFactor = 1.; // ScaleFactor(7.5); CB now printing bare counts
    
    std::vector<RateAlgo*> muTrigs;

    muTrigs.push_back(new RateAlgoDT(1,plotterDtPt,0,0.85));
    muTrigs.push_back(new RateAlgoDT(1,plotterDtEtaPhi16,16,0.85));
    muTrigs.push_back(new RateAlgoDT(1,plotterDtEtaPhi20,20,0.85));
    muTrigs.push_back(new RateAlgoDT(1,plotterDtEtaPhi25,25,0.85));

    muTrigs.push_back(new RateAlgoGMT(1,plotterGMTPt,0,2.10));
    muTrigs.push_back(new RateAlgoGMT(1,plotterGMTPt16,16,2.10));
    muTrigs.push_back(new RateAlgoGMT(1,plotterGMTPt20,20,2.10));
    muTrigs.push_back(new RateAlgoGMT(1,plotterGMTPt25,25,2.10));
    muTrigs.push_back(new RateAlgoGMT(1,plotterGMTBarrelPt,0,0.85));

    int nevents = nEvents == 0 ? GetEntries() : nEvents;
    int counts = 0;

    std::cout << "Running on " << nevents << " events." << std::endl;
    for (Long64_t event=0; event<nevents; ++event)
    { 
        Long64_t eventEntry = LoadTree(event); 
        if (eventEntry < 0) break;
        GetEntry(event);

        if (event%200000 == 0) {
	  std::cout << "Processed " << event << " events." << std::endl;
        }

        FillBits();
	
        bool hasTrigger = false;
	
        std::vector<RateAlgo*>::iterator muTrigIt  = muTrigs.begin();
        std::vector<RateAlgo*>::iterator muTrigEnd = muTrigs.end();
	
        for (;muTrigIt!=muTrigEnd;++muTrigIt) {
	  hasTrigger |= (*muTrigIt)->runAlgo(PhysicsBits,recoMuon_,gmt_,dttf_);
        }
	
        if (hasTrigger) counts ++;
    } // end event loop
    
    
    std::cout << "Overall counts : " << counts << std::endl; 
    
    std::vector<RateAlgo*>::iterator muTrigIt  = muTrigs.begin();
    std::vector<RateAlgo*>::iterator muTrigEnd = muTrigs.end();
    
    std::cout << std::endl;
    
    for (;muTrigIt!=muTrigEnd;++muTrigIt) {
      cout << (*muTrigIt)->trigger() << " :\t" << (*muTrigIt)->getRate(scaleFactor)  << endl;
    }
    
    plotterDtPt->plot();

    plotterDtEtaPhi16->plot();
    plotterDtEtaPhi20->plot();
    plotterDtEtaPhi25->plot();
    
    plotterGMTPt->plot();
    plotterGMTBarrelPt->plot();

    plotterGMTPt16->plot();
    plotterGMTPt20->plot();
    plotterGMTPt25->plot();

    outFile->Write();

}

// --------------------------------------------------------------------

void goRates(int nEvents = 0) 
{

  //DTRate dtRate("/data2/battilan/L1Trigger/L1Tree_Run2012C_Commissioning_SingleMuOpenSingleMu12_NewDTTFLUTs_RPC_v10.root");
  DTRate dtRate("/data2/battilan/L1Trigger/L1Tree_Run2012C_Commissioning_SingleMuOpenSingleMu12_RPC_v7.root");
  dtRate.runRates(nEvents);

}

