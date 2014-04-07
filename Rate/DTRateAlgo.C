#include "DTRatePlotter.C"
#include "../../UserCode/L1TriggerDPG/interface/L1AnalysisGMTDataFormat.h" // CB FIX path
#include "../../UserCode/L1TriggerDPG/interface/L1AnalysisDTTFDataFormat.h" // CB FIX path
#include "../../UserCode/L1TriggerDPG/interface/L1AnalysisRecoMuonDataFormat.h" // CB FIX path


// --------------------------------
// Generic Muon Trigger class
// (virtual, not possible to make
// an instance of this class)
// --------------------------------

class RateAlgo
{

    public :
        RateAlgo();
        RateAlgo(const RateAlgo & l1RateMu);
        ~RateAlgo() {};

        bool increment() {counts_++; return  (counts_%prescale_) == 0; };
        int getPrescaledCounts() { return counts_/prescale_; };

        // virtual function, can be changes by the specific classes  
        virtual bool runAlgo(bool PhysicsBits[128], 
			     L1Analysis::L1AnalysisRecoMuonDataFormat * mu,
			     L1Analysis::L1AnalysisGMTDataFormat      * gmt,
			     L1Analysis::L1AnalysisDTTFDataFormat     * dttf) = 0;
        virtual std::string trigger() { return ""; };
        virtual float getRate(float scaleFactor) ; 

    protected :
        Plotter* plotter_;
        int prescale_;
        float trigPt_;
        float trigEta_;
        int counts_;

};

RateAlgo::RateAlgo() :
    plotter_(0) , prescale_(0) , trigPt_(0) , trigEta_(0) , counts_(0) {};

RateAlgo::RateAlgo(const RateAlgo & l1RateMu) 
{

    plotter_  = l1RateMu.plotter_;
    prescale_ = l1RateMu.prescale_;
    trigPt_  = l1RateMu.trigPt_;
    trigEta_  = l1RateMu.trigEta_;
    counts_   = l1RateMu.counts_;

}

float RateAlgo::getRate(float scaleFactor) 
{ 

    if (plotter_) 
        plotter_->scale(scaleFactor) ;

    return scaleFactor*getPrescaledCounts(); 

}


// --------------------------------
// DT Muon rate class
// (perform detailed DTTF rate
// studies based on pt/eta cuts)
// --------------------------------

class RateAlgoDT : public RateAlgo {

    public :
        RateAlgoDT(int prescale, Plotter *plotter, float trigPt, float trigEta);
        RateAlgoDT(const RateAlgoDT & l1RateMu);

        std::string trigger();
        bool runAlgo(bool PhysicsBits[128], 
		     L1Analysis::L1AnalysisRecoMuonDataFormat * mu,
		     L1Analysis::L1AnalysisGMTDataFormat      * gmt,
		     L1Analysis::L1AnalysisDTTFDataFormat     * dttf);
        virtual float getRate(float scaleFactor);

  ControlPlotter * controlPlotter_;
  

};

RateAlgoDT::RateAlgoDT(int prescale, Plotter *plotter, float trigPt, float trigEta)
{

    plotter_  = plotter;
    prescale_ = prescale;
    trigPt_  = trigPt;
    trigEta_  = trigEta;
    counts_   = 0;

    system("mkdir -p results");
    TFile *outFile = new TFile(("results/controlPlots_" + trigger() + ".root").c_str(),"recreate");
    controlPlotter_ = new ControlPlotter(outFile,trigger());

}

RateAlgoDT::RateAlgoDT(const RateAlgoDT & l1RateMu) : RateAlgo(l1RateMu) 
{

    plotter_  = l1RateMu.plotter_;
    controlPlotter_  = l1RateMu.controlPlotter_;
    prescale_ = l1RateMu.prescale_;
    trigPt_  = l1RateMu.trigPt_;
    trigEta_  = l1RateMu.trigEta_;
    counts_   = l1RateMu.counts_;

}

std::string RateAlgoDT::trigger() {

    stringstream triggerName;
    triggerName << "L1DTTrigger_pt" << trigPt_;
    return triggerName.str();

}

float RateAlgoDT::getRate(float scaleFactor) 
{ 

  controlPlotter_->plot();
  return RateAlgo::getRate(scaleFactor) ;

}

bool RateAlgoDT::runAlgo(bool PhysicsBits[128],
			 L1Analysis::L1AnalysisRecoMuonDataFormat * mu,
			 L1Analysis::L1AnalysisGMTDataFormat      * gmt,
			 L1Analysis::L1AnalysisDTTFDataFormat     * dttf) 
{

    float PTmax  = -1;
    TriggeredMuon tMuPtMax(0,-1,0);
    TriggeredMuons tMu(mu,gmt,dttf);

    if ( PhysicsBits[55] ) // CB single mu open 
    {	
      tMu.findMuons(controlPlotter_) && tMu.runGmtDttfMatching(controlPlotter_) && tMu.runDTTriggerMatching(controlPlotter_);
      
      triggeredMuonsIt tMuIt  = tMu.dt_trigger_muons.begin();
      triggeredMuonsIt tMuEnd = tMu.dt_trigger_muons.end();
  
      for (; tMuIt != tMuEnd; ++tMuIt) 
        {
	  
	  int bx = tMuIt->gmt()->Bxdt[tMuIt->igmt()];
            if (bx != 0) continue;

	    float eta= tMuIt->gmt()->Etadt[tMuIt->igmt()];
	    if (fabs(eta) > trigEta_) continue;

            float pt = tMuIt->gmt()->Ptdt[tMuIt->igmt()];
            if ( pt > PTmax) 
            {
                PTmax = pt;
                tMuPtMax = (*tMuIt);
            }
        }	
    }    

    if (plotter_ && PTmax >= trigPt_)
      plotter_->fill(gmt,-1,tMuPtMax.igmt(),-1,-1,-1,tMuPtMax);

    return (PTmax >= trigPt_) ? increment() : false ;

}


// bool RateAlgoDT::runAlgo(bool PhysicsBits[128],
// 			 L1Analysis::L1AnalysisRecoMuonDataFormat * mu,
// 			 L1Analysis::L1AnalysisGMTDataFormat      * gmt,
// 			 L1Analysis::L1AnalysisDTTFDataFormat     * dttf) 
// {

//     float PTmax  = -1;
//     int   iPTmax = -1;
//     if ( PhysicsBits[55] ) // CB single mu open 
//     {	
//         int NmuDt = gmt -> Ndt;
//         for (int imu=0; imu < NmuDt; imu++) 
//         {
//             int bx = gmt -> Bxdt[imu];
//             if (bx != 0) continue;

// 	    float eta= gmt -> Etadt[imu];
// 	    if (fabs(eta) > trigEta_) continue;

//             float pt = gmt -> Ptdt[imu];
//             if ( pt > PTmax) 
//             {
//                 PTmax = pt;
//                 iPTmax = imu;
//             }
//         }	
//     }

//     if (plotter_ && PTmax >= trigPt_)
//       plotter_->fill(gmt,-1,iPTmax,-1,-1,-1);

//     return (PTmax >= trigPt_) ? increment() : false ;

// }


//  --------------------------------
// Global muon trigger rate class
// (triggers when there is at least
// one muon candidate from DT, RPC,
// CSC or a GMT muon with 
// sufficient pt)
// --------------------------------

class RateAlgoGMT : public RateAlgo {

    public :
        RateAlgoGMT(int prescale, Plotter *plotter, float trigPt, float trigEta);
        RateAlgoGMT(const RateAlgoGMT & l1RateMu);

        std::string trigger();
        bool runAlgo(bool PhysicsBits[128],
			     L1Analysis::L1AnalysisRecoMuonDataFormat * mu,
			     L1Analysis::L1AnalysisGMTDataFormat      * gmt,
			     L1Analysis::L1AnalysisDTTFDataFormat     * dttf);
        virtual float getRate(float scaleFactor);

};

RateAlgoGMT::RateAlgoGMT(int prescale, Plotter *plotter, 
			 float trigPt, float trigEta   )
{

    plotter_  = plotter;
    prescale_ = prescale;
    trigPt_  = trigPt;
    trigEta_  = trigEta;
    counts_   = 0;

}

RateAlgoGMT::RateAlgoGMT(const RateAlgoGMT & l1RateMu) : RateAlgo(l1RateMu) 
{

    plotter_  = l1RateMu.plotter_;
    prescale_ = l1RateMu.prescale_;
    trigPt_  = l1RateMu.trigPt_;
    trigEta_  = l1RateMu.trigEta_;
    counts_   = l1RateMu.counts_;

}

std::string RateAlgoGMT::trigger() {

    stringstream triggerName;
    triggerName << "L1GmtMu_pt" << trigPt_;

    return triggerName.str();

}

float RateAlgoGMT::getRate(float scaleFactor) 
{ 

    return RateAlgo::getRate(scaleFactor) ;

}

bool RateAlgoGMT::runAlgo(bool PhysicsBits[128],
			  L1Analysis::L1AnalysisRecoMuonDataFormat * mu,
			  L1Analysis::L1AnalysisGMTDataFormat      * gmt,
			  L1Analysis::L1AnalysisDTTFDataFormat     * dttf) 
{

    float gmtPTmax  = -1;
    int   iGmtPTmax = -1;

    float dtPTmax  = -1;
    int   iDtPTmax = -1;

    float rpcbPTmax  = -1;
    int   iRpcbPTmax = -1;

    float cscPTmax  = -1;
    int   iCscPTmax = -1;

    float rpcfPTmax  = -1;
    int   iRpcfPTmax = -1;

    if ( PhysicsBits[55] ) // CB single mu open
    {	
        int Nmu = gmt -> N;
        for (int imu=0; imu < Nmu; imu++) 
        {
            int bx = gmt -> CandBx[imu];
            if (bx != 0) continue;

            int qual = gmt -> Qual[imu]; // CB single mu quality cut
            if (qual < 4) continue;

	    float eta = gmt -> Eta[imu];
	    if (fabs(eta) > trigEta_) continue;

            float pt = gmt -> Pt[imu];
            if ( pt > gmtPTmax) 
            {
                gmtPTmax = pt;
                iGmtPTmax = imu;
            }
        }

        int NmuDt = gmt -> Ndt;
        for (int imu=0; imu < NmuDt; imu++) 
        {
            int bx = gmt -> Bxdt[imu];
            if (bx != 0) continue;

	    float eta = gmt -> Etadt[imu];
	    if (fabs(eta) > trigEta_) continue;

            float pt = gmt -> Ptdt[imu];
            if ( pt > dtPTmax) 
            {
                dtPTmax = pt;
                iDtPTmax = imu;
            }
        }	

        int NmuRpcb = gmt -> Nrpcb;
        for (int imu=0; imu < NmuRpcb; imu++) 
        {
            int bx = gmt -> Bxrpcb[imu];
            if (bx != 0) continue;

	    float eta = gmt -> Etarpcb[imu];
	    if (fabs(eta) > trigEta_) continue;

	    float pt = gmt -> Ptrpcb[imu];
            if ( pt > rpcbPTmax) 
            {
                rpcbPTmax = pt;
                iRpcbPTmax = imu;
            }
        }	

        int NmuCsc = gmt -> Ncsc;
        for (int imu=0; imu < NmuCsc; imu++) 
        {
            int bx = gmt -> Bxcsc[imu];
            if (bx != 0) continue;

	    float eta = gmt -> Etacsc[imu];
	    if (fabs(eta) > trigEta_) continue;

	    float pt = gmt -> Ptcsc[imu];
            if ( pt > cscPTmax) 
            {
                cscPTmax = pt;
                iCscPTmax = imu;
            }
        }	

        int NmuRpcf = gmt -> Nrpcf;
        for (int imu=0; imu < NmuRpcf; imu++) 
        {
            int bx = gmt -> Bxrpcf[imu];
            if (bx != 0) continue;

	    float eta = gmt -> Etarpcf[imu];
	    if (fabs(eta) > trigEta_) continue;

	    float pt = gmt -> Ptrpcf[imu];
            if ( pt > rpcfPTmax) 
            {
                rpcfPTmax = pt;
                iRpcfPTmax = imu;
            }
        }	

    }


    bool hasFiredGmt = (gmtPTmax  >= 0);

    if (gmtPTmax  < trigPt_) iGmtPTmax  = -1;
    if (dtPTmax   < trigPt_) iDtPTmax   = -1;
    if (rpcbPTmax < trigPt_) iRpcbPTmax = -1;
    if (cscPTmax  < trigPt_) iCscPTmax  = -1;
    if (rpcfPTmax < trigPt_) iRpcfPTmax = -1; 

    TriggeredMuon dummy(0,-1,0);

    if (plotter_ && ( gmtPTmax   >= trigPt_ ||
		      ( hasFiredGmt           &&   // CB emulates GMT quality cut    
			( dtPTmax    >= trigPt_ || // a trigger is accepred only if 
			  rpcbPTmax  >= trigPt_ || // it has fired GMT
			  cscPTmax   >= trigPt_ ||
			  rpcfPTmax  >= trigPt_
			)
		      )
		    ) 
	)
      plotter_->fill(gmt,iGmtPTmax,iDtPTmax,iRpcbPTmax,iCscPTmax,iRpcfPTmax,dummy);

    return ( gmtPTmax   >= trigPt_  ||
	     ( hasFiredGmt           &&   // CB emulates GMT quality cut    
	       ( dtPTmax    >= trigPt_ || // a trigger is accepred only if 
		 rpcbPTmax  >= trigPt_ || // it has fired GMT
		 cscPTmax   >= trigPt_ ||
		 rpcfPTmax  >= trigPt_
	       )
	     )
	   ) ? increment() : false ;

}
