// ***********************************
// Triggered object collection class
// ( Provide collections of triggered
//   muon objects, apply selections
//   for tight/hlt/probe matching
// )
// ***********************************

class TriggeredMuons
{

public :

  TriggeredMuons(L1Analysis::L1AnalysisRecoMuonDataFormat * mu,
		 L1Analysis::L1AnalysisGMTDataFormat      * gmt,
		 L1Analysis::L1AnalysisDTTFDataFormat     * dttf) : 
    mu_(mu), gmt_(gmt) , dttf_(dttf) { };   

  ~TriggeredMuons() 
  { 
    global_muons.clear();
    dt_trigger_muons.clear();
    gmt_dttf_muons.clear();
  }

  bool findMuons(ControlPlotter* plot); // find exactely 2 tight muons get 0 otherwise
  bool findProbes();     // find all valid probes out of tight muons (true if >=1)

  bool runGmtDttfMatching(ControlPlotter* plot); // attempt matching between GMT DTTF Trigger
  bool runDTTriggerMatching(ControlPlotter* plot); // attempt matching with DT Trigger

public :

  L1Analysis::L1AnalysisRecoMuonDataFormat * mu_;  
  L1Analysis::L1AnalysisGMTDataFormat      * gmt_;
  L1Analysis::L1AnalysisDTTFDataFormat     * dttf_;
  

  triggeredMuons global_muons;
  triggeredMuons dt_trigger_muons;
  gmtDttfMuons   gmt_dttf_muons;


};


bool TriggeredMuons::findMuons(ControlPlotter* plot)
{
      
  for(int iMu=0; iMu<mu_->nMuons; ++iMu)
    {
      
      if ( mu_->howmanytypes.at(iMu) & 0x01  )
	{
	  TriggeredMuon triggeredMu(mu_,iMu,0);
	  global_muons.push_back(triggeredMu);
	  if (triggeredMu.hasTightMuon()) plot->fillTight(triggeredMu);
	}
    }

  return true; // CB it is always true for rates
  
}

bool TriggeredMuons::runGmtDttfMatching(ControlPlotter *plot)
{

  for (int iDt=0; iDt<gmt_->Ndt; ++iDt)
    {

      if (abs(gmt_->Bxdt.at(iDt))>1)
	continue;

      GmtDttfMuon gmtDttf(0,-1,gmt_,iDt);

      for (int iDttf=0; iDttf<dttf_->trSize; ++iDttf)
	{
	  GmtDttfMuon tmpGmtDttf(dttf_,iDttf,gmt_,iDt);
	  if ( tmpGmtDttf.hasDtTriggerMatch() )
	    {
	      gmtDttf = tmpGmtDttf;
	      break;
	    }
	}
      
      if (gmtDttf.gmt_->Ptdt.at(gmtDttf.igmt_) >= 16) plot->fillGmtDttf(gmtDttf);
      gmt_dttf_muons.push_back(gmtDttf);

    }

  return gmt_dttf_muons.size() > 0;

}

bool TriggeredMuons::runDTTriggerMatching(ControlPlotter *plot)
{

  gmtDttfMuonsIt gmtDttfIt  = gmt_dttf_muons.begin();
  gmtDttfMuonsIt gmtDttfEnd = gmt_dttf_muons.end();
  
  for (; gmtDttfIt!=gmtDttfEnd; ++gmtDttfIt)
    {
      
      TriggeredMuon bestCand(0,-1,&(*gmtDttfIt));

      triggeredMuonsIt globalMuonIt  = global_muons.begin();
      triggeredMuonsIt globalMuonEnd = global_muons.end();

      for (;globalMuonIt!=globalMuonEnd;++globalMuonIt)
	{

	  TriggeredMuon cand(globalMuonIt->mu_,globalMuonIt->imu_,&(*gmtDttfIt));
	  
	  if (cand.hasDtTriggerMatch(0)) // CB event set to 0
	    {
	      if (bestCand.hasMuon()) 
		{
		  int dPhiCand     = cand.deltaPhi();
		  int dPhiBestCand = bestCand.deltaPhi();
		  
		  if (dPhiCand > dPhiBestCand)
		    bestCand = cand;

		} 
	      else
		{
		  bestCand = cand;
		}
		  
		  
	    }
	}

      plot->fillTrigger(bestCand);
      dt_trigger_muons.push_back(bestCand);

    }

  return true;

}
