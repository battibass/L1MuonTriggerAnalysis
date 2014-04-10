// ***********************************
// Triggered object collection class
// ( Provide collections of triggered muon objects, apply selections
// for tight/hlt/probe matching )
// ***********************************

class TriggeredMuons
{

public :

  TriggeredMuons(L1Analysis::L1AnalysisRecoMuonDataFormat * mu,
		 L1Analysis::L1AnalysisGMTDataFormat      * gmt) : 
    my_mu(mu), my_gmt(gmt) { };   

  ~TriggeredMuons() 
  { 
    my_tight_muons.clear();
    my_probe_muons.clear();
  }

  bool runTriggerMatching(ControlPlotter* plots); // attempt matching between RECO muon and GMT

  bool findTightMuons(ControlPlotter* plots);       // find exactely 2 tight muons far in eta and phi return 0 otherwise
  bool findProbes();                                // find all valid probes out of tight muons (true if >=1)

public :

  L1Analysis::L1AnalysisRecoMuonDataFormat * my_mu;  
  L1Analysis::L1AnalysisGMTDataFormat      * my_gmt;  

  triggeredMuons my_tight_muons;
  triggeredMuons my_probe_muons;
  triggeredMuons my_triggered_muons;

};


bool TriggeredMuons::findTightMuons(ControlPlotter* plots)
{
      
  for(int iMu=0; iMu<my_mu->nMuons; ++iMu)
    {
      
      if ( (my_mu->howmanytypes.at(iMu) >> 4) & 0x01 &&
	   (my_mu->howmanytypes.at(iMu) >> 5) & 0x01   )
	{
	  my_tight_muons.push_back(TriggeredMuon(my_mu,iMu,0,0));
	}
    }

  int nTightMu = my_tight_muons.size();  

  bool isGoodEvent = false;

  if (nTightMu == 2 && (effCompType == TWO_MUON_NO_TRIG || effCompType == COUNT_TNP) )
    {
      TriggeredMuon & mu1 = my_tight_muons.at(0);
      TriggeredMuon & mu2 = my_tight_muons.at(1);
      
      plots->fillTight(mu1);
      plots->fillTight(mu2);
      
      float fDPhi = fabs(acos(cos( mu1.my_mu->sa_phi_mb2[mu1.my_imu] - 
			           mu2.my_mu->sa_phi_mb2[mu2.my_imu] ) ) );
      
      float fDEta = fabs(acos(cos( mu1.my_mu->eta[mu1.my_imu] - 
			           mu2.my_mu->eta[mu2.my_imu] ) ) );

      isGoodEvent = fDPhi > MAX_MU_MU_DPHI &&
	            fDEta > MAX_MU_MU_DETA ;
    }
  else if (nTightMu == 1 && effCompType == ONE_MUON_NO_TRIG )
    {
      isGoodEvent = true;
    }

  return isGoodEvent;

}

bool TriggeredMuons::findProbes()
{

  if (effCompType == ONE_MUON_NO_TRIG)
    {
      my_probe_muons = my_tight_muons;
    }
  else
    {
      triggeredMuonsIt probeIt  = my_tight_muons.begin();
      triggeredMuonsIt probeEnd = my_tight_muons.end();

      for (;probeIt!=probeEnd;++probeIt)
	{
	  triggeredMuonsIt tagIt  = my_tight_muons.begin();
	  triggeredMuonsIt tagEnd = my_tight_muons.end();
	  
	  for (;tagIt!=tagEnd;++tagIt)
	    {
	      
	      bool isGoodTag = effCompType == COUNT_TNP ? tagIt->my_mu->hlt_isomu.at(tagIt->my_imu) : true;
	      
	      if (tagIt!=probeIt && isGoodTag)
		{
		  my_probe_muons.push_back(*probeIt);
		}
	      
	    }
	}
    }

  return my_probe_muons.size() > 0;

}

bool TriggeredMuons::runTriggerMatching(ControlPlotter* plots)
{

  triggeredMuonsIt probeIt  = my_probe_muons.begin();
  triggeredMuonsIt probeEnd = my_probe_muons.end();

  for (;probeIt!=probeEnd;++probeIt)
    {

      TriggeredMuon bestCand(probeIt->my_mu,probeIt->my_imu,0,0);

      //       if (abs(my_gmt->CandBx.at(iDt))>1) CB what to do with BX
      // 	continue;
      
      int nMatchedCands = 0;

      for (int iGmt=0; iGmt<my_gmt->N; ++iGmt)
	{
	  TriggeredMuon cand(probeIt->my_mu,probeIt->my_imu,my_gmt,iGmt);
	  
	  plots->fillTrigger(cand);
	  
	  if (cand.hasTriggerMatch())
	    {
	      if (bestCand.hasGmt()) 
		{
		  int qualCand     = cand.my_gmt->Qual.at(cand.my_igmt);
		  int qualBestCand = bestCand.my_gmt->Qual.at(bestCand.my_igmt);
		  
		  if (qualCand > qualBestCand)
		    bestCand = cand;
		} 
	      else
		{
		  bestCand = cand;
		}

	      nMatchedCands++;
	      
	    }
	}

      plots->fillTrigger(nMatchedCands);

      my_triggered_muons.push_back(bestCand);

    }

  return true;

}
