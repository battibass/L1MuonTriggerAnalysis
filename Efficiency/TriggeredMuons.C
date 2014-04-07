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
	  my_tight_muons.push_back(TriggeredMuon(my_mu,iMu,0));
	}
    }

  int nTightMu = my_tight_muons.size();  

  if (nTightMu == 2)
    {
      TriggeredMuon & mu1 = my_tight_muons.at(0);
      TriggeredMuon & mu2 = my_tight_muons.at(1);

      plots->fillTight(mu1);
      plots->fillTight(mu2);

      float fDPhi = fabs(acos(cos( mu1.my_mu->sa_phi_mb2[mu1.imy_mu] - 
			           mu2.my_mu->sa_phi_mb2[mu2.imy_mu] ) ) );

      float fDEta = fabs(acos(cos( mu1.my_mu->eta[mu1.imy_mu] - 
			           mu2.my_mu->eta[mu2.imy_mu] ) ) );

      return fDPhi > MAX_MU_MU_DPHI &&
	     fDEta > MAX_MU_MU_DETA ;
    }

  return false;

}

bool TriggeredMuons::findProbes()
{

  triggeredMuonsIt probeIt  = my_tight_muons.begin();
  triggeredMuonsIt probeEnd = my_tight_muons.end();

  for (;probeIt!=probeEnd;++probeIt)
    {
      triggeredMuonsIt tagIt  = my_tight_muons.begin();
      triggeredMuonsIt tagEnd = my_tight_muons.end();

      for (;tagIt!=tagEnd;++tagIt)
	{
	  
	  if (tagIt!=probeIt && tagIt->my_mu->hlt_isomu.at(tagIt->imy_mu))
	    {
	      my_probe_muons.push_back(*probeIt);
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

      TriggeredMuon bestCand(probeIt->my_mu,probeIt->imy_mu,0,0);

      //       if (abs(my_gmt->CandBx.at(iDt))>1) CB what to do with BX
      // 	continue;
      
      for (int iGmt=0; iGmt<my_gmt->Ncand; ++iGmt)
	{
	  TriggeredMuon cand(probeIt->my_mu,probeIt->my_imu,my_gmt,iGmt);
	  
	  plots->fillTrigger(cand);
	  
	  if (cand.hasTriggerMatch())
	    {
	      if (bestCand.hasGmt()) 
		{
		  int ptCand     = cand.gmt()->Pt.at(cand.igmt());
		  int ptBestCand = bestCand.gmt()->Pt.at(bestCand.igmt());
		  
		  if (ptCand > ptBestCand)
		    bestCand = cand;
		} 
	      else
		{
		  bestCand = cand;
		}  
	    }
	}

      my_triggered_muons.push_back(bestCand);

    }

  return true;

}