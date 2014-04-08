// ******************************
// GMT/Muon Matched object
// ( basic GMT/Muon matching object, applies and returns
// deltaPhi,deltaEta cuts with values from CommonUtils )
// ******************************

class TriggeredMuon
{

public :

  TriggeredMuon(L1Analysis::L1AnalysisRecoMuonDataFormat * mu,  int imu,
		L1Analysis::L1AnalysisGMTDataFormat      * gmt, int igmt) : 
    my_mu(mu), my_imu(imu), my_gmt(gmt), my_igmt(igmt) { };   

  TriggeredMuon(TriggeredMuon const & trigMu) 
  {
    my_mu  = trigMu.my_mu;  my_imu  = trigMu.my_imu;
    my_gmt = trigMu.my_gmt; my_igmt = trigMu.my_igmt;
  }

  ~TriggeredMuon() {};


  bool  hasTriggerMatch() const  
  { 
    return fabs(deltaPhi()) < MAX_MU_GMT_DPHI &&
           fabs(deltaEta()) < MAX_MU_GMT_DETA ; 
  };

  bool  hasHLTTriggerMatch() const 
  { 
    return  my_mu != 0 ? my_mu->hlt_isomu.at(my_imu) : false; 
  };

  float deltaPhi() const 
  { 
    if (my_gmt == 0 || my_mu == 0)
      return 999. ;

    float dPhiBarrel = acos(cos(my_gmt->Phi.at(my_igmt) - my_mu->sa_phi_mb2.at(my_imu)));
    float dPhiEndcap = acos(cos(my_mu->eta.at(my_imu) > 0 ? 
				my_gmt->Phi.at(my_igmt) - my_mu->sa_phi_me2_p.at(my_imu) :
				my_gmt->Phi.at(my_igmt) - my_mu->sa_phi_me2_n.at(my_imu) 
				));

    float dPhi = fabs(dPhiBarrel) <= fabs(dPhiEndcap) ?
                 dPhiBarrel : dPhiEndcap;

    return dPhi;
  }
  
  float deltaEta() const
  { 
    return my_gmt != 0 && my_mu != 0 ? 
           my_gmt->Eta.at(my_igmt) - my_mu->eta.at(my_imu) : 999. ;
  }

  bool hasGmt() const
  {
    return my_gmt != 0;
  }

  bool hasMuon() const
  {
    return my_mu != 0;
  }

  bool hasTightMuon() const
  {
    return my_mu != 0 && 
      (my_mu->howmanytypes.at(my_imu) >> 4) & 0x01 &&
      (my_mu->howmanytypes.at(my_imu) >> 5) & 0x01    ;
  }

public :

  L1Analysis::L1AnalysisRecoMuonDataFormat * my_mu;
  int my_imu;

  L1Analysis::L1AnalysisGMTDataFormat * my_gmt;
  int my_igmt;

};
