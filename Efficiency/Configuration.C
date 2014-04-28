// ******************************
// Generic definitions
// ******************************

// dPhi / dEta rectangle definition 
// for GMT - RECO matching
#define MAX_MU_GMT_DPHI .2
#define MAX_MU_GMT_DETA .5

// dPhi / dEta rectangle definition
// to veto close by RECO muons 
#define MAX_MU_MU_DPHI .4
#define MAX_MU_MU_DETA 1.

#define MAX_MU_ETA 2.4

// T&P selection logic
enum EffCompType { COUNT_TNP, TWO_MUON_NO_TRIG, ONE_MUON_NO_TRIG };

// EffCompType effCompType = ONE_MUON_NO_TRIG; // one muon samples (e.g. Ws MC)
EffCompType effCompType = TWO_MUON_NO_TRIG; // two muon samples (e.g. DY MC)
// EffCompType effCompType = COUNT_TNP; // two muon samples with trigger bias (e.g. Z skims from Sinfle Mu dataset on data) 

// GMT quality cuts
int gmtQualityMask[7] = {0, 0, 0, 1, 1, 1, 1}; // single mu qualities
//int gmtQualityMask[7] = {0, 0, 1, 0, 1, 1, 1}; // double mu qualities
