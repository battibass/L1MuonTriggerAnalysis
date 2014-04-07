#include <map>
#include <vector>

#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"

#include "TMath.h"
#include "THStack.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TColor.h"

#include "../CommonUtils/tdrstyle.C"

// ******************************
// Generic definitions
// ******************************

#define MAX_MU_GMT_DPHI .2
#define MAX_MU_GMT_DETA .5

#define MAX_MU_MU_DPHI .4
#define MAX_MU_MU_DETA 1.

#define MAX_MU_ETA 2.4

int colorMap[7] = {kYellow,kOrange+7,kRed,kBlue+1,
		   kAzure+1,kCyan+1,kGreen+2}; // CB the color map for DTTF qualities

int colorMapFine[11] = {kYellow, kOrange+7,kOrange+1,
			kRed+1,kRed-4,kRed-7,
			kBlue+1, kAzure+1, kCyan+1,kCyan,
			kGreen+2}; // CB the color map for chamber based qualities

typedef std::pair<int,int> qualPair;
typedef std::map<qualPair,int> qualToFine;

typedef std::map<std::string,TH1*> hTH1Map;
typedef hTH1Map::const_iterator hTH1MapIt; 

typedef std::map<std::string,TEfficiency*> histoMap;
typedef histoMap::const_iterator histoMapIt; 

class TriggeredMuon;

typedef std::vector<TriggeredMuon> triggeredMuons;
typedef triggeredMuons::const_iterator triggeredMuonsIt; 