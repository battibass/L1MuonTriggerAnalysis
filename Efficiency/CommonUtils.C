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

// CB the color map for GMT qualities
int colorMap[7] = {kYellow,kOrange+7,kRed,kBlue+1,
		   kAzure+1,kCyan+1,kGreen+2}; 

typedef std::pair<int,int> qualPair;
typedef std::map<qualPair,int> qualToFine;

typedef std::map<std::string,TH1*> hTH1Map;
typedef hTH1Map::const_iterator hTH1MapIt; 

typedef std::map<std::string,TEfficiency*> histoMap;
typedef histoMap::const_iterator histoMapIt; 

class TriggeredMuon;

typedef std::vector<TriggeredMuon> triggeredMuons;
typedef triggeredMuons::const_iterator triggeredMuonsIt; 
