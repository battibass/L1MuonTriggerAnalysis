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

#define MAX_MU_GMT_DPHI .2
#define MAX_MU_GMT_DETA .5

#define MAX_DTTF_GMT_DPHI .1
#define MAX_DTTF_GMT_DETA .5

#define MAX_MU_ETA 1.05


int colorMap[7] = {kYellow,kOrange+7,kRed,kBlue+1,
		   kAzure+1,kCyan+1,kGreen+2}; // CB the color map for DTTF qualities

int colorMapFine[11] = {kYellow, 
			kOrange+7,kOrange+1,
			kRed+1,kRed-4,kRed-7,
			kBlue+1, 
			kAzure+1, 
			kCyan+1,kCyan,
			kGreen+2}; // CB the color map for chamber based qualities


int chambToIndex(int qual, int outerCh)
{

  int index = 0;

  switch (qual) 
    {

    case 1:
      index = outerCh == 4 ? 1 : 0;
      break;

    case 2:
      index = outerCh == 4 ? 2 : outerCh == 3 ? 3 : 0;
      break;

    case 3:
      index = outerCh == 4 ? 4 : outerCh == 3 ? 5 : outerCh == 2 ? 6 :  0;
      break;

    case 4:
      index = outerCh == 4 ? 7 : 0;
      break;

    case 5:
      index = outerCh == 4 ? 8 : 0;
      break;

    case 6:
      index = outerCh == 4 ? 9 : outerCh == 3 ? 10 : 0;
      break;

    case 7:
      index = outerCh == 4 ? 11 : 0 ;
      break;
    }

  if (index == 0)
    {
      std::cout << "index is 0! qual : " << qual 
		<< " outerCh : " <<outerCh << std::endl;
    }

  return index;

}

typedef std::map<std::string,TH1*> hTH1Map;
typedef hTH1Map::const_iterator hTH1MapIt; 

typedef std::map<std::string,TEfficiency*> histoMap;
typedef histoMap::const_iterator histoMapIt; 

class TriggeredMuon;

typedef std::vector<TriggeredMuon> triggeredMuons;
typedef triggeredMuons::const_iterator triggeredMuonsIt; 

class GmtDttfMuon;

typedef std::vector<GmtDttfMuon> gmtDttfMuons;
typedef gmtDttfMuons::iterator gmtDttfMuonsIt; 

