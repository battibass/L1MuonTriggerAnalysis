#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TObject.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TList.h"
#include "THStack.h"
#include "TIterator.h"
#include "TObject.h"
#include "TClass.h"
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>

#include "../CommonUtils/tdrstyle.C"

void getPlotList( const std::string & fileName,
		  std::map<std::string, std::vector<std::string> > & plotList )
{
  TH1::AddDirectory(kFALSE);

  // build list of histogram names
  TFile rootfile( fileName.c_str() );

  std::vector< std::string > dirList;
	  
  TList *dirs = rootfile.GetListOfKeys();
  TIterator *itdir = dirs->MakeIterator();
  TObject *nextdir;

  while ( (nextdir = itdir->Next()) ) {
    
    if( nextdir->IsFolder())
      dirList.push_back( nextdir->GetName() );
    else 
      plotList[""].push_back( ( nextdir->GetName() ) );

    
  }
  
  std::vector<std::string>::const_iterator dirIt  = dirList.begin();
  std::vector<std::string>::const_iterator dirEnd = dirList.end();
  
  for (;dirIt!=dirEnd;++dirIt){
    
    TDirectory * thisdir = (TDirectory*)( rootfile.Get( dirIt->c_str() ) );
    
    TList * dircontent = thisdir->GetListOfKeys();
    TIterator * thisplot = dircontent->MakeIterator();
    TObject * nextplot;
 
    const std::string & dirName = (*dirIt); 
    
    while ( (nextplot = thisplot->Next()) ) {
      plotList[dirName].push_back( (  dirName + "/" + nextplot->GetName() ) );
    }

  }

  rootfile.Close();

}


void plot( std::vector<TEfficiency*> plots,
	   std::string &baseDir, std::string outputDir ) 
{
  
  if (plots.at(0) && dynamic_cast<TEfficiency*>(plots.at(0)) && plots.at(0)->GetDimension() == 1)
    {
      std::cout << "Plotting : " << plots.at(0)->GetName() << std::endl;
  
      // plot everything
      TCanvas *c = new TCanvas();

      c->cd();

      c->Draw();
      c->SetGrid();

      for (size_t iPlot=0; iPlot<plots.size(); ++iPlot) 
	{

	  c->cd();
	  
	  plots.at(iPlot)->SetLineColor( iPlot+1 );
	  plots.at(iPlot)->SetFillColor( iPlot+1 );
	  //plots.at(iPlot)->SetMarkerColor( iPlot+1 );
	  //plots.at(iPlot)->SetMarkerStyle( 21 + iPlot );
	  
	  plots.at(iPlot)->Draw( iPlot ? "same" : "" );
	  
	  c->Update();

	  std::string plotName = plots.at(iPlot)->GetName();
	  float minY = plotName.find("EffVsPt") == std::string::npos ? 0.5 : 0.;
	  float maxY = 1.1;

	  std::cout << plotName << " " <<  minY << std::endl;
	  
	  plots.at(iPlot)->GetPaintedGraph()->GetYaxis()->SetRangeUser( minY, maxY );
	  
	}

      std::string path = baseDir + "/" + outputDir;
      system( (std::string("mkdir -p ") + path).c_str() );
      
      c->Update();
      std::string printname = path + "/" + plots.at(0)->GetName();
      c->Print ( ( printname + ".gif" ).c_str() ); 
      c->Print ( ( printname + ".C" ).c_str() ); 
    }
  
}

void plotAll(std::vector<std::string> &files,
	     std::string &baseDir) 
{

  
  size_t nFiles = 0;
  std::vector<TFile*> filesRoot;

  for (size_t iFile=0; iFile<files.size(); ++iFile) {
    filesRoot.push_back(new TFile(files.at(iFile).c_str(),"READONLY"));
  }

  system( (std::string("mkdir -p ") + baseDir).c_str() );
  
  std::map<std::string, std::vector<std::string> > plotNames;

  getPlotList(files.at(0),plotNames);

  std::map<std::string, std::vector<std::string> >::const_iterator plotDirIt  = plotNames.begin();
  std::map<std::string, std::vector<std::string> >::const_iterator plotDirEnd = plotNames.end();

  for(;plotDirIt!=plotDirEnd;++plotDirIt) {

    std::vector<std::string>::const_iterator plotIt  = plotDirIt->second.begin();
    std::vector<std::string>::const_iterator plotEnd = plotDirIt->second.end();

    for(;plotIt!=plotEnd;++plotIt) {

      std::vector<TEfficiency*> plots;

      for (size_t iFile=0; iFile<filesRoot.size(); ++iFile) {
	TObject * effPlot = filesRoot.at(iFile)->Get( plotIt->c_str() );
	if (effPlot->InheritsFrom("TEfficiency")){
	  plots.push_back(static_cast<TEfficiency *>(effPlot));
	}
      }
      
      if (plots.size() > 0)
	plot(plots,baseDir,plotDirIt->first);
      
    }
  }
  
}

int main(int argc, char* argv[]) 
{  

  setTDRStyle();
  
  if ( argc<3 ) {
    std::cout << "Error in number of arguments: " << argc << std::endl;
    std::cout << "Passed args: " << argc << std::endl;
    for ( int i = 1; i < argc; ++i ) {
      std::cout << "\t" << argv[i] << std::endl;
    }
    std::cout << "Usage: \n\t\t " <<  argv[0] << " <first inputfile> <second inputfile> ... "
	      << std::endl << std::endl;
    return -1;
  }

  
  std::vector<std::string> files;
  for (int iArg=1; iArg<argc; ++iArg) {
    files.push_back(argv[iArg]);
  }

  std::string baseDir = "results/comparePlots";

  plotAll(files,baseDir);

}
