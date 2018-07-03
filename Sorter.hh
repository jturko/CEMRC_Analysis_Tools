
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TRandom3.h"

class Sorter
{
  public:
    Sorter(std::string filename);
    Sorter();
    ~Sorter();
  
    // These functions should do most of the heavy lifting...
    bool Sort(double fraction);
    void InitializeHistograms();
    void DeleteHistograms();

    // To get the histograms
    TH1F * GetHistogram(std::string name = "lo", int det = -1);    
    TH1F * GetHistogramRes(std::string name = "lo", int det = -1);    
    TH1F * GetHistogramSup(std::string name = "lo", int det = -1);    
    TH1F * GetHistogramResSup(std::string name = "lo", int det = -1);    

  private:
    void FillAllHistograms();  
    void FitResolution();
    void ClearSums();
    
    bool fVerbose;
    
    std::string fFileName;
    TFile * fRawFile;
    
    TTree * fTree;
    double fEdep;
    int fDetID;
    std::vector<double> * fEdep_vector;
    std::vector<double> * fTime_vector;
    std::vector<int> * fDetID_vector;

    TGraph * fResGraph[12];
    TF1 * fResFit[12];
    TRandom3 fRandom;   
 
    double fwhm2sigma;

    double fEdep_Lung_sum[4];
    double fEdep_WB_sum[4];
    double fEdep_BGO_sum;
    double fTimeWindow;

    // Pointers for LoLung detector histograms and binning information
    TH1F * fLoLung[4];
    TH1F * fLoLung_res[4];
    TH1F * fLoLung_sum;
    TH1F * fLoLung_sum_res;
    TH1F * fLoLung_sup[4];
    TH1F * fLoLung_res_sup[4];
    TH1F * fLoLung_sum_sup;
    TH1F * fLoLung_sum_res_sup;
    int fLoLung_nbins;
    double fLoLung_binw;
    double fLoLung_lowbin;
    double fLoLung_highbin;
    
    // Pointers for HiLung detector histograms and binning information
    TH1F * fHiLung[4];
    TH1F * fHiLung_res[4];
    TH1F * fHiLung_sum;
    TH1F * fHiLung_sum_res;
    TH1F * fHiLung_sup[4];
    TH1F * fHiLung_res_sup[4];
    TH1F * fHiLung_sum_sup;
    TH1F * fHiLung_sum_res_sup;
    int fHiLung_nbins;
    double fHiLung_lowbin;
    double fHiLung_highbin;
    
    // Pointers for WB detector histograms and binning information
    TH1F * fWB[4];
    TH1F * fWB_res[4];
    TH1F * fWB_sum;
    TH1F * fWB_sum_res;
    TH1F * fWB_sup[4];
    TH1F * fWB_res_sup[4];
    TH1F * fWB_sum_sup;
    TH1F * fWB_sum_res_sup;
    int fWB_nbins;
    double fWB_lowbin;
    double fWB_highbin;

};
