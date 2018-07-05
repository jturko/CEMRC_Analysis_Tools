
#ifndef ANALYZER
#define ANALYZER 1

#include "Sorter.cc"

#include "TGraph2D.h"

class Analyzer
{
  public:
    // Constructors and destructors
    Analyzer(int plate, double radon=0., double k=0., double cosmic=0.);
    ~Analyzer();

    // Setters and getters
    void Set222RnActivity(double radon) { f222RnActivity = radon; }
    double Get222RnActivity() { return f222RnActivity; }
    void Set222RnMaxActivity(double radon) { f222RnMaxActivity = radon; }
    double Get222RnMaxActivity() { return f222RnMaxActivity; }
    
    void Set40KActivity(double k) { f40KActivity = k; }
    double Get40KActivity() { return f40KActivity; }
    void Set40KMaxActivity(double k) { f40KMaxActivity = k; }
    double Get40KMaxActivity() { return f40KMaxActivity; }

    void SetCosmicActivity(double cosmic) { fCosmicActivity = cosmic; }
    double GetCosmicActivity() { return fCosmicActivity; }
    void SetCosmicMaxActivity(double cosmic) { fCosmicMaxActivity = cosmic; }
    double GetCosmicMaxActivity() { return fCosmicMaxActivity; }

    void SetActivities(double radon, double k, double cosmic);
    void SetMaxActivities(double radon, double k, double cosmic);

    Sorter * Get222RnSorter() { return f222RnSorter; }
    Sorter * Get40KSorter() { return f40KSorter; }
    Sorter * GetCosmicSorter(int i) { return fCosmicSorter[i]; }

    void Sort();
    void PlotAll(std::string det="hi");
    void GenerateBkgd();
    void DeleteHistograms();
    TH1F * GetLoLungHist(int i) { return fLoLung_bkgd[i]; }
    TH1F * GetHiLungHist(int i) { return fHiLung_bkgd[i]; }
    TH1F * GetWBHist(int i) { return fWB_bkgd[i]; }

    // MDA calculating functions
    void SetMDAEfficiency(double eff) { fMDAEfficiency = eff; }
    void SetMDAIGamma(double igam) { fMDAIGamma = igam; }
    void SetMDATime(double time) { fMDATime = time; }
    void SetMDALowBin(int lowbin) { fMDALowBin = lowbin; }
    void SetMDAHighBin(int highbin) { fMDAHighBin = highbin; }
    const double * CalculateMDA(int lowbin, int highbin, double eff, double I_gamma, double time);
    const double * CalculateMDA();
    void GenerateMDATable(double cosmic=7e-3, double radon_low=0., double radon_high=4., double radon_inc=0.5, double k_low=0., double k_high=160., double k_inc=20.);
    TGraph2D * GetMDAGraph2D(int i) const { return fMDAGraph2D[i]; }

  private:
    bool fVerbose;
    bool fHalfHour;    

    double f222RnActivity;
    double f222RnMaxActivity;

    double f40KActivity;
    double f40KMaxActivity;

    double fCosmicActivity;
    double fCosmicMaxActivity;

    double fCWT; 
   
    double fMDA[3];
    int fMDALowBin;
    int fMDAHighBin;
    double fMDAEfficiency;
    double fMDAIGamma;
    double fMDATime;

    TGraph2D * fMDAGraph2D[3];

    Sorter * f222RnSorter;
    Sorter * f40KSorter;
    Sorter * fCosmicSorter[3];
    
    TH1F * fLoLung_bkgd[3];
    TH1F * fHiLung_bkgd[3];
    TH1F * fWB_bkgd[3];
};

#endif
