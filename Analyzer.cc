
#include "Analyzer.hh"

Analyzer::Analyzer(int plate, double radon, double k, double cosmic) :
    fVerbose(true),
    f222RnActivity(radon),
    f40KActivity(k),
    fCosmicActivity(cosmic),
    f222RnSorter(NULL),
    f40KSorter(NULL),
    fHalfHour(true)
{
    // Default maximum file sizes
    f40KMaxActivity = 160.0; // 160 [nCi] 40K
    f222RnMaxActivity = 4.0; // 4 [pCi/L] 222Rn
    fCosmicMaxActivity = 1e-2;   // 1e-2 [sr-1 s-1 cm-2] cosmic muons
    
    double cwt;
    if     (plate == 0) cwt = 0.00;
    else if(plate == 1) cwt = 1.60;
    else if(plate == 2) cwt = 2.22;
    else if(plate == 3) cwt = 3.01;
    else if(plate == 4) cwt = 3.33;
    else if(plate == 5) cwt = 4.18;
    else if(plate == 6) cwt = 5.10;
    else if(plate == 7) cwt = 6.00;
    else                cwt = 0.00;

    f222RnSorter =     new Sorter(Form("/home/jturko/CEMRC/Contract1/Simulation_data/CEMRC_bkgd_CWT/CWT_%.2fcm/Rn222_4pCi_CEMRC_bkgd/yesBOMAB/output.root",cwt));
    f40KSorter =       new Sorter(Form("/home/jturko/CEMRC/Contract1/Simulation_data/CEMRC_bkgd_CWT/CWT_%.2fcm/K40_160nCi_CEMRC_bkgd/output.root",cwt));
    fCosmicSorter[0] = new Sorter(Form("/home/jturko/CEMRC/Contract2/simulation_data/CEMRC_cosmic_bkgd/CWT_%.2fcm/CEMRC_flux_Grieder327EnergyDist/output.root",cwt));
    fCosmicSorter[1] = new Sorter(Form("/home/jturko/CEMRC/Contract2/simulation_data/CEMRC_cosmic_bkgd/CWT_%.2fcm/CEMRC_flux_Grieder328EnergyDist/output.root",cwt));
    fCosmicSorter[2] = new Sorter(Form("/home/jturko/CEMRC/Contract2/simulation_data/CEMRC_cosmic_bkgd/CWT_%.2fcm/CEMRC_flux_MeiHimeEnergyDist/output.root",cwt));
}

Analyzer::~Analyzer() {}

void Analyzer::Sort() 
{
    double fraction;

    fraction = f222RnActivity / f222RnMaxActivity;
    if(fHalfHour) fraction /= 48.;
    if(fraction > 1.) { std::cout << " ---> The 222Rn activity provided is too large! max activity = " << f222RnMaxActivity << "pCi/L" << std::endl; return; }
    f222RnSorter->Sort(fraction);

    fraction = f40KActivity / f40KMaxActivity;
    if(fHalfHour) fraction /= 48.;
    if(fraction > 1.) { std::cout << " ---> The 40K activity provided is too large! max activity = " << f40KMaxActivity << "nCi" << std::endl; return; }
    f40KSorter->Sort(fraction);

    fraction = fCosmicActivity / fCosmicMaxActivity;
    if(fraction > 1.) { std::cout << " ---> The cosmic activity provided is too large! max activity = " << fCosmicMaxActivity << "/(sr s cm^2)" << std::endl; return; }
    fCosmicSorter[0]->Sort(fraction);
    fCosmicSorter[1]->Sort(fraction);
    fCosmicSorter[2]->Sort(fraction);

}

void Analyzer::SetActivities(double radon, double k, double cosmic)
{
    f222RnActivity = radon;
    f40KActivity = k;
    fCosmicActivity = cosmic;
}

void Analyzer::SetMaxActivities(double radon, double k, double cosmic)
{
    f222RnMaxActivity = radon;
    f40KMaxActivity = k;
    fCosmicMaxActivity = cosmic;
}

void Analyzer::PlotAll(std::string det)
{
    new TCanvas();
    f222RnSorter->GetHistogramRes(det)->SetTitle(Form("^{222}Rn: %.1fpCi/L",f222RnActivity));
    f222RnSorter->GetHistogramRes(det)->GetXaxis()->SetTitle("Energy [keV]");
    f222RnSorter->GetHistogramRes(det)->GetYaxis()->SetTitle("Counts");
    f222RnSorter->GetHistogramRes(det)->Draw();
    
    new TCanvas();
    f40KSorter->GetHistogramRes(det)->SetTitle(Form("^{222}Rn: %.1fnCi",f40KActivity));
    f40KSorter->GetHistogramRes(det)->GetXaxis()->SetTitle("Energy [keV]");
    f40KSorter->GetHistogramRes(det)->GetYaxis()->SetTitle("Counts");
    f40KSorter->GetHistogramRes(det)->Draw();
    
    new TCanvas();
    fCosmicSorter[0]->GetHistogramRes(det)->SetTitle(Form("Grieder 3.27: %.2e/(s sr cm2)",fCosmicActivity));
    fCosmicSorter[0]->GetHistogramRes(det)->GetXaxis()->SetTitle("Energy [keV]");
    fCosmicSorter[0]->GetHistogramRes(det)->GetYaxis()->SetTitle("Counts");
    fCosmicSorter[0]->GetHistogramRes(det)->Draw();
    
    new TCanvas();
    fCosmicSorter[1]->GetHistogramRes(det)->SetTitle(Form("Grieder 3.28: %.2e/(s sr cm2)",fCosmicActivity));
    fCosmicSorter[1]->GetHistogramRes(det)->GetXaxis()->SetTitle("Energy [keV]");
    fCosmicSorter[1]->GetHistogramRes(det)->GetYaxis()->SetTitle("Counts");
    fCosmicSorter[1]->GetHistogramRes(det)->Draw();
    
    new TCanvas();
    fCosmicSorter[2]->GetHistogramRes(det)->SetTitle(Form("Mei/Hime: %.2e/(s sr cm2)",fCosmicActivity));
    fCosmicSorter[2]->GetHistogramRes(det)->GetXaxis()->SetTitle("Energy [keV]");
    fCosmicSorter[2]->GetHistogramRes(det)->GetYaxis()->SetTitle("Counts");
    fCosmicSorter[2]->GetHistogramRes(det)->Draw();
}
