
#include "Analyzer.hh"

Analyzer::Analyzer(double radon, double k, double cosmic, int plate) :
    f222RnActivity(radon),
    f40KActivity(k),
    fCosmicActivity(cosmic)
{
    // Default maximum file sizes
    f40KMaxActivity = 160.0; // 160 [nCi] 40K
    f222RnMaxActivity = 4.0; // 4 [pCi/L] 222Rn
    fCosmicMaxActivity = 1e-2;   // 1e-2 [sr-1 s-1 cm-2] cosmic muons
    
    double cwt;
    switch(plate) {
        case 0:
            cwt = 0.00;
        case 1:
            cwt = 1.60;
        case 2:
            cwt = 2.22;
        case 3:
            cwt = 3.01;
        case 4: 
            cwt = 3.33;
        case 5:
            cwt = 4.18;
        case 6:
            cwt = 5.10;
        case 7:
            cwt = 6.00;
        default:
            cwt = 0.00;
    }
    f222RnSorter = new Sorter(Form("/home/jturko/CEMRC/Contract1/Simulation_data/CEMRC_bkgd_CWT/CWT_%.2fcm/Rn222_4pCi_CEMRC_bkgd/yesBOMAB/output.root",cwt));
    f40KSorter = new Sorter(  Form("/home/jturko/CEMRC/Contract1/Simulation_data/CEMRC_bkgd_CWT/CWT_%.2fcm/K40_160nCi_CEMRC_bkgd/output.root",cwt));
}

Analyzer::~Analyzer() {}

void Analyzer::Sort() 
{
    double fraction;

    fraction = f222RnActivity / f222RnMaxActivity;
    if(fraction > 1.) { std::cout << " ---> The 222Rn activity provided is too large!" << std::endl; return; }
    f222RnSorter->Sort(fraction);

    fraction = f40KActivity / f40KMaxActivity;
    if(fraction > 1.) { std::cout << " ---> The 40K activity provided is too large!" << std::endl; return; }
    f40KSorter->Sort(fraction);
}

