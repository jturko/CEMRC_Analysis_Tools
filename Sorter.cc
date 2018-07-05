
#include "Sorter.hh"
#include "Utilities.hh"

#include "TMath.h"

Sorter::Sorter() {}

Sorter::Sorter(std::string filename) :
    fFileName(filename),
    fRawFile(NULL),
    fTree(NULL),
    fVerbose(false),
    fEdep_vector(NULL),
    fTime_vector(NULL),
    fDetID_vector(NULL)
{
    // Initialize some variables
    fwhm2sigma = 1.0 / ( 2.0*TMath::Sqrt(2.0*TMath::Log(2.)) );
    fTimeWindow = 300; // in nanoseconds
    //fRandom = TRandom3(0);

    // LoLung binning
    fLoLung_nbins = 4097;
    fLoLung_binw   = 0.06109179;
    fLoLung_lowbin = 0.000 + fLoLung_binw/2.;
    fLoLung_highbin = fLoLung_lowbin + fLoLung_nbins*fLoLung_binw;

    // HiLung binning
    fHiLung_nbins = 2000;
    fHiLung_lowbin = 0.;
    fHiLung_highbin = 2000.;

    // Whole-body binning
    fWB_nbins = 2000;
    fWB_lowbin = 0.;
    fWB_highbin = 2000.;

    // Check if the raw data file provided exists
    if(Utilities::Exists(filename.c_str())) {
        if(fVerbose) std::cout << " ---> Opening file: " << filename;
        fRawFile = TFile::Open(filename.c_str());
        fTree = (TTree*)fRawFile->Get("tree");    
        if(fVerbose) std::cout << "; Found tree w/ " << fTree->GetEntries() << " entries" << std::endl;
        fTree->SetBranchAddress("edep",&fEdep);
        fTree->SetBranchAddress("detectorID",&fDetID);
        fTree->SetBranchAddress("energy_vector",&fEdep_vector);
        fTree->SetBranchAddress("detectorID_vector",&fDetID_vector);
        fTree->SetBranchAddress("time_vector",&fTime_vector);    
        
        //fTree->SetDirectory(0);
        //fRawFile->Close();
    }
    else {
        std::cout << " ---> File doesn't exists: " << filename << std::endl;
    }
 
    // Initialize the histograms and sums 
    InitializeHistograms();
    ClearSums();    
 
    // Fit the resolution data
    FitResolution();
}

Sorter::~Sorter() {}

bool Sorter::Sort(double fraction)
{
    if( (fLoLung_sum->GetEntries() > 0) || (fWB_sum->GetEntries() > 0) ) {
        DeleteHistograms();
        InitializeHistograms();
    }

    int nentries = int(fraction*(double)fTree->GetEntries());
    int ninteractions;    
    double dt;

    // Loop over all entries in the tree and bin them in the histograms created above
    for(int i=0; i<nentries; i++) {
        if(i%1000==0 && fVerbose) std::cout << " " << double(i)/double(nentries)*100. << " % complete   \r";
        
        // Clear all the edep counters
        ClearSums();

        // Get entry (event)
        fTree->GetEntry(i);
        // Get the number of interactions for the entry (event)
        ninteractions = (int)fEdep_vector->size();
     
        // Loop over all interactions
        for(int j=0; j<ninteractions; ++j) {
            // Calculate time difference between interactions
            if(j==0) dt = 0;
            else     dt = TMath::Abs( fTime_vector->at(j)-fTime_vector->at(j-1) );

            if(dt > fTimeWindow) {
                FillAllHistograms();
                ClearSums();
            }

            if(fDetID_vector->at(j) >= 1 && fDetID_vector->at(j) <= 4) {
                fEdep_Lung_sum[fDetID_vector->at(j)-1] += fEdep_vector->at(j);
            }
            if(fDetID_vector->at(j) >= 9 && fDetID_vector->at(j) <= 12) {
                fEdep_WB_sum[fDetID_vector->at(j)-9] += fEdep_vector->at(j);
            }
            if(fDetID_vector->at(j) == 20) fEdep_BGO_sum += fEdep_vector->at(j);
        } // interaction loop
        FillAllHistograms();       
        ClearSums();

    } // event loop
    if(fVerbose) std::cout << " 100% complete!       " << std::endl;
    
    return true;
}

void Sorter::FitResolution()
{
    TF1 * tmpfunc = NULL;
    //if(fVerbose) std::cout << " ---> Fitting the detector resolution data w/:  a1 + a2*sqrt(E)" << std::endl;
    for(int i=0; i<12; i++) {
        // access the information stored in the resolution_data directory
        // for the ith detector
        fResGraph[i] = new TGraph(Form("~/GEANT4_sims/CEMRC_ID/resolution_data/det%d.dat",i+1));
        
        // create a new TF1 for the ith detector, with the same range as the
        // detector (same low and high bin positions), stored as a temporary
        // function
        if(i<=3)             tmpfunc = new TF1(Form("ResFit_LoLung_det%d",i+1),"[0]+[1]*TMath::Sqrt(x)",fLoLung_lowbin,fLoLung_highbin);
        if(i>=4 && i<=7)     tmpfunc = new TF1(Form("ResFit_HiLung_det%d",i+1),"[0]+[1]*TMath::Sqrt(x)",fHiLung_lowbin,fHiLung_highbin);
        if(i>=8 && i<=11)    tmpfunc = new TF1(Form("ResFit_WB_det%d",i+1),"[0]+[1]*TMath::Sqrt(x)",fWB_lowbin,fWB_highbin);
        
        // fit the resolution data stored in the TGraph and output fit
        // parameters
        fResGraph[i]->Fit(tmpfunc,"Q");
        //if(fVerbose) std::cout << " det" << i+1 << ": a1 = " << tmpfunc->GetParameter(0) << "\ta2 = " << tmpfunc->GetParameter(1) << std::endl;

        // save the temporary function in the res_fit pointer array
        fResFit[i] = tmpfunc;
    }
}

void Sorter::FillAllHistograms()
{
    int ndet = 4;
    double edep_res;

    // fill lolung and hilung detectors
    for(int k=0; k<ndet; k++) {
        if(fEdep_Lung_sum[k] > 0.) {
            fLoLung[k]->Fill(fEdep_Lung_sum[k]);
            fLoLung_sum->Fill(fEdep_Lung_sum[k]);
            if(fEdep_BGO_sum <= 0.) {
                fLoLung_sup[k]->Fill(fEdep_Lung_sum[k]);
                fLoLung_sum_sup->Fill(fEdep_Lung_sum[k]);
            }

            fHiLung[k]->Fill(fEdep_Lung_sum[k]);
            fHiLung_sum->Fill(fEdep_Lung_sum[k]);
            if(fEdep_BGO_sum <= 0.) {
                fHiLung_sup[k]->Fill(fEdep_Lung_sum[k]);
                fHiLung_sum_sup->Fill(fEdep_Lung_sum[k]);
            }

            edep_res = fRandom.Gaus(fEdep_Lung_sum[k],fwhm2sigma*fResFit[k]->Eval(fEdep_Lung_sum[k]));
            fLoLung_res[k]->Fill(edep_res);
            fLoLung_sum_res->Fill(edep_res);
            if(fEdep_BGO_sum <= 0.) {
                fLoLung_res_sup[k]->Fill( edep_res);
                fLoLung_sum_res_sup->Fill(edep_res);
            }
            
            edep_res = fRandom.Gaus(fEdep_Lung_sum[k],fwhm2sigma*fResFit[k+4]->Eval(fEdep_Lung_sum[k]));
            fHiLung_res[k]->Fill(edep_res);
            fHiLung_sum_res->Fill(edep_res);
            if(fEdep_BGO_sum <= 0.) {
                fHiLung_res_sup[k]->Fill( edep_res);
                fHiLung_sum_res_sup->Fill(edep_res);
            }

        }
    }
    // fill wbc detectors
    for(int k=0; k<ndet; k++) {
        if(fEdep_WB_sum[k] > 0.) {
            fWB[k]->Fill(fEdep_WB_sum[k]);
            fWB_sum->Fill(fEdep_WB_sum[k]);
            if(fEdep_BGO_sum <= 0.) {
                fWB_sup[k]->Fill(fEdep_WB_sum[k]);
                fWB_sum_sup->Fill(fEdep_WB_sum[k]);
            }
            
            edep_res = fRandom.Gaus(fEdep_WB_sum[k],fwhm2sigma*fResFit[k+8]->Eval(fEdep_WB_sum[k]));
            fWB_res[k]->Fill(edep_res);
            fWB_sum_res->Fill(edep_res);
            if(fEdep_BGO_sum <= 0.) {
                fWB_res_sup[k]->Fill( edep_res);
                fWB_sum_res_sup->Fill(edep_res);
            }
        }
    }
    ClearSums();
}

void Sorter::ClearSums()
{
    // reset detector sums to zero
    int ndet = 4;
    for(int k=0; k<ndet; k++) {
        fEdep_Lung_sum[k] = 0;
        fEdep_WB_sum[k] = 0;
    }
    fEdep_BGO_sum = 0;
}

TH1F * Sorter::GetHistogram(std::string name, int det)
{
    if(det == -1) {
        if(name=="hi") return fHiLung_sum;
        else if(name=="lo") return fLoLung_sum;
        else if(name=="wb") return fWB_sum;
        else std::cout << " ---> Unknown histogram type: " << name << std::endl;
    }
    else if( (det >= 0) && (det <= 3) ) {
        if(name=="hi") return fHiLung[det];
        else if(name=="lo") return fLoLung[det];
        else if(name=="wb") return fWB[det];
        else std::cout << " ---> Unknown histogram type: " << name << std::endl;
    }
    else std::cout << " ---> Can't return this histogram... (?!)" << std::endl;
    return NULL;
}

TH1F * Sorter::GetHistogramRes(std::string name, int det)
{
    if(det == -1) {
        if(name=="hi") return fHiLung_sum_res;
        else if(name=="lo") return fLoLung_sum_res;
        else if(name=="wb") return fWB_sum_res;
        else std::cout << " ---> Unknown histogram type: " << name << std::endl;
    }
    else if( (det >= 0) && (det <= 3) ) {
        if(name=="hi") return fHiLung_res[det];
        else if(name=="lo") return fLoLung_res[det];
        else if(name=="wb") return fWB_res[det];
        else std::cout << " ---> Unknown histogram type: " << name << std::endl;
    }
    else std::cout << " ---> Can't return this histogram... (?!)" << std::endl;
    return NULL;
}

TH1F * Sorter::GetHistogramSup(std::string name, int det)
{
    if(det == -1) {
        if(name=="hi") return fHiLung_sum_sup;
        else if(name=="lo") return fLoLung_sum_sup;
        else if(name=="wb") return fWB_sum_sup;
        else std::cout << " ---> Unknown histogram type: " << name << std::endl;
    }
    else if( (det >= 0) && (det <= 3) ) {
        if(name=="hi") return fHiLung_sup[det];
        else if(name=="lo") return fLoLung_sup[det];
        else if(name=="wb") return fWB_sup[det];
        else std::cout << " ---> Unknown histogram type: " << name << std::endl;
    }
    else std::cout << " ---> Can't return this histogram... (?!)" << std::endl;
    return NULL;
}

TH1F * Sorter::GetHistogramResSup(std::string name, int det)
{
    if(det == -1) {
        if(name=="hi") return fHiLung_sum_res_sup;
        else if(name=="lo") return fLoLung_sum_res_sup;
        else if(name=="wb") return fWB_sum_res_sup;
        else std::cout << " ---> Unknown histogram type: " << name << std::endl;
    }
    else if( (det >= 0) && (det <= 3) ) {
        if(name=="hi") return fHiLung_res_sup[det];
        else if(name=="lo") return fLoLung_res_sup[det];
        else if(name=="wb") return fWB_res_sup[det];
        else std::cout << " ---> Unknown histogram type: " << name << std::endl;
    }
    else std::cout << " ---> Can't return this histogram... (?!)" << std::endl;
    return NULL;
}

void Sorter::InitializeHistograms()
{
    std::string append = fFileName.substr(0,fFileName.size()-5);

    // create histogram objects for each individual detector
    for(int i=0; i<4; i++) {
        // the unresolved histograms
        fLoLung[i] = new TH1F(Form("LoLung_det%d_%s",i+1,append.c_str()),Form("LoLung_det%d",i+1),fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
        fHiLung[i] = new TH1F(Form("HiLung_det%d_%s",i+5,append.c_str()),Form("HiLung_det%d",i+5),fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
        fWB[i] =     new TH1F(Form("WB_det%d_%s",i+9,append.c_str())    ,Form("WB_det%d",i+9)    ,fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
        // the resolved histograms
        fLoLung_res[i] = new TH1F(Form("LoLung_det%d_res_%s",i+1,append.c_str()),Form("LoLung_det%d_res",i+1),fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
        fHiLung_res[i] = new TH1F(Form("HiLung_det%d_res_%s",i+5,append.c_str()),Form("HiLung_det%d_res",i+5),fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
        fWB_res[i] =     new TH1F(Form("WB_det%d_res_%s",i+9,append.c_str())    ,Form("WB_det%d_res",i+9)    ,fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
        
        // the unresolved histograms (compton suppressed)
        fLoLung_sup[i] = new TH1F(Form("LoLung_det%d_sup_%s",i+1,append.c_str()),Form("LoLung_det%d_sup",i+1),fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
        fHiLung_sup[i] = new TH1F(Form("HiLung_det%d_sup_%s",i+5,append.c_str()),Form("HiLung_det%d_sup",i+5),fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
        fWB_sup[i] =     new TH1F(Form("WB_det%d_sup_%s",i+9,append.c_str())    ,Form("WB_det%d_sup",i+9)    ,fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
        // the resolved histograms (compton supressed)
        fLoLung_res_sup[i] = new TH1F(Form("LoLung_det%d_res_sup_%s",i+1,append.c_str()),Form("LoLung_det%d_res_sup",i+1),fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
        fHiLung_res_sup[i] = new TH1F(Form("HiLung_det%d_res_sup_%s",i+5,append.c_str()),Form("HiLung_det%d_res_sup",i+5),fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
        fWB_res_sup[i] =     new TH1F(Form("WB_det%d_res_sup_%s",i+9,append.c_str())    ,Form("WB_det%d_res_sup",i+9)    ,fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    }
    // create histogram objects for each summed geometry (lolung, hilung, wbc)
    // - both resolved and unresolved
    fLoLung_sum =     new TH1F(Form("LoLung_sum_%s",append.c_str()),     "LoLung_sum",      fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    fHiLung_sum =     new TH1F(Form("HiLung_sum_%s",append.c_str()),     "HiLung_sum",      fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    fWB_sum     =     new TH1F(Form("WB_sum_%s",append.c_str()),         "WB_sum",          fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    fLoLung_sum_res = new TH1F(Form("LoLung_sum_res_%s",append.c_str()), "LoLung_sum_res",  fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    fHiLung_sum_res = new TH1F(Form("HiLung_sum_res_%s",append.c_str()), "HiLung_sum_res",  fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    fWB_sum_res     = new TH1F(Form("WB_sum_res_%s",append.c_str()),     "WB_sum_res",      fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    // compton supressed
    fLoLung_sum_sup =     new TH1F(Form("LoLung_sum_sup_%s",append.c_str()),     "LoLung_sum_sup",      fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    fHiLung_sum_sup =     new TH1F(Form("HiLung_sum_sup_%s",append.c_str()),     "HiLung_sum_sup",      fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    fWB_sum_sup     =     new TH1F(Form("WB_sum_sup_%s",append.c_str()),         "WB_sum_sup",          fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    fLoLung_sum_res_sup = new TH1F(Form("LoLung_sum_res_sup_%s",append.c_str()), "LoLung_sum_res_sup",  fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    fHiLung_sum_res_sup = new TH1F(Form("HiLung_sum_res_sup_%s",append.c_str()), "HiLung_sum_res_sup",  fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    fWB_sum_res_sup     = new TH1F(Form("WB_sum_res_sup_%s",append.c_str()),     "WB_sum_res_sup",      fWB_nbins    ,fWB_lowbin    ,fWB_highbin);

    //// create histogram objects for each individual detector
    //for(int i=0; i<4; i++) {
    //    // the unresolved histograms
    //    fLoLung[i] = new TH1F(Form("LoLung_det%d",i+1),Form("LoLung_det%d",i+1),fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    //    fHiLung[i] = new TH1F(Form("HiLung_det%d",i+5),Form("HiLung_det%d",i+5),fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    //    fWB[i] =     new TH1F(Form("WB_det%d",i+9)    ,Form("WB_det%d",i+9)    ,fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    //    // the resolved histograms
    //    fLoLung_res[i] = new TH1F(Form("LoLung_det%d_res",i+1),Form("LoLung_det%d_res",i+1),fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    //    fHiLung_res[i] = new TH1F(Form("HiLung_det%d_res",i+5),Form("HiLung_det%d_res",i+5),fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    //    fWB_res[i] =     new TH1F(Form("WB_det%d_res",i+9)    ,Form("WB_det%d_res",i+9)    ,fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    //    
    //    // the unresolved histograms (compton suppressed)
    //    fLoLung_sup[i] = new TH1F(Form("LoLung_det%d_sup",i+1),Form("LoLung_det%d_sup",i+1),fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    //    fHiLung_sup[i] = new TH1F(Form("HiLung_det%d_sup",i+5),Form("HiLung_det%d_sup",i+5),fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    //    fWB_sup[i] =     new TH1F(Form("WB_det%d_sup",i+9)    ,Form("WB_det%d_sup",i+9)    ,fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    //    // the resolved histograms (compton supressed)
    //    fLoLung_res_sup[i] = new TH1F(Form("LoLung_det%d_res_sup",i+1),Form("LoLung_det%d_res_sup",i+1),fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    //    fHiLung_res_sup[i] = new TH1F(Form("HiLung_det%d_res_sup",i+5),Form("HiLung_det%d_res_sup",i+5),fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    //    fWB_res_sup[i] =     new TH1F(Form("WB_det%d_res_sup",i+9)    ,Form("WB_det%d_res_sup",i+9)    ,fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    //}
    //// create histogram objects for each summed geometry (lolung, hilung, wbc)
    //// - both resolved and unresolved
    //fLoLung_sum =     new TH1F("LoLung_sum",     "LoLung_sum",      fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    //fHiLung_sum =     new TH1F("HiLung_sum",     "HiLung_sum",      fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    //fWB_sum     =     new TH1F("WB_sum",         "WB_sum",          fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    //fLoLung_sum_res = new TH1F("LoLung_sum_res", "LoLung_sum_res",  fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    //fHiLung_sum_res = new TH1F("HiLung_sum_res", "HiLung_sum_res",  fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    //fWB_sum_res     = new TH1F("WB_sum_res",     "WB_sum_res",      fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    //// compton supressed
    //fLoLung_sum_sup =     new TH1F("LoLung_sum_sup",     "LoLung_sum_sup",      fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    //fHiLung_sum_sup =     new TH1F("HiLung_sum_sup",     "HiLung_sum_sup",      fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    //fWB_sum_sup     =     new TH1F("WB_sum_sup",         "WB_sum_sup",          fWB_nbins    ,fWB_lowbin    ,fWB_highbin);
    //fLoLung_sum_res_sup = new TH1F("LoLung_sum_res_sup", "LoLung_sum_res_sup",  fLoLung_nbins,fLoLung_lowbin,fLoLung_highbin);
    //fHiLung_sum_res_sup = new TH1F("HiLung_sum_res_sup", "HiLung_sum_res_sup",  fHiLung_nbins,fHiLung_lowbin,fHiLung_highbin);
    //fWB_sum_res_sup     = new TH1F("WB_sum_res_sup",     "WB_sum_res_sup",      fWB_nbins    ,fWB_lowbin    ,fWB_highbin);

}

void Sorter::DeleteHistograms()
{
    for(int i=0; i<4; i++) {
        // the unresolved histograms
        delete fLoLung[i];
        delete fHiLung[i];
        delete fWB[i];
        // the resolved histograms
        delete fLoLung_res[i];
        delete fHiLung_res[i];
        delete fWB_res[i];
        
        // the unresolved histograms (compton suppressed)
        delete fLoLung_sup[i];
        delete fHiLung_sup[i];
        delete fWB_sup[i];
        // the resolved histograms (compton supressed)
        delete fLoLung_res_sup[i];
        delete fHiLung_res_sup[i];
        delete fWB_res_sup[i];
    }
    // - both resolved and unresolved
    delete fLoLung_sum;
    delete fHiLung_sum; 
    delete fWB_sum;  
    delete fLoLung_sum_res;
    delete fHiLung_sum_res;
    delete fWB_sum_res;    
    // compton supressed
    delete fLoLung_sum_sup;
    delete fHiLung_sum_sup;  
    delete fWB_sum_sup;  
    delete fLoLung_sum_res_sup;
    delete fHiLung_sum_res_sup;
    delete fWB_sum_res_sup;    
}



