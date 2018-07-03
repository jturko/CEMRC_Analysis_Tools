
#include "Sorter.cc"

class Analyzer
{
  public:
    // Constructors and destructors
    Analyzer(double radon, double k, double cosmic, int plate);
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

    Sorter * Get222RnSorter() { return f222RnSorter; }
    Sorter * Get40KSorter() { return f40KSorter; }
    Sorter * GetCosmicSorter() { return fCosmicSorter; }

    void Sort();

  private:
    double f222RnActivity;
    double f222RnMaxActivity;

    double f40KActivity;
    double f40KMaxActivity;

    double fCosmicActivity;
    double fCosmicMaxActivity;
    
    Sorter * f222RnSorter;
    Sorter * f40KSorter;
    Sorter * fCosmicSorter;

};
