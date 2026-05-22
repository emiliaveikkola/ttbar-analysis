#include <TChain.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

void SumGenWeight(const char* fileList="xsec_input_files/mcFiles_TTtoLNu2Q_MCSummer24_skims.txt") {
    TChain runs("Runs");

    std::ifstream fin(fileList);
    std::string filename;
    int nfiles = 0;

    while (fin >> filename) {
        runs.Add(filename.c_str());
        nfiles++;
    }

    double genEventSumw = 0.0;
    runs.SetBranchAddress("genEventSumw", &genEventSumw);

    double sumGenWeight = 0.0;
    Long64_t nRuns = runs.GetEntries();

    for (Long64_t i = 0; i < nRuns; i++) {
        runs.GetEntry(i);
        sumGenWeight += genEventSumw;
    }

    std::cout << std::setprecision(17);
    std::cout << "Files          = " << nfiles << std::endl;
    std::cout << "Runs entries   = " << nRuns << std::endl;
    std::cout << "sumGenWeight   = " << sumGenWeight << std::endl;
}
