#include <iostream>
#include <fstream>
#include <vector>
#include <TStyle.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TVirtualFFT.h>
#include <string>
#include "runner.h"
#include <TFile.h>

using namespace std;

void dupli3(TTree* tree) {

	vector<double>* w0 = new vector<double>;
	vector<double>* w1 = new vector<double>;
	vector<double>* w2 = new vector<double>;

	vector<double> wave0;
	vector<double> wave1;
	vector<double> wave2; 
	vector<double> times;

	unsigned int initTime = 0;
	unsigned int nSamples = 0;
	double dt = 0.0;

	int dim = 0;
	int nEntries = tree->GetEntries();

	string filename;
	string subfilename;
	int end;

	tree->SetBranchAddress("Waveform000", &w0);
	tree->SetBranchAddress("Waveform001", &w1);
	tree->SetBranchAddress("Waveform002", &w2);
	tree->SetBranchAddress("NumberOfSamples", &nSamples);
	tree->SetBranchAddress("SamplingWidth_s", &dt);

	tree->GetEntry(0);
	tree->SetBranchAddress("Timestamp_s", &initTime);

	for (int j = 0; j < nEntries; j++) {

		tree->GetEntry(j);

		wave0.insert(wave0.end(), w0->begin(), w0->end());
		wave1.insert(wave1.end(), w1->begin(), w1->end());
		wave2.insert(wave2.end(), w2->begin(), w2->end());

		for (int i = 0; i < nSamples; i++) {
			times.push_back(j + i * dt);
		}

	}

	if (wave0.size() == wave1.size() && wave1.size() == wave2.size()) {
		dim = wave1.size();
	}
	else {
		cout << "Problemo!" << endl;
	}

	filename = tree->GetCurrentFile()->GetName();
	end = filename.find(".");
	subfilename = filename.substr(0,end) + ".csv";
	ofstream myfile;
	myfile.open(subfilename.c_str());
	for (int i = 0; i < dim; i++) {
		myfile << times[i] << "," << wave0[i] << "," << wave1[i] << "," << wave2[i] << endl;
	}
	myfile.close();
}


void root2csv(const char *dirname="/home/Kenny/Analysis/pythonAnalysis/", const char *ext=".root") {

	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();

	if (files) {
		TSystemFile *file; TString fname; TIter next(files);

		while ((file=(TSystemFile*)next())) {
			fname = file->GetName();

			if (!file->IsDirectory() && fname.EndsWith(ext)) {
				TFile *hfile = new TFile(fname);
				TTree* tree = nullptr;
				hfile->GetObject("data_tree",tree);				
				dupli3(tree);

				cout << fname.Data() <<  endl;
			}

		}
	}
}
