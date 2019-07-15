#include <iostream>
#include <vector>
#include <TStyle.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TVirtualFFT.h>
#include <string>
using namespace std;

void plotTimeXYZ(TString path, TString prefix, int nFiles, double start, double end, TString title){

/*	TString titlepiece;
	cout << "Enter Title Suffix: " << endl;
	getline(cin, titlepiece);
*/

	TChain *tree = new TChain("data_tree");
	for(int filenum = 1; filenum <= nFiles; filenum++)
		tree->Add(Form("%s/%s_p%05d.root",path.Data(),prefix.Data(),filenum));	

	vector<double> *w0 = new vector<double>;
	vector<double> *w1 = new vector<double>;
	vector<double> *w2 = new vector<double>;

	vector<double> wave0;
	vector<double> wave1;
	vector<double> wave2;
	vector<double> times;

	unsigned int initTime = 0;
	unsigned int nSamples = 0;
	double dt = 0.0;

	int dim = 0;
	int nEntries = tree->GetEntries();

	tree->SetBranchAddress("Waveform000",&w0);
	tree->SetBranchAddress("Waveform001",&w1);
	tree->SetBranchAddress("Waveform002",&w2);
	tree->SetBranchAddress("NumberOfSamples", &nSamples);
	tree->SetBranchAddress("SamplingWidth_s", &dt);

	tree->GetEntry(0);
	tree->SetBranchAddress("Timestamp_s", &initTime);

	for(int j = 0; j < nEntries; j++){

		tree->GetEntry(j);

		wave0.insert(wave0.end(), w0->begin(), w0->end());
		wave1.insert(wave1.end(), w1->begin(), w1->end());
		wave2.insert(wave2.end(), w2->begin(), w2->end());

		for(int i = 0; i < nSamples; i++){
			times.push_back(j + i*dt);
		}

	}

	if (wave0.size() == wave1.size() && wave1.size() == wave2.size()){ 
		dim = wave1.size();
	}
	else{ 
		cout << "Problemo!" << endl;
	}

	TGraph* grx = new TGraph(dim, &times[0], &wave0[0]);	
	TGraph* gry = new TGraph(dim, &times[0], &wave1[0]);
	TGraph* grz = new TGraph(dim, &times[0], &wave2[0]);

	auto c0 = new TCanvas("Timestream","Multigraph");
	auto mg = new TMultiGraph();
	mg->SetTitle(title);

	grx->SetTitle("X");
	grx->SetLineColor(kRed);
	grx->SetLineWidth(1);
	gry->SetTitle("Y");
	gry->SetLineColor(kGreen);
	gry->SetLineWidth(1);
	grz->SetTitle("Z");
	grz->SetLineColor(kBlue);
	grz->SetLineWidth(1);
	mg->Add( grx );
	mg->Add( gry );
	mg->Add( grz );
	mg->GetXaxis()->SetLimits(start,end);
	mg->SetMinimum(-0.1);
	mg->SetMaximum(0.05);
	mg->GetXaxis()->SetTitle("time (s)");
	mg->GetYaxis()->SetTitle("Accelerometer Response (V)");	
	mg->Draw("AL");

	c0->BuildLegend();


} 
