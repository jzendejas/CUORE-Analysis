
#include <iostream>
#include <vector>
#include <TStyle.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TVirtualFFT.h>



using namespace std;

void plotTimestream(TTree *tree){

vector<double> *w0 = new vector<double>;
vector<double> *w1 = new vector<double>;
vector<double> *w2 = new vector<double>;

vector<double> wave0;
vector<double> wave1;
vector<double> wave2; // = new vector<double>;
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

	//cout << j << endl;
	//cout << "The vectors are now length: " << wave0.size() << endl;
	
	
	for(int i = 0; i < nSamples; i++){
		times.push_back(j + i*dt);
	}

}
cout << times.size() << endl;
cout << wave0.size() << endl;

if (wave0.size() == wave1.size() && wave1.size() == wave2.size()){ 
dim = wave1.size();
}
else{ 
cout << "Problemo!" << endl;
}
 
cout << wave0[0] << endl;
TCanvas *c1 = new TCanvas("c1", "Timestream"); 
TGraph* gr = new TGraph(dim, &times[0], &wave0[0]);
gr->Draw();
}
