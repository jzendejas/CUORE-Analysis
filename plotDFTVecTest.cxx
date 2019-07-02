#include <iostream>
#include <vector>
#include <TStyle.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TVirtualFFT.h>

using namespace std;

void plotDFTVecTest(TString path, TString prefix, int numFiles){

TChain *tree = new TChain("data_tree");
for(int filenum = 1; filenum <= numFiles; filenum++)
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
 
TCanvas *c1 = new TCanvas("c1", "Timestream"); 
TGraph* gr = new TGraph(dim, &times[0], &wave0[0]);
gr->Draw();

cout << "Creating the DFT..." << endl;
Int_t n_size = dim + 1;
TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &n_size, "R2C M K");
if (!fft_own) return;
fft_own->SetPoints(&wave0[0]);
fft_own->Transform();

fft_own->GetPoints(&wave0[0]);

//Draw the graphs of the Fourier transform (both Real and Im parts)
TCanvas *c2 = new TCanvas("FFT Canvas", "Fast Fourier Transform", 800, 600);
c2->Divide(2,1);
c2->Draw();

c2->cd(1);
TH1 *hr = 0;
hr = TH1::TransformHisto(fft_own, hr, "RE");
hr->SetTitle("Real part of array transform");
hr->Draw();
hr->SetStats(kFALSE);

c2->cd(2);
TH1 *him = 0;
him = TH1::TransformHisto(fft_own, him, "IM");
him->SetTitle("Im. part of transform");
him->Draw();
him->SetStats(kFALSE);
him->GetXaxis()->SetLabelSize(0.05);
him->GetYaxis()->SetLabelSize(0.05);

int N = dim;
cout << "N = " << N << endl;
int f_samp = 1/dt;
int f_nyq = f_samp/2;
double df = double(f_samp)/N;
double freqs[N/2+1];

for (int k = 0; k < N/2+1; k++){
	freqs[k] = k*df;
}

cout << "Sampling Frequency is: " << f_samp << endl;
cout << "Nyquist frequency is: " << f_nyq << endl;
cout << "The frequency spacing is: " << df << endl;
cout << "The max frequency plotted is: " << freqs[N/2] << endl;

Double_t *re_full = new Double_t[N];
Double_t *im_full = new Double_t[N];

for(int g = 0; g < N; g++)
fft_own->GetPointComplex(g, re_full[g], im_full[g]);

vector<double> power_spectrum; 
for(int i = 0; i < N/2+1; i++){
	power_spectrum.push_back(pow(re_full[i],2) + pow(im_full[i],2));
}

TCanvas *c3 = new TCanvas("c3", "Power Spectrum");
c3->SetLogy();
c3->SetLogx();
TGraph* gr3 = new TGraph(N/2+1, &freqs[0], &power_spectrum[0]);
gr3->SetTitle("Power Spectrum"); 
gr3->GetXaxis()->SetRangeUser(0.01,10000);
gr3->Draw();

} 
