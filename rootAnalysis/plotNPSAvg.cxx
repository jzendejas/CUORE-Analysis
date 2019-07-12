#include <iostream>
#include <vector>
#include <TStyle.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TVirtualFFT.h>

using namespace std;

int N;
int f_samp;
int f_nyq;
double df;


vector<double> singleNPS(TTree *tree){

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
		cout << "Error. Waveform arrays are not the same length!" << endl;
	}
	
	Int_t n_size = dim + 1;
	TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &n_size, "R2C M K");
	fft_own->SetPoints(&wave0[0]);
	fft_own->Transform();
	fft_own->GetPoints(&wave0[0]);

	N = dim;
	f_samp = 1/dt;
	f_nyq = f_samp/2;
	df = double(f_samp)/N;

cout << "N = " << N << endl;
cout << "Sampling Frequency is: " << f_samp << endl;
cout << "Nyquist frequency is: " << f_nyq << endl;
cout << "The frequency spacing is: " << df << endl;

Double_t *re_full = new Double_t[N];
Double_t *im_full = new Double_t[N];

for(int g = 0; g < N; g++)
fft_own->GetPointComplex(g, re_full[g], im_full[g]);

vector<double> power_spectrum;
for(int i = 0; i < N/2+1; i++){
	power_spectrum.push_back(pow(re_full[i],2) + pow(im_full[i],2));
}
w0->clear(); w1->clear(); w2->clear();
return power_spectrum;

}

void plotNPSAvg(TString path, TString prefix, const char *title, int numFiles)
{
	vector<double> sum_ps;
	vector<double> avg_ps;
	vector<double> temp_ps;
	vector<double> freqs;

	for(int filenum = 1; filenum <= numFiles; filenum++){
		temp_ps = {};
		TFile *hfile = new TFile(Form("%s/%s_p%05d.root",path.Data(),prefix.Data(),filenum));	
		TTree* tree = nullptr;
		hfile->GetObject("data_tree",tree);
		cout << "\n" << "Creating Power Spectrum number " << filenum << endl;
		temp_ps = singleNPS(tree);
		if (filenum == 1){
			sum_ps = temp_ps;
		}else{
			for(int i = 0; i < sum_ps.size(); i++){
				sum_ps[i] += temp_ps[i];
			}
		}	

	}
	for (int k = 0; k < N/2+1; k++){
		freqs.push_back(k*df);
	}

		int len = sum_ps.size();
		for(int i = 0; i < len; i++){
			avg_ps.push_back(sum_ps[i]/numFiles/N);
		}	

		TCanvas *c3 = new TCanvas("c3", "Power Spectrum");
		c3->SetLogy();
		c3->SetLogx();
		TGraph* gr3 = new TGraph(len, &freqs[0], &avg_ps[0]);
		gr3->SetTitle(title); 
		gr3->GetXaxis()->SetRangeUser(0.01,10000);
		gr3->GetXaxis()->SetTitle("Frequency (Hz)");
		gr3->GetYaxis()->SetTitle("Spectral Power (V^2/Hz)");
		gr3->Draw();
	}

