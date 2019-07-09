#include <iostream>
#include <vector>
#include <TStyle.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include "TVirtualFFT.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>
#include <set>
#include "TMultiGraph.h"


using namespace std;

int N;
int f_samp;
int f_nyq;
double df;


void singleNPS(TTree *tree, vector<double>& psx, vector<double>& psy, vector<double>& psz, int filenum){

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

	N = dim;
	f_samp = 1/dt;
	f_nyq = f_samp/2;
	df = double(f_samp)/N;


	cout << "Transforming X..." << endl;
	TVirtualFFT *fft_ownx = TVirtualFFT::FFT(1, &n_size, "R2C M K");
	fft_ownx->SetPoints(&wave0[0]);
	fft_ownx->Transform();
	fft_ownx->GetPoints(&wave0[0]);

	cout << "Transforming Y..." << endl;	
	TVirtualFFT *fft_owny = TVirtualFFT::FFT(1, &n_size, "R2C M K");
	fft_owny->SetPoints(&wave1[0]);
	fft_owny->Transform();
	fft_owny->GetPoints(&wave1[0]);
	
	cout << "Transforming Z..." << endl;
	TVirtualFFT *fft_ownz = TVirtualFFT::FFT(1, &n_size, "R2C M K");
	fft_ownz->SetPoints(&wave2[0]);
	fft_ownz->Transform();
	fft_ownz->GetPoints(&wave2[0]);

	Double_t *re_full_x = new Double_t[N];
	Double_t *im_full_x = new Double_t[N];
	Double_t *re_full_y = new Double_t[N];
	Double_t *im_full_y = new Double_t[N];
	Double_t *re_full_z = new Double_t[N];
	Double_t *im_full_z = new Double_t[N];


	for(int g = 0; g < N; g++){
		fft_ownx->GetPointComplex(g, re_full_x[g], im_full_x[g]);
		fft_owny->GetPointComplex(g, re_full_y[g], im_full_y[g]);
		fft_ownz->GetPointComplex(g, re_full_z[g], im_full_z[g]);
	}
	if(filenum == 0){
		for(int i = 0; i < (N/2); i++){
			psx.push_back(pow(re_full_x[i],2) + pow(im_full_x[i],2));
			psy.push_back(pow(re_full_y[i],2) + pow(im_full_y[i],2));
			psz.push_back(pow(re_full_z[i],2) + pow(im_full_z[i],2));
		}
	}else{
		for(int i = 0; i < (N/2); i++){
			psx[i] += (pow(re_full_x[i],2) + pow(im_full_x[i],2));
			psy[i] += (pow(re_full_y[i],2) + pow(im_full_y[i],2));
			psz[i] += (pow(re_full_z[i],2) + pow(im_full_z[i],2));
		}
	}
	w0->clear(); w1->clear(); w2->clear();
}

void plotNPSAvgXYZ(TString path, TString prefix, int nFiles){

	vector<double> avg_psx;
	vector<double> avg_psy;
	vector<double> avg_psz;
	vector<double> temp_psx;
	vector<double> temp_psy;
	vector<double> temp_psz;
	vector<double> freqs;

	for(int filenum = 0; filenum < nFiles; filenum++){
		TString inFile = Form("%s/%s_p%05d.root",path.Data(),prefix.Data(),filenum+1);
		std::cout << inFile.Data() << std::endl;
		TFile *hfile = new TFile(inFile);
		TTree* tree = nullptr;
		hfile->GetObject("data_tree",tree);

		cout << "\n" << "Creating Power Spectrum number " << filenum + 1 << endl;
		singleNPS(tree,temp_psx,temp_psy,temp_psz,filenum);
		//temp_psx.clear(); temp_psy.clear(); temp_psz.clear();    

	}


	for(int i = 0; i < temp_psx.size(); i++){
		temp_psx[i] /= (nFiles*N);
		temp_psy[i] /= (nFiles*N);
		temp_psz[i] /= (nFiles*N);
	}

	for (int k = 0; k < N/2; k++){
		freqs.push_back(k*df);
	}
	cout << "Frequency axis created..." << endl;
	cout << "Avg vectors made." << endl;
	cout << "Size of ps x,y,z = " << temp_psx.size() << " " << temp_psy.size() << " " << temp_psz.size() << endl;

	cout << "N = " << N << endl;

	TCanvas* c0 = new TCanvas("nps_xyz", "NPS", 1600, 1200);

	TGraph* grx = new TGraph(N/2, &freqs[0], &temp_psx[0]);
	TGraph* gry = new TGraph(N/2, &freqs[0], &temp_psy[0]);
	TGraph* grz = new TGraph(N/2, &freqs[0], &temp_psz[0]);

	TMultiGraph* mg = new TMultiGraph();

	c0->SetLogx();
	c0->SetLogy();

	grx->SetTitle("X");
	grx->SetLineColor(kRed);
	grx->SetLineWidth(2);
	gry->SetTitle("Y");
	gry->SetLineColor(kGreen);
	gry->SetLineWidth(2);
	grz->SetTitle("Z");
	grz->SetLineColor(kBlue);
	grz->SetLineWidth(2); 
	mg->Add( grx );
	mg->Add( gry );
	mg->Add( grz );
	mg->GetXaxis()->SetLimits(0.01,f_nyq);
	mg->SetMinimum(1E-08);
	mg->SetMaximum(100);

	mg->SetTitle(Form("NPS Averaged over %d runs: Speaker Tests; frequency (Hz); Noise Power Density (V^2/Hz)", nFiles));

	mg->Draw("AL");
	c0->BuildLegend();
	c0->Modified();
	c0->Update();
	c0->SaveAs(Form("SpeakerTests_%d.png", nFiles));
	TFile* file1 = new TFile(Form("SpeakerTests_%d.root",nFiles),"recreate");
	file1->cd();
	c0->Write("c0");
	file1->Close();
}

//#ifndef __CINT__
int main(int argc, char **argv) {
	//char* workingDir = getenv(“CUORE_EXT_INSTALL”);

	//gErrorIgnoreLevel = kWarning;

	if(argc!=4)
	{
		cout<<"Run as ./"<<argv[0]<<" <path> <prefix> <nFiles>"<<endl;
		return 1;
	}

	TString path = Form("%s",argv[1]);
	TString prefix = TString(argv[2]);
	int nFiles = atoi(argv[3]);

	plotNPSAvgXYZ(path, prefix, nFiles);

	cout << "N = " << N << endl;
	cout << "Sampling Frequency is: " << f_samp << endl;
	cout << "Nyquist frequency is: " << f_nyq << endl;
	cout << "The frequency spacing is: " << df << endl;
	return 0;
}
//#endif /* __CINT __ */
