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
#include <TLegend.h>
#include <numeric>
#include <fftw3.h>


using namespace std;

int N;
int f_samp;
int f_nyq;
double df;


std::vector<double> singleNPS(TTree *tree, TString branchAddress, int filenum){
	vector<double> *w0 = new vector<double>;
	vector<double> wave0;
	std::vector<double> psx;
	unsigned int initTime = 0;
	unsigned int nSamples = 0;
	double dt = 0.0;

	int nEntries = tree->GetEntries();

	tree->SetBranchAddress(branchAddress,&w0);
	tree->SetBranchAddress("NumberOfSamples", &nSamples);
	tree->SetBranchAddress("SamplingWidth_s", &dt);

	tree->GetEntry(0);
	tree->SetBranchAddress("Timestamp_s", &initTime);

	for(int j = 0; j < nEntries; j++){
		tree->GetEntry(j);
		wave0.insert(wave0.end(), w0->begin(), w0->end());
	}
	N = wave0.size();
	Int_t n_size = N + 1;
	f_samp = 1/dt;
	f_nyq = f_samp/2;
	df = double(f_samp)/N;
	
	TVirtualFFT *fft_ownx = TVirtualFFT::FFT(1, &n_size, "R2C ES");
	fft_ownx->SetPoints(&wave0[0]);
	fft_ownx->Transform();
	fft_ownx->GetPoints(&wave0[0]);

	Double_t *re_full_x = new Double_t[N];
	Double_t *im_full_x = new Double_t[N];

	for(int g = 0; g < N; g++){
		fft_ownx->GetPointComplex(g, re_full_x[g], im_full_x[g]);
	}
	for(int i = 0; i < (N/2); i++){
		psx.push_back(pow(re_full_x[i],2) + pow(im_full_x[i],2));
	}
	w0->clear();
	delete re_full_x;
	delete im_full_x;
	delete fft_ownx;
	delete w0;
	wave0.clear(); 
	
	return psx; 
}

void plot3dNPS(TString path, TString prefix, int nFiles, int initfile = 0){

	vector<double> avg_psx; vector<double> avg_psy; vector<double> avg_psz;
	vector<double> freqs;
	for(int filenum = initfile; filenum < nFiles + initfile; filenum++){
		TString inFile = Form("%s/%s_p%05d.root",path.Data(),prefix.Data(),filenum+1);
		std::cout << "\n" << inFile.Data() << std::endl;
		TFile *hfile = new TFile(inFile);
		TTree* tree = nullptr;
		hfile->GetObject("data_tree",tree);

		cout << "Creating Power Spectrum number " << filenum + 1 << endl;

		cout << "Transforming X..." << endl;
		TString xname("Waveform000");
		std::vector<double> temp_psx = singleNPS(tree,xname,filenum);
			
		cout << "Transforming Y ..." << endl;
		TString yname("Waveform001");
		std::vector<double> temp_psy = singleNPS(tree,yname,filenum);

		cout << "Transforming Z ..." << endl;
		TString zname("Waveform002");
		std::vector<double> temp_psz = singleNPS(tree,zname,filenum);
		
		if (filenum == initfile){
			for(int i = 0; i < temp_psx.size(); i++){
				avg_psx.push_back(temp_psx[i]/(nFiles*N));
				avg_psy.push_back(temp_psy[i]/(nFiles*N));
				avg_psz.push_back(temp_psz[i]/(nFiles*N));
			}
		}
		else {
			for(int i = 0; i < temp_psx.size(); i++){
				avg_psx[i] += temp_psx[i]/(nFiles*N);
				avg_psy[i] += temp_psy[i]/(nFiles*N);
				avg_psz[i] += temp_psz[i]/(nFiles*N);
			}
		}
/*		cout << "psx: ";
		for (int k = 1; k < avg_psx.size(); k = k + 10000){
			cout << temp_psx[k] << " ";
		}
		cout << "\npsy: "; 
		for (int k = 1; k < avg_psy.size(); k = k + 10000){
			cout << temp_psy[k] << " ";
		}
		cout << "\n";
*/
		double avgX = accumulate( temp_psx.begin(),temp_psx.end(), 0.0)/temp_psx.size();
		cout << "The average of psx is " << avgX << endl;	
		
		double avgY = accumulate( temp_psy.begin(),temp_psy.end(), 0.0)/ temp_psy.size();
		cout << "The average of psy is " << avgY << endl;	
	

	temp_psx.clear(); temp_psy.clear(); temp_psz.clear();
	}
	for (int k = 0; k < N/2; k++){
		freqs.push_back(k*df);
	}
	
	TCanvas* c0 = new TCanvas("nps_xyz", "NPS", 1600, 1200);

	TGraph* grx = new TGraph(N/2, &freqs[0], &avg_psx[0]);
	TGraph* gry = new TGraph(N/2, &freqs[0], &avg_psy[0]);
	TGraph* grz = new TGraph(N/2, &freqs[0], &avg_psz[0]);

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
//	mg->Add( grz );
	//mg->GetXaxis()->SetLimits(1,1000);
	//mg->SetMinimum(1E-08);
	//mg->SetMaximum(100);

	mg->SetTitle(Form("NPS Averaged over %d runs: Speaker Tests; frequency (Hz); Noise Power Density (V^2/Hz)", nFiles));

	mg->Draw("AL");
	TLegend *l = new TLegend(1,2,5,10);
	c0->BuildLegend();
	c0->Modified();
	c0->Update();
	gPad->SetGridx(); gPad->SetGridy();
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

		if(argc!=5)
		{
			cout<<"Run as "<<argv[0]<<" <path> <prefix> <nFiles> <initfilenumber>"<<endl;
			return 1;
		}

		TString path = Form("%s",argv[1]);
		TString prefix = TString(argv[2]);
		int nFiles = atoi(argv[3]);
		int initfile = atoi(argv[4]); 
		plot3dNPS(path, prefix, nFiles, initfile);

		cout << "N = " << N << endl;
		cout << "Sampling Frequency is: " << f_samp << endl;
		cout << "Nyquist frequency is: " << f_nyq << endl;
		cout << "The frequency spacing is: " << df << endl;
		return 0;
	}
	//#endif /* __CINT __ */
