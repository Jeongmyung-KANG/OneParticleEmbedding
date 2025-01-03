#include <iostream>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TParticle.h" 
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom.h"
#include <vector> 
#include "TROOT.h" 
//#include "/alice/data/kgm8091/gradTask_star/jetFinde_mc/jmJetVec.h"  
#include "jmJetVec.h"


const int TrueJetIndex = -99999; 

void FillPhiHistogram(double phi, vector<double> v, TH1F *h){
	for (int z = 0; z < v.size(); z++){
		double data = TVector2::Phi_0_2pi(v[z]) - TVector2::Phi_0_2pi(phi); 
		h->Fill(TVector2::Phi_mpi_pi(data)); 
	}
}

double GetDeltaR (double th, double ch, double tf,double  cf){
	double dEta = th - ch;
	double dPhi = TVector2::Phi_mpi_pi(TVector2::Phi_mpi_pi(tf) - TVector2::Phi_mpi_pi(cf));  
	return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

int GetDeltaRCase (double dR) {
	int val = -999; 
	if (dR < 0.4) val = 0; 
	if (dR >= 0.4 && dR < 0.8) val = 1; 
	if (dR >= 0.8 && dR < 1.2) val = 2;
	if (dR >= 1.2 && dR < 1.6) val = 3;
	if (dR >= 1.6) val = 4; 

	return val;
}

bool IsRecoil(double x){
	bool out = false; 
	double phi = TVector2::Phi_0_2pi(x); 
	if (phi < TMath::Pi() / 2) out = true; 

	return out;
}

bool IsNear(double x){

	double phi = TVector2::Phi_0_2pi(x);  
	bool out = false;

	if (phi > TMath::Pi()/4 - 0.4 && phi < TMath::Pi() /4 + 0.4) out = true;
	return out; 
}

bool IsOut (double x){

	double phi = TVector2::Phi_0_2pi(x); 
	bool out = false; 
	if (phi > TMath::Pi()) out = true;
	return out;
}

bool IsMatchedJet (jmJetVec *j1, jmJetVec *j2){
	bool Val = false; 
	return Val;
}

void DeltaR(){
	
	TString NameHist[5] = {"core", "near", "middle", "outter", "out"}; 

	TH2F *hPhiEta_SE[5];  
	TH2F *hPhiEta_SE_Area[5];  
	TH1F *hPt_SE[5]; 
	TH1F *hPt_SE_Area[5]; 
	TH1F *hPt_SE_Miss = new TH1F("hPt_SE_Miss", "hPt_SE_Miss", 140, -20, 50); 
	TH1F *hPt_SE_Matched = new TH1F("hPt_SE_Matched", "hPt_SE_Matched", 140, -20, 50); 

	TH1F *hDeltaPhi_SE[5];		 
	TH1F *hDeltaPhi_SE_Area[5];		 
	TH1F *hPt_SE_Inc = new TH1F("hPt_SE_Inc", "hPt_SE_Inc", 140, -20, 50); 
	TH1F *hPt_SE_Inc_Area = new  TH1F("hPt_SE_Inc_Area", "hPt_SE_Inc_Area", 140, -20, 50); 

	TH1F *hDeltaPhi_ME = new TH1F("hDeltaPhi_ME", "hDeltaPhi_ME", 200, -4, 4);   
	TH1F *hDeltaPhi_ME_Area = new TH1F("hDeltaPhi_ME_Area", "hDeltaPhi_ME", 200, -4, 4); 

	TH2F *hPhiEta_ME = new TH2F("hPhiEta_ME", "hPhiEta_ME", 100, -4, 4, 50, -1, 1);      
	TH2F *hPhiEta_ME_Area = new TH2F("hPhiEta_ME_Area", "hPhiEta_ME_Area", 100, -4, 4, 50, -1, 1); 
	TH1F *hPt_ME = new TH1F("hPt_ME", "hPt_ME", 140, -20, 50); 
	TH1F *hPt_ME_Area = new TH1F("hPt_ME_Area", "hPt_ME_Area", 140, -20, 50); 
	TH1F *hCombinatorialPt = new TH1F("hCombonatorialPt", "hCombinatorialPt", 140, -20, 50); 
	TH1F *hPt_SE_Miss_IncludeTrue = new TH1F("hPt_SE_Miss_IncludeTrue", "hPt_SE_Miss_IncludeTrue", 140, -20, 50); 

	hDeltaPhi_ME->Sumw2(); 
	hDeltaPhi_ME_Area->Sumw2(); 
	hPt_SE_Miss->Sumw2(); 
	hPt_SE_Matched->Sumw2();
	hPt_ME_Area->Sumw2(); 
	hPt_ME->Sumw2(); 
	hPt_SE_Inc->Sumw2();   
	hPt_SE_Miss_IncludeTrue->Sumw2(); 
	hPt_SE_Inc_Area->Sumw2();   
	hCombinatorialPt->Sumw2(); 

	TH2F *hM_SE[5]; 
	for (int i = 0; i < 5; i ++){
		hPhiEta_SE[i] = new TH2F(Form("hPhiEta_SE_%s", NameHist[i].Data()), Form("hPhiEta_SE_%s", NameHist[i].Data()), 100, -4, 4, 50, -1, 1); 
		hPhiEta_SE_Area[i] = new TH2F(Form("hPhiEta_SE_Area_%s", NameHist[i].Data()), Form("hPhiEta_SE_Area_%s", NameHist[i].Data()), 100, -4, 4, 50, -1, 1); 
		hDeltaPhi_SE[i] = new TH1F(Form("hDeltaPhi_SE_%s", NameHist[i].Data()),Form("hDeltaPhi_SE_%s", NameHist[i].Data()) , 200, -4, 4);   
		hDeltaPhi_SE_Area[i] = new TH1F(Form("hDeltaPhi_SE_Area_%s", NameHist[i].Data()),Form("hDeltaPhi_SE_Area_%s", NameHist[i].Data()) , 200, -4, 4);   
		hM_SE[i]  = new TH2F(Form("hM_SE_%s", NameHist[i].Data()), Form("hM_SE_%s", NameHist[i].Data()), 140, -20, 50, 90, 0, 30);   
		
		hPt_SE[i] = new TH1F(Form("hPt_SE_%s", NameHist[i].Data()), Form("hPt_SE_%s", NameHist[i].Data()), 140, -20, 50); 
		hPt_SE_Area[i] = new TH1F(Form("hPt_SE_Area_%s", NameHist[i].Data()), Form("hPt_SE_Area_%s", NameHist[i].Data()), 140, -20, 50); 

		hDeltaPhi_SE[i]->Sumw2();  
		hDeltaPhi_SE_Area[i]->Sumw2(); 
		hM_SE[i]->Sumw2();
		hDeltaPhi_SE[i]->Sumw2();  
		hPt_SE[i]->Sumw2(); 
		hPt_SE_Area[i]->Sumw2(); 
	}

	gROOT->ProcessLine(".L ./jmJetVec.h"); 
	TFile *fSE = TFile::Open("SE.root"); 
	TFile *fME = TFile::Open("ME.root"); 

	TTree *SEjetTree = (TTree*)fSE->Get("JetTree"); 
	TTree *MEjetTree = (TTree*)fME->Get("JetTree"); 
	
	TClonesArray *SE_Tracks = new TClonesArray("jmJetVec", 10000); 
	TClonesArray *ME_Tracks = new TClonesArray("jmJetVec", 10000);  
	
	SEjetTree->SetBranchAddress("Jets", &SE_Tracks); 
	MEjetTree->SetBranchAddress("Jets", &ME_Tracks); 

	int ME_N_Events = SEjetTree->GetEntries(); 
	int SE_N_Events = MEjetTree->GetEntries(); 
	
	int Tot_num_ME = 0;  
	for (int i = 0; i < SE_N_Events; i++){  

		cout << i << endl; 

		SE_Tracks->Clear(); 
		ME_Tracks->Clear(); 

		SEjetTree->GetEntry(i); 
		MEjetTree->GetEntry(i);  
		
		int Size_SE = SE_Tracks->GetEntriesFast(); 
		int Size_ME = ME_Tracks->GetEntriesFast(); 
		cout << " " << Size_SE << " " << Size_ME << endl; 
		double TrueJetPhi = 0.0; 
		double TrueJetEta = 0.0; 
		
		for (int j = 0; j < Size_ME; j++){  
			jmJetVec *jet3 = (jmJetVec*)ME_Tracks->At(j);  
			double pt3 = jet3->pt();  
			double rho3 = jet3->rho(); 
			double a3 = jet3->A(); 
			double eta3 = jet3->eta();  
			double phi3 = jet3->phi();
			double ptc3 = pt3 - a3 * rho3; 
			vector<double> PhiVector3 = jet3->constiVectorPhi();   
			vector<double> EtaVector3 = jet3->constiVectorEta();  

			
			Tot_num_ME += PhiVector3.size();
			if (TMath::Abs(eta3) > 1.0 - 0.4) continue;
			if (pt3 < 0.2) continue;  

			//cout << "ME Jet Index : " << j << " " << pt3 << " " << eta3 << " " << phi3 << endl; 
			hPhiEta_ME->Fill(phi3, eta3); 
			FillPhiHistogram(phi3, PhiVector3, hDeltaPhi_ME); 
			hPt_ME->Fill(pt3); 
			if (a3 > 0.35) hPt_ME_Area->Fill(pt3);  
			hPhiEta_ME_Area->Fill(phi3, eta3);  
			FillPhiHistogram(phi3, PhiVector3, hDeltaPhi_ME_Area); 
		}
		
		int NumUnmatched = 0; 

		for (int j = 0; j < Size_SE; j++){
			jmJetVec *jet1 = (jmJetVec*)SE_Tracks->At(j); 
			bool isMatched = false;
			double pt1 = jet1->pt();   
			if (pt1 < 0.2) continue;  
			
			double phi1 = jet1->phi(); 
			double eta1 = jet1->eta();  

			if (TMath::Abs(eta1) > 1.0 - 0.4) continue; 
			double a1 = jet1->A(); 
			double m1 = jet1->m();  
			double rho1 = jet1->rho(); 
			double ptc1 = pt1  -rho1 * a1;
			int jetIndex = jet1->jetIndex();
			if (jetIndex != TrueJetIndex){
				//cout << "ss" << endl; 
				hCombinatorialPt->Fill(pt1); 
			} 
			vector <double> ConstiPhi = jet1->constiVectorPhi();   

			//cout << "SE Jet Index : " << j << " " << pt1 << " " << eta1 << " " << phi1<< endl; 
			hPt_SE_Inc->Fill(pt1);  

			if (a1 > 0.35){
				hPt_SE_Inc_Area->Fill(pt1); 
			}

			if (jetIndex == TrueJetIndex) {
				TrueJetPhi = phi1; 
				TrueJetEta = eta1;  
				hPt_SE_Miss_IncludeTrue->Fill(pt1); 
			}
			
			for (int k = 0; k < Size_ME; k++){  
			
				jmJetVec *jet2 = (jmJetVec*) ME_Tracks->At(k); 
				double pt2 = jet2->pt(); 
				if (pt2 < 0.2) continue;
				double phi2 = jet2->phi();
				double eta2 = jet2->eta(); 
				if (TMath::Abs(eta2) < 1.0 - 0.4) { 
					double a2 = jet2->A(); 
					double m2 = jet2->m();
					double rho2 = jet2->rho(); 
					vector <double> MEConstiPhi = jet2->constiVectorPhi();  
					if (pt1 == pt2) {
						isMatched = true;
					}
				}
			}

			if (isMatched){
				if (TMath::Abs(eta1) < 1.0 - 0.4){
					hPt_SE_Matched->Fill(pt1); 
				}
			}
			
			if (!isMatched){ 

				double DeltaR = GetDeltaR(TrueJetEta, eta1, TrueJetPhi, phi1); 
				int FillHistIndex = GetDeltaRCase(DeltaR);   
				//cout << Form("%.6f", TrueJetPhi) << " " << Form("%.6f",TrueJetEta) << " " << Form("%.6f", DeltaR) << " " << FillHistIndex << " " << phi1 << " " << eta1 << endl; 
				
				//cout << i << " " << NumUnmatched << " " << j << " " << pt1  << " " << ptc1 << " " << DeltaR <<  " " << eta1 << " " << phi1 << endl;  
				if (TMath::Abs(eta1) < 1.0 - 0.4) {
					hPt_SE_Miss->Fill(pt1);
				}

				FillPhiHistogram(phi1, ConstiPhi, hDeltaPhi_SE[FillHistIndex]);
				hPhiEta_SE[FillHistIndex]->Fill(TVector2::Phi_mpi_pi(phi1), eta1);  
				hM_SE[FillHistIndex]->Fill(pt1, m1);  
				hPt_SE[FillHistIndex]->Fill(pt1); 

				if (a1 > 0.35) {  
					FillPhiHistogram(phi1, ConstiPhi, hDeltaPhi_SE_Area[FillHistIndex]);
					hPhiEta_SE_Area[FillHistIndex]->Fill(TVector2::Phi_mpi_pi(phi1), eta1);    
					hPt_SE_Area[FillHistIndex]->Fill(pt1); 
				}
				 
				NumUnmatched += 1; 
			}

			
		}
	}


	cout << Tot_num_ME << endl; 

	TFile *fout = new TFile("dR.root", "RECREATE"); 

	fout->cd();  
	for (int i = 0; i < 5; i++){
		hDeltaPhi_SE[i]->Write(); 
		hDeltaPhi_SE_Area[i]->Write(); 
		hPhiEta_SE[i]->Write(); 
		hPhiEta_SE_Area[i]->Write();
		hM_SE[i]->Write();  
		hPt_SE[i]->Write(); 
		hPt_SE_Area[i]->Write(); 
	} 
	
	hPt_SE_Inc->Write(); 
	hPt_SE_Inc_Area->Write(); 
	hDeltaPhi_ME->Write(); 
	hDeltaPhi_ME_Area->Write(); 
	hPhiEta_ME->Write(); 
	hPhiEta_ME_Area->Write(); 
	hPt_ME->Write(); 
	hPt_ME_Area->Write();  
	hPt_SE_Miss->Write(); 
	hPt_SE_Matched->Write(); 
	hCombinatorialPt->Write();   
	hPt_SE_Miss_IncludeTrue->Write(); 
	fout->Close();
 }


