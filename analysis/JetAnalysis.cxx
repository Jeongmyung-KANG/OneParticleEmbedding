#include <iostream>
#include "TParticle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom.h"
#include <vector>
#include "TRandom3.h" 
#include "/Users/gangjeongmyeong/flat_study/jmJetVec.h"
#include "fastjet/ClusterSequenceArea.hh"  // use this instead of the "usual" ClusterSequence to get area support
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfold.h"



using namespace std;
using namespace fastjet; 

Double_t get_ME_Rho(Double_t R, vector<PseudoJet> v, int NumEx){
    //cout<<"JET RADIUS : " << R <<endl;
    Double_t ghost_maxrap = 1.0;
    Double_t TrackEtaMax = 1.0;
    Int_t repeat = 1;
    Double_t ghost_area = 0.01;
    JetDefinition jet_def_bkgd(kt_algorithm, R);
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap,repeat,ghost_area));
    Selector selector = SelectorAbsEtaMax(TrackEtaMax - R) * !SelectorNHardest(NumEx);
    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);   
    vector <int> seed;
    seed.push_back(1); 
    seed.push_back(2); 

    bkgd_estimator.set_particles_with_seed(v, seed); 
    Double_t rho = bkgd_estimator.rho();
    return rho;    
}

JetMedianBackgroundEstimator GetBackgroundEstimator(Double_t R, vector<PseudoJet> Particles, Int_t NumEx){
    Double_t ghost_maxrap = 1.0;
    Double_t TrackEtaMax = 1.0;
    Int_t repeat = 1;
    Double_t ghost_area = 0.01;
    JetDefinition jet_def_bkgd(kt_algorithm, R);
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap,repeat,ghost_area));
    Selector selector = SelectorAbsEtaMax(TrackEtaMax - R) * !SelectorNHardest(NumEx) * SelectorPtMin (0.001);
    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);   
    vector <int> seed;
    seed.push_back(1); 
    seed.push_back(2); 

   // bkgd_estimator.set_particles_with_seed(Particles, seed);
    bkgd_estimator.set_particles(Particles);

    return bkgd_estimator;
}

PseudoJet SetTrueParticlePtEtaPhiM(Double_t TruePt, Double_t TrueEta, Double_t TruePhi, Double_t TrueM){
    TLorentzVector TrueParticle;
    TrueParticle.SetPtEtaPhiM(TruePt, TrueEta, TruePhi, TrueM);
    double TruePx = TrueParticle.Px(); 
    double TruePy = TrueParticle.Py(); 
    double TruePz = TrueParticle.Pz(); 
    double TrueE = TrueParticle.E(); 
    
    PseudoJet Particle(TruePx, TruePy, TruePz, TrueE); 

    return Particle;
}

vector<PseudoJet> SetEmbeddedParticlesVector(Int_t NumBkgParticles, TF1* fPt, vector<PseudoJet> BkgJets){
        
    Double_t Bg_pT, Bg_Px, Bg_Py, Bg_Pz, Bg_E, Bg_Eta, Bg_Phi, Bg_P;
    Double_t pion_mass = 0.13975; 

    for(int i_em = 0; i_em < NumBkgParticles; i_em++){  //Generate Thermal BackGround

        Bg_pT = fPt->GetRandom();
        Bg_Eta = gRandom -> Uniform(2) - 1;
        Bg_Phi = gRandom -> Uniform(2*TMath::Pi());

        Bg_Px = Bg_pT*TMath::Cos(Bg_Phi);
        Bg_Py = Bg_pT*TMath::Sin(Bg_Phi);
        Bg_P = Bg_pT*TMath::CosH(Bg_Eta);
		
        Bg_Pz = Bg_pT*TMath::SinH(Bg_Eta);
        Bg_E = TMath::Sqrt(Bg_Px*Bg_Px + Bg_Py*Bg_Py + Bg_Pz*Bg_Pz + pion_mass*pion_mass);
        Double_t E_test = TMath::Sqrt(Bg_P*Bg_P + pion_mass*pion_mass); 

        BkgJets.push_back(PseudoJet(Bg_Px, Bg_Py, Bg_Pz, Bg_E));
    }

    return BkgJets;
}

Int_t IsThisJetTrueJet(PseudoJet Jet, Int_t TrueIndex) { 
    
    Bool_t Is = false;
    vector<PseudoJet> Consties = Jet.constituents(); 
    for (Int_t IndexTrueFind = 0; IndexTrueFind < Consties.size(); IndexTrueFind ++){
        PseudoJet Consti = Consties[IndexTrueFind];
        Int_t JetIndex = Consti.user_index(); 
        if (JetIndex == TrueIndex)
	  return IndexTrueFind; 
	  }

    return -1;
}

void FillRecoilHistogram(vector<PseudoJet> JetVector,
			 TH1F *hRecoilJetPt, TH2F *hRecoilJetPtM,
			 Double_t Rho, Double_t RhoM,
			 Double_t AreaCut, Double_t R,
			 Bool_t IsSE, Double_t TmpRho, TH1F *hCombinatorialJetPt, TH1F *hCombinatorialJetPt_RhoME,
			 TH1F *hTrueJetPt = nullptr,
			 RooUnfoldResponse *ResponseMatrixJetPt = nullptr, 
             RooUnfoldResponse *ResponseMatrixJetPtMJet = nullptr, 
             TH2F *hTrueJetPtJetM = nullptr,
             TH2F *hMatchedRecoJetPtJetM = nullptr, 
             TH2F *hMatchedTrueJetPtJetM = nullptr,
             TH1F *hRecordSequence = nullptr) {

    Int_t FillRecoilJetIndex = 0;
    for (Int_t FillRecoilJetIndex = 0; FillRecoilJetIndex < JetVector.size(); FillRecoilJetIndex++){
                
                Double_t A = JetVector[FillRecoilJetIndex].area();
                Double_t H = JetVector[FillRecoilJetIndex].eta(); 
                Double_t F = TVector2::Phi_0_2pi(JetVector[FillRecoilJetIndex].phi()); 
                Double_t pt = JetVector[FillRecoilJetIndex].pt();
		        Double_t M = JetVector[FillRecoilJetIndex].m(); 
                Double_t Ptc = pt - Rho * A; 
				Double_t TmpPtc = pt - TmpRho * A;
		        Double_t Mc = M - RhoM * A; 
		
				Double_t Ptc_shift = pt - Rho * A; 

                Bool_t Hcut = (TMath::Abs(H) < 1. - R);
                Bool_t Acut = (A > AreaCut);
                Bool_t isRecoil = (F < TMath::Pi() / 2);
                Bool_t isInAcceptance = (Hcut && Acut && isRecoil);
                    
                Int_t JetIndex = IsThisJetTrueJet(JetVector[FillRecoilJetIndex], 9999);
                if (JetIndex != -1) {
                    hTrueJetPt->Fill(Ptc);
                    hRecordSequence->Fill(FillRecoilJetIndex);
                }
                if (Hcut && Acut && isRecoil) {
                    hRecoilJetPt->Fill(Ptc);
			        hRecoilJetPtM->Fill(Ptc, Mc);
					if (JetIndex == -1 && IsSE == true) {
						hCombinatorialJetPt->Fill(Ptc); 
						hCombinatorialJetPt_RhoME->Fill(TmpPtc);
					}
                if (JetIndex != -1) { 
					/*
                PseudoJet MatchedTrueJet = JetVector[FillRecoilJetIndex].constituents()[JetIndex];
                Double_t TrueJetPt = MatchedTrueJet.pt();
                Double_t TrueJetM = MatchedTrueJet.m(); 
                ResponseMatrixJetPt->Fill(Ptc, TrueJetPt);
                ResponseMatrixJetPtMJet->Fill(Ptc, Mc, TrueJetPt, TrueJetM); 
                hMatchedRecoJetPtJetM->Fill(Ptc, Mc);
                hMatchedTrueJetPtJetM->Fill(TrueJetPt, TrueJetM); 

				*/
                //cout << "Matching Result : " << Ptc << " " << TrueJetPt << " " << " Rho M : " << " " <<Mc << " " << TrueJetM << " " << JetIndex <<  endl; 
			}
        } 
    }
}
