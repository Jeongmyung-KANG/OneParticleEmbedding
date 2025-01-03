#include <iostream>
#include "include/jmJetVec.h" 
//#include "JetAnalysis.cxx" 
//#include "RooUnfoldResponse.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfold.h"
#include <thread>
#include "fastjet/tools/Subtractor.hh"
#include <cstdlib> // system() �Լ� ����� ���� ���



const Int_t NumBkgParticles = 800;
const Int_t TrueParticleDummyIndex = -99999;

class BackgroundParticles {
	private : 
		vector<fastjet::PseudoJet> Particles;  	
		TF1* fPt; // = new TF1("fPt", Form("x*exp(-2.0*x/%s)", ParmBackgroundParticlesPt.Data()), 0.2, 30.0);   
	public :   
		BackgroundParticles() : Particles(0), fPt(0){}; 
		BackgroundParticles(double _Particles, TF1 *_fPt) : Particles(_Particles), fPt(_fPt){};
		~BackgroundParticles() {}; 
															
		//BackgroundParticles () : PseudoJet BackgroundParticles(0);  
		void ClearVector ();
		void SetBackgroundPt       (TString ParmName); 
		void SetBackgroundParticles(int Seed);   
		void SetTrueParticlePtPhiEtaM(Double_t TruePt, Double_t TruePhi, Double_t TrueEta, Double_t TrueM); 
		vector<fastjet::PseudoJet>GetBackgroundParticles() const {return Particles;} }; 


/*
class SameEvent : public TObject { 
		private : PseudoJet TrueParticle; vector <PseudoJet> Jets;  
		public : //SameEvent() : TrueParticle(0), Jets(0) {};  
		void SetTrueParticlePtPhiEtaM(Double_t TruePt, Double_t TruePhi, Double_t TrueEta, Double_t TrueM);
};
*/ 

class MixedEvent : public TObject {
	private : 
		vector <fastjet::PseudoJet> Jets; 
		TClonesArray *JetsClonesArray;

	public : 
		MixedEvent() : Jets(0), JetsClonesArray(0){};  
		MixedEvent(vector<fastjet::PseudoJet> _Jets, TClonesArray *_JetsClonesArray) : Jets(_Jets), JetsClonesArray(_JetsClonesArray) {}; 
		~MixedEvent() {}; 

		void SetEmbedding(BackgroundParticles *BkgParticles, double JetRadius, Int_t NumEx);
		void SetJetsClonesArray();  
		void ClearVector();

		TClonesArray *GetJetsClonseArray() const {return JetsClonesArray;}

};


inline void BackgroundParticles::SetTrueParticlePtPhiEtaM(Double_t TruePt, Double_t TruePhi, Double_t TrueEta, Double_t TrueM){
	TLorentzVector TLV_TrueParticle; 
	fastjet::PseudoJet TrueParticle;
	TLV_TrueParticle.SetPtEtaPhiM(TruePt, TrueEta, TruePhi, TrueM);  

	Double_t TmpTruePx = TLV_TrueParticle.Px();
	Double_t TmpTruePy = TLV_TrueParticle.Py();
	Double_t TmpTruePz = TLV_TrueParticle.Pz();
	Double_t TmpTrueE =  TLV_TrueParticle.E();  

	TrueParticle = fastjet::PseudoJet(TmpTruePx, TmpTruePy, TmpTruePz, TmpTrueE);
	TrueParticle.set_user_index(TrueParticleDummyIndex);
	Particles.push_back(TrueParticle); 

}

inline void BackgroundParticles::ClearVector() {
	Particles.clear(); 
} 

inline void MixedEvent::ClearVector() { 
	Jets.clear(); 
}

inline void BackgroundParticles::SetBackgroundPt(TString ParmName){ 
	fPt = new TF1("fPt", Form("x*exp(-2.0*x/%s)", ParmName.Data()), 0.2, 30.0);   
}

inline void BackgroundParticles::SetBackgroundParticles(int Seed){ 
	Particles  = SetEmbeddedParticlesVector(NumBkgParticles, fPt, Particles); 
}

inline void MixedEvent::SetJetsClonesArray(){
	JetsClonesArray = new TClonesArray("jmJetVec", 1000000); 
}

inline void MixedEvent::SetEmbedding(BackgroundParticles *BkgParticles, double JetRadius, Int_t NumEx){
	
	fastjet::JetDefinition   jet_def(antikt_algorithm, JetRadius); 
	Double_t ghost_maxrap = 1.0; 
	Double_t TrackEtaMax = 1.0;
	
	vector <int> seed;
    seed.push_back(1); 
    seed.push_back(2); 


	GhostedAreaSpec  area_spec(ghost_maxrap); 
	AreaDefinition  area_def(active_area, area_spec);
	//ClusterSequence testJets(BkgParticles->GetBackgroundParticles(), jet_def);
	fastjt::ClusterSequenceArea testJets(BkgParticles->GetBackgroundParticles(), jet_def, area_def);
	vector<fastjet::PseudoJet> sort1 = sorted_by_pt(testJets.inclusive_jets());   
    Int_t repeat = 1;
    Double_t ghost_area = 0.01;
    JetDefinition jet_def_bkgd(kt_algorithm, JetRadius);
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap,repeat,ghost_area));
    Selector selector = SelectorAbsEtaMax(TrackEtaMax - JetRadius) * !SelectorNHardest(NumEx) * SelectorPtMin (0.001);
    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);  

	vector <int> vec_seed;
	vec_seed.push_back(1); 
	vec_seed.push_back(2); 

    bkgd_estimator.set_particles_with_seed(BkgParticles->GetBackgroundParticles(), vec_seed); 
	bkgd_estimator.set_particles(BkgParticles->GetBackgroundParticles()); 
	
	int TestNum = 0;
	Double_t rho = bkgd_estimator.rho(); 
	Double_t rho_m = bkgd_estimator.rho_m(); 
	//cout << sort1.size() << " " << BkgParticles->GetBackgroundParticles().size() << endl; 
	for (int tmp_particles_index = 0; tmp_particles_index < 2; tmp_particles_index ++){
		//cout << BkgParticles->GetBackgroundParticles()[tmp_particles_index].pt() << " " << BkgParticles->GetBackgroundParticles()[tmp_particles_index].eta() << BkgParticles->GetBackgroundParticles()[tmp_particles_index].phi() << endl; 
	}
	for (int tmp_jet_index = 0; tmp_jet_index< sort1.size(); tmp_jet_index ++){
		jmJetVec *JetsVec = new((*JetsClonesArray)[tmp_jet_index]) jmJetVec(); 
		Double_t JetArea = sort1[tmp_jet_index].area(); 
		//Double_t JetArea = 0.;	
		Double_t JetPx = sort1[tmp_jet_index].px();
		Double_t JetPy = sort1[tmp_jet_index].py();
		Double_t JetPz = sort1[tmp_jet_index].pz();
		Double_t JetE = sort1[tmp_jet_index].e();
		Double_t JetM = sort1[tmp_jet_index].m(); 
		Double_t JetPhi = sort1[tmp_jet_index].phi();
		Double_t JetEta = sort1[tmp_jet_index].eta();
		vector<PseudoJet> JetConsti = sort1[tmp_jet_index].constituents();
		vector<PseudoJet> consti = sort1[tmp_jet_index].constituents(); 
		vector<double> vec_consti_phi;
		vector<double> vec_consti_eta;  
		vector<double> vec_consti_pt; 
		JetsVec->setA(JetArea); 
		JetsVec->setJetKin(JetPx, JetPy, JetPz, JetE, JetM);  
		JetsVec->setJetPhiEta(JetPhi, JetEta); 
		JetsVec->setRho(rho); 
		JetsVec->setRhoM(rho_m);
		JetsVec->setConstiSize(consti.size()); 
		bool IsTrueJet = false;

		for (int tmp_consti_index = 0; tmp_consti_index < consti.size(); tmp_consti_index++){
			PseudoJet consti_jet = consti[tmp_consti_index];   
			vec_consti_phi.push_back(consti_jet.phi()); 
			vec_consti_eta.push_back(consti_jet.eta()); 
			vec_consti_pt.push_back(consti_jet.pt()); 

			int tmp_consti_jet_true_index = consti_jet.user_index();  
			if (tmp_consti_jet_true_index == TrueParticleDummyIndex) IsTrueJet = true; 
			TestNum += 1;
		}
		
		JetsVec->setConstiVectorPhi(vec_consti_phi); 
		JetsVec->setConstiVectorEta(vec_consti_eta);  
		JetsVec->setConstiVectorPt(vec_consti_pt); 

		if (IsTrueJet) JetsVec->setJetIndex(TrueParticleDummyIndex);
	}

	cout << TestNum << endl; 
	
}


