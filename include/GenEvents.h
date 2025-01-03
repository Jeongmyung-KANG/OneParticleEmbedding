#include <iostream>
#include "jmJetVec.h" 
//#include "JetAnalysis.cxx" 
#include <thread>
#include "fastjet/tools/Subtractor.hh"
#include <cstdlib> 
#include "TF1.h"
#include "TClonesArray.h"
#include "gROOT.h"

const Int_t NumBkgParticles = 800;
const Int_t TrueParticleDummyIndex = -99999;

class BackgroundParticles {
	private : 
		std::vector<fastjet::PseudoJet> Particles;  	
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
		std::vector<fastjet::PseudoJet>GetBackgroundParticles() const {return Particles;} }; 


class MixedEvent : public TObject {
	private : 
		std::vector <fastjet::PseudoJet> Jets; 
		TClonesArray *JetsClonesArray;

	public : 
		MixedEvent() : Jets(0), JetsClonesArray(0){};  
		MixedEvent(std::vector<fastjet::PseudoJet> _Jets, TClonesArray *_JetsClonesArray) : Jets(_Jets), JetsClonesArray(_JetsClonesArray) {}; 
		~MixedEvent() {}; 

		void SetEmbedding(BackgroundParticles *BkgParticles, double JetRadius, Int_t NumEx);
		void SetJetsClonesArray();  
		void ClearVector();

		TClonesArray *GetJetsClonseArray() const {return JetsClonesArray;}

};





