#include <iostream>
#include "jmJetVec.h"
//#include "JetAnalysis.cxx"
#include <thread>
#include "fastjet/tools/Subtractor.hh"
#include <cstdlib> // system() �Լ� ����� ���� ���
#include <cassert> 
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
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
#include "fastjet/ClusterSequenceArea.hh"  // use this instead of the "usual" ClusterSequence to get area support
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfold.h"

using namespace fastjet;
using namespace std; 

using namespace fastjet;
using namespace std;
const Int_t NumBkgParticles = 800;
const Int_t TrueParticleDummyIndex = -99999;

class BackgroundParticles {
    private :
        vector<PseudoJet> Particles;
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
        vector<PseudoJet>GetBackgroundParticles() const {return Particles;} 
        
        };


/*
class SameEvent : public TObject {
        private : PseudoJet TrueParticle; vector <PseudoJet> Jets;
        public : //SameEvent() : TrueParticle(0), Jets(0) {};
        void SetTrueParticlePtPhiEtaM(Double_t TruePt, Double_t TruePhi, Double_t TrueEta, Double_t TrueM);
};
*/

class MixedEvent : public TObject {
    private :
        vector <PseudoJet> Jets;
        TClonesArray *JetsClonesArray;
        Double_t event_rho;
        Double_t event_rho_m;
    public :
        MixedEvent() : Jets(0), JetsClonesArray(0), event_rho(0), event_rho_m(0){};
        MixedEvent(vector<PseudoJet> _Jets, TClonesArray *_JetsClonesArray) : Jets(_Jets), JetsClonesArray(_JetsClonesArray), event_rho(0), event_rho_m(0) {};
        ~MixedEvent() {};

        void SetEmbedding(BackgroundParticles *BkgParticles, double JetRadius, Int_t NumEx);
        void SetJetsClonesArray();
        void ClearVector();

        TClonesArray *GetJetsClonseArray() const {return JetsClonesArray;}
        Double_t getRho() const {return event_rho;}
        Double_t getRhoM() const {return event_rho_m;}
};

inline void BackgroundParticles::SetTrueParticlePtPhiEtaM(Double_t TruePt, Double_t TruePhi, Double_t TrueEta, Double_t TrueM){
    TLorentzVector TLV_TrueParticle;
    PseudoJet TrueParticle;
    TLV_TrueParticle.SetPtEtaPhiM(TruePt, TrueEta, TruePhi, TrueM);

    Double_t TmpTruePx = TLV_TrueParticle.Px();
    Double_t TmpTruePy = TLV_TrueParticle.Py();
    Double_t TmpTruePz = TLV_TrueParticle.Pz();
    Double_t TmpTrueE =  TLV_TrueParticle.E();

    TrueParticle = PseudoJet(TmpTruePx, TmpTruePy, TmpTruePz, TmpTrueE);
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

inline void BackgroundParticles::SetBackgroundParticles(int BkgSize){
    //Particles  = SetEmbeddedParticlesVector(BkgSize, fPt, Particles);
}

inline void MixedEvent::SetJetsClonesArray(){
    JetsClonesArray = new TClonesArray("jmJetVec", 1000000);
}

inline void MixedEvent::SetEmbedding(BackgroundParticles *BkgParticles, double JetRadius, Int_t NumEx){
    
    JetDefinition   jet_def(antikt_algorithm, JetRadius);
    Double_t ghost_maxrap = 1.0;
    Double_t TrackEtaMax = 1.0;
    
    vector <int> seed;
    seed.push_back(1);
    seed.push_back(2);


    GhostedAreaSpec  area_spec(ghost_maxrap);
    AreaDefinition  area_def(active_area, area_spec);
    //ClusterSequence testJets(BkgParticles->GetBackgroundParticles(), jet_def);
    ClusterSequenceArea testJets(BkgParticles->GetBackgroundParticles(), jet_def, area_def);
    vector<PseudoJet> sort1 = sorted_by_pt(testJets.inclusive_jets());
    Int_t repeat = 1;
    Double_t ghost_area = 0.01;
    JetDefinition jet_def_bkgd(kt_algorithm, JetRadius);
    AreaDefinition area_def_bkgd(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap,repeat,ghost_area));
    Selector selector = SelectorAbsEtaMax(TrackEtaMax - JetRadius) * !SelectorNHardest(NumEx) * SelectorPtMin (0.001);
    Selector selector_for_rhoM = SelectorAbsEtaMax(TrackEtaMax - JetRadius) * !SelectorNHardest(NumEx) * SelectorPtMin (0.001);

    JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    JetMedianBackgroundEstimator bkgd_estimator_for_rhoM(selector_for_rhoM, jet_def_bkgd, area_def_bkgd); 

    vector <int> vec_seed;
    vec_seed.push_back(1);
    vec_seed.push_back(2);
    Subtractor subtractor(&bkgd_estimator);


#if FASTJET_VERSION_NUMBER >= 30100
    subtractor.set_use_rho_m(true); 
    //subtractor.set_safe_mass(true); 
#endif 

    bkgd_estimator.set_particles(BkgParticles->GetBackgroundParticles());
    bkgd_estimator_for_rhoM.set_particles(BkgParticles->GetBackgroundParticles()); 

    vector<PseudoJet> subtracted_jets = subtractor(sort1);

    //cout << subtracted_jets.size() << " " << sort1.size() << endl; 
    int TestNum = 0;
    Double_t rho = bkgd_estimator.rho();
    Double_t rho_m = bkgd_estimator_for_rhoM.rho_m();
    event_rho = rho; 
    event_rho_m = rho_m; 
    //cout << sort1.size() << " " << BkgParticles->GetBackgroundParticles().size() << endl;
    for (int tmp_particles_index = 0; tmp_particles_index < 2; tmp_particles_index ++){
        //cout << BkgParticles->GetBackgroundParticles()[tmp_particles_index].pt() << " " << BkgParticles->GetBackgroundParticles()[tmp_particles_index].eta() << BkgParticles->GetBackgroundParticles()[tmp_particles_index].phi() << endl;
    }

    assert (subtracted_jets.size() == sort1.size());

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
        Double_t JetSubM = subtracted_jets[tmp_jet_index].m(); 
        Double_t JetEta = sort1[tmp_jet_index].eta();
        PseudoJet a4 = sort1[tmp_jet_index].area_4vector();
        Double_t recon_JetArea = TMath::Sqrt(a4.px()*a4.px() + a4.py()*a4.py());
       // cout << "Area Check " << "Quadratic sum : " << recon_JetArea << " jet area : " << JetArea << " area error : " << sort1[tmp_jet_index].area_error() << endl; 
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
        JetsVec->setJetSubtractedMass(JetSubM);
        JetsVec->setConstiSize(consti.size());
        JetsVec->setAreaFourVector(a4.e(), a4.px(), a4.py(), a4.pz()); 

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
        
        //JetsVec->setConstiVectorPhi(vec_consti_phi);
        //JetsVec->setConstiVectorEta(vec_consti_eta);
        //JetsVec->setConstiVectorPt(vec_consti_pt);

        //cout << tmp_jet_index << " " << subtracted_jets[tmp_jet_index].pt() << " SubTracted Jet Mass : " << JetSubM << " JetMass : " << JetM << " IsTrueJet : " << IsTrueJet << endl; 

        if (IsTrueJet) JetsVec->setJetIndex(TrueParticleDummyIndex);
    }

    //cout << TestNum << endl;
    
}


