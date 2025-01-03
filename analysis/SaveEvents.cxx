#include "/Users/gangjeongmyeong/OneParticleEmbedding/include/Test.cc"
#include "/Users/gangjeongmyeong/OneParticleEmbedding/include/jmJetVec.h"

void JetFind() {
    
    gRandom->SetSeed(2);
    TTree *tree = new TTree("JetTree", "JetTree");
    gROOT->ProcessLine(".L /Users/gangjeongmyeong/OneParticleEmbedding/include/jmJetVec.h");
    
    Double_t TruePt = 20.0;
    Double_t TruePhi = TMath::Pi() / 4;
    Double_t TrueEta = 0.0;
    //Double_t TrueM = 0.13795;
    Double_t TrueM = 15; 

    MixedEvent *ME = new MixedEvent();
    ME->SetJetsClonesArray();
    
    BackgroundParticles *BkgParticle = new BackgroundParticles();

    tree->Branch("Jets", "TClonesArray", ME->GetJetsClonseArray());
    Double_t rho;
    Double_t rho_m; 
    tree->Branch("rho", &rho); 
    tree->Branch("rho_m", &rho_m); 

    BkgParticle->SetBackgroundPt("0.45");
    int NumEvents = 100000;

    for (int e = 0; e < NumEvents; e++){
        if (e % 1000 ==0) cout << e << endl;
        BkgParticle->SetBackgroundParticles(800);
        BkgParticle->SetTrueParticlePtPhiEtaM(TruePt, TruePhi, TrueEta, TrueM);
        //BkgParticle->SetTrueParticlePtPhiEtaM(50, TMath::Pi(), 0.5, TrueM);
        ME->SetEmbedding(BkgParticle, 0.4, 1); 		
        rho = ME->getRho(); 
        rho_m = ME->getRhoM(); 
        tree->Fill();
        ME->GetJetsClonseArray()->Clear();
        BkgParticle->ClearVector();
        ME->ClearVector();
    }

    TFile *fout = new TFile("SE_0.root", "RECREATE");
    fout->cd();
    tree->Write();
    fout->Close();
}

int main(){
        JetFind();
        return 0;
}
