//#include "Test.cc"


void JetFind() {
	
	gRandom->SetSeed(1);
	TFile *fout = new TFile("SE.root", "RECREATE"); 
	TTree *tree = new TTree("JetTree", "JetTree"); 

    gROOT->ProcessLine(".L ./jmJetVec.h");	
	
	Double_t TruePt = 20.0; 
	Double_t TruePhi = TMath::Pi() / 4; 
	Double_t TrueEta = 0.0; 
	Double_t TrueM = 0.13579; 
	
	MixedEvent *ME = new MixedEvent();
	ME->SetJetsClonesArray(); 
	
	BackgroundParticles *BkgParticle = new BackgroundParticles();    

	tree->Branch("Jets", "TClonesArray", ME->GetJetsClonseArray());

	BkgParticle->SetBackgroundPt("0.45"); 
	int NumEvents = 10000;

	for (int e = 0; e < NumEvents; e++){  
		//cout << e << endl;
		BkgParticle->SetBackgroundParticles(1); 
		BkgParticle->SetTrueParticlePtPhiEtaM(TruePt, TruePhi, TrueEta, TrueM); 
		//BkgParticle->SetTrueParticlePtPhiEtaM(50, TMath::Pi(), 0.5, TrueM); 
		ME->SetEmbedding(BkgParticle,  0.4,  0); 

		tree->Fill(); 
		ME->GetJetsClonseArray()->Clear(); 
		BkgParticle->ClearVector();
		ME->ClearVector(); 
	}
	
	fout->cd(); 
	tree->Write(); 
	fout->Close();
}

int main(){
		JetFind();
		return 0;
}
