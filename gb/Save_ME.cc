#include "include/GenEvents.h" 
void JetFind() {
	
	gRandom->SetSeed(1);
	TFile *fout = new TFile("ME.root", "RECREATE"); 
	TTree *tree = new TTree("JetTree", "JetTree"); 

    gROOT->ProcessLine(".L ./jmJetVec.h");	
	
	MixedEvent *ME = new MixedEvent();
	ME->SetJetsClonesArray(); 
	
	BackgroundParticles *BkgParticle = new BackgroundParticles();    
	TH1F *hTrackEta = new TH1F("hTrackEta", "hTrackEta", 200, -2, 2); 
	tree->Branch("Jets", "TClonesArray", ME->GetJetsClonseArray());

	BkgParticle->SetBackgroundPt("0.45"); 
	int NumEvents = 10000;

	for (int e = 0; e < NumEvents; e++){  
		//cout << e << endl;

		BkgParticle->SetBackgroundParticles(1); 
		ME->SetEmbedding(BkgParticle,  0.4, 1); 
		/*
		cout << BkgParticle->GetBackgroundParticles().size() << endl;
		for (int t = 0; t < BkgParticle->GetBackgroundParticles().size(); t ++){ 
		if (t < 10) cout << e << " " << t << " " << BkgParticle->GetBackgroundParticles()[t].pt() << " " << BkgParticle->GetBackgroundParticles()[t].eta() << BkgParticle->GetBackgroundParticles()[t].eta() << endl; 
		}

		*/ 
		tree->Fill(); 
		ME->GetJetsClonseArray()->Clear(); 
		BkgParticle->ClearVector();
		ME->ClearVector(); 
	}
	
	fout->cd(); 
	tree->Write(); 
	hTrackEta->Write(); 
	fout->Close();
}

int main(){
		JetFind();
		return 0;
}
