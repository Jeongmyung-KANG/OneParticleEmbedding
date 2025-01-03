void SetHist(TH1F *h, int MarkerStyle, int Color){
	h->SetMarkerStyle(MarkerStyle); 
	h->SetLineColor(Color);
	h->SetMarkerColor(Color); 
}
	

void AreaStudy(){


	
	TFile *fin = TFile::Open("dR.root"); 
	TH1F *SE = (TH1F*)fin->Get("hPt_SE_Inc");
	TH1F *ME = (TH1F*)fin->Get("hPt_ME");  
	TH1F *SE_Area = (TH1F*)fin->Get("hPt_SE_Inc_Area"); 
	TH1F *ME_Area = (TH1F*)fin->Get("hPt_ME_Area");  
	TH1F *SE_Miss = (TH1F*)fin->Get("hPt_SE_Miss"); 
	TH1F *SE_Matched = (TH1F*)fin->Get("hPt_SE_Matched"); 
	TH1F *Combinatorial = (TH1F*)fin->Get("hCombonatorialPt"); 
	TH1F *SE_Miss_IncludeTrue = (TH1F*)fin->Get("hPt_SE_Miss_IncludeTrue"); 

	cout << SE->GetEntries() << " " << ME->GetEntries() << " " << SE_Area->GetEntries() << " " << ME_Area->GetEntries() << endl;
	SetHist(SE, 20, 1); 
	SetHist(ME, 24, 4); 
	SetHist(SE_Area, 20, 1); 
	SetHist(ME_Area, 24, 4); 
	SetHist(SE_Miss, 21, kGreen + 2); 
	SetHist(SE_Matched, 22, kMagenta); 
	SetHist(Combinatorial, 23, 4); 
	SetHist(SE_Miss_IncludeTrue, 24, 5); 

	TCanvas *c1 = new TCanvas("c1", "c1", 500, 500, 1000, 1500) ;
	TPad *p1 = new TPad("p1", "p1", 0.0, 0.0, 1.0, 0.3); 
	TPad *p2 = new TPad("p2", "p2", 0.0, 0.3, 1.0, 1.0); 
	TH1F *hRatio = (TH1F*)Combinatorial->Clone("r1"); 
	hRatio->Divide(ME);  
	hRatio->GetYaxis()->SetRangeUser(0, 1.2); 
	c1->cd(); 
	p1->Draw();
	p1->cd(); 
	hRatio->Draw(); 
	c1->cd();	
	p2->Draw(); 
	p2->cd(); 
	SE->Draw("same");   
	ME->Draw("same"); 
	SE_Miss->Draw("same"); 
	Combinatorial->Draw("same"); 
	SE_Matched->Draw("same");  
	SE_Miss_IncludeTrue->Draw("same"); 
	TCanvas *c2 = new TCanvas("c2", "c2", 500, 1000, 1000, 1500);    
	TPad *p11 = new TPad("p11", "p11", 0.0, 0.0, 1.0, 0.3); 
	TPad *p22 = new TPad("p22", "p22", 0.0, 0.3, 1.0, 1.0); 

	c2->cd(); 
	p22->Draw(); 
	p22->cd(); 
	SE_Area->Draw(); 
	ME_Area->Draw("Same");  
	TH1F *hRatio2 = (TH1F*)SE_Area->Clone("r2"); 
	hRatio2->Divide(ME_Area);  
	hRatio2->GetYaxis()->SetRangeUser(0, 1.2); 
	c2->cd(); 
	p11->Draw(); 
	p11->cd();
	hRatio2->Draw(); 
}
