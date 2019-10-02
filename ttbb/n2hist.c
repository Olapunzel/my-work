void n2hist(){

	TFile *file = new TFile("/home/juhee5819/public/data/tight_default_0911/4f_63.root", "read");
	TTree *tree = (TTree*) file->Get("dnn_input");
//Variables	
	double addbjet1_pt;
	tree->SetBranchAddress("addbjet1_pt", &addbjet1_pt);
	double addbjet1_eta;
	tree->SetBranchAddress("addbjet1_eta", &addbjet1_eta);
	double addbjet1_phi;
	tree->SetBranchAddress("addbjet1_phi", &addbjet1_phi);
	double addbjet1_e;
	tree->SetBranchAddress("addbjet1_e", &addbjet1_e);
	double addbjet2_pt;
	tree->SetBranchAddress("addbjet2_pt", &addbjet2_pt);
	double addbjet2_eta;
	tree->SetBranchAddress("addbjet2_eta", &addbjet2_eta);
	double addbjet2_phi;
	tree->SetBranchAddress("addbjet2_phi", &addbjet2_phi);
	double addbjet2_e;
	tree->SetBranchAddress("addbjet2_e", &addbjet2_e);
	double pt1;
	
	UInt_t signal;
	tree->SetBranchAddress("signal", &signal);
	tree->SetBranchAddress("pt1", &pt1);
	double eta1;
	tree->SetBranchAddress("eta1", &eta1);
	double phi1;
	tree->SetBranchAddress("phi1", &phi1);
	double e1;
	tree->SetBranchAddress("e1", &e1);
	double pt2;
	tree->SetBranchAddress("pt2", &pt2);
	double eta2;
	tree->SetBranchAddress("eta2", &eta2);
	double phi2;
	tree->SetBranchAddress("phi2", &phi2);
	double e2;
	tree->SetBranchAddress("e2", &e2);

//Book histogrmas
 	TH1 *hist_addbjet1_pt = new TH1F("h_addbjet1_pt", "addbjet_pt", 40, 0, 260);
 	TH1 *hist_addbjet1_eta = new TH1F("h_addbjet1_eta", "addbjet_eta", 80, -4, 4);
 	TH1 *hist_addbjet1_phi = new TH1F("h_addbjet1_phi", "addbjet_phi", 80, -4, 4);
 	TH1 *hist_addbjet1_e = new TH1F("h_addbjet1_en", "addbjet_energy", 40, 0, 260);
 	TH1 *hist_addbjet2_pt = new TH1F("h_addbjet2_pt", "addbjet_pt", 40, 0, 260);
 	TH1 *hist_addbjet2_eta = new TH1F("h_addbjet2_eta", "addbjet_eta", 80, -4, 4);
 	TH1 *hist_addbjet2_phi = new TH1F("h_addbjet2_phi", "addbjet_phi", 80, -4, 4);
 	TH1 *hist_addbjet2_e = new TH1F("h_addbjet2_en", "addbjet_energy", 40, 0, 260);

	TH1 *hist_s1_pt1 = new TH1F("h_s1_pt1", "signal_PT1", 40, 0, 260);
 	TH1 *hist_s1_eta1 = new TH1F("h_s1_eta1", "signal_Eta1", 80, -4, 4);
 	TH1 *hist_s1_phi1 = new TH1F("h_s1_phi1", "signal_Phi1", 80, -4, 4);
 	TH1 *hist_s1_e1 = new TH1F("h_s1_en1", "signal_Energy1",40, 0, 300); 
 	TH1 *hist_s1_pt2 = new TH1F("h_s1_pt2", "signal_PT2", 40, 0, 200);
 	TH1 *hist_s1_eta2 = new TH1F("h_s1_eta2", "signal_Eta2", 80, -4, 4);
 	TH1 *hist_s1_phi2 = new TH1F("h_s1_phi2", "signal_Phi2", 80, -4, 4);
 	TH1 *hist_s1_e2 = new TH1F("h_s1_en2", "signal_Energy2",40, 0, 240);

 	TH1 *hist_s0_pt1 = new TH1F("h_s0_pt1", "signal_PT1", 40, 0, 260);
 	TH1 *hist_s0_eta1 = new TH1F("h_s0_eta1", "signal_Eta1", 80, -4, 4);
 	TH1 *hist_s0_phi1 = new TH1F("h_s0_phi1", "signal_Phi1", 80, -4, 4);
 	TH1 *hist_s0_e1 = new TH1F("h_s0_en1", "signal_Energy1",40, 0, 300); 
 	TH1 *hist_s0_pt2 = new TH1F("h_s0_pt2", "signal_PT2", 40, 0, 200);
 	TH1 *hist_s0_eta2 = new TH1F("h_s0_eta2", "signal_Eta2", 80, -4, 4);
 	TH1 *hist_s0_phi2 = new TH1F("h_s0_phi2", "signal_Phi2", 80, -4, 4);
 	TH1 *hist_s0_e2 = new TH1F("h_s0_en2", "signal_Energy2",40, 0, 240);

	for(int irow=0;irow<tree->GetEntries();++irow){
		tree->GetEntry(irow);
		hist_addbjet1_pt->Fill(addbjet1_pt);
		hist_addbjet1_eta->Fill(addbjet1_eta);
		hist_addbjet1_phi->Fill(addbjet1_phi);
		hist_addbjet1_e->Fill(addbjet1_e);
		hist_addbjet2_pt->Fill(addbjet2_pt);
		hist_addbjet2_eta->Fill(addbjet2_eta);
		hist_addbjet2_phi->Fill(addbjet2_phi);
		hist_addbjet2_e->Fill(addbjet2_e);
		if(signal==1){
			hist_s1_pt1->Fill(pt1);
			hist_s1_eta1->Fill(eta1);
			hist_s1_phi1->Fill(phi1);
			hist_s1_e1->Fill(e1);
			hist_s1_pt2->Fill(pt2);
			hist_s1_eta2->Fill(eta2);
			hist_s1_phi2->Fill(phi2);
			hist_s1_e2->Fill(e2);
		}
		else{
			hist_s0_pt1->Fill(pt1);
			hist_s0_eta1->Fill(eta1);
			hist_s0_phi1->Fill(phi1);
			hist_s0_e1->Fill(e1);
			hist_s0_pt2->Fill(pt2);
			hist_s0_eta2->Fill(eta2);
			hist_s0_phi2->Fill(phi2);
			hist_s0_e2->Fill(e2);
		}
	}
	TList *list1 = new TList;
	list1->Add(hist_addbjet1_pt);
	list1->Add(hist_addbjet1_eta);
	list1->Add(hist_addbjet1_phi);
	list1->Add(hist_addbjet1_e);
	list1->Add(hist_addbjet2_pt);
	list1->Add(hist_addbjet2_eta);
	list1->Add(hist_addbjet2_phi);
	list1->Add(hist_addbjet2_e);

	TList *list2 = new TList;	
	list2->Add(hist_s1_pt1);
	list2->Add(hist_s1_eta1);
	list2->Add(hist_s1_phi1);
	list2->Add(hist_s1_e1);
	list2->Add(hist_s1_pt2);
	list2->Add(hist_s1_eta2);
	list2->Add(hist_s1_phi2);
	list2->Add(hist_s1_e2);
	list2->Add(hist_s0_pt1);
	list2->Add(hist_s0_eta1);
	list2->Add(hist_s0_phi1);
	list2->Add(hist_s0_e1);
	list2->Add(hist_s0_pt2);
	list2->Add(hist_s0_eta2);
	list2->Add(hist_s0_phi2);
	list2->Add(hist_s0_e2);

	for(auto&& KeyAddb1 : *list1){
		auto key1 = (TKey*)KeyAddb1;
		const char *name1 = key1->GetName();

		if(strstr(name1, "1")!=NULL){
			for(auto&& KeyAddb2 : *list1){
				auto key2 = (TKey*)KeyAddb2;
				const char *name2 = key2->GetName();
				string str = name1;
				const char *sub = str.substr(10).c_str();
//				cout << sub << endl;
				if(strstr(name2, "2")!=NULL && strstr(name2, sub)!=NULL){
					TH1F *hist1 = (TH1F*) key1;
					TH1F *hist2 = (TH1F*) key2;
				
					TCanvas *c = new TCanvas;
					gStyle->SetOptStat(0);

					hist1->SetYTitle("Normalized Entries");
					hist1->SetLineWidth(2);
				//	hist1->SetFillColor(4);
					hist1->SetLineColor(6);
					hist1->SetLineWidth(2);
					hist2->SetLineColor(8);
					hist2->SetLineWidth(2);
					hist1->Scale(1/hist1->Integral());
					hist2->Scale(1/hist2->Integral());
					
					Double_t maxi = max(hist1->GetMaximum(), hist2->GetMaximum());
					hist1->SetMaximum(maxi*1.2);
					hist1->Draw("hist");
					hist2->Draw("hist Same");

    				auto leg=new TLegend(.6,.8,.9,.9,"");
	    			leg->AddEntry(hist1, "4f_63_addbjet1");
   				 	leg->AddEntry(hist2, "4f_63_addbjet2");
  			  		leg->SetTextSize(.04);
		    		leg->Draw();
	//				c->Print(Form("4f_63_addbjet%s.pdf", str.substr(10).c_str()));
				}
			}
		}
	}

	for(auto&& KeyAsObj1 : *list2){
		auto key1 = (TKey*)KeyAsObj1;
		const char *name1 = key1->GetName();

		if(strstr(name1, "s1")!=NULL){
			for(auto&& KeyAsObj2 : *list2){
				auto key2 = (TKey*)KeyAsObj2;
				const char *name2 = key2->GetName();
				string str = name1;
				const char *sub = str.substr(4).c_str();
				
				if(strstr(name2, "s0")!=NULL && strstr(name2, sub)) {

					TH1F *hist1 = (TH1F*) key1;
					TH1F *hist2 = (TH1F*) key2;
				
					TCanvas *c = new TCanvas;
		
					gStyle->SetOptStat(0);

					hist1->SetYTitle("Normalized Entries");
					hist1->SetLineWidth(2);
					hist1->SetFillColor(kRed-10);
					hist1->SetLineColor(kRed+1);
					hist1->SetLineWidth(2);
					hist2->SetLineColor(kBlue+1);
					hist2->SetLineWidth(2);
					hist1->Scale(1/hist1->Integral());
					hist2->Scale(1/hist2->Integral());
					
					Double_t maxi = max(hist1->GetMaximum(), hist2->GetMaximum());
					hist1->SetMaximum(maxi*1.2);
					hist1->Draw("hist");
					hist2->Draw("hist Same");

    				auto leg=new TLegend(.725,.8,.9,.9,"");
	    			leg->AddEntry(hist1, "4f_63_s1");
   				 	leg->AddEntry(hist2, "4f_63_s0");
  			  		leg->SetTextSize(.04);
		    		leg->Draw();
			//		c->Print(Form("4f_63_sig%s.pdf", str.substr(4).c_str()));


					}
				else continue;
		}
	
			}
		
		else continue;
}
}
