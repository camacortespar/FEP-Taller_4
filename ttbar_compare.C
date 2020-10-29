// TTBar selection with data to compare with MC sample.

void ttbar_compare(){

  // Let's open the TFile
  
  TFile *file = TFile::Open("Data_8TeV.root");
  TTree *tree = (TTree*) file->Get("mini");
  tree->Print();

  // Extract and declare variables

  Bool_t e_trig;    //Cut 2
  Bool_t mu_trig;
  Bool_t good_vtx;  //Cut 1
  UInt_t lep_n;     //Cut 3 - Leptons (good)
  UInt_t jet_n;     //Cut 4,5,6 - Jets (good), b-jets
  Float_t MET;      //Cut 7,8 - About MET
  Float_t MET_phi;

  Float_t lep_pt[10];
  Float_t lep_eta[10];
  Float_t lep_phi[10];
  Float_t lep_E[10];
  Int_t lep_type[10];
  Float_t lep_ptcone30[10];
  Float_t lep_etcone20[10];

  Float_t jet_pt[10];
  Float_t jet_eta[10];
  Float_t jet_phi[10];
  Float_t jet_jvf[10];
  Float_t jet_MV1[10];
  
  tree->SetBranchAddress("hasGoodVertex", &good_vtx);
  tree->SetBranchAddress("trigE", &e_trig);
  tree->SetBranchAddress("trigM", &mu_trig);

  //Leptons
  tree->SetBranchAddress("lep_n", &lep_n);
  tree->SetBranchAddress("lep_pt", &lep_pt);
  tree->SetBranchAddress("lep_eta", &lep_eta);
  tree->SetBranchAddress("lep_phi", &lep_phi);
  tree->SetBranchAddress("lep_E", &lep_E);
  tree->SetBranchAddress("lep_type", &lep_type);
  tree->SetBranchAddress("lep_ptcone30", &lep_ptcone30);
  tree->SetBranchAddress("lep_etcone20", &lep_etcone20);

  //Jets
  tree->SetBranchAddress("jet_n", &jet_n);
  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_eta", &jet_eta);
  tree->SetBranchAddress("jet_jvf", &jet_jvf);
  tree->SetBranchAddress("jet_MV1", &jet_MV1);

  //MET
  tree->SetBranchAddress("met_et", &MET);
  tree->SetBranchAddress("met_phi", &MET_phi);

  
  // Create a Canvas 800x600

  TCanvas *canvas = new TCanvas("Canvas","",800,600);

  // Create a histograms

  TH1F *cutflow = new TH1F("Cutflow","Cutflow; Cut; Events",10,0,10);
  TH1F *hist_lep_pt = new TH1F("Leptons p_{T}","Leptons p_{T}; p_{T} (MeV); Events",20,0,500000);
  TH1F *hist_lep_pt30 = new TH1F("Track isolation","Track isolation; Lep_ptcone30/Lep_pt; Events",20,0,10);
  TH1F *hist_lep_et20 = new TH1F("Calorimeter isolation","Calorimeter isolation; Lep_etcone20/Lep_pt; Events",20,-2,2);
  TH1F *hist_lep_eta = new TH1F("Leptons #eta","Leptons #eta; #eta; Events",20,-5,5);
  TH1F *hist_njets = new TH1F("n-jets","Number of jets; Jet multiplicity; Events",10,0,10);
  TH1F *hist_jets_pt = new TH1F("Jets p_{T}","Jets p_{T}; p_{T} (MeV); Events",20,0,500000);
  TH1F *hist_jets_eta = new TH1F("Jets #eta","Jets #eta; #eta; Events",20,-5,5);
  TH1F *hist_jets_JVF = new TH1F("Jets JVF","JVF; Jets JVF; Events",10,-1,1);
  TH1F *hist_jets_MV1 = new TH1F("Jets MV1","MV1; Jets MV1; Events",10,-1,1);
  TH1F *hist_bjets = new TH1F("b-jets","Number of bjets; Jet multiplicity; Events",10,0,10);
  TH1F *hist_MET = new TH1F("MET","MET; MET; Events",200,0,100000);
  TH1F *hist_MTW = new TH1F("MTW","MTW; MTW; Events",200,0,100000);
  
  // Loop and fills histograms

  int nentries, nbytes, i;
  nentries = (Int_t)tree->GetEntries();

  int cut1 = 0;
  int cut2 = 0;
  int cut3 = 0;
  int cut4 = 0;
  int cut5 = 0;
  int cut6 = 0;
  int cut7 = 0;
  int cut8 = 0;

  for (i = 0; i < nentries; i++)
    {
      nbytes = tree->GetEntry(i);
    
      //First cut: Good vertex
      if(!good_vtx) continue;
      cut1++;
      cutflow->Fill(1);

      //Second cut: Trigger
      if(!e_trig && !mu_trig) continue;
      cut2++;
      cutflow->Fill(2);

      //Preselection of good leptons
      int n_mu=0;
      int n_el=0;
      int n_lep=0;
      int good_lep=0;
      
      //Loop over leptons
      for(unsigned int i=0; i<lep_n; i++){
	//Cut
	if( lep_pt[i] < 25000.) continue;
	if( lep_ptcone30[i]/lep_pt[i] > 0.15 ) continue;          //Track isolation
	if( lep_etcone20[i]/lep_pt[i] > 0.15 ) continue;          //Calorimeter isolation
	if( lep_type [i]==13 && TMath::Abs(lep_eta[i]) < 2.5 )    //Muons
	  {n_mu++;
	   good_lep=i;}
	if( lep_type [i]==11 && TMath::Abs(lep_eta[i]) < 2.47 ){
	  if( TMath::Abs(lep_eta[i]) < 1.37 || TMath::Abs(lep_eta[i]) > 1.52 )      //Electrons
	    {n_el++;
	     good_lep=i;}
          }
      }
     

      n_lep=n_mu+n_el;
      
      //Select events with only 1 good lepton and fill the cutflow histogram
      //Third cut (one good lepton):
      if(n_lep!=1) continue;
      cutflow->Fill(3);
      cut3++;

      //Go for cut 4
      int n_jets=0;
      int n_bjets=0;

      //Fourth cut: At least 4 jets
      if(jet_n<4) continue;
      cutflow->Fill(4);
      cut4++;

      //Number of jets distribution
      hist_njets->Fill(jet_n);

      for(unsigned int j=0; j<jet_n; j++){
	//Cut
	if(jet_pt[j] < 25000.) continue;
	if(TMath::Abs(jet_eta[j]) > 2.5) continue;
	if(jet_pt[j] < 50000. && TMath::Abs(jet_eta[j]) < 2.4 && jet_jvf[j] < 0.5) continue;
	  n_jets++;
	if(TMath::Abs(jet_MV1[j]) < 0.7892) continue;
	  n_bjets++;
	}

      //Fifth cut: At least 4 good jets
      if(n_jets<4) continue;
      cutflow->Fill(5);
      cut5++;

      //Sixth cut: at least one b-jet
      if(n_bjets<2) continue;
      cutflow->Fill(6);
      cut6++;

      //Number of bjets distribution
      hist_bjets->Fill(n_bjets);

      //Seventh cut: MET > 30 GeV
      if(MET<30000.) continue;
      cutflow->Fill(7);
      cut7++;

      //MET Histogram
      hist_MET->Fill(MET);

      //TLorentz definitions
      TLorentzVector Lepton = TLorentzVector();
      TLorentzVector MeT = TLorentzVector();

      Lepton.SetPtEtaPhiE(lep_pt[good_lep],lep_eta[good_lep],lep_phi[good_lep],lep_E[good_lep]);
      MeT.SetPtEtaPhiE(MET,0,MET_phi,MET);

      //mTW
      float mTW = sqrt(2*Lepton.Pt()*MeT.Et()*(1-cos(Lepton.DeltaPhi(MeT))));
     
      //Eight cut: mTW > 30 GeV
      if(mTW>30000.)
	{cutflow->Fill(8);
	  cut8++;}

      //mTW histogram
      hist_MTW->Fill(mTW);

      for(unsigned int i=0; i<lep_n; i++){
	//Leptons histograms
	hist_lep_pt->Fill(lep_pt[i]);
	hist_lep_pt30->Fill(lep_ptcone30[i]/lep_pt[i]);
	hist_lep_et20->Fill(lep_etcone20[i]/lep_pt[i]);
	hist_lep_eta->Fill(lep_eta[i]);
      }
      for(unsigned int j=0; j<jet_n; j++){
	//Jets histograms
	hist_jets_pt->Fill(jet_pt[j]);
	hist_jets_eta->Fill(jet_eta[j]);
	hist_jets_JVF->Fill(jet_jvf[j]);
	hist_jets_MV1->Fill(jet_MV1[j]);
      }
    }


  std::cout << "Done!" << std::endl;
  std::cout << "All events:" << nentries << std::endl;
  std::cout << "Cut1:" << cut1 << std::endl;
  std::cout << "Cut2:" << cut2 << std::endl;
  std::cout << "Cut3:" << cut3 << std::endl;
  std::cout << "Cut4:" << cut4 << std::endl;
  std::cout << "Cut5:" << cut5 << std::endl;
  std::cout << "Cut6:" << cut6 << std::endl;
  std::cout << "Cut7:" << cut7 << std::endl;
  std::cout << "Cut8:" << cut8 << std::endl;

  // Let's plot, just choose one
  
  //cutflow->Draw("");
  //hist_lep_pt->Draw();
  //hist_lep_pt30->Draw();
  //hist_lep_et20->Draw();
  //hist_lep_eta->Draw();
  //hist_njets->Draw();
  hist_jets_pt->Draw();
  //hist_jets_eta->Draw();
  //hist_jets_JVF->Draw();
  //hist_jets_MV1->Draw();
  //hist_bjets->Draw();
  //hist_MET->Draw();
  //hist_MTW->Draw();
  //canvas->SetLogy();
  canvas->Draw();

}
