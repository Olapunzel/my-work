/*
root -l -b -q ana.C'("step_1.root", "step_1_plots.root")'
*/

//------------------------------------------------------------------------------
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include <vector>

GenParticle* MotherParticle(const GenParticle* p, const TClonesArray* branchParticle, int i = 0){

  if( p == 0) return 0;

  int PID = p->PID;
  int index[2];
  index[0] = p->M1;
  index[1] = p->M2;

  GenParticle * m = 0;

  if( index[i] > 0){
    m = (GenParticle *) branchParticle->At( index[i] );
    double mPID = m->PID;
    cout << "particle " << PID << " (" << p->P4().Pt() << ") " << " is from " << " mother (M" << i << ") " << mPID << " ( " << m->P4().Pt() << " ) " << endl;
  }else{
    cout << "No mother particle " << "M" << i << " from  " << PID << "  !!!" << endl;
  }

  return m; 

}

GenParticle* DaughterParticle(const GenParticle* p, const TClonesArray* branchParticle, int i = 0){

  if( p == 0) return 0;

  int PID = p->PID;
  int index[2];
  index[i] = p->D1;
  index[i] = p->D2;

  GenParticle * m = 0;

  if( index[i] > 0){
    m = (GenParticle *) branchParticle->At( index[i] );
    double mPID = m->PID;
    cout << "particle " << PID << " (" << p->P4().Pt() << ") " << " decays to " << " particle (D" << i << ") " << mPID << " ( " << m->P4().Pt() << " ) " << endl;
  }else{
    cout << "No decaying particle " << "D" << i << " from  " << PID << "  !!!" << endl;
  }

  return m;

}

bool isFromTop(const GenParticle* p, const TClonesArray* branchParticle){

  bool output = false;
  bool debug = false;
  bool printout = false;
  double m1 = p->M1;
  double m2 = p->M2;
  int PID = p->PID;
  int mPID1 = -1;
  int mPID2 = -1;

  while( m1 >= 0 || m2 >= 0){
    if( m1 >= 0 && m2 < 0){
      GenParticle *mother = (GenParticle *) branchParticle->At(m1);
      m1 = mother->M1;
      m2 = mother->M2;
      mPID1 = mother->PID;
      if(debug) cout << "case1 : mother = " << mother->PID << endl;
      if( abs(mother->PID) == 6 ) { 
        output = true;
        break;
      }
    }else if( m1 < 0 && m2 >= 0){
      GenParticle *mother = (GenParticle *) branchParticle->At(m2);
      m1 = mother->M1;
      m2 = mother->M2;
      mPID2 = mother->PID;
      if(debug) cout << "case2 : mother = " << mother->PID << endl;
      if( abs(mother->PID) == 6 ) {
        output = true;
        break;
      }
    }else if( m1 == m2 ){
      GenParticle *mother = (GenParticle *) branchParticle->At(m1);
      m1 = mother->M1;
      m2 = mother->M2;
      mPID1 = mother->PID;
      mPID2 = mother->PID;
      if(debug) cout << "case 3 : mother = " << mother->PID << endl;
      if( abs(mother->PID) == 6 ) {
        output = true;
        break;
      }
    }else{
      GenParticle *mother1 = (GenParticle *) branchParticle->At(m1);
      GenParticle *mother2 = (GenParticle *) branchParticle->At(m2);

      mPID1 = mother1->PID;
      mPID2 = mother2->PID;
      if(debug) cout << "case4 : mother = " << mother1->PID << " " << mother2->PID << endl;

      double m11 = mother1->M1;
      double m12 = mother1->M2;
      double m21 = mother2->M1;
      double m22 = mother2->M2;
      double mo1_m1_PID = -1;
      double mo1_m2_PID = -1;
      double mo2_m1_PID = -1;
      double mo2_m2_PID = -1;
      if( m11 >= 0 ) {
        GenParticle * mo1_m1 = (GenParticle*) branchParticle->At(m11);
        mo1_m1_PID = mo1_m1->PID;
      }
      if( m12 >= 0 ) {
        GenParticle * mo1_m2 = (GenParticle*) branchParticle->At(m12);
        mo1_m2_PID =  mo1_m2->PID;
      }
      if( m21 >= 0 ) {
        GenParticle * mo2_m1 = (GenParticle*) branchParticle->At(m21);
        mo2_m1_PID = mo2_m1->PID;
      }
      if( m22 >= 0 ) {
        GenParticle * mo2_m2 = (GenParticle*) branchParticle->At(m22);
        mo2_m2_PID = mo2_m2->PID;
      }

      bool fromtop = abs(mother1->PID) == 6 || abs(mother2->PID) == 6; 
      bool fromtop2 = abs(mo1_m1_PID) == 6 || abs(mo1_m2_PID) == 6;
      bool fromtop3 = abs(mo2_m1_PID) == 6 || abs(mo2_m2_PID) == 6;
      if( fromtop ) {
        mPID1 = mother1->PID;
        mPID2 = mother2->PID;
        output = true;
        break;
      }
      if( fromtop2 ) {
        mPID1 = mo1_m1_PID;
        mPID2 = mo1_m2_PID;
        output = true;
        break;
      }
      if( fromtop3 ) {
        mPID1 = mo2_m1_PID;
        mPID2 = mo2_m2_PID;
        output = true;
        break;
      }

  if( (m11 >= 0 && m12 >= 0) && (m21 < 0 && m22 < 0) ) {
        m1 = m11;
        m2 = m12;
      }else if( (m11 < 0 && m12 < 0) && (m21 >= 0 && m22 >= 0) ){
        m1 = m21;
        m2 = m22;
      }else if( (m11 >= 0 && m12 < 0) && (m21 >= 0 && m22 < 0) ){
        m1 = m11;
        m2 = m21;
      }else if( (m11 < 0 && m12 >= 0) && (m21 < 0 && m22 >= 0) ){
        m1 = m12;
        m2 = m22;
      }else if( (m11 >= 0 && m12 < 0) && (m21 < 0 && m22 >= 0) ){
        m1 = m11;
        m2 = m22;
      }else if( (m11 < 0 && m12 >= 0) && (m21 >= 0 && m22 < 0) ){ 
        m1 = m12;
        m2 = m21; 
      }else if( (m11 < 0 && m12 < 0) && (m21 < 0 && m22 < 0)  ){
        break;
      }else{
        if( abs(mo1_m1_PID) == 5 || abs(mo1_m2_PID) == 5){
          m1 = m11;
          m2 = m12;
        }else if( abs(mo2_m1_PID) == 5 || abs(mo2_m2_PID) == 5)  {
          m1 = m21;
          m2 = m22;
        }else{
          if(debug){
            cout << "===== DEBUG =====" << endl;
            cout << m11 << " " << m12 << " " << m21 << " " << m22 << endl;
            cout << mo1_m1_PID << " " << mo1_m2_PID << " " << mo2_m1_PID << " " << mo2_m2_PID << endl;
            cout << "===== BREAK =====" << endl;
          }
          break;
        }
      }

    }
  }
  if(printout) cout << "Ancestor PID = " << mPID1 << " " << mPID2 << " test particle ID = " << PID << endl;

  return output;
}

void ana_dilep(const char *inputFile, const char *outputFile, int jcut, int bcut)
{
 //jcut is for the number of jets cut and bcut is for the number of b-jets cut
 
 gSystem->Load("libDelphes");

 // Create chain of root trees
 TChain chain("Delphes");
 chain.Add(inputFile);

 TString filename = inputFile;
 
 // Switch for single lepton or dilepton
 bool isdilepton = false;
 if( filename.Contains("di") == true ){
   isdilepton = true;
   cout<<"Dilepton Channel"<<endl;
 }
 else cout<<"Single Lepton Channel"<<endl;

 //Output
 TFile *fout = TFile::Open(outputFile,"RECREATE");
 fout->cd();

 //DNN variables
 int signal, event;
 int njets, nbjets;
 double addbjet1_pt, addbjet1_eta, addbjet1_phi, addbjet1_e;
 double addbjet2_pt, addbjet2_eta, addbjet2_phi, addbjet2_e;
 double pt1, eta1, phi1, e1;
 double pt2, eta2, phi2, e2;
 double bbdR, bbdEta, bbdPhi, bbPt, bbEta, bbPhi, bbMass, bbHt, bbMt;
 double nubbdR, nubbdEta, nubbdPhi, nubbPt, nubbEta, nubbPhi, nubbMass, nubbHt, nubbMt;
 double lbbdR, lbbdEta, lbbdPhi, lbbPt, lbbEta, lbbPhi, lbbMass, lbbHt, lbbMt;
 double Wjb1dR, Wjb1dEta, Wjb1dPhi, Wjb1Pt, Wjb1Eta, Wjb1Phi, Wjb1Mass, Wjb1Ht, Wjb1Mt;
 double Wjb2dR, Wjb2dEta, Wjb2dPhi, Wjb2Pt, Wjb2Eta, Wjb2Phi, Wjb2Mass, Wjb2Ht, Wjb2Mt;
 double nub1dR, nub1dEta, nub1dPhi, nub1Pt, nub1Eta, nub1Phi, nub1Mass, nub1Ht, nub1Mt;
 double nub2dR, nub2dEta, nub2dPhi, nub2Pt, nub2Eta, nub2Phi, nub2Mass, nub2Ht, nub2Mt;
 double lb1dR, lb1dEta, lb1dPhi, lb1Pt, lb1Eta, lb1Phi, lb1Mass, lb1Ht, lb1Mt;
 double lb2dR, lb2dEta, lb2dPhi, lb2Pt, lb2Eta, lb2Phi, lb2Mass, lb2Ht, lb2Mt;
 
 //Tree variables ( minimum deltaR )
 double bjet1_pt;
 double bjet1_eta;
 double bjet1_phi;
 double bjet1_e;
 double bjet2_pt;
 double bjet2_eta;
 double bjet2_phi;
 double bjet2_e;
 double bjet3_pt;
 double bjet3_eta;
 double bjet3_phi;
 double bjet3_e;

 unsigned short nMuon, nElectron, nLepton;
 double Jet_pt, Jet_eta, Jet_phi, Jet_e;
 double Electron1_pt, Electron1_eta, Electron1_phi, Electron1_e;
 double Electron2_pt, Electron2_eta, Electron2_phi, Electron2_e;
 double Muon1_pt, Muon1_eta, Muon1_phi, Muon1_e;
 double Muon2_pt, Muon2_eta, Muon2_phi, Muon2_e;
 double MET_px, MET_py;
 int MATCHED;
 // Selected Events (Cut Flow)
 int s1 = 0;
 int s2 = 0;
 int s3 = 0;

 //Tree for Deep learning input 
 TTree * dnn_tree = new TTree( "dnn_input", "tree for dnn");
 dnn_tree->Branch("signal",&signal,"signal/i");
 dnn_tree->Branch("event",&event,"event/i");
 dnn_tree->Branch("njets",&njets,"njets/i");
 dnn_tree->Branch("nbjets",&nbjets,"nbjets/i"); 
 dnn_tree->Branch("addbjet1_pt",&addbjet1_pt,"addbjet1_pt/d");
 dnn_tree->Branch("addbjet1_eta",&addbjet1_eta,"addbjet1_eta/d");
 dnn_tree->Branch("addbjet1_phi",&addbjet1_phi,"addbjet1_phi/d");
 dnn_tree->Branch("addbjet1_e",&addbjet1_e,"addbjet1_e/d");
 dnn_tree->Branch("addbjet2_pt",&addbjet2_pt,"addbjet2_pt/d");
 dnn_tree->Branch("addbjet2_eta",&addbjet2_eta,"addbjet2_eta/d");
 dnn_tree->Branch("addbjet2_phi",&addbjet2_phi,"addbjet2_phi/d");
 dnn_tree->Branch("addbjet2_e",&addbjet2_e,"addbjet2_e/d");
 
 dnn_tree->Branch("pt1",&pt1,"pt1/d");
 dnn_tree->Branch("eta1",&eta1,"eta1/d");
 dnn_tree->Branch("phi1",&phi1,"phi1/d");
 dnn_tree->Branch("e1",&e1,"e1/d");
 dnn_tree->Branch("pt2",&pt2,"pt2/d");
 dnn_tree->Branch("eta2",&eta2,"eta2/d");
 dnn_tree->Branch("phi2",&phi2,"phi2/d");
 dnn_tree->Branch("e2",&e2,"e2/d");

 dnn_tree->Branch("bbdR",&bbdR,"bbdR/d");
 dnn_tree->Branch("bbdEta",&bbdEta,"bbdEta/d");
 dnn_tree->Branch("bbdPhi",&bbdPhi,"bbdPhi/d");
 dnn_tree->Branch("bbPt",&bbPt,"bbPt/d");
 dnn_tree->Branch("bbEta",&bbEta,"bbEta/d");
 dnn_tree->Branch("bbPhi",&bbPhi,"bbPhi/d");
 dnn_tree->Branch("bbMass",&bbMass,"bbMass/d");
 dnn_tree->Branch("bbHt",&bbHt,"bbHt/d");
 dnn_tree->Branch("bbMt",&bbMt,"bbMt/d");

 dnn_tree->Branch("nubbdR",&nubbdR,"nubbdR/d");
 dnn_tree->Branch("nubbdEta",&nubbdEta,"nubbdEta/d");
 dnn_tree->Branch("nubbdPhi",&nubbdPhi,"nubbdPhi/d");
 dnn_tree->Branch("nubbPt",&nubbPt,"nubbPt/d");
 dnn_tree->Branch("nubbEta",&nubbEta,"nubbEta/d");
 dnn_tree->Branch("nubbPhi",&nubbPhi,"nubbPhi/d");
 dnn_tree->Branch("nubbMass",&nubbMass,"nubbMass/d");
 dnn_tree->Branch("nubbHt",&nubbHt,"nubbHt/d");
 dnn_tree->Branch("nubbMt",&nubbMt,"nubbMt/d");

 dnn_tree->Branch("nub1dR",&nub1dR,"nub1dR/d");
 dnn_tree->Branch("nub1dEta",&nub1dEta,"nub1dEta/d");
 dnn_tree->Branch("nub1dPhi",&nub1dPhi,"nub1dPhi/d");
 dnn_tree->Branch("nub1Pt",&nub1Pt,"nub1Pt/d");
 dnn_tree->Branch("nub1Eta",&nub1Eta,"nub1Eta/d");
 dnn_tree->Branch("nub1Phi",&nub1Phi,"nub1Phi/d");
 dnn_tree->Branch("nub1Mass",&nub1Mass,"nub1Mass/d");
 dnn_tree->Branch("nub1Ht",&nub1Ht,"nub1Ht/d");
 dnn_tree->Branch("nub1Mt",&nub1Mt,"nub1Mt/d");
 dnn_tree->Branch("nub2dR",&nub2dR,"nub2dR/d");
 dnn_tree->Branch("nub2dEta",&nub2dEta,"nub2dEta/d");
 dnn_tree->Branch("nub2dPhi",&nub2dPhi,"nub2dPhi/d");
 dnn_tree->Branch("nub2Pt",&nub2Pt,"nub2Pt/d");
 dnn_tree->Branch("nub2Eta",&nub2Eta,"nub2Eta/d");
 dnn_tree->Branch("nub2Phi",&nub2Phi,"nub2Phi/d");
 dnn_tree->Branch("nub2Mass",&nub2Mass,"nub2Mass/d");
 dnn_tree->Branch("nub2Ht",&nub2Ht,"nub2Ht/d");
 dnn_tree->Branch("nub2Mt",&nub2Mt,"nub2Mt/d");

 dnn_tree->Branch("lbbdR",&lbbdR,"lbbdR/d");
 dnn_tree->Branch("lbbdEta",&lbbdEta,"lbbdEta/d");
 dnn_tree->Branch("lbbdPhi",&lbbdPhi,"lbbdPhi/d");
 dnn_tree->Branch("lbbPt",&lbbPt,"lbbPt/d");
 dnn_tree->Branch("lbbEta",&lbbEta,"lbbEta/d");
 dnn_tree->Branch("lbbPhi",&lbbPhi,"lbbPhi/d");
 dnn_tree->Branch("lbbMass",&lbbMass,"lbbMass/d");
 dnn_tree->Branch("lbbHt",&lbbHt,"lbbHt/d");
 dnn_tree->Branch("lbbMt",&lbbMt,"lbbMt/d");

 dnn_tree->Branch("lb1dR",&lb1dR,"lb1dR/d");
 dnn_tree->Branch("lb1dEta",&lb1dEta,"lb1dEta/d");
 dnn_tree->Branch("lb1dPhi",&lb1dPhi,"lb1dPhi/d");
 dnn_tree->Branch("lb1Pt",&lb1Pt,"lb1Pt/d");
 dnn_tree->Branch("lb1Eta",&lb1Eta,"lb1Eta/d");
 dnn_tree->Branch("lb1Phi",&lb1Phi,"lb1Phi/d");
 dnn_tree->Branch("lb1Mass",&lb1Mass,"lb1Mass/d");
 dnn_tree->Branch("lb1Ht",&lb1Ht,"lb1Ht/d");
 dnn_tree->Branch("lb1Mt",&lb1Mt,"lb1Mt/d");
 dnn_tree->Branch("lb2dR",&lb2dR,"lb2dR/d");
 dnn_tree->Branch("lb2dEta",&lb2dEta,"lb2dEta/d");
 dnn_tree->Branch("lb2dPhi",&lb2dPhi,"lb2dPhi/d");
 dnn_tree->Branch("lb2Pt",&lb2Pt,"lb2Pt/d");
 dnn_tree->Branch("lb2Eta",&lb2Eta,"lb2Eta/d");
 dnn_tree->Branch("lb2Phi",&lb2Phi,"lb2Phi/d");
 dnn_tree->Branch("lb2Mass",&lb2Mass,"lb2Mass/d");
 dnn_tree->Branch("lb2Ht",&lb2Ht,"lb2Ht/d");
 dnn_tree->Branch("lb2Mt",&lb2Mt,"lb2Mt/d");

 dnn_tree->Branch("Wjb1dR",&Wjb1dR,"Wjb1dR/d");
 dnn_tree->Branch("Wjb1dEta",&Wjb1dEta,"Wjb1dEta/d");
 dnn_tree->Branch("Wjb1dPhi",&Wjb1dPhi,"Wjb1dPhi/d");
 dnn_tree->Branch("Wjb1Pt",&Wjb1Pt,"Wjb1Pt/d");
 dnn_tree->Branch("Wjb1Eta",&Wjb1Eta,"Wjb1Eta/d");
 dnn_tree->Branch("Wjb1Phi",&Wjb1Phi,"Wjb1Phi/d");
 dnn_tree->Branch("Wjb1Mass",&Wjb1Mass,"Wjb1Mass/d");
 dnn_tree->Branch("Wjb1Ht",&Wjb1Ht,"Wjb1Ht/d");
 dnn_tree->Branch("Wjb1Mt",&Wjb1Mt,"Wjb1Mt/d");
 dnn_tree->Branch("Wjb2dR",&Wjb2dR,"Wjb2dR/d");
 dnn_tree->Branch("Wjb2dEta",&Wjb2dEta,"Wjb2dEta/d");
 dnn_tree->Branch("Wjb2dPhi",&Wjb2dPhi,"Wjb2dPhi/d");
 dnn_tree->Branch("Wjb2Pt",&Wjb2Pt,"Wjb2Pt/d");
 dnn_tree->Branch("Wjb2Eta",&Wjb2Eta,"Wjb2Eta/d");
 dnn_tree->Branch("Wjb2Phi",&Wjb2Phi,"Wjb2Phi/d");
 dnn_tree->Branch("Wjb2Mass",&Wjb2Mass,"Wjb2Mass/d");
 dnn_tree->Branch("Wjb2Ht",&Wjb2Ht,"Wjb2Ht/d");
 dnn_tree->Branch("Wjb2Mt",&Wjb2Mt,"Wjb2Mt/d");

 //Tree for minimum dR analysis
 TTree * tree = new TTree( "tree", "tree for ttbb");
 tree->Branch("nbjets",&nbjets,"nbjets/s");
 tree->Branch("bjet1_pt",&bjet1_pt,"bjet1_pt/d");
 tree->Branch("bjet1_eta",&bjet1_eta,"bjet1_eta/d");
 tree->Branch("bjet1_phi",&bjet1_phi,"bjet1_phi/d");
 tree->Branch("bjet1_e",&bjet2_e,"bjet1_e/d");
 tree->Branch("bjet2_pt",&bjet2_pt,"bjet2_pt/d");
 tree->Branch("bjet2_eta",&bjet2_eta,"bjet2_eta/d");
 tree->Branch("bjet2_phi",&bjet2_phi,"bjet2_phi/d");
 tree->Branch("bjet2_e",&bjet2_e,"bjet2_e/d");
 tree->Branch("bjet3_pt",&bjet3_pt,"bjet3_pt/d");
 tree->Branch("bjet3_eta",&bjet3_eta,"bjet3_eta/d");
 tree->Branch("bjet3_phi",&bjet3_phi,"bjet3_phi/d");
 tree->Branch("bjet3_e",&bjet3_e,"bjet3_e/d");

 tree->Branch("njets",&njets,"njets/s");
 tree->Branch("Jet_pt",&Jet_pt,"Jet_pt/d");
 tree->Branch("Jet_eta",&Jet_eta,"Jet_eta/d");
 tree->Branch("Jet_phi",&Jet_phi,"Jet_phi/d");
 tree->Branch("Jet_e",&Jet_e,"Jet_e/d");

 tree->Branch("nElectron",&nElectron,"nElectron/s");
 tree->Branch("Electron1_pt",&Electron1_pt,"Electron1_pt/d");
 tree->Branch("Electron1_eta",&Electron1_eta,"Electron1_eta/d");
 tree->Branch("Electron1_phi",&Electron1_phi,"Electron1_phi/d");
 tree->Branch("Electron1_e",&Electron1_e,"Electron1_e/d");
 tree->Branch("Electron2_pt",&Electron2_pt,"Electron2_pt/d");
 tree->Branch("Electron2_eta",&Electron2_eta,"Electron2_eta/d");
 tree->Branch("Electron2_phi",&Electron2_phi,"Electron2_phi/d");
 tree->Branch("Electron2_e",&Electron2_e,"Electron2_e/d");
 
 tree->Branch("nMuon",&nMuon,"nMuon/s");
 tree->Branch("Muon1_pt",&Muon1_pt,"Muon1_pt/d");
 tree->Branch("Muon1_eta",&Muon1_eta,"Muon1_eta/d");
 tree->Branch("Muon1_phi",&Muon1_phi,"Muon1_phi/d");
 tree->Branch("Muon1_e",&Muon1_e,"Muon1_e/d");
 tree->Branch("Muon2_pt",&Muon2_pt,"Muon2_pt/d");
 tree->Branch("Muon2_eta",&Muon2_eta,"Muon2_eta/d");
 tree->Branch("Muon2_phi",&Muon2_phi,"Muon2_phi/d");
 tree->Branch("Muon2_e",&Muon2_e,"Muon2_e/d");

 tree->Branch("MET_px",&MET_px,"MET_px/d");
 tree->Branch("MET_py",&MET_py,"MET_py/d");
 tree->Branch("nLepton",&nLepton,"nLepton/s");
 tree->Branch("MATCHED",&MATCHED,"MATCHED/i");
 
 // Create object of class ExRootTreeReader
 ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
 Long64_t numberOfEntries = treeReader->GetEntries();

 // Get pointers to branches used in this analysis
 TClonesArray *branchGenJet  = treeReader->UseBranch("GenJet");
 TClonesArray *branchJet  = treeReader->UseBranch("Jet");
 TClonesArray *branchParticle  = treeReader->UseBranch("Particle");
 TClonesArray *branchElectron = treeReader->UseBranch("Electron");
 TClonesArray *branchMuon = treeReader->UseBranch("Muon");
 TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

 // Book histograms
 TH1 *histnjet = new TH1F("h_njet", "Number of jets", 14, 0.0, 14.0);
 TH1 *histnElectron = new TH1F("h_nElectron", "Number of electrons", 3, 0.0, 3.0);
 TH1 *histnMuon = new TH1F("h_nMuon", "Number of muons", 3, 0.0, 3.0);
 
 TH1 *histnbjet = new TH1F("h_nbjet", "Number of b-jets", 10, 0.0, 10.0);
 TH1 *histMbb = new TH1F("h_mbb", "M_{inv}(b, b)", 20, 20.0, 180.0);
 TH1 *histdRbb = new TH1F("h_dRbb", "dR(b, b)", 20, 0, 4.0);

 TH1 *hist_gennbjet = new TH1F("h_gennbjet", "Number of b-jets", 6, 0.0, 6.0);
 TH1 *hist_genMbb = new TH1F("h_genmbb", "M_{inv}(b, b)", 20, 20.0, 180.0);
 TH1 *hist_gendRbb = new TH1F("h_gendRbb", "dR(b, b)", 20, 0, 4.0);

 TH1 *hist_matchednbjet = new TH1F("h_matchednbjet", "Number of b-jets", 6, 0.0, 6.0);
 TH1 *hist_matchedMbb = new TH1F("h_matchedmbb", "M_{inv}(b, b)", 20, 20.0, 180.0);
 TH1 *hist_matcheddRbb = new TH1F("h_matcheddRbb", "dR(b, b)", 20, 0, 4.0);

 TH1 *hist_selection = new TH1F("h_selection","Selection",3,0.0,3.0);

 Int_t numberOfSelectedEvents = 0;
 Int_t numberOfMatchedEvents = 0;
 Int_t numberOfNewMatchedEvents = 0;

 vector<Jet *> bJets;
 vector<Jet *> Jets;
 vector<Electron *> Electrons;
 vector<Muon *> Muons;

 TLorentzVector p4[2];
 Jet *jet;
 Electron *electron;
 Muon *muon;
 MissingET *met;

 int entry, i;
 bool leptonSelection = false;

 // Loop over all events
 for(entry = 0; entry < numberOfEntries; ++entry)
 {
   if(entry%1000 == 0 && entry < 10000) cout << "event number: " << entry << endl;
   else if(entry%10000 == 0) cout<< "event number: " << entry << endl;
   // Load selected branches with data from specified event
   treeReader->ReadEntry(entry);
   
   // Initialization
   Jets.clear();
   Electrons.clear();
   Muons.clear(); 
   bJets.clear();
   
   // dnn variables
   signal = -1;
   event = -1;
   njets = 999;
   nbjets = 999;
   addbjet1_pt = 999;
   addbjet1_eta = 999;
   addbjet1_phi = 999;
   addbjet1_e = 999;
   addbjet2_pt = 999;
   addbjet2_eta = 999;
   addbjet2_phi = 999;
   addbjet2_e = 999;
   pt1 = 999;
   eta1 = 999;
   phi1 = 999;
   e1 = 999;
   pt2 = 999;
   eta2 = 999;
   phi2 = 999;
   e2 = 999;
   bbdR = 999;
   bbdEta = 999;
   bbdPhi = 999;
   bbPt = 999;
   bbEta = 999;
   bbPhi = 999;
   bbMass = 999;
   bbHt = 999;
   bbMt = 999;

   nubbdR = 999;
   nubbdEta = 999;
   nubbdPhi = 999;
   nubbPt = 999;
   nubbEta = 999;
   nubbPhi = 999;
   nubbMass = 999;
   nubbHt = 999;
   nubbMt = 999;
   nub1dR = 999;
   nub1dEta = 999;
   nub1dPhi = 999;
   nub1Pt = 999;
   nub1Eta = 999;
   nub1Phi = 999;
   nub1Mass = 999;
   nub1Ht = 999;
   nub1Mt = 999;
   nub2dR = 999;
   nub2dEta = 999;
   nub2dPhi = 999;
   nub2Pt = 999;
   nub2Eta = 999;
   nub2Phi = 999;
   nub2Mass = 999;
   nub2Ht = 999;
   nub2Mt = 999;

   lbbdR = 999;
   lbbdEta = 999;
   lbbdPhi = 999;
   lbbPt = 999;
   lbbEta = 999;
   lbbPhi = 999;
   lbbMass = 999;
   lbbHt = 999;
   lbbMt = 999;
   lb1dR = 999;
   lb1dEta = 999;
   lb1dPhi = 999;
   lb1Pt = 999;
   lb1Eta = 999;
   lb1Phi = 999;
   lb1Mass = 999;
   lb1Ht = 999;
   lb1Mt = 999;
   lb2dR = 999;
   lb2dEta = 999;
   lb2dPhi = 999;
   lb2Pt = 999;
   lb2Eta = 999;
   lb2Phi = 999;
   lb2Mass = 999;
   lb2Ht = 999;
   lb2Mt = 999;

   Wjb1dR = 999;
   Wjb1dEta = 999;
   Wjb1dPhi = 999;
   Wjb1Pt = 999;
   Wjb1Eta = 999;
   Wjb1Phi = 999;
   Wjb1Mass = 999;
   Wjb1Ht = 999;
   Wjb1Mt = 999;
   Wjb2dR = 999;
   Wjb2dEta = 999;
   Wjb2dPhi = 999;
   Wjb2Pt = 999;
   Wjb2Eta = 999;
   Wjb2Phi = 999;
   Wjb2Mass = 999;
   Wjb2Ht = 999;
   Wjb2Mt = 999;

   // tree variables
   Jet_pt = 999;
   Jet_eta = 999;
   Jet_phi = 999;
   Jet_e = 999;
   Electron1_pt = 999;
   Electron1_eta = 999;
   Electron1_phi = 999;
   Electron1_e = 999;
   Electron2_pt = 999;
   Electron2_eta = 999;
   Electron2_phi = 999;
   Electron2_e = 999;
   Muon1_pt = 999;
   Muon1_eta = 999;
   Muon1_phi = 999;
   Muon1_e = 999;
   Muon2_pt = 999;
   Muon2_eta = 999;
   Muon2_phi = 999;
   Muon2_e = 999;
   bjet1_pt = 999;
   bjet1_eta = 999;
   bjet1_phi = 999;
   bjet1_e = 999;
   bjet2_pt = 999 ;
   bjet2_eta = 999;
   bjet2_phi = 999;
   bjet2_e = 999;
   bjet3_pt = 999 ;
   bjet3_eta = 999;
   bjet3_phi = 999;
   bjet3_e = 999;
   MET_px = 999;
   MET_py = 999;
   MATCHED = -1;

   //Genaddbjet Selection (S1)
   int nb = 0;
   int nbFromTop = 0;
   int nb_status3 = 0; 
   int ntop = 0;
   vector<GenParticle*> GenAddbQuarks;
   vector<Jet*> GenAddbJets;
   for(i = 0; i < branchParticle->GetEntriesFast(); ++i){
     GenParticle *genP = (GenParticle*) branchParticle->At(i);
     if( abs(genP->PID) == 6 ) {
       GenParticle *dauP1 = (GenParticle *) branchParticle->At(genP->D1);
       GenParticle *dauP2 = (GenParticle *) branchParticle->At(genP->D2);
       double dauPID1 = dauP1->PID;
       double dauPID2 = dauP2->PID;
       if( abs(dauPID1) != 6 && abs(dauPID2) != 6){
         ntop++;
       }
     }
     if( abs(genP->PID) == 5){
       if( genP->Status == 2) nb_status3++;
       //check if this is the last b quark 
       GenParticle *dauP1 = (GenParticle *) branchParticle->At(genP->D1);
       GenParticle *dauP2 = (GenParticle *) branchParticle->At(genP->D2);
       double dauPID1 = dauP1->PID;
       double dauPID2 = dauP2->PID;
       if( abs(dauPID1) == 5 || abs(dauPID2) == 5) continue; 
       //cout << "test b quark = " << genP->P4().Pt() << " decays to " << dauPID1 << " and " << dauPID2 << endl;
       nb++;
       bool fromTop = isFromTop(genP, branchParticle);
       if(fromTop) {
         nbFromTop++;
       }else{
         GenAddbQuarks.push_back(genP);
       }
     }
   } 

   for(i = 0; i < branchGenJet->GetEntriesFast(); ++i){
     Jet *genjet = (Jet*) branchGenJet->At(i);
     for( int j = 0 ; j < GenAddbQuarks.size() ; j++){
       double dR = GenAddbQuarks[j]->P4().DeltaR( genjet->P4() );
       if( dR < 0.5 ) GenAddbJets.push_back(genjet);
     }
   }

   //cout << "=========" << " Number of top = " << ntop << " number of b = " << nb << " (from top = " << nbFromTop << " )" << "=========" << endl;
   
   //Select 2 addbjets in descending order of Pt
   TLorentzVector addbjet[2];
   //cout<<"gen addbjet selection"<<endl;
   for( int j = 0; j < GenAddbJets.size(); j++){
     if( GenAddbJets[j]->P4().Pt() < 20 || fabs(GenAddbJets[j]->P4().Eta()) > 2.5 ) continue;
     if( GenAddbJets[j]->P4().Pt() > addbjet[0].Pt() ){
       addbjet[1] = addbjet[0];
       addbjet[0] = GenAddbJets[j]->P4();
     }
     else if( (GenAddbJets[j]->P4().Pt() < addbjet[0].Pt()) && (GenAddbJets[j]->P4().Pt() > addbjet[1].Pt() )){
       addbjet[1] = GenAddbJets[j]->P4();
     }
   }
   if( addbjet[0].Pt() == 0 || addbjet[1].Pt() == 0 ) continue;   // S1
   s1++;

   //Lepton Selection (S2)
   //Electron Selection
   nElectron = 0; //unsigned short
	for(i = 0; i < branchElectron->GetEntries(); ++i){
		electron = (Electron*) branchElectron->At(i);
		if(Electrons.size() == 0 && (electron->PT > 30) && (fabs(electron->Eta < 2.5))){
				Electrons.push_back(electron);
      	}
		else if(Electrons.size() > 0 && (electron->PT > 25) && (fabs(electron->Eta < 2.5))){
				Electrons.push_back(electron);
		}
		else continue;
	}
	nElectron = Electrons.size();
	//cout << "# of e = " << nElectron << endl;

	//Muon Selection
	nMuon = 0;
	for(i = 0; i < branchMuon->GetEntries(); ++i){
		muon = (Muon*) branchMuon->At(i);
		if(Muons.size() == 0 && (muon->PT > 30) && (fabs(muon->Eta < 2.5))){
			Muons.push_back(muon);
		}
		else if(Muons.size() > 0 && (muon->PT > 25) && (fabs(muon->Eta < 2.5))){
			Muons.push_back(muon);
		}
		else continue;
	}
	nMuon = Muons.size();
  // 	cout << "# of m = " << nMuon <<endl;
   
   nLepton = nElectron + nMuon;
   //cout << "# of Lepton = " << nLepton << endl;
   if( isdilepton ){
     leptonSelection = (nLepton >= 2);
   }
   else{
     leptonSelection = (nLepton == 1);
   }
   if (!leptonSelection) continue;   // S2
   s2++;

   // Jet and b-tag Selections (S3)
   njets = 0;
   nbjets = 0;
   for(i = 0; i < branchJet->GetEntriesFast(); ++i){
     jet = (Jet*) branchJet->At(i);
     if( (jet->PT < 30) || (fabs(jet->Eta) > 2.5) ) continue;
     njets++;
     Jets.push_back(jet);
     if( jet->BTag ) {
       nbjets++;
       bJets.push_back(jet);
     }
   }
   if ( njets < jcut || nbjets < bcut ) continue;   // S3
   s3++;
   //cout<<"njet "<<njets<<" nbjet "<<nbjets<<" nlepton "<<nLepton<<endl;

   TLorentzVector nu;
   if( branchMissingET->GetEntriesFast() > 0){
     met = (MissingET*) branchMissingET->At(0);
     MET_px = met->MET * cos( met->Phi);
     MET_py = met->MET * sin( met->Phi);
     //cout << "px " << MET_px <<" py "<< MET_py<<" MET "<<met->MET<<endl;
     nu.SetPxPyPzE( MET_px, MET_py, 0, met->MET );
   }

   // MissingET 4-vector
   //TLorentzVector nu;
   //if( branchMissingET->GetEntriesFast() > 0){
   //  met = (MissingET*) branchMissingET->At(0);
     //cout << "px " << MET_px <<" py "<< MET_py<<" MET "<<met->MET<<endl;
     
     //nu = met->P4();
   //}

   //Lepton 4-vector ( only for lep+jet )
   TLorentzVector lep;  //lep[2]
   if ( nElectron == 1)
     lep = Electrons[0]->P4();
   else if ( nMuon == 1)
     lep = Muons[0]->P4();
// e1m1  e2  m2
   // Fill the tree ntuples (minimum dR)
   Jet_pt = Jets[0]->P4().Pt();
   Jet_eta = Jets[0]->P4().Eta();
   Jet_phi = Jets[0]->P4().Phi();
   Jet_e = Jets[0]->P4().E();
   if( nElectron == 1){
     Electron1_pt = Electrons[0]->P4().Pt();
     Electron1_eta = Electrons[0]->P4().Eta();
     Electron1_phi = Electrons[0]->P4().Phi();
     Electron1_e = Electrons[0]->P4().E();
   }
   if( nMuon == 1 ){
     Muon1_pt = Muons[0]->P4().Pt();
     Muon1_eta = Muons[0]->P4().Eta();
     Muon1_phi = Muons[0]->P4().Phi();
     Muon1_e = Muons[0]->P4().E();
   }
   if(isdilepton){
     if(nElectron == 2){
       Electron2_pt = Electrons[1]->P4().Pt();
       Electron2_eta = Electrons[1]->P4().Eta();
       Electron2_phi = Electrons[1]->P4().Phi();
       Electron2_e = Electrons[1]->P4().E();
     }
     if(nMuon == 2){
       Muon2_pt = Muons[1]->P4().Pt();
       Muon2_eta = Muons[1]->P4().Eta();
       Muon2_phi = Muons[1]->P4().Phi();
       Muon2_e = Muons[1]->P4().E();
     }
   } 
   bjet1_pt = bJets[0]->P4().Pt();
   bjet1_eta = bJets[0]->P4().Eta();
   bjet1_phi = bJets[0]->P4().Phi();
   bjet1_e = bJets[0]->P4().E();
   bjet2_pt = bJets[1]->P4().Pt();
   bjet2_eta = bJets[1]->P4().Eta();
   bjet2_phi = bJets[1]->P4().Phi();
   bjet2_e = bJets[1]->P4().E();

   if(bJets.size() >=3){
     bjet3_pt = bJets[2]->P4().Pt();
     bjet3_eta = bJets[2]->P4().Eta();
     bjet3_phi = bJets[2]->P4().Phi();
     bjet3_e = bJets[2]->P4().E();
   }
   
   histnjet->Fill( njets );
   histnElectron->Fill( nElectron );
   histnMuon->Fill( nMuon );

   //Matched b-jets
   /*vector<Jet*> MatchedbJets;
   for( int j = 0; j < bJets.size(); j++){
     TLorentzVector recobjet = bJets[j]->P4();
     for(int k = 0 ; k < GenAddbJets.size(); k++){
       TLorentzVector genbjet = GenAddbJets[k]->P4();
       double dR = recobjet.DeltaR( genbjet );
       //cout <<" test dR = " << dR << endl;
       if( dR < 0.4 ) {
         MatchedbJets.push_back( bJets[j] ) ;
         //cout<<entry<<" "<<j<<k<<"MatchedbJet pt : "<< MatchedbJets[MatchedbJets.size()-1]->P4().Pt() << endl;
       }
     }
   }*/
   //for(int a=0; a < MatchedbJets.size(); a++) cout<<MatchedbJets[a]->P4().Pt()<<endl;
   histnbjet->Fill(bJets.size());
   hist_gennbjet->Fill(GenAddbJets.size());
   if( GenAddbJets.size() > 1){
     double gen_Mbb = ( GenAddbJets[0]->P4() + GenAddbJets[1]->P4() ).M();
     double gen_dRbb = GenAddbJets[0]->P4().DeltaR( GenAddbJets[1]->P4() );
     hist_genMbb->Fill( gen_Mbb );
     hist_gendRbb->Fill( gen_dRbb );
   }
   /*hist_matchednbjet->Fill( MatchedbJets.size() );
   if( MatchedbJets.size() > 1){
     double matched_mbb = ( MatchedbJets[0]->P4() + MatchedbJets[1]->P4() ).M();
     double matched_dRbb = MatchedbJets[0]->P4().DeltaR( MatchedbJets[1]->P4() ); 
     hist_matchedMbb->Fill(matched_mbb);
     hist_matcheddRbb->Fill(matched_dRbb);
   }*/
   
   //Jet Combination
   TLorentzVector Wj;
   TLorentzVector tmpWj;
   double tmpWjM = -1;
   for(int j1 = 0; j1 < njets-1; j1++) {
     if(Jets[j1]->BTag) continue;
     TLorentzVector jet1 = Jets[j1]->P4();
     for(int j2 = j1+1; j2 < njets; j2++) {
       if(Jets[j2]->BTag) continue;
       TLorentzVector jet2 = Jets[j2]->P4();

       tmpWj = jet1 + jet2;
       if( abs(tmpWj.M() - 80.2 ) < tmpWjM ) {
         Wj = tmpWj;
         tmpWjM = tmpWj.M();
       }
     }
   }//Jet combi

   float mbb = 999;
   float dRbb = 999;

   // Select two bjets with minimum dR and fill the dnn ntuples
   TLorentzVector RecoAddJets[2];
   for(int b1 = 0; b1 < bJets.size()-1; b1++){
     for(int b2 = b1+1; b2 < bJets.size(); b2++){
       p4[0] = bJets[b1]->P4();
       p4[1] = bJets[b2]->P4();

       float tmp_mbb = ((p4[0]) + (p4[1])).M();
       float tmp_dRbb = p4[0].DeltaR(p4[1]);
       if(tmp_dRbb < dRbb) {
          dRbb = tmp_dRbb;
          mbb = tmp_mbb; 
          RecoAddJets[0] = p4[0];
          RecoAddJets[1] = p4[1];
       }
 
       //bool p4_1_matched = std::find(MatchedbJets.begin(), MatchedbJets.end(), bJets[b1]) != MatchedbJets.end(); 
       //bool p4_2_matched = std::find(MatchedbJets.begin(), MatchedbJets.end(), bJets[b2]) != MatchedbJets.end(); 
       //bool matched_signal = p4_1_matched && p4_2_matched;
       bool matched_signal1 = ( p4[0].DeltaR( addbjet[0] ) < 0.4 ) && ( p4[1].DeltaR( addbjet[1] ) < 0.4 );
       bool matched_signal2 = ( p4[0].DeltaR( addbjet[1] ) < 0.4 ) && ( p4[1].DeltaR( addbjet[0] ) < 0.4 );
       
       if (matched_signal1 || matched_signal2) signal = 1;
       else signal = 0;
       //cout<<signal<<endl;

       event = entry;

       pt1 = p4[0].Pt();
       eta1 = p4[0].Eta();
       phi1 = p4[0].Phi();
       e1 = p4[0].E();
       pt2 = p4[1].Pt();
       eta2 = p4[1].Eta();
       phi2 = p4[1].Phi();
       e2 = p4[1].E();

       addbjet1_pt = addbjet[0].Pt();
       addbjet1_eta = addbjet[0].Eta();
       addbjet1_phi = addbjet[0].Phi();
       addbjet1_e = addbjet[0].E();
       addbjet2_pt = addbjet[1].Pt();
       addbjet2_eta = addbjet[1].Eta();
       addbjet2_phi = addbjet[1].Phi();
       addbjet2_e = addbjet[1].E();

       bbdR = tmp_dRbb;
       bbdEta = abs(p4[0].Eta()-p4[1].Eta());
       bbdPhi = p4[0].DeltaPhi(p4[1]);
       bbPt = (p4[0]+p4[1]).Pt();
       bbEta = (p4[0]+p4[1]).Eta();
       bbPhi = (p4[0]+p4[1]).Phi();
       bbMass = tmp_mbb;
       bbHt = p4[0].Pt()+p4[1].Pt();
       bbMt = (p4[0]+p4[1]).Mt();

       nubbdR = (p4[0]+p4[1]).DeltaR(nu);
       nubbdEta = abs((p4[0]+p4[1]).Eta()-nu.Eta());
       nubbdPhi = (p4[0]+p4[1]).DeltaPhi(nu);
       nubbPt = (p4[0]+p4[1]+nu).Pt();
       nubbEta = (p4[0]+p4[1]+nu).Eta();
       nubbPhi = (p4[0]+p4[1]+nu).Phi();
       nubbMass = (p4[0]+p4[1]+nu).M();
       nubbHt = (p4[0]+p4[1]).Pt()+nu.Pt();
       nubbMt = (p4[0]+p4[1]+nu).Mt();

       nub1dR = p4[0].DeltaR(nu);
       nub1dEta = abs(p4[0].Eta()-nu.Eta());
       nub1dPhi = p4[0].DeltaPhi(nu);
       nub1Pt = (p4[0]+nu).Pt();
       nub1Eta = (p4[0]+nu).Eta();
       nub1Phi = (p4[0]+nu).Phi();
       nub1Mass = (p4[0]+nu).M();
       nub1Ht = p4[0].Pt()+nu.Pt();
       nub1Mt = (p4[0]+nu).Mt();

       nub2dR = p4[1].DeltaR(nu);
       nub2dEta = abs(p4[1].Eta()-nu.Eta());
       nub2dPhi = p4[1].DeltaPhi(nu);
       nub2Pt = (p4[1]+nu).Pt();
       nub2Eta = (p4[1]+nu).Eta();
       nub2Phi = (p4[1]+nu).Phi();
       nub2Mass = (p4[1]+nu).M();
       nub2Ht = p4[1].Pt()+nu.Pt();
       nub2Mt = (p4[1]+nu).Mt();

       lbbdR = (p4[0]+p4[1]).DeltaR(lep);
       lbbdEta = abs((p4[0]+p4[1]).Eta()-lep.Eta());
       lbbdPhi = (p4[0]+p4[1]).DeltaPhi(lep);
       lbbPt = (p4[0]+p4[1]+lep).Pt();
       lbbEta = (p4[0]+p4[1]+lep).Eta();
       lbbPhi = (p4[0]+p4[1]+lep).Phi();
       lbbMass = (p4[0]+p4[1]+lep).M();
       lbbHt = (p4[0]+p4[1]).Pt()+lep.Pt();
       lbbMt = (p4[0]+p4[1]+lep).Mt();

       lb1dR = p4[0].DeltaR(lep);
       lb1dEta = abs(p4[0].Eta()-lep.Eta());
       lb1dPhi = p4[0].DeltaPhi(lep);
       lb1Pt = (p4[0]+lep).Pt();
       lb1Eta = (p4[0]+lep).Eta();
       lb1Phi = (p4[0]+lep).Phi();
       lb1Mass = (p4[0]+lep).M();
       lb1Ht = p4[0].Pt()+lep.Pt();
       lb1Mt = (p4[0]+lep).Mt();
       lb2dR = p4[1].DeltaR(lep);
       lb2dEta = abs(p4[1].Eta()-lep.Eta());
       lb2dPhi = p4[1].DeltaPhi(lep);
       lb2Pt = (p4[1]+lep).Pt();
       lb2Eta = (p4[1]+lep).Eta();
       lb2Phi = (p4[1]+lep).Phi();
       lb2Mass = (p4[1]+lep).M();
       lb2Ht = p4[1].Pt()+lep.Pt();
       lb2Mt = (p4[1]+lep).Mt();

       Wjb1dR = p4[0].DeltaR(Wj);
       Wjb1dEta = abs(p4[0].Eta()-Wj.Eta());
       Wjb1dPhi = p4[0].DeltaPhi(Wj);
       Wjb1Pt = (p4[0]+Wj).Pt();
       Wjb1Eta = (p4[0]+Wj).Eta();
       Wjb1Phi = (p4[0]+Wj).Phi();
       Wjb1Mass = (p4[0]+Wj).M();
       Wjb1Ht = p4[0].Pt()+Wj.Pt();
       Wjb1Mt = (p4[0]+Wj).Mt();
       Wjb2dR = p4[1].DeltaR(Wj);
       Wjb2dEta = abs(p4[1].Eta()-Wj.Eta());
       Wjb2dPhi = p4[1].DeltaPhi(Wj);
       Wjb2Pt = (p4[1]+Wj).Pt();
       Wjb2Eta = (p4[1]+Wj).Eta();
       Wjb2Phi = (p4[1]+Wj).Phi();
       Wjb2Mass = (p4[1]+Wj).M();
       Wjb2Ht = p4[1].Pt()+Wj.Pt();
       Wjb2Mt = (p4[1]+Wj).Mt();

       dnn_tree->Fill();
     }
   }
   hist_selection->SetBinContent(1,s1);
   hist_selection->SetBinContent(2,s2);
   hist_selection->SetBinContent(3,s3);

   bool matched = false;
   bool matched1 = ( RecoAddJets[0].DeltaR( addbjet[0] ) < 0.4 ) && ( RecoAddJets[1].DeltaR( addbjet[1] ) < 0.4 );
   bool matched2 = ( RecoAddJets[0].DeltaR( addbjet[1] ) < 0.4 ) && ( RecoAddJets[1].DeltaR( addbjet[0] ) < 0.4 );
   if ( matched1 || matched2 ) matched = true;
   //cout<<"nElectron = "<< nElectron <<" nMuon = "<< nMuon <<" njets = "<< njets <<" nbjets = "<< nbjets << endl;

   ++numberOfSelectedEvents;
   if(matched) {
     numberOfMatchedEvents++;
     MATCHED = 1;
   }
   else MATCHED = 0;

   histMbb->Fill(mbb);
   histdRbb->Fill(dRbb);
   tree->Fill();

 }

 cout<<"Event Info : jet >= "<<jcut<<" bjet >= "<<bcut<<endl;
 //cout << "Total number of events = " << numberOfSelectedEvents << endl;
 //cout << "Total number of matched events = " << numberOfMatchedEvents << endl;
 double eff = (double) numberOfMatchedEvents/ (double) numberOfSelectedEvents;
 cout << "Matching eff. = " << eff << " ( "<<numberOfMatchedEvents<<" / "<<numberOfSelectedEvents<<" )"<< endl;
 double accept1 = (double) s1 / (double) numberOfEntries; 
 double accept2 = (double) s2 / (double) numberOfEntries;
 double accept3 = (double) s3 / (double) numberOfEntries;
 cout << "Entries "<<numberOfEntries<<endl;
 cout << "Acceptance1 (S1/Entry) = "<<accept1<<" ( "<<s1<<" )"<<endl;
 cout << "Acceptance2 (S2/Entry) = "<<accept2<<" ( "<<s2<<" )"<<endl;
 cout << "Acceptance3 (S3/Entry) = "<<accept3<<" ( "<<s3<<" )"<<endl;

 fout->Write();
 fout->Close();

}
