#include <iostream>
#include <fstream>
#include <string>

#include <TStyle.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TPaveStats.h>
#include <TFile.h>

//#include "/gpfs/home/youhei/atlasstyle-00-03-05/AtlasUtils.C"
//#include "/gpfs/home/yfukuhar/RootUtils/rootlogon.C"
//#include "/gpfs/home/youhei/atlasstyle-00-03-05/AtlasStyle.C"

void macro_closeby_yfukuhara(string TrigChain="mu4", TString epsdir="") {
  //TString pdf = "temp424100.pdf";
  //TString filename = "/gpfs/fs7001/yfukuhar/workspace/TestCalcEffTool/CalcEffTool/run/FTK_424100_100_.root";

  TString pdf = "temp3p5.pdf";
  TString filename = "/gpfs/fs7001/yfukuhar/workspace/TestCalcEffTool/CalcEffTool/run/Jpsi_mu3p5mu3p5_10.root";

  // TString filename = "run/test.root";
  // TString filename = "~/L2MuonSA/dataset/calceff_ntuple/Jpsi_mu3p5mu3p5_NoTag_noDR.root";
  //TString filename = "~youhei/L2MuonSA/dataset/calceff_ntuple/Jpsi_mu3p5mu3p5_NoTag_noDR_RoInumCut.root";
  //TString filename = "/gpfs/fs7001/yfukuhar/workspace/TestCalcEffTool/CalcEffTool/run/FTK_424103_50.root";
  //TString filename = "/gpfs/fs7001/yfukuhar/workspace/TestCalcEffTool/CalcEffTool/run/FTK_424100_200.root";
  // TString filename = "~/L2MuonSA/dataset/calceff_ntuple/Jpsi_default1k_NoTag_roiWord.root";

  rootlogon();
  TColor::InvertPalette();
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetEndErrorSize(2);


  Color_t col0_dR = kAzure+1;
  Color_t col1_dR = kViolet+1;
  Color_t col2_dR = kPink+1;
  Color_t col3_dR = kOrange+1;

  Color_t col01_dR = kViolet-6;
  Color_t col12_dR = kPink+9;
  Color_t col23_dR = kOrange+7;
  Color_t col14_dR = kAzure-3;
  Color_t col05_dR = kGreen+1;
  Color_t col56_dR = kOrange;

  TFile* ifile = TFile::Open(filename);
  TTree* tree = (TTree*)ifile->Get("t_tap");

  int EventNumber;
  int mes_n;
  vector < string >* mes_name;
  int tag_proc;
  double tag_dRL1;
  double tag_pt;
  double tag_eta;
  double tag_exteta;
  double tag_extinneta;
  double tag_phi;
  double tag_extphi;
  double tag_extinnphi;
  double tag_charge;
  double tag_L1_pt;
  double tag_L1_eta;
  double tag_L1_phi;
  double tag_SA_pt;
  double tag_SA_eta;
  double tag_SA_phi;
  double tag_CB_pt;
  double tag_CB_eta;
  double tag_CB_phi;
  double probe_pt;
  double probe_eta;
  double probe_phi;
  double probe_charge;
  double tp_dR;
  double tp_deta;
  double tp_dphi;
  double tp_extdR;
  double tp_extdeta;
  double tp_extdphi;
  double tp_mass;
  vector < int >* probe_mesL1_pass;
  vector < double >* probe_mesL1_dR;
  vector < double >* probe_mesL1_tpdR;
  vector < double >* probe_mesL1_pt;
  vector < double >* probe_mesL1_eta;
  vector < double >* probe_mesL1_phi;
  vector < int >* probe_mesL1_MoreCandInRoI;
  vector < int >* probe_mesL1_MoreCandInSector;
  vector < int >* probe_mesL1_charge;
  vector < int >* probe_mesSA_pass;
  vector < double >* probe_mesSA_sAddress;
  vector < double >* probe_mesSA_dR;
  vector < double >* probe_mesSA_tpdR;
  vector < double >* probe_mesSA_pt;
  vector < double >* probe_mesSA_eta;
  vector < double >* probe_mesSA_phi;
  vector < double >* probe_mesSA_etams;
  vector < double >* probe_mesSA_phims;
  vector < double >* probe_mesSA_etabe;
  vector < double >* probe_mesSA_phibe;
  vector < double >* probe_mesSA_ptFtk;
  vector < double >* probe_mesSA_etaFtk;
  vector < double >* probe_mesSA_phiFtk;
  vector < int >* probe_mesSA_roadAlgo;
  vector < int >* probe_mesCB_pass;
  vector < double >* probe_mesCB_dR;
  vector < double >* probe_mesCB_tpdR;
  vector < double >* probe_mesCB_pt;
  vector < double >* probe_mesCB_eta;
  vector < double >* probe_mesCB_phi;

  EventNumber = 0;
  mes_n             = 0;
  mes_name          = new vector < string > ();
  tag_proc          = -99999.;
  tag_dRL1          = -99999.;
  tag_pt            = -99999.;
  tag_eta           = -99999.;
  tag_exteta        = -99999.;
  tag_extinneta     = -99999.;
  tag_phi           = -99999.;
  tag_extphi        = -99999.;
  tag_extinnphi     = -99999.;
  tag_charge        = -99999.;
  tag_L1_pt         = -99999.;
  tag_L1_eta        = -99999.;
  tag_L1_phi        = -99999.;
  tag_SA_pt         = -99999.;
  tag_SA_eta        = -99999.;
  tag_SA_phi        = -99999.;
  tag_CB_pt         = -99999.;
  tag_CB_eta        = -99999.;
  tag_CB_phi        = -99999.;
  probe_pt          = -99999.;
  probe_eta         = -99999.;
  probe_phi         = -99999.;
  probe_charge      = -99999.;
  tp_dR             = -99999.;
  tp_deta           = -99999.;
  tp_dphi           = -99999.;
  tp_extdR          = -99999.;
  tp_extdeta        = -99999.;
  tp_extdphi        = -99999.;
  tp_mass           = -99999.;
  probe_mesL1_pass  = new vector < int > ();
  probe_mesL1_dR    = new vector < double > ();
  probe_mesL1_tpdR  = new vector < double > ();
  probe_mesL1_pt    = new vector < double > ();
  probe_mesL1_eta   = new vector < double > ();
  probe_mesL1_phi   = new vector < double > ();
  probe_mesL1_MoreCandInRoI    = new vector < int > ();
  probe_mesL1_MoreCandInSector = new vector < int > ();
  probe_mesL1_charge           = new vector < int > ();
  probe_mesSA_pass  = new vector < int > ();
  probe_mesSA_sAddress    = new vector < double > ();
  probe_mesSA_dR    = new vector < double > ();
  probe_mesSA_tpdR  = new vector < double > ();
  probe_mesSA_pt    = new vector < double > ();
  probe_mesSA_eta   = new vector < double > ();
  probe_mesSA_phi   = new vector < double > ();
  probe_mesSA_etams = new vector < double > ();
  probe_mesSA_phims = new vector < double > ();
  probe_mesSA_etabe = new vector < double > ();
  probe_mesSA_phibe = new vector < double > ();
  probe_mesSA_ptFtk = new vector < double > ();
  probe_mesSA_etaFtk = new vector < double > ();
  probe_mesSA_phiFtk = new vector < double > ();
  probe_mesSA_roadAlgo = new vector < int > ();
  probe_mesCB_pass  = new vector < int > ();
  probe_mesCB_dR    = new vector < double > ();
  probe_mesCB_tpdR  = new vector < double > ();
  probe_mesCB_pt    = new vector < double > ();
  probe_mesCB_eta   = new vector < double > ();
  probe_mesCB_phi   = new vector < double > ();

  //--------------------------------------------------
  // BRANCH SET UP
  //--------------------------------------------------

  tree->SetBranchAddress( "EventNumber",        &EventNumber );
  tree->SetBranchAddress( "mes_n",        &mes_n );
  tree->SetBranchAddress( "mes_name",     &mes_name );
  tree->SetBranchAddress( "tag_proc",     &tag_proc );
  tree->SetBranchAddress( "tag_dRL1",     &tag_dRL1 );
  tree->SetBranchAddress( "tag_pt",       &tag_pt );
  tree->SetBranchAddress( "tag_eta",      &tag_eta );
  tree->SetBranchAddress( "tag_phi",      &tag_phi );
  tree->SetBranchAddress( "tag_charge",   &tag_charge );
  tree->SetBranchAddress( "tag_L1_pt",    &tag_L1_pt );
  tree->SetBranchAddress( "tag_L1_eta",   &tag_L1_eta );
  tree->SetBranchAddress( "tag_L1_phi",   &tag_L1_phi );
  tree->SetBranchAddress( "tag_SA_pt",    &tag_SA_pt );
  tree->SetBranchAddress( "tag_SA_eta",   &tag_SA_eta );
  tree->SetBranchAddress( "tag_SA_phi",   &tag_SA_phi );
  tree->SetBranchAddress( "tag_CB_pt",    &tag_CB_pt );
  tree->SetBranchAddress( "tag_CB_eta",   &tag_CB_eta );
  tree->SetBranchAddress( "tag_CB_phi",   &tag_CB_phi );
  tree->SetBranchAddress( "probe_pt",           &probe_pt );
  tree->SetBranchAddress( "probe_eta",          &probe_eta );
  tree->SetBranchAddress( "probe_phi",          &probe_phi );
  tree->SetBranchAddress( "probe_charge",       &probe_charge );
  tree->SetBranchAddress( "tp_dR",              &tp_dR );
  tree->SetBranchAddress( "tp_deta",            &tp_deta );
  tree->SetBranchAddress( "tp_dphi",            &tp_dphi );
  tree->SetBranchAddress( "tp_extdR",           &tp_extdR );
  tree->SetBranchAddress( "tp_extdeta",         &tp_extdeta );
  tree->SetBranchAddress( "tp_extdphi",         &tp_extdphi );
  tree->SetBranchAddress( "tp_mass",            &tp_mass );
  tree->SetBranchAddress( "probe_mesL1_pass",   &probe_mesL1_pass );
  tree->SetBranchAddress( "probe_mesL1_dR",     &probe_mesL1_dR );
  tree->SetBranchAddress( "probe_mesL1_tpdR",   &probe_mesL1_tpdR );
  tree->SetBranchAddress( "probe_mesL1_pt",     &probe_mesL1_pt );
  tree->SetBranchAddress( "probe_mesL1_eta",    &probe_mesL1_eta );
  tree->SetBranchAddress( "probe_mesL1_phi",    &probe_mesL1_phi );
  tree->SetBranchAddress( "probe_mesL1_MoreCandInRoI",    &probe_mesL1_MoreCandInRoI );
  tree->SetBranchAddress( "probe_mesL1_MoreCandInSector", &probe_mesL1_MoreCandInSector );
  tree->SetBranchAddress( "probe_mesL1_charge",           &probe_mesL1_charge );
  tree->SetBranchAddress( "probe_mesSA_pass",   &probe_mesSA_pass );
  tree->SetBranchAddress( "probe_mesSA_sAddress",     &probe_mesSA_sAddress );
  tree->SetBranchAddress( "probe_mesSA_dR",     &probe_mesSA_dR );
  tree->SetBranchAddress( "probe_mesSA_tpdR",   &probe_mesSA_tpdR );
  tree->SetBranchAddress( "probe_mesSA_pt",     &probe_mesSA_pt );
  tree->SetBranchAddress( "probe_mesSA_eta",    &probe_mesSA_eta );
  tree->SetBranchAddress( "probe_mesSA_phi",    &probe_mesSA_phi );
  tree->SetBranchAddress( "probe_mesSA_etams",  &probe_mesSA_etams );
  tree->SetBranchAddress( "probe_mesSA_phims",  &probe_mesSA_phims );
  tree->SetBranchAddress( "probe_mesSA_etabe",  &probe_mesSA_etabe );
  tree->SetBranchAddress( "probe_mesSA_phibe",  &probe_mesSA_phibe );
  tree->SetBranchAddress( "probe_mesSA_ptFtk",  &probe_mesSA_ptFtk );
  tree->SetBranchAddress( "probe_mesSA_etaFtk",  &probe_mesSA_etaFtk );
  tree->SetBranchAddress( "probe_mesSA_phiFtk",  &probe_mesSA_phiFtk );
  tree->SetBranchAddress( "probe_mesSA_roadAlgo",  &probe_mesSA_roadAlgo );
  tree->SetBranchAddress( "probe_mesCB_pass",   &probe_mesCB_pass );
  tree->SetBranchAddress( "probe_mesCB_dR",     &probe_mesCB_dR );
  tree->SetBranchAddress( "probe_mesCB_tpdR",   &probe_mesCB_tpdR );
  tree->SetBranchAddress( "probe_mesCB_pt",     &probe_mesCB_pt );
  tree->SetBranchAddress( "probe_mesCB_eta",    &probe_mesCB_eta );
  tree->SetBranchAddress( "probe_mesCB_phi",    &probe_mesCB_phi );


  Long64_t nentries_d = tree->GetEntriesFast();
  Int_t nbytes = 0, nb = 0;

  cout << "nentries_d: " << nentries_d << endl;

  TH1D* h_EventNumber = new TH1D("h_EventNumber","",10000,5000,10000);
  TH1D* h_mass = new TH1D("h_mass","",100,2600,3600);
  TH1D* h_tag_pt   = new TH1D("h_tag_pt",  "",100,0,500);
  TH1D* h_probe_pt = new TH1D("h_probe_pt","",100,0,500);

  TH2F* hh_pt  = new TH2F("hh_pt", "hh_pt;Offline pt [GeV];FTK pT [GeV]", 100, 0,    20,  100, 0,    20);
  TH2F* hh_eta = new TH2F("hh_eta", "hh_eta;Offline eta;FTK eta", 100, -3,   3,   100, -3,   3);
  TH2F* hh_phi = new TH2F("hh_phi", "hh_phi;Offline phi;FTK phi", 100, -3.5, 3.5, 100, -3.5, 3.5);

  TH2F* hh_ptl2  = new TH2F("hh_ptl2", "hh_ptl2;Offline pT [GeV];L2MuonSA pT [GeV]", 100, 0,    20,  100, 0,    20);
  TH2F* hh_etal2 = new TH2F("hh_etal2", "hh_etal2;Offline eta;L2MuonSA eta", 100, -3,   3,   100, -3,   3);
  TH2F* hh_phil2 = new TH2F("hh_phil2", "hh_phil2;Offline phi;L2MuonSA phi", 100, -3.5, 3.5, 100, -3.5, 3.5);

  TH2F* hh_ptftk  = new TH2F("hh_ptftk", "hh_ptftk;FTK pT [GeV];L2MuonSA pT [GeV]", 100, 0,    20,  100, 0,    20);
  TH2F* hh_etaftk = new TH2F("hh_etaftk", "hh_etaftk;FTK eta;L2MuonSA eta", 100, -3,   3,   100, -3,   3);
  TH2F* hh_phiftk = new TH2F("hh_phiftk", "hh_phiftk;FTK phi;L2MuonSA phi", 100, -3.5, 3.5, 100, -3.5, 3.5);

  TH2F* hh_pt_12   = new TH2F("hh_pt_12",   "hh_pt_12;mu1 pt [GeV];mu2 pT [GeV]",       100, 0, 100, 100, 0, 100);
  TH2F* hh_pt_12SA = new TH2F("hh_pt_12SA", "hh_pt_12;mu1 SA pt [GeV];mu2 SA pT [GeV]", 100, 0, 100, 100, 0, 100);

  const int nbin_dR = 50;
  const double xmin_dR = 0;
  const double xmax_dR = 1;

  TH1D* h0_br_dR = new TH1D("h0_br_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h1_br_dR = new TH1D("h1_br_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h2_br_dR = new TH1D("h2_br_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h3_br_dR = new TH1D("h3_br_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h4_br_dR = new TH1D("h4_br_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h5_br_dR = new TH1D("h5_br_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h6_br_dR = new TH1D("h6_br_dR","",nbin_dR,xmin_dR,xmax_dR);

  TH1D* h_br_dR_ftk = new TH1D("h_br_dR_ftk",";#DeltaR_{tag,probe}^{ext};",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h_br_dR_ftkpass = new TH1D("h_br_dR_ftkpass","",nbin_dR,xmin_dR,xmax_dR);

  TH1D* h0_ec_dR = new TH1D("h0_ec_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h1_ec_dR = new TH1D("h1_ec_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h2_ec_dR = new TH1D("h2_ec_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h3_ec_dR = new TH1D("h3_ec_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h4_ec_dR = new TH1D("h4_ec_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h5_ec_dR = new TH1D("h5_ec_dR","",nbin_dR,xmin_dR,xmax_dR);
  TH1D* h6_ec_dR = new TH1D("h6_ec_dR","",nbin_dR,xmin_dR,xmax_dR);

  h0_br_dR->SetFillColor(col0_dR);
  h1_br_dR->SetFillColor(col1_dR);
  h2_br_dR->SetFillColor(col2_dR);
  h3_br_dR->SetFillColor(col3_dR);

  h0_ec_dR->SetFillColor(col0_dR);
  h1_ec_dR->SetFillColor(col1_dR);
  h2_ec_dR->SetFillColor(col2_dR);
  h3_ec_dR->SetFillColor(col3_dR);

  int index = -1;
  if(TrigChain=="mu4") index = 0;
  else if(TrigChain=="mu6") index = 1;

  if(index <0) {
    cout << "unknown trigger chain: " << TrigChain << endl;
    return;
  }


  double current_mass = -999;
  double previous_mass = 999;

  int evt_count = 0;

  double ZERO = 1E-10;

  double even_L1_eta = -999;
  double even_L1_phi = -999;
  double even_SA_pt = -999;
  double even_SA_eta = -999;
  double even_SA_phi = -999;
  double odd_L1_eta = -999;
  double odd_L1_phi = -999;
  double odd_SA_pt = -999;
  double odd_SA_eta = -999;
  double odd_SA_phi = -999;

  double even_SA_ptFtk = -999;
  double even_SA_etaFtk = -999;
  double even_SA_phiFtk = -999;
  double odd_SA_ptFtk = -999;
  double odd_SA_etaFtk = -999;
  double odd_SA_phiFtk = -999;

  bool sameSA = false;
  bool sameL1 = false;
  bool sameFtk = false;

  bool even_pass_L1 = false;
  bool even_pass_SA = false;
  bool even_pass_CB = false;
  bool odd_pass_L1 = false;
  bool odd_pass_SA = false;
  bool odd_pass_CB = false;

  int count;
  for (Long64_t jentry=0; jentry<nentries_d;jentry++) {
    count++;
    nb = tree->GetEntry(jentry);   nbytes += nb;

    if(mes_name->at(index) != TrigChain) {
      cout << "inconsistent trigger chain " << mes_name->at(index) << " " << TrigChain << endl;
      return;
    }

    cout << "" << endl;
    cout << "" << endl;
    cout << "" << endl;
    cout << "EventNumber1: " << EventNumber << "/" << count << "/" << evt_count << endl;
    cout << "probe_pt/eta/phi: " << probe_pt/1000. << "/" << probe_eta << "/" << probe_phi << endl;
    cout << "probe_mesSA_pt/eta/phi: " << probe_mesSA_pt->at(index) << "/" << probe_mesSA_eta->at(index) << "/" << probe_mesSA_phi->at(index) << endl;

    h_EventNumber->Fill(EventNumber);
    h_mass->Fill(tp_mass);

    cout << "pt size: " << probe_mesSA_pt->size() << endl;
    cout << "ptFtk size: " << probe_mesSA_ptFtk->size() << endl;

    //if (probe_mesSA_ptFtk->size()>0 && probe_mesSA_sAddress->at(index)>-1 && probe_mesSA_roadAlgo->at(index)==2)
    if (probe_mesSA_ptFtk->size()>0 && probe_mesSA_sAddress->at(index)>-1 ){
      cout << "roadAlgo: " << probe_mesSA_roadAlgo->at(index) << endl;
      cout << "ptFtk/etaFtk/phiFtk: " << probe_mesSA_ptFtk->at(index)/1000. << "/" << probe_mesSA_etaFtk->at(index) << "/" << probe_mesSA_phiFtk->at(index) << endl;
      hh_pt->Fill(probe_pt/1000., probe_mesSA_ptFtk->at(index)/1000.);
      hh_eta->Fill(probe_eta,probe_mesSA_etaFtk->at(index));
      hh_phi->Fill(probe_phi,probe_mesSA_phiFtk->at(index));
      hh_ptl2->Fill(probe_pt/1000., abs(probe_mesSA_pt->at(index)));
      hh_etal2->Fill(probe_eta,probe_mesSA_eta->at(index));
      hh_phil2->Fill(probe_phi,probe_mesSA_phi->at(index));
      hh_ptftk->Fill(probe_mesSA_ptFtk->at(index)/1000., abs(probe_mesSA_pt->at(index)));
      hh_etaftk->Fill(probe_mesSA_etaFtk->at(index),probe_mesSA_eta->at(index));
      hh_phiftk->Fill(probe_mesSA_phiFtk->at(index),probe_mesSA_phi->at(index));

    }



    if(3099. - 300. > tp_mass || tp_mass > 3099. + 300) continue;

    //cout << "EventNumber2: " << EventNumber << "/" << count << "/" << evt_count << endl;

    h_tag_pt  ->Fill(tag_pt*0.001);
    h_probe_pt->Fill(probe_pt*0.001);
    hh_pt_12->Fill(tag_pt/1000.,probe_pt/1000);

    // pt cut
    if(tag_pt*0.001<4. || probe_pt*0.001<4) continue;

    //cout << "EventNumber3: " << EventNumber << "/" << count << "/" << evt_count << endl;


    current_mass = tp_mass;
    if(fabs(current_mass-previous_mass) > ZERO &&
        evt_count%2) {
      cout << "odd event should have same mass as even event: " << current_mass << " " << previous_mass << endl;
      return;
    }

    if(evt_count%2 == 0) {
      even_L1_eta = probe_mesL1_eta->at(index);
      even_L1_phi = probe_mesL1_phi->at(index);
      even_SA_pt  = probe_mesSA_pt->at(index);
      even_SA_eta = probe_mesSA_eta->at(index);
      even_SA_phi = probe_mesSA_phi->at(index);

      even_pass_L1 = (probe_mesL1_pass->at(index) > 0);
      even_pass_SA = (probe_mesSA_pass->at(index) > 0);
      even_pass_CB = (probe_mesCB_pass->at(index) > 0);

      even_SA_ptFtk = probe_mesSA_ptFtk->at(index);
      even_SA_etaFtk = probe_mesSA_etaFtk->at(index);
      even_SA_phiFtk = probe_mesSA_phiFtk->at(index);

      // cout << "even" << evt_count << " " << even_SA_eta << " " << even_SA_phi << endl;
      hh_pt_12SA->Fill(odd_SA_pt,even_SA_pt);
    }

    if(evt_count%2 == 1) {
      odd_L1_eta = probe_mesL1_eta->at(index);
      odd_L1_phi = probe_mesL1_phi->at(index);
      odd_SA_pt  = probe_mesSA_pt->at(index);
      odd_SA_eta = probe_mesSA_eta->at(index);
      odd_SA_phi = probe_mesSA_phi->at(index);

      odd_pass_L1 = (probe_mesL1_pass->at(index) > 0);
      odd_pass_SA = (probe_mesSA_pass->at(index) > 0);
      odd_pass_CB = (probe_mesCB_pass->at(index) > 0);

      odd_SA_ptFtk = probe_mesSA_ptFtk->at(index);
      odd_SA_etaFtk = probe_mesSA_etaFtk->at(index);
      odd_SA_phiFtk = probe_mesSA_phiFtk->at(index);

      // cout << "odd " << evt_count << " " << odd_SA_eta << " " << odd_SA_phi << endl;
      hh_pt_12SA->Fill(odd_SA_pt,even_SA_pt);
    }

    if(evt_count%2 == 1) {

      sameL1 = false;
      if(fabs(even_L1_eta - odd_L1_eta) < ZERO &&
          fabs(even_L1_phi - odd_L1_phi) < ZERO)
        sameL1 = true;

      sameSA = false;
      if(fabs(even_SA_pt  - odd_SA_pt) < ZERO &&
          fabs(even_SA_eta - odd_SA_eta) < ZERO &&
          fabs(even_SA_phi - odd_SA_phi) < ZERO)
        sameSA = true;

      sameFtk = false;
      if(fabs(even_SA_ptFtk  - odd_SA_ptFtk) < ZERO &&
          fabs(even_SA_etaFtk - odd_SA_etaFtk) < ZERO &&
          fabs(even_SA_phiFtk - odd_SA_phiFtk) < ZERO)
        sameFtk = true;


      if(fabs(even_L1_eta)<1.05) {

        h0_br_dR->Fill(tp_extdR);
        if(even_pass_L1) h1_br_dR->Fill(tp_extdR);
        if(even_pass_L1 && !sameL1) h5_br_dR->Fill(tp_extdR);
        if(even_pass_L1 && even_pass_SA) h2_br_dR->Fill(tp_extdR);
        if(even_pass_L1 && even_pass_SA && !sameSA) h4_br_dR->Fill(tp_extdR);
        if(even_pass_L1 && !sameL1 && even_pass_SA && !sameSA) h6_br_dR->Fill(tp_extdR);
        if(even_pass_L1 && even_pass_SA && even_pass_CB) h3_br_dR->Fill(tp_extdR);

        if(even_pass_L1 && !sameFtk) h_br_dR_ftk->Fill(tp_extdR);
        if(even_pass_L1 && even_pass_SA && !sameFtk) h_br_dR_ftkpass->Fill(tp_extdR);

      }
      else {

        h0_ec_dR->Fill(tp_extdR);
        if(even_pass_L1) h1_ec_dR->Fill(tp_extdR);
        if(even_pass_L1 && !sameL1) h5_ec_dR->Fill(tp_extdR);
        if(even_pass_L1 && even_pass_SA) h2_ec_dR->Fill(tp_extdR);
        if(even_pass_L1 && even_pass_SA && !sameSA) h4_ec_dR->Fill(tp_extdR);
        if(even_pass_L1 && !sameL1 && even_pass_SA && !sameSA) h6_ec_dR->Fill(tp_extdR);
        if(even_pass_L1 && even_pass_SA && even_pass_CB) h3_ec_dR->Fill(tp_extdR);

      }

      if(fabs(odd_L1_eta)<1.05) {

        h0_br_dR->Fill(tp_extdR);
        if(odd_pass_L1) h1_br_dR->Fill(tp_extdR);
        if(odd_pass_L1 && !sameL1) h5_br_dR->Fill(tp_extdR);
        if(odd_pass_L1 && odd_pass_SA) h2_br_dR->Fill(tp_extdR);
        if(odd_pass_L1 && odd_pass_SA && !sameSA) h4_br_dR->Fill(tp_extdR);
        if(odd_pass_L1 && !sameL1 && odd_pass_SA && !sameSA) h6_br_dR->Fill(tp_extdR);
        if(odd_pass_L1 && odd_pass_SA && odd_pass_CB) h3_br_dR->Fill(tp_extdR);

        if(odd_pass_L1 && !sameFtk) h_br_dR_ftk->Fill(tp_extdR);
        if(odd_pass_L1 && odd_pass_SA && !sameFtk) h_br_dR_ftkpass->Fill(tp_extdR);

      }
      else {

        h0_ec_dR->Fill(tp_extdR);
        if(odd_pass_L1) h1_ec_dR->Fill(tp_extdR);
        if(odd_pass_L1 && !sameL1) h5_ec_dR->Fill(tp_extdR);
        if(odd_pass_L1 && odd_pass_SA) h2_ec_dR->Fill(tp_extdR);
        if(odd_pass_L1 && odd_pass_SA && !sameSA) h4_ec_dR->Fill(tp_extdR);
        if(odd_pass_L1 && !sameL1 && odd_pass_SA && !sameSA) h6_ec_dR->Fill(tp_extdR);
        if(odd_pass_L1 && odd_pass_SA && odd_pass_CB) h3_ec_dR->Fill(tp_extdR);

      }

    }

    previous_mass = current_mass;
    evt_count++;
  }


  TEfficiency* eff01_br_dR = new TEfficiency(*h1_br_dR,*h0_br_dR);
  TEfficiency* eff12_br_dR = new TEfficiency(*h2_br_dR,*h1_br_dR);
  TEfficiency* eff23_br_dR = new TEfficiency(*h3_br_dR,*h2_br_dR);
  TEfficiency* eff14_br_dR = new TEfficiency(*h4_br_dR,*h1_br_dR);
  TEfficiency* eff05_br_dR = new TEfficiency(*h5_br_dR,*h0_br_dR);
  TEfficiency* eff56_br_dR = new TEfficiency(*h6_br_dR,*h5_br_dR);

  TEfficiency* eff_br_dR_ftk = new TEfficiency(*h_br_dR_ftk,*h1_br_dR);
  TEfficiency* eff_br_dR_ftkpass = new TEfficiency(*h_br_dR_ftkpass,*h1_br_dR);

  TEfficiency* eff01_ec_dR = new TEfficiency(*h1_ec_dR,*h0_ec_dR);
  TEfficiency* eff12_ec_dR = new TEfficiency(*h2_ec_dR,*h1_ec_dR);
  TEfficiency* eff23_ec_dR = new TEfficiency(*h3_ec_dR,*h2_ec_dR);
  TEfficiency* eff14_ec_dR = new TEfficiency(*h4_ec_dR,*h1_ec_dR);
  TEfficiency* eff05_ec_dR = new TEfficiency(*h5_ec_dR,*h0_ec_dR);
  TEfficiency* eff56_ec_dR = new TEfficiency(*h6_ec_dR,*h5_ec_dR);

  eff01_br_dR->SetLineColor(col01_dR);
  eff12_br_dR->SetLineColor(col12_dR);
  eff23_br_dR->SetLineColor(col23_dR);
  eff14_br_dR->SetLineColor(col14_dR);
  eff05_br_dR->SetLineColor(col05_dR);
  eff56_br_dR->SetLineColor(col56_dR);

  eff_br_dR_ftk->SetLineColor(kGreen+2);
  eff_br_dR_ftkpass->SetLineColor(kRed);
  eff_br_dR_ftk->SetMarkerColor(kGreen+2);
  eff_br_dR_ftkpass->SetMarkerColor(kRed);

  eff01_ec_dR->SetLineColor(col01_dR);
  eff12_ec_dR->SetLineColor(col12_dR);
  eff23_ec_dR->SetLineColor(col23_dR);
  eff14_ec_dR->SetLineColor(col14_dR);
  eff05_ec_dR->SetLineColor(col05_dR);
  eff56_ec_dR->SetLineColor(col56_dR);

  eff01_br_dR->SetMarkerColor(col01_dR);
  eff12_br_dR->SetMarkerColor(col12_dR);
  eff23_br_dR->SetMarkerColor(col23_dR);
  eff14_br_dR->SetMarkerColor(col14_dR);
  eff05_br_dR->SetMarkerColor(col05_dR);
  eff56_br_dR->SetMarkerColor(col56_dR);

  eff01_ec_dR->SetMarkerColor(col01_dR);
  eff12_ec_dR->SetMarkerColor(col12_dR);
  eff23_ec_dR->SetMarkerColor(col23_dR);
  eff14_ec_dR->SetMarkerColor(col14_dR);
  eff05_ec_dR->SetMarkerColor(col05_dR);
  eff56_ec_dR->SetMarkerColor(col56_dR);

  eff05_br_dR->SetMarkerSize(1.1);
  eff05_ec_dR->SetMarkerSize(1.1);

  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1020, 700);
  c1->SetGrid();
  c1->SetRightMargin(0.20);
  c1->SetLeftMargin(0.23);
  c1->SetBottomMargin(0.20);
  c1 -> Print( pdf + "[", "pdf" );

  h_EventNumber->Draw();
  c1->Print(pdf , "pdf");
  c1->Clear();

  h_mass->Draw();
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_pt_12->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_pt_12SA->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_pt->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_eta->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_phi->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_ptl2->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_etal2->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_phil2->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_ptftk->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_etaftk->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  hh_phiftk->Draw("colz");
  c1->Print(pdf , "pdf");
  c1->Clear();

  h1_br_dR->Draw("");
  h1_br_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  h_br_dR_ftk->SetFillColor(kRed);
  h_br_dR_ftk->Draw("same");
  h_br_dR_ftkpass->SetFillColor(kGreen+2);
  h_br_dR_ftkpass->Draw("same");
  c1->Print(pdf , "pdf");
  c1->Clear();

  //h_probe_pt->SetFillColor(col0_dR);
  //h_probe_pt->Draw();
  //h_tag_pt->Draw("same");
  //c1->Print(pdf , "pdf");
  //c1->Clear();

  TCanvas* c3 = new TCanvas("c3");
  TH1F* frame3 = gPad->DrawFrame(xmin_dR,0,xmax_dR,1);
  frame3->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  frame3->GetYaxis()->SetTitle("efficiency");
  eff_br_dR_ftk->Draw("psame");
  eff_br_dR_ftkpass->Draw("psame");
  //eff14_br_dR->Draw("psame");
  //eff12_br_dR->Draw("psame");
  myText(0.65,0.55,1, "J/#psi MC");
  myText(0.65,0.49,1, "MU4");
  // myText(0.65,0.37,1, "\"JpsitapMC\"");
  myText(0.65,0.43,1, "\"NoTag\"");
  myText(0.65,0.37,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.31,1, "probe |#it{#eta}^{L1}| < 1.05");

  myText(0.23,0.37,kRed, "(FTK separated && SA)/L1");
  myText(0.23,0.31,kGreen+2, "(FTK separated)/L1");
  //myText(0.23,0.25,col12_dR, "(SA)/L1");
  //myText(0.23,0.19,col14_dR, "(SA && separated)/L1");
  // myText(0.23,0.25,col56_dR, "(SA && SA_separated)/(L1 && L1_separated)");
  c3->Print(pdf , "pdf");


  TCanvas* c1_br_dR = new TCanvas("c1_br_dR");
  h0_br_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  h0_br_dR->GetYaxis()->SetTitle("");
  h0_br_dR->Draw("hist");
  h1_br_dR->Draw("histsame");
  h2_br_dR->Draw("histsame");
  h3_br_dR->Draw("histsame");
  myText(0.65,0.86,1, "J/#psi MC");
  myText(0.65,0.80,1, "MU4");
  // myText(0.65,0.80,1, "\"JpsitapMC\"");
  myText(0.65,0.74,1, "\"NoTag\"");
  myText(0.65,0.68,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.62,1, "probe |#it{#eta}^{L1}| < 1.05");

  myText(0.23,0.86,col0_dR, "Probe");
  myText(0.23,0.80,col1_dR, "L1");
  myText(0.23,0.74,col2_dR, "SA");
  myText(0.23,0.68,col3_dR, "CB");

  if(epsdir!=""){
    //c1_br_dR->SaveAs(epsdir+"/dR_barrel.eps");
    c1_br_dR->Print(pdf , "pdf");
  }

  TCanvas* c2_br_dR = new TCanvas("c2_br_dR");
  TH1F* frame_br_dR = gPad->DrawFrame(xmin_dR,0,xmax_dR,1);
  frame_br_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  frame_br_dR->GetYaxis()->SetTitle("efficiency");
  eff01_br_dR->Draw("psame");
  eff12_br_dR->Draw("psame");
  eff23_br_dR->Draw("psame");
  myText(0.65,0.49,1, "J/#psi MC");
  myText(0.65,0.43,1, "MU4");
  // myText(0.65,0.37,1, "\"JpsitapMC\"");
  myText(0.65,0.37,1, "\"NoTag\"");
  myText(0.65,0.31,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.25,1, "probe |#it{#eta}^{L1}| < 1.05");

  myText(0.23,0.37,col01_dR, "L1");
  myText(0.23,0.31,col12_dR, "SA/L1");
  myText(0.23,0.25,col23_dR, "CB/SA");

  if(epsdir!=""){
    //c2_br_dR->SaveAs(epsdir+"/dR_barrel_ratio.eps");
    c2_br_dR->Print(pdf , "pdf");
  }

  TCanvas* c3_br_dR = new TCanvas("c3_br_dR");
  TH1F* frame2_br_dR = gPad->DrawFrame(xmin_dR,0,xmax_dR,1);
  frame2_br_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  frame2_br_dR->GetYaxis()->SetTitle("efficiency");
  eff12_br_dR->Draw("psame");
  eff14_br_dR->Draw("psame");
  // eff56_br_dR->Draw("psame");
  myText(0.65,0.55,1, "J/#psi MC");
  myText(0.65,0.49,1, "MU4");
  // myText(0.65,0.37,1, "\"JpsitapMC\"");
  myText(0.65,0.43,1, "\"NoTag\"");
  myText(0.65,0.37,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.31,1, "probe |#it{#eta}^{L1}| < 1.05");

  myText(0.23,0.37,col12_dR, "SA/L1");
  myText(0.23,0.31,col14_dR, "(SA && separated)/L1");
  // myText(0.23,0.25,col56_dR, "(SA && SA_separated)/(L1 && L1_separated)");

  if(epsdir!=""){
    //c3_br_dR->SaveAs(epsdir+"/dR_barrel_ratio2.eps");
    c3_br_dR->Print(pdf , "pdf");
  }


  TCanvas* c4_br_dR = new TCanvas("c4_br_dR");
  TH1F* frame3_br_dR = gPad->DrawFrame(xmin_dR,0,xmax_dR,1);
  frame3_br_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  frame3_br_dR->GetYaxis()->SetTitle("efficiency");
  eff01_br_dR->Draw("psame");
  eff05_br_dR->Draw("psame");
  myText(0.65,0.49,1, "J/#psi MC");
  myText(0.65,0.43,1, "MU4");
  // myText(0.65,0.37,1, "\"JpsitapMC\"");
  myText(0.65,0.37,1, "\"NoTag\"");
  myText(0.65,0.31,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.25,1, "probe |#it{#eta}^{L1}| < 1.05");

  myText(0.23,0.31,col01_dR, "L1");
  myText(0.23,0.25,col05_dR, "L1 && L1_separated");

  if(epsdir!=""){
    //c4_br_dR->SaveAs(epsdir+"/dR_barrel_ratio3.eps");
    c4_br_dR->Print(pdf , "pdf");
  }



  TCanvas* c1_ec_dR = new TCanvas("c1_ec_dR");
  h0_ec_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  h0_ec_dR->GetYaxis()->SetTitle("");
  h0_ec_dR->Draw("hist");
  h1_ec_dR->Draw("histsame");
  h2_ec_dR->Draw("histsame");
  h3_ec_dR->Draw("histsame");
  myText(0.65,0.86,1, "J/#psi MC");
  myText(0.65,0.80,1, "MU4");
  // myText(0.65,0.80,1, "\"JpsitapMC\"");
  myText(0.65,0.74,1, "\"NoTag\"");
  myText(0.65,0.68,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.62,1, "probe |#it{#eta}^{L1}| > 1.05");

  myText(0.23,0.86,col0_dR, "Probe");
  myText(0.23,0.80,col1_dR, "L1");
  myText(0.23,0.74,col2_dR, "SA");
  myText(0.23,0.68,col3_dR, "CB");

  if(epsdir!=""){
    c1_ec_dR->Print(pdf , "pdf");
    //c1_ec_dR->SaveAs(epsdir+"/dR_endcap.eps");
  }

  TCanvas* c2_ec_dR = new TCanvas("c2_ec_dR");
  TH1F* frame_ec_dR = gPad->DrawFrame(xmin_dR,0,xmax_dR,1);
  frame_ec_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  frame_ec_dR->GetYaxis()->SetTitle("efficiency");
  eff01_ec_dR->Draw("psame");
  eff12_ec_dR->Draw("psame");
  eff23_ec_dR->Draw("psame");
  myText(0.65,0.49,1, "J/#psi MC");
  myText(0.65,0.43,1, "MU4");
  // myText(0.65,0.37,1, "\"JpsitapMC\"");
  myText(0.65,0.37,1, "\"NoTag\"");
  myText(0.65,0.31,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.25,1, "probe |#it{#eta}^{L1}| > 1.05");

  myText(0.23,0.37,col01_dR, "L1");
  myText(0.23,0.31,col12_dR, "SA/L1");
  myText(0.23,0.25,col23_dR, "CB/SA");

  if(epsdir!=""){
    //c2_ec_dR->SaveAs(epsdir+"/dR_endcap_ratio.eps");
    c2_ec_dR->Print(pdf , "pdf");
  }

  TCanvas* c3_ec_dR = new TCanvas("c3_ec_dR");
  TH1F* frame2_ec_dR = gPad->DrawFrame(xmin_dR,0,xmax_dR,1);
  frame2_ec_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  frame2_ec_dR->GetYaxis()->SetTitle("efficiency");
  eff12_ec_dR->Draw("psame");
  eff14_ec_dR->Draw("psame");
  // eff56_ec_dR->Draw("psame");
  myText(0.65,0.55,1, "J/#psi MC");
  myText(0.65,0.49,1, "MU4");
  // myText(0.65,0.37,1, "\"JpsitapMC\"");
  myText(0.65,0.43,1, "\"NoTag\"");
  myText(0.65,0.37,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.31,1, "probe |#it{#eta}^{L1}| > 1.05");

  myText(0.23,0.37,col12_dR, "SA/L1");
  myText(0.23,0.31,col14_dR, "(SA && separated)/L1");
  // myText(0.23,0.25,col56_dR, "(SA && SA_separated)/(L1 && L1_separated)");

  if(epsdir!=""){
    //c3_ec_dR->SaveAs(epsdir+"/dR_endcap_ratio2.eps");
    c3_ec_dR->Print(pdf , "pdf");
  }

  TCanvas* c4_ec_dR = new TCanvas("c4_ec_dR");
  TH1F* frame3_ec_dR = gPad->DrawFrame(xmin_dR,0,xmax_dR,1);
  frame3_ec_dR->GetXaxis()->SetTitle("#DeltaR_{tag,probe}^{ext}");
  frame3_ec_dR->GetYaxis()->SetTitle("efficiency");
  eff01_ec_dR->Draw("psame");
  eff05_ec_dR->Draw("psame");
  myText(0.65,0.49,1, "J/#psi MC");
  myText(0.65,0.43,1, "MU4");
  // myText(0.65,0.37,1, "\"JpsitapMC\"");
  myText(0.65,0.37,1, "\"NoTag\"");
  myText(0.65,0.31,1, "#it{p}_{T}^{offline} > 4 GeV");
  myText(0.65,0.25,1, "probe |#it{#eta}^{L1}| > 1.05");

  myText(0.23,0.31,col01_dR, "L1");
  myText(0.23,0.25,col05_dR, "L1 && L1_separated");

  if(epsdir!=""){
    //c4_ec_dR->SaveAs(epsdir+"/dR_endcap_ratio3.eps");
    c4_ec_dR->Print(pdf , "pdf");
    c4_ec_dR -> Print( pdf + "]", "pdf" );
  }


  return;
}
