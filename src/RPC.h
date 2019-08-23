//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 25 13:30:35 2018 by ROOT version 6.12/04
// from TTree t_tap/TrigMuonTagAndProbe
// found on file: /home/yfukuhar/gpfs/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0.root
//////////////////////////////////////////////////////////

#ifndef RPC_h
#define RPC_h

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include <vector>
#include <iostream>


// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class RPC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   const int N4 = 1;
   const int N50 = 0;
   const int NTagProc = 3;
   // Declaration of leaf types
   Int_t EventNumber;
   Int_t RunNumber;
   Int_t LumiBlock;
   Float_t AverageInteractionsPerCrossing;
   Int_t           mes_n;
   vector<string>  *mes_name;
   Double_t        sumReqdRL1;
   Double_t        sumReqdREF;
   Int_t           tag_proc;
   Double_t        tag_ReqdRL1;
   Double_t        tag_ReqdREF;
   Double_t        tag_dRL1;
   Double_t        tag_dREF;
   Double_t        tag_pt;
   Double_t        tag_eta;
   Double_t        tag_exteta;
   Double_t        tag_extinneta;
   Double_t        tag_phi;
   Double_t        tag_extphi;
   Double_t        tag_extinnphi;
   Double_t        tag_d0;
   Double_t        tag_z0;
   Double_t        tag_charge;
   Double_t        tag_L1_pt;
   Double_t        tag_L1_eta;
   Double_t        tag_L1_phi;
   Double_t        tag_SA_pt;
   Double_t        tag_SA_eta;
   Double_t        tag_SA_phi;
   Double_t        tag_CB_pt;
   Double_t        tag_CB_eta;
   Double_t        tag_CB_phi;
   Double_t        tag_EF_pt;
   Double_t        tag_EF_eta;
   Double_t        tag_EF_phi;
   Double_t        probe_pt;
   Double_t        probe_eta;
   Double_t        probe_exteta;
   Double_t        probe_extinneta;
   Double_t        probe_phi;
   Double_t        probe_extphi;
   Double_t        probe_extinnphi;
   Double_t        probe_d0;
   Double_t        probe_z0;
   Double_t        probe_charge;
   Double_t        probe_segment_n;
   Double_t        probe_segment_x[10];
   Double_t        probe_segment_y[10];
   Double_t        probe_segment_z[10];
   Double_t        probe_segment_px[10];
   Double_t        probe_segment_py[10];
   Double_t        probe_segment_pz[10];
   Double_t        probe_segment_chiSquared[10];
   Double_t        probe_segment_numberDoF[10];
   Double_t        probe_segment_sector[10];
   Double_t        probe_segment_chamberIndex[10];
   Double_t        probe_segment_etaIndex[10];
   Double_t        probe_segment_nPrecisionHits[10];
   Double_t        probe_segment_nPhiLayers[10];
   Double_t        probe_segment_nTrigEtaLayers[10];
   Double_t        tp_dR;
   Double_t        tp_deta;
   Double_t        tp_dphi;
   Double_t        tp_extdR;
   Double_t        tp_extdeta;
   Double_t        tp_extdphi;
   Double_t        tp_extinndR;
   Double_t        tp_extinndeta;
   Double_t        tp_extinndphi;
   Double_t        tp_mass;
   Double_t        tp_vftlxy;
   Double_t        tp_vftchi2;
   Int_t           tp_vftndof;
   vector<int>     *probe_mesEFTAG_pass;
   vector<double>  *probe_mesEFTAG_dR;
   vector<double>  *probe_mesEFTAG_tpdR;
   vector<double>  *probe_mesEFTAG_pt;
   vector<double>  *probe_mesEFTAG_eta;
   vector<double>  *probe_mesEFTAG_phi;
   vector<int>     *probe_mesL1_pass;
   vector<double>  *probe_mesL1_dR;
   vector<double>  *probe_mesL1_tpdR;
   vector<double>  *probe_mesL1_pt;
   vector<double>  *probe_mesL1_eta;
   vector<double>  *probe_mesL1_phi;
   vector<int>     *probe_mesSA_pass;
   vector<double>  *probe_mesSA_dR;
   vector<double>  *probe_mesSA_tpdR;
   vector<double>  *probe_mesSA_pt;
   vector<double>  *probe_mesSA_eta;
   vector<double>  *probe_mesSA_phi;
   vector<double>  *probe_mesSA_etams;
   vector<double>  *probe_mesSA_phims;
   vector<double>  *probe_mesSA_etabe;
   vector<double>  *probe_mesSA_phibe;
   vector<double>  *probe_mesSA_tgcpt;
   vector<double>  *probe_mesSA_ptBarrelRadius;
   vector<double>  *probe_mesSA_ptBarrelSagitta;
   vector<double>  *probe_mesSA_ptEndcapAlpha;
   vector<double>  *probe_mesSA_ptEndcapBeta;
   vector<double>  *probe_mesSA_ptEndcapRadius;
   vector<double>  *probe_mesSA_ptCSC;
   vector<double>  *probe_mesSA_sAddress;
   vector<double>  *probe_mesSA_superPointR_BI;
   vector<double>  *probe_mesSA_superPointR_BM;
   vector<double>  *probe_mesSA_superPointR_BO;
   vector<double>  *probe_mesSA_superPointR_EI;
   vector<double>  *probe_mesSA_superPointR_EM;
   vector<double>  *probe_mesSA_superPointR_EO;
   vector<double>  *probe_mesSA_superPointR_EE;
   vector<double>  *probe_mesSA_superPointR_CSC;
   vector<double>  *probe_mesSA_superPointR_BEE;
   vector<double>  *probe_mesSA_superPointR_BME;
   vector<double>  *probe_mesSA_superPointZ_BI;
   vector<double>  *probe_mesSA_superPointZ_BM;
   vector<double>  *probe_mesSA_superPointZ_BO;
   vector<double>  *probe_mesSA_superPointZ_EI;
   vector<double>  *probe_mesSA_superPointZ_EM;
   vector<double>  *probe_mesSA_superPointZ_EO;
   vector<double>  *probe_mesSA_superPointZ_EE;
   vector<double>  *probe_mesSA_superPointZ_CSC;
   vector<double>  *probe_mesSA_superPointZ_BEE;
   vector<double>  *probe_mesSA_superPointZ_BME;
   vector<double>  *probe_mesSA_superPointSlope_BI;
   vector<double>  *probe_mesSA_superPointSlope_BM;
   vector<double>  *probe_mesSA_superPointSlope_BO;
   vector<double>  *probe_mesSA_superPointSlope_EI;
   vector<double>  *probe_mesSA_superPointSlope_EM;
   vector<double>  *probe_mesSA_superPointSlope_EO;
   vector<double>  *probe_mesSA_superPointSlope_EE;
   vector<double>  *probe_mesSA_superPointSlope_CSC;
   vector<double>  *probe_mesSA_superPointSlope_BEE;
   vector<double>  *probe_mesSA_superPointSlope_BME;
   vector<double>  *probe_mesSA_superPointIntercept_BI;
   vector<double>  *probe_mesSA_superPointIntercept_BM;
   vector<double>  *probe_mesSA_superPointIntercept_BO;
   vector<double>  *probe_mesSA_superPointIntercept_EI;
   vector<double>  *probe_mesSA_superPointIntercept_EM;
   vector<double>  *probe_mesSA_superPointIntercept_EO;
   vector<double>  *probe_mesSA_superPointIntercept_EE;
   vector<double>  *probe_mesSA_superPointIntercept_CSC;
   vector<double>  *probe_mesSA_superPointIntercept_BEE;
   vector<double>  *probe_mesSA_superPointIntercept_BME;
   vector<double>  *probe_mesSA_superPointChi2_BI;
   vector<double>  *probe_mesSA_superPointChi2_BM;
   vector<double>  *probe_mesSA_superPointChi2_BO;
   vector<double>  *probe_mesSA_superPointChi2_EI;
   vector<double>  *probe_mesSA_superPointChi2_EM;
   vector<double>  *probe_mesSA_superPointChi2_EO;
   vector<double>  *probe_mesSA_superPointChi2_EE;
   vector<double>  *probe_mesSA_superPointChi2_CSC;
   vector<double>  *probe_mesSA_superPointChi2_BEE;
   vector<double>  *probe_mesSA_superPointChi2_BME;
   vector<int>  *probe_mesSA_isRpcFailure;
   vector<int>  *probe_mesSA_isTgcFailure;
   vector<vector<float> > *probe_mesSA_rpcHitX;
   vector<vector<float> > *probe_mesSA_rpcHitY;
   vector<vector<float> > *probe_mesSA_rpcHitZ;
   vector<vector<double> > *probe_mesSA_rpcHitR;
   vector<vector<double> > *probe_mesSA_rpcHitEta;
   vector<vector<double> > *probe_mesSA_rpcHitPhi;
   vector<vector<double> > *probe_mesSA_rpcHitStationNumber;
   vector<vector<string> > *probe_mesSA_rpcHitStationName;
   vector<vector<int> > *probe_mesSA_mdtHitIsOutlier;
   vector<vector<int> > *probe_mesSA_mdtHitChamber;
   vector<vector<float> > *probe_mesSA_mdtHitR;
   vector<vector<float> > *probe_mesSA_mdtHitZ;
   vector<vector<float> > *probe_mesSA_mdtHitPhi;
   vector<vector<float> > *probe_mesSA_mdtHitResidual;

   vector < vector < float > > *probe_mesSA_roadAw;
   vector < vector < float > > *probe_mesSA_roadBw;
   vector < vector < float > > *probe_mesSA_zMin;
   vector < vector < float > > *probe_mesSA_zMax;
   vector < vector < float > > *probe_mesSA_rMin;
   vector < vector < float > > *probe_mesSA_rMax;
   vector < vector < float > > *probe_mesSA_etaMin;
   vector < vector < float > > *probe_mesSA_etaMax;

   vector < vector < double > >* probe_mesL2_pass;
   vector < vector < double > >* probe_mesL2_pt;
   vector < vector < double > >* probe_mesL2_eta;
   vector < vector < double > >* probe_mesL2_phi;

   vector < vector < double > >* probe_mesL2_ptFtk;
   vector < vector < double > >* probe_mesL2_etaFtk;
   vector < vector < double > >* probe_mesL2_phiFtk;
   vector < vector < double > >* probe_mesL2_roadAlgo;

   vector < vector < double > >* probe_mesL2_passSA;
   vector < vector < double > >* probe_mesL2_numSpSA;
   vector < vector < double > >* probe_mesL2_ptSA;
   vector < vector < double > >* probe_mesL2_etaSA;
   vector < vector < double > >* probe_mesL2_phiSA;
   vector < vector < double > >* probe_mesL2_sAddressSA;


   vector < vector < float > >* probe_mesSA_roadFtkAw;
   vector < vector < float > >* probe_mesSA_roadFtkBw;
   vector < vector < float > >* probe_mesSA_extFtkEta;
   vector < vector < float > >* probe_mesSA_extFtkPhi;
   vector < vector < float > >* probe_mesSA_extFtkR;
   vector < vector < float > >* probe_mesSA_extFtkZ;
   vector < vector < float > >* probe_mesSA_roadRpcAw;
   vector < vector < float > >* probe_mesSA_roadRpcBw;
   vector < int >* probe_mesSA_roadAlgo;
   vector < double >* probe_mesSA_ptFtk;
   vector < double >* probe_mesSA_etaFtk;
   vector < double >* probe_mesSA_phiFtk;

   vector <float> *probe_mesSA_roiEta;
   vector <float> *probe_mesSA_roiPhi;

   vector<int>     *probe_mesCB_pass;
   vector<double>  *probe_mesCB_dR;
   vector<double>  *probe_mesCB_tpdR;
   vector<double>  *probe_mesCB_pt;
   vector<double>  *probe_mesCB_eta;
   vector<double>  *probe_mesCB_phi;
   vector<int>     *probe_mesFTF_pass;
   vector<double>  *probe_mesFTF_dR;
   vector<double>  *probe_mesFTF_pt;
   vector<double>  *probe_mesFTF_eta;
   vector<double>  *probe_mesFTF_phi;
   vector<int>     *probe_mesEF_pass;
   vector<double>  *probe_mesEF_dR;
   vector<double>  *probe_mesEF_tpdR;
   vector<double>  *probe_mesEF_pt;
   vector<double>  *probe_mesEF_eta;
   vector<double>  *probe_mesEF_phi;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_AverageInteractionsPerCrossing;   //!
   TBranch        *b_mes_n;   //!
   TBranch        *b_mes_name;   //!
   TBranch        *b_sumReqdRL1;   //!
   TBranch        *b_sumReqdREF;   //!
   TBranch        *b_tag_proc;   //!
   TBranch        *b_tag_ReqdRL1;   //!
   TBranch        *b_tag_ReqdREF;   //!
   TBranch        *b_tag_dRL1;   //!
   TBranch        *b_tag_dREF;   //!
   TBranch        *b_tag_pt;   //!
   TBranch        *b_tag_eta;   //!
   TBranch        *b_tag_exteta;   //!
   TBranch        *b_tag_extinneta;   //!
   TBranch        *b_tag_phi;   //!
   TBranch        *b_tag_extphi;   //!
   TBranch        *b_tag_extinnphi;   //!
   TBranch        *b_tag_d0;   //!
   TBranch        *b_tag_z0;   //!
   TBranch        *b_tag_charge;   //!
   TBranch        *b_tag_L1_pt;   //!
   TBranch        *b_tag_L1_eta;   //!
   TBranch        *b_tag_L1_phi;   //!
   TBranch        *b_tag_SA_pt;   //!
   TBranch        *b_tag_SA_eta;   //!
   TBranch        *b_tag_SA_phi;   //!
   TBranch        *b_tag_CB_pt;   //!
   TBranch        *b_tag_CB_eta;   //!
   TBranch        *b_tag_CB_phi;   //!
   TBranch        *b_tag_EF_pt;   //!
   TBranch        *b_tag_EF_eta;   //!
   TBranch        *b_tag_EF_phi;   //!
   TBranch        *b_probe_pt;   //!
   TBranch        *b_probe_eta;   //!
   TBranch        *b_probe_exteta;   //!
   TBranch        *b_probe_extinneta;   //!
   TBranch        *b_probe_phi;   //!
   TBranch        *b_probe_extphi;   //!
   TBranch        *b_probe_extinnphi;   //!
   TBranch        *b_probe_d0;   //!
   TBranch        *b_probe_z0;   //!
   TBranch        *b_probe_charge;   //!
   TBranch        *b_probe_segment_n;   //!
   TBranch        *b_probe_segment_x;   //!
   TBranch        *b_probe_segment_y;   //!
   TBranch        *b_probe_segment_z;   //!
   TBranch        *b_probe_segment_px;   //!
   TBranch        *b_probe_segment_py;   //!
   TBranch        *b_probe_segment_pz;   //!
   TBranch        *b_probe_segment_chiSquared;   //!
   TBranch        *b_probe_segment_numberDoF;   //!
   TBranch        *b_probe_segment_sector;   //!
   TBranch        *b_probe_segment_chamberIndex;   //!
   TBranch        *b_probe_segment_etaIndex;   //!
   TBranch        *b_probe_segment_nPrecisionHits;   //!
   TBranch        *b_probe_segment_nPhiLayers;   //!
   TBranch        *b_probe_segment_nTrigEtaLayers;   //!
   TBranch        *b_tp_dR;   //!
   TBranch        *b_tp_deta;   //!
   TBranch        *b_tp_dphi;   //!
   TBranch        *b_tp_extdR;   //!
   TBranch        *b_tp_extdeta;   //!
   TBranch        *b_tp_extdphi;   //!
   TBranch        *b_tp_extinndR;   //!
   TBranch        *b_tp_extinndeta;   //!
   TBranch        *b_tp_extinndphi;   //!
   TBranch        *b_tp_mass;   //!
   TBranch        *b_tp_vftlxy;   //!
   TBranch        *b_tp_vftchi2;   //!
   TBranch        *b_tp_vftndof;   //!
   TBranch        *b_probe_mesEFTAG_pass;   //!
   TBranch        *b_probe_mesEFTAG_dR;   //!
   TBranch        *b_probe_mesEFTAG_tpdR;   //!
   TBranch        *b_probe_mesEFTAG_pt;   //!
   TBranch        *b_probe_mesEFTAG_eta;   //!
   TBranch        *b_probe_mesEFTAG_phi;   //!
   TBranch        *b_probe_mesL1_pass;   //!
   TBranch        *b_probe_mesL1_dR;   //!
   TBranch        *b_probe_mesL1_tpdR;   //!
   TBranch        *b_probe_mesL1_pt;   //!
   TBranch        *b_probe_mesL1_eta;   //!
   TBranch        *b_probe_mesL1_phi;   //!
   TBranch        *b_probe_mesSA_pass;   //!
   TBranch        *b_probe_mesSA_dR;   //!
   TBranch        *b_probe_mesSA_tpdR;   //!
   TBranch        *b_probe_mesSA_pt;   //!
   TBranch        *b_probe_mesSA_eta;   //!
   TBranch        *b_probe_mesSA_phi;   //!
   TBranch        *b_probe_mesSA_etams;   //!
   TBranch        *b_probe_mesSA_phims;   //!
   TBranch        *b_probe_mesSA_etabe;   //!
   TBranch        *b_probe_mesSA_phibe;   //!
   TBranch        *b_probe_mesSA_tgcpt;   //!
   TBranch        *b_probe_mesSA_ptBarrelRadius;   //!
   TBranch        *b_probe_mesSA_ptBarrelSagitta;   //!
   TBranch        *b_probe_mesSA_ptEndcapAlpha;   //!
   TBranch        *b_probe_mesSA_ptEndcapBeta;   //!
   TBranch        *b_probe_mesSA_ptEndcapRadius;   //!
   TBranch        *b_probe_mesSA_ptCSC;   //!
   TBranch        *b_probe_mesSA_sAddress;   //!
   TBranch        *b_probe_mesSA_superPointR_BI;   //!
   TBranch        *b_probe_mesSA_superPointR_BM;   //!
   TBranch        *b_probe_mesSA_superPointR_BO;   //!
   TBranch        *b_probe_mesSA_superPointR_EI;   //!
   TBranch        *b_probe_mesSA_superPointR_EM;   //!
   TBranch        *b_probe_mesSA_superPointR_EO;   //!
   TBranch        *b_probe_mesSA_superPointR_EE;   //!
   TBranch        *b_probe_mesSA_superPointR_CSC;   //!
   TBranch        *b_probe_mesSA_superPointR_BEE;   //!
   TBranch        *b_probe_mesSA_superPointR_BME;   //!
   TBranch        *b_probe_mesSA_superPointZ_BI;   //!
   TBranch        *b_probe_mesSA_superPointZ_BM;   //!
   TBranch        *b_probe_mesSA_superPointZ_BO;   //!
   TBranch        *b_probe_mesSA_superPointZ_EI;   //!
   TBranch        *b_probe_mesSA_superPointZ_EM;   //!
   TBranch        *b_probe_mesSA_superPointZ_EO;   //!
   TBranch        *b_probe_mesSA_superPointZ_EE;   //!
   TBranch        *b_probe_mesSA_superPointZ_CSC;   //!
   TBranch        *b_probe_mesSA_superPointZ_BEE;   //!
   TBranch        *b_probe_mesSA_superPointZ_BME;   //!
   TBranch        *b_probe_mesSA_superPointSlope_BI;   //!
   TBranch        *b_probe_mesSA_superPointSlope_BM;   //!
   TBranch        *b_probe_mesSA_superPointSlope_BO;   //!
   TBranch        *b_probe_mesSA_superPointSlope_EI;   //!
   TBranch        *b_probe_mesSA_superPointSlope_EM;   //!
   TBranch        *b_probe_mesSA_superPointSlope_EO;   //!
   TBranch        *b_probe_mesSA_superPointSlope_EE;   //!
   TBranch        *b_probe_mesSA_superPointSlope_CSC;   //!
   TBranch        *b_probe_mesSA_superPointSlope_BEE;   //!
   TBranch        *b_probe_mesSA_superPointSlope_BME;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_BI;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_BM;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_BO;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_EI;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_EM;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_EO;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_EE;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_CSC;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_BEE;   //!
   TBranch        *b_probe_mesSA_superPointIntercept_BME;   //!
   TBranch        *b_probe_mesSA_superPointChi2_BI;   //!
   TBranch        *b_probe_mesSA_superPointChi2_BM;   //!
   TBranch        *b_probe_mesSA_superPointChi2_BO;   //!
   TBranch        *b_probe_mesSA_superPointChi2_EI;   //!
   TBranch        *b_probe_mesSA_superPointChi2_EM;   //!
   TBranch        *b_probe_mesSA_superPointChi2_EO;   //!
   TBranch        *b_probe_mesSA_superPointChi2_EE;   //!
   TBranch        *b_probe_mesSA_superPointChi2_CSC;   //!
   TBranch        *b_probe_mesSA_superPointChi2_BEE;   //!
   TBranch        *b_probe_mesSA_superPointChi2_BME;   //!
   TBranch        *b_probe_mesSA_isRpcFailure;   //!
   TBranch        *b_probe_mesSA_isTgcFailure;   //!
   TBranch        *b_probe_mesSA_rpcHitX;   //!
   TBranch        *b_probe_mesSA_rpcHitY;   //!
   TBranch        *b_probe_mesSA_rpcHitZ;   //!
   TBranch        *b_probe_mesSA_rpcHitR;   //!
   TBranch        *b_probe_mesSA_rpcHitEta;   //!
   TBranch        *b_probe_mesSA_rpcHitPhi;   //!
   TBranch        *b_probe_mesSA_rpcHitStationNumber;   //!
   TBranch        *b_probe_mesSA_rpcHitStationName;   //!
   TBranch *b_probe_mesSA_mdtHitIsOutlier; //!
   TBranch *b_probe_mesSA_mdtHitChamber; //!
   TBranch *b_probe_mesSA_mdtHitR; //!
   TBranch *b_probe_mesSA_mdtHitZ; //!
   TBranch *b_probe_mesSA_mdtHitPhi; //!
   TBranch *b_probe_mesSA_mdtHitResidual; //!
   TBranch *b_probe_mesSA_roadAw; //!
   TBranch *b_probe_mesSA_roadBw; //!
   TBranch *b_probe_mesSA_zMin; //!
   TBranch *b_probe_mesSA_zMax; //!
   TBranch *b_probe_mesSA_rMin; //!
   TBranch *b_probe_mesSA_rMax; //!
   TBranch *b_probe_mesSA_etaMin; //!
   TBranch *b_probe_mesSA_etaMax; //!

   TBranch* b_probe_mesSA_roadFtkAw; //!
   TBranch* b_probe_mesSA_roadFtkBw; //!
   TBranch* b_probe_mesSA_extFtkEta; //!
   TBranch* b_probe_mesSA_extFtkPhi; //!
   TBranch* b_probe_mesSA_extFtkR; //!
   TBranch* b_probe_mesSA_extFtkZ; //!
   TBranch* b_probe_mesSA_roadRpcAw; //!
   TBranch* b_probe_mesSA_roadRpcBw; //!
   TBranch* b_probe_mesSA_roadAlgo; //!
   TBranch* b_probe_mesSA_ptFtk; //!
   TBranch* b_probe_mesSA_etaFtk; //!
   TBranch* b_probe_mesSA_phiFtk; //!

   TBranch* b_probe_mesL2_pass; //!
   TBranch* b_probe_mesL2_pt; //!
   TBranch* b_probe_mesL2_eta; //!
   TBranch* b_probe_mesL2_phi; //!

   TBranch* b_probe_mesL2_ptFtk; //!
   TBranch* b_probe_mesL2_etaFtk; //!
   TBranch* b_probe_mesL2_phiFtk; //!
   TBranch* b_probe_mesL2_roadAlgo; //!

   TBranch* b_probe_mesL2_passSA; //!
   TBranch* b_probe_mesL2_numSpSA; //!
   TBranch* b_probe_mesL2_ptSA; //!
   TBranch* b_probe_mesL2_etaSA; //!
   TBranch* b_probe_mesL2_phiSA; //!
   TBranch* b_probe_mesL2_sAddressSA; //!

   TBranch *b_probe_mesSA_roiEta; //!
   TBranch *b_probe_mesSA_roiPhi; //!

   TBranch        *b_probe_mesCB_pass;   //!
   TBranch        *b_probe_mesCB_dR;   //!
   TBranch        *b_probe_mesCB_tpdR;   //!
   TBranch        *b_probe_mesCB_pt;   //!
   TBranch        *b_probe_mesCB_eta;   //!
   TBranch        *b_probe_mesCB_phi;   //!
   TBranch        *b_probe_mesFTF_pass;   //!
   TBranch        *b_probe_mesFTF_dR;   //!
   TBranch        *b_probe_mesFTF_pt;   //!
   TBranch        *b_probe_mesFTF_eta;   //!
   TBranch        *b_probe_mesFTF_phi;   //!
   TBranch        *b_probe_mesEF_pass;   //!
   TBranch        *b_probe_mesEF_dR;   //!
   TBranch        *b_probe_mesEF_tpdR;   //!
   TBranch        *b_probe_mesEF_pt;   //!
   TBranch        *b_probe_mesEF_eta;   //!
   TBranch        *b_probe_mesEF_phi;   //!

   // Histgrams
   vector< string > m_chain; //!
   vector< TH1F* > m_probe_pt_mu4; //!
   vector< TH1F* > m_probe_mu_mu4; //!
   vector< TH1F* > m_probe_phi_mu4; //!
   vector< TH1F* > m_probe_eta_mu4; //!
   vector< TH2F* > m_probe_qetapt_mu4; //!
   vector< TH2F* > m_probe_etaphi_mu4; //!

   TH2F *m_ftk_size_pt_mu4;       //!
   TH2F *m_ftk_size_eta_mu4;       //!

   TH2F *m_ftk_size_pt_mu50;       //!
   TH2F *m_ftk_size_eta_mu50;       //!

   TH1F *m_eff_pt_mu4_L1;       //!
   TH1F *m_eff_pt_mu4_L1SA;       //!
   TH1F *m_eff_pt_mu4_SACB;       //!
   TH1F *m_eff_pt_mu4_L1L2;       //!
   TH1F *m_eff_pt_mu4_L1CB;       //!

   TH1F *m_eff_mu_mu4_L1;       //!
   TH1F *m_eff_mu_mu4_L1SA;       //!
   TH1F *m_eff_mu_mu4_SACB;       //!
   TH1F *m_eff_mu_mu4_L1L2;       //!
   TH1F *m_eff_mu_mu4_L1CB;       //!

   TH1F *m_eff_eta_mu4_L1;       //!
   TH1F *m_eff_eta_mu4_L1SA;       //!
   TH1F *m_eff_eta_mu4_SACB;       //!
   TH1F *m_eff_eta_mu4_L1L2;       //!
   TH1F *m_eff_eta_mu4_L1CB;       //!

   TH1F *m_eff_phi_mu4_L1;       //!
   TH1F *m_eff_phi_mu4_L1SA;       //!
   TH1F *m_eff_phi_mu4_SACB;       //!
   TH1F *m_eff_phi_mu4_L1L2;       //!
   TH1F *m_eff_phi_mu4_L1CB;       //!

   TH2F *m_eff_qetapt_mu4_L1L2;       //!
   TH2F *m_eff_etaphi_mu4_L1L2;       //!

   vector< TH1F* > h_probe_pt_mu6; //!
   vector< TH1F* > h_probe_mu_mu6; //!
   vector< TH1F* > h_probe_phi_mu6; //!
   vector< TH1F* > h_probe_eta_mu6; //!

   vector< TH1F* > h_probe_pt_mu50; //!
   vector< TH1F* > h_probe_mu_mu50; //!
   vector< TH1F* > h_probe_phi_mu50; //!
   vector< TH1F* > h_probe_eta_mu50; //!

   TH1F *h_probe_mu_mu4_L1;       //!
   TH1F *h_probe_mu_mu4_SA;       //!
   TH1F *h_eff_mu_mu4_L1;         //!
   TH1F *h_eff_mu_mu4_L1SA;       //!

   TH1F *h_probe_pt_mu4_PROBE;       //!
   TH1F *h_probe_pt_mu4_L1;       //!
   TH1F *h_probe_pt_mu4_SA;       //!
   TH1F *h_probe_pt_mu4_CB;       //!
   TH1F *h_probe_pt_mu4_EF;       //!
   TH1F *h_eff_pt_mu4_L1;         //!
   TH1F *h_eff_pt_mu4_L1SA;       //!
   TH1F *h_eff_pt_mu4_SACB;       //!
   TH1F *h_eff_pt_mu4_CBEF;       //!

   TH1F *h_probe_phi_mu4_L1;      //!
   TH1F *h_probe_phi_mu4_SA;      //!
   TH1F *h_eff_phi_mu4_L1;        //!
   TH1F *h_eff_phi_mu4_L1SA;      //!

   TH1F *h_probe_eta_mu4_L1;      //!
   TH1F *h_probe_eta_mu4_SA;      //!
   TH1F *h_eff_eta_mu4_L1;        //!
   TH1F *h_eff_eta_mu4_L1SA;      //!

   TH2F *hh_probe_etaphi_mu4_L1;  //!
   TH2F *hh_probe_etaphi_mu4_SA;  //!
   TH2F *hh_eff_etaphi_mu4_L1;    //!
   TH2F *hh_eff_etaphi_mu4_L1SA;  //!

   TH2F *hh_probe_qetapt_mu4_L1;  //!
   TH2F *hh_probe_qetapt_mu4_SA;  //!
   TH2F *hh_eff_qetapt_mu4_L1;    //!
   TH2F *hh_eff_qetapt_mu4_L1SA;  //!

   TH1F *h_probe_mu_mu50_L1;     //!
   TH1F *h_probe_mu_mu50_SA;     //!
   TH1F *h_probe_mu_mu50_CB;     //!
   TH1F *h_eff_mu_mu50_L1;       //!
   TH1F *h_eff_mu_mu50_L1SA;     //!
   TH1F *h_eff_mu_mu50_SACB;     //!

   TH1F *h_probe_pt_mu50_L1;     //!
   TH1F *h_probe_pt_mu50_SA;     //!
   TH1F *h_probe_pt_mu50_CB;     //!
   TH1F *h_eff_pt_mu50_L1;       //!
   TH1F *h_eff_pt_mu50_L1SA;     //!
   TH1F *h_eff_pt_mu50_SACB;     //!

   TH1F *h_probe_phi_mu50_L1;     //!
   TH1F *h_probe_phi_mu50_SA;     //!
   TH1F *h_probe_phi_mu50_CB;     //!
   TH1F *h_eff_phi_mu50_L1;       //!
   TH1F *h_eff_phi_mu50_L1SA;     //!
   TH1F *h_eff_phi_mu50_SACB;     //!

   TH1F *h_probe_eta_mu50_L1;     //!
   TH1F *h_probe_eta_mu50_SA;     //!
   TH1F *h_probe_eta_mu50_CB;     //!
   TH1F *h_eff_eta_mu50_L1;       //!
   TH1F *h_eff_eta_mu50_L1SA;     //!
   TH1F *h_eff_eta_mu50_SACB;     //!

   TH2F *hh_probe_etaphi_mu50_L1; //!
   TH2F *hh_probe_etaphi_mu50_SA; //!
   TH2F *hh_eff_etaphi_mu50_L1;   //!
   TH2F *hh_eff_etaphi_mu50_L1SA; //!

   TH2F *hh_probe_qetapt_mu50_L1; //!
   TH2F *hh_probe_qetapt_mu50_SA; //!
   TH2F *hh_eff_qetapt_mu50_L1;   //!
   TH2F *hh_eff_qetapt_mu50_L1SA; //!

   TH2F *h_superPointRZ_BIS; //!
   TH2F *h_superPointRZ_BIL; //!
   TH2F *h_superPointRZ_BMS; //!
   TH2F *h_superPointRZ_BML; //!
   TH2F *h_superPointRZ_BOS; //!
   TH2F *h_superPointRZ_BOL; //!

   TH2F *h_segmentRZ_BIS; //!
   TH2F *h_segmentRZ_BIL; //!
   TH2F *h_segmentRZ_BMS; //!
   TH2F *h_segmentRZ_BML; //!
   TH2F *h_segmentRZ_BOS; //!
   TH2F *h_segmentRZ_BOL; //!


   TH2F *h_residualRZ_BIS; //!
   TH2F *h_residualRZ_BIL; //!
   TH2F *h_residualRZ_BMS; //!
   TH2F *h_residualRZ_BML; //!
   TH2F *h_residualRZ_BOS; //!
   TH2F *h_residualRZ_BOL; //!

   TH2F *h_NumberOfSP_LumiBlock; //!
   TH2F *h_NumberOfSP_eta; //!
   TH2F *h_NumberOfSP_qeta; //!
   TH2F *h_NumberOfSP_pt_barrel; //!
   TH2F *h_NumberOfSP_pass_barrel; //!

   TH2F *h_NumberOfMdt_LumiBlock; //!
   TH2F *h_NumberOfMdt_LumiBlock_BI; //!
   TH2F *h_NumberOfMdt_LumiBlock_BM; //!
   TH2F *h_NumberOfMdt_LumiBlock_BO; //!

   TH2F *h_NumberOfMdt_eta; //!
   TH2F *h_NumberOfMdt_eta_BI; //!
   TH2F *h_NumberOfMdt_eta_BM; //!
   TH2F *h_NumberOfMdt_eta_BO; //!

   TH2F *h_NumberOfMdt_qeta; //!
   TH2F *h_NumberOfMdt_qeta_BI; //!
   TH2F *h_NumberOfMdt_qeta_BM; //!
   TH2F *h_NumberOfMdt_qeta_BO; //!

   TH2F *h_NumberOfMdt_pt_barrel; //!
   TH2F *h_NumberOfMdt_pt_barrel_BI; //!
   TH2F *h_NumberOfMdt_pt_barrel_BM; //!
   TH2F *h_NumberOfMdt_pt_barrel_BO; //!

   TH2F *h_ResidualSegment_eta; //!
   TH2F *h_ResidualSegment_eta_BI; //!
   TH2F *h_ResidualSegment_eta_BI_2d; //!
   TH2F *h_ResidualSegment_eta_BM; //!
   TH2F *h_ResidualSegment_eta_BO; //!

   TH1F *h_ResidualSegment_eta_rpc_BI; //!
   TH1F *h_ResidualSegment_eta_rpc_BM; //!
   TH1F *h_ResidualSegment_eta_rpc_BO; //!

   TH1F *h_ResidualSegment_eta_ftk_BI; //!
   TH1F *h_ResidualSegment_eta_ftk_BM; //!
   TH1F *h_ResidualSegment_eta_ftk_BO; //!

   TH2F *h_ResidualMdt_Outlier_eta; //!
   TH2F *h_ResidualMdt_Outlier_eta_BI; //!
   TH2F *h_ResidualMdt_Outlier_eta_BM; //!
   TH2F *h_ResidualMdt_Outlier_eta_BO; //!

   TH2F *h_ResidualMdt_Inlier_eta; //!
   TH2F *h_ResidualMdt_Inlier_eta_BI; //!
   TH2F *h_ResidualMdt_Inlier_eta_BM; //!
   TH2F *h_ResidualMdt_Inlier_eta_BO; //!

   TH2F *h_ResidualMdt_Outlier_pt_barrel; //!
   TH2F *h_ResidualMdt_Outlier_pt_barrel_BI; //!
   TH2F *h_ResidualMdt_Outlier_pt_barrel_BM; //!
   TH2F *h_ResidualMdt_Outlier_pt_barrel_BO; //!

   TH2F *h_ResidualMdt_Inlier_pt_barrel; //!
   TH2F *h_ResidualMdt_Inlier_pt_barrel_BI; //!
   TH2F *h_ResidualMdt_Inlier_pt_barrel_BM; //!
   TH2F *h_ResidualMdt_Inlier_pt_barrel_BO; //!

   TH2F *h_PtResidual_pt1; //!
   TH2F *h_PtResidual_pt2; //!
   TH2F *h_PtResidual_pt3; //!
   TH2F *h_PtResidual_ptFtk; //!
   TH2F *h_PtResidual_ptLut; //!
   TH2F *h_PtResidual_eta; //!
   TH2F *h_pt_vs_pt; //!
   TH2F *h_ptl2_vs_ptFtk; //!
   TH2F *h_pt_vs_ptFtk; //!

   TH2F* h_InEff_pt;
   TH2F* h_InEff_eta;
   TH2F* h_InEff_qeta;
   TH3F* h_InEff_qetapt;


   RPC(TChain *tree);
   virtual ~RPC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Loop( int Nevents, int DisplayNumber );
   virtual void     DrawHist(TString pdf);
   virtual void     End();
   virtual Bool_t   Notify();
   virtual void     Display(int tap_type, int trig_chain, Long64_t begin_entry, Long64_t limit_entry, TString pdf);
   virtual double rWidthToBw( double aw, double rwidth);
   virtual void     Show(Long64_t entry = -1);
   virtual double   calc_residual(double aw, double bw, double Z, double R);

   // Initialize
   virtual void Init(TTree *tree);
   virtual void InitHist();
   virtual void InitInEffHist();
   virtual void InitEffHist();
   virtual void InitPtResidualHist();
   virtual void InitMdtHist();
   virtual void InitSPHist();

   // SuperPoints
   virtual void FillSPHist();
   virtual int  NumberOfSP();
   virtual int  NumberOfSP(int NTrigChain);
   virtual void DrawFractionOfnSPs(TH2F* h_NumberOfSP, TCanvas* c1, TString pdf);

   // MDT hits
   virtual void FillMdtHist();
   virtual void DrawFractionOfnMDTs(TH2F* h_NumberOfMdt, TCanvas* c1, TString pdf);

   // Efficiency
   virtual void FillProbeHist();
   virtual void CalcEff();
   virtual void DrawEffHist(TString pdf);
   void CalcHistToHist( TH1F* h1, TH1F* h2, TH1F* hout );
   void CalcHistToHist( TH2F* h1, TH2F* h2, TH2F* hout );

   // pT residual
   virtual double calc_pTresidual( double offline_pt, double trig_pt);
   virtual void FillPtResidualHist();
   virtual void DrawPtResidualHist(TCanvas* c1, TString pdf);

   // InEff
   virtual void FillInEffHist(int tap_type, int trig_chain);
   virtual void DrawInEffHist(TString pdf);

};

#endif

#ifdef RPC_cxx
RPC::RPC(TChain *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/yfukuhar/gpfs/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/yfukuhar/gpfs/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0.root");
      }
      f->GetObject("t_tap",tree);

   }
   Init(tree);
}

RPC::~RPC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RPC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RPC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void RPC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mes_name = 0;
   probe_mesEFTAG_pass = 0;
   probe_mesEFTAG_dR = 0;
   probe_mesEFTAG_tpdR = 0;
   probe_mesEFTAG_pt = 0;
   probe_mesEFTAG_eta = 0;
   probe_mesEFTAG_phi = 0;
   probe_mesL1_pass = 0;
   probe_mesL1_dR = 0;
   probe_mesL1_tpdR = 0;
   probe_mesL1_pt = 0;
   probe_mesL1_eta = 0;
   probe_mesL1_phi = 0;
   probe_mesSA_pass = 0;
   probe_mesSA_dR = 0;
   probe_mesSA_tpdR = 0;
   probe_mesSA_pt = 0;
   probe_mesSA_eta = 0;
   probe_mesSA_phi = 0;
   probe_mesSA_etams = 0;
   probe_mesSA_phims = 0;
   probe_mesSA_etabe = 0;
   probe_mesSA_phibe = 0;
   probe_mesSA_tgcpt = 0;
   probe_mesSA_ptBarrelRadius = 0;
   probe_mesSA_ptBarrelSagitta = 0;
   probe_mesSA_ptEndcapAlpha = 0;
   probe_mesSA_ptEndcapBeta = 0;
   probe_mesSA_ptEndcapRadius = 0;
   probe_mesSA_ptCSC = 0;
   probe_mesSA_sAddress = 0;
   probe_mesSA_superPointR_BI = 0;
   probe_mesSA_superPointR_BM = 0;
   probe_mesSA_superPointR_BO = 0;
   probe_mesSA_superPointR_EI = 0;
   probe_mesSA_superPointR_EM = 0;
   probe_mesSA_superPointR_EO = 0;
   probe_mesSA_superPointR_EE = 0;
   probe_mesSA_superPointR_CSC = 0;
   probe_mesSA_superPointR_BEE = 0;
   probe_mesSA_superPointR_BME = 0;
   probe_mesSA_superPointZ_BI = 0;
   probe_mesSA_superPointZ_BM = 0;
   probe_mesSA_superPointZ_BO = 0;
   probe_mesSA_superPointZ_EI = 0;
   probe_mesSA_superPointZ_EM = 0;
   probe_mesSA_superPointZ_EO = 0;
   probe_mesSA_superPointZ_EE = 0;
   probe_mesSA_superPointZ_CSC = 0;
   probe_mesSA_superPointZ_BEE = 0;
   probe_mesSA_superPointZ_BME = 0;
   probe_mesSA_superPointSlope_BI = 0;
   probe_mesSA_superPointSlope_BM = 0;
   probe_mesSA_superPointSlope_BO = 0;
   probe_mesSA_superPointSlope_EI = 0;
   probe_mesSA_superPointSlope_EM = 0;
   probe_mesSA_superPointSlope_EO = 0;
   probe_mesSA_superPointSlope_EE = 0;
   probe_mesSA_superPointSlope_CSC = 0;
   probe_mesSA_superPointSlope_BEE = 0;
   probe_mesSA_superPointSlope_BME = 0;
   probe_mesSA_superPointIntercept_BI = 0;
   probe_mesSA_superPointIntercept_BM = 0;
   probe_mesSA_superPointIntercept_BO = 0;
   probe_mesSA_superPointIntercept_EI = 0;
   probe_mesSA_superPointIntercept_EM = 0;
   probe_mesSA_superPointIntercept_EO = 0;
   probe_mesSA_superPointIntercept_EE = 0;
   probe_mesSA_superPointIntercept_CSC = 0;
   probe_mesSA_superPointIntercept_BEE = 0;
   probe_mesSA_superPointIntercept_BME = 0;
   probe_mesSA_superPointChi2_BI = 0;
   probe_mesSA_superPointChi2_BM = 0;
   probe_mesSA_superPointChi2_BO = 0;
   probe_mesSA_superPointChi2_EI = 0;
   probe_mesSA_superPointChi2_EM = 0;
   probe_mesSA_superPointChi2_EO = 0;
   probe_mesSA_superPointChi2_EE = 0;
   probe_mesSA_superPointChi2_CSC = 0;
   probe_mesSA_superPointChi2_BEE = 0;
   probe_mesSA_superPointChi2_BME = 0;
   probe_mesSA_isRpcFailure = 0;
   probe_mesSA_isTgcFailure = 0;
   probe_mesSA_rpcHitX = 0;
   probe_mesSA_rpcHitY = 0;
   probe_mesSA_rpcHitZ = 0;
   probe_mesSA_rpcHitR = 0;
   probe_mesSA_rpcHitEta = 0;
   probe_mesSA_rpcHitPhi = 0;
   probe_mesSA_rpcHitStationNumber = 0;
   probe_mesSA_rpcHitStationName = 0;
   probe_mesSA_mdtHitIsOutlier = 0;
   probe_mesSA_mdtHitChamber = 0;
   probe_mesSA_mdtHitR = 0;
   probe_mesSA_mdtHitZ = 0;
   probe_mesSA_mdtHitPhi = 0;
   probe_mesSA_mdtHitResidual = 0;
   probe_mesSA_roadAw = 0;
   probe_mesSA_roadBw = 0;
   probe_mesSA_zMin   = 0;
   probe_mesSA_zMax   = 0;
   probe_mesSA_rMin   = 0;
   probe_mesSA_rMax   = 0;
   probe_mesSA_etaMin = 0;
   probe_mesSA_etaMax = 0;

   probe_mesSA_roadFtkAw = 0;
   probe_mesSA_roadFtkBw = 0;
   probe_mesSA_extFtkEta = 0;
   probe_mesSA_extFtkPhi = 0;
   probe_mesSA_extFtkR = 0;
   probe_mesSA_extFtkZ = 0;
   probe_mesSA_roadRpcAw = 0;
   probe_mesSA_roadRpcBw = 0;
   probe_mesSA_roadAlgo = 0;
   probe_mesSA_ptFtk = 0;
   probe_mesSA_etaFtk = 0;
   probe_mesSA_phiFtk = 0;

   probe_mesSA_roiEta = 0;
   probe_mesSA_roiPhi = 0;

   probe_mesL2_pass       = 0;
   probe_mesL2_pt         = 0;
   probe_mesL2_eta        = 0;
   probe_mesL2_phi        = 0;

   probe_mesL2_ptFtk      = 0;
   probe_mesL2_etaFtk     = 0;
   probe_mesL2_phiFtk     = 0;
   probe_mesL2_roadAlgo   = 0;

   probe_mesL2_passSA     = 0;
   probe_mesL2_numSpSA    = 0;
   probe_mesL2_ptSA       = 0;
   probe_mesL2_etaSA      = 0;
   probe_mesL2_phiSA      = 0;
   probe_mesL2_sAddressSA = 0;

   probe_mesCB_pass = 0;
   probe_mesCB_dR = 0;
   probe_mesCB_tpdR = 0;
   probe_mesCB_pt = 0;
   probe_mesCB_eta = 0;
   probe_mesCB_phi = 0;
   probe_mesFTF_pass = 0;
   probe_mesFTF_dR = 0;
   probe_mesFTF_pt = 0;
   probe_mesFTF_eta = 0;
   probe_mesFTF_phi = 0;
   probe_mesEF_pass = 0;
   probe_mesEF_dR = 0;
   probe_mesEF_tpdR = 0;
   probe_mesEF_pt = 0;
   probe_mesEF_eta = 0;
   probe_mesEF_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("AverageInteractionsPerCrossing", &AverageInteractionsPerCrossing, &b_AverageInteractionsPerCrossing);
   fChain->SetBranchAddress("mes_n", &mes_n, &b_mes_n);
   fChain->SetBranchAddress("mes_name", &mes_name, &b_mes_name);
   fChain->SetBranchAddress("sumReqdRL1", &sumReqdRL1, &b_sumReqdRL1);
   fChain->SetBranchAddress("sumReqdREF", &sumReqdREF, &b_sumReqdREF);
   fChain->SetBranchAddress("tag_proc", &tag_proc, &b_tag_proc);
   fChain->SetBranchAddress("tag_ReqdRL1", &tag_ReqdRL1, &b_tag_ReqdRL1);
   fChain->SetBranchAddress("tag_ReqdREF", &tag_ReqdREF, &b_tag_ReqdREF);
   fChain->SetBranchAddress("tag_dRL1", &tag_dRL1, &b_tag_dRL1);
   fChain->SetBranchAddress("tag_dREF", &tag_dREF, &b_tag_dREF);
   fChain->SetBranchAddress("tag_pt", &tag_pt, &b_tag_pt);
   fChain->SetBranchAddress("tag_eta", &tag_eta, &b_tag_eta);
   fChain->SetBranchAddress("tag_exteta", &tag_exteta, &b_tag_exteta);
   fChain->SetBranchAddress("tag_extinneta", &tag_extinneta, &b_tag_extinneta);
   fChain->SetBranchAddress("tag_phi", &tag_phi, &b_tag_phi);
   fChain->SetBranchAddress("tag_extphi", &tag_extphi, &b_tag_extphi);
   fChain->SetBranchAddress("tag_extinnphi", &tag_extinnphi, &b_tag_extinnphi);
   fChain->SetBranchAddress("tag_d0", &tag_d0, &b_tag_d0);
   fChain->SetBranchAddress("tag_z0", &tag_z0, &b_tag_z0);
   fChain->SetBranchAddress("tag_charge", &tag_charge, &b_tag_charge);
   fChain->SetBranchAddress("tag_L1_pt", &tag_L1_pt, &b_tag_L1_pt);
   fChain->SetBranchAddress("tag_L1_eta", &tag_L1_eta, &b_tag_L1_eta);
   fChain->SetBranchAddress("tag_L1_phi", &tag_L1_phi, &b_tag_L1_phi);
   fChain->SetBranchAddress("tag_SA_pt", &tag_SA_pt, &b_tag_SA_pt);
   fChain->SetBranchAddress("tag_SA_eta", &tag_SA_eta, &b_tag_SA_eta);
   fChain->SetBranchAddress("tag_SA_phi", &tag_SA_phi, &b_tag_SA_phi);
   fChain->SetBranchAddress("tag_CB_pt", &tag_CB_pt, &b_tag_CB_pt);
   fChain->SetBranchAddress("tag_CB_eta", &tag_CB_eta, &b_tag_CB_eta);
   fChain->SetBranchAddress("tag_CB_phi", &tag_CB_phi, &b_tag_CB_phi);
   fChain->SetBranchAddress("tag_EF_pt", &tag_EF_pt, &b_tag_EF_pt);
   fChain->SetBranchAddress("tag_EF_eta", &tag_EF_eta, &b_tag_EF_eta);
   fChain->SetBranchAddress("tag_EF_phi", &tag_EF_phi, &b_tag_EF_phi);
   fChain->SetBranchAddress("probe_pt", &probe_pt, &b_probe_pt);
   fChain->SetBranchAddress("probe_eta", &probe_eta, &b_probe_eta);
   fChain->SetBranchAddress("probe_exteta", &probe_exteta, &b_probe_exteta);
   fChain->SetBranchAddress("probe_extinneta", &probe_extinneta, &b_probe_extinneta);
   fChain->SetBranchAddress("probe_phi", &probe_phi, &b_probe_phi);
   fChain->SetBranchAddress("probe_extphi", &probe_extphi, &b_probe_extphi);
   fChain->SetBranchAddress("probe_extinnphi", &probe_extinnphi, &b_probe_extinnphi);
   fChain->SetBranchAddress("probe_d0", &probe_d0, &b_probe_d0);
   fChain->SetBranchAddress("probe_z0", &probe_z0, &b_probe_z0);
   fChain->SetBranchAddress("probe_charge", &probe_charge, &b_probe_charge);
   fChain->SetBranchAddress("probe_segment_n", &probe_segment_n, &b_probe_segment_n);
   fChain->SetBranchAddress("probe_segment_x", probe_segment_x, &b_probe_segment_x);
   fChain->SetBranchAddress("probe_segment_y", probe_segment_y, &b_probe_segment_y);
   fChain->SetBranchAddress("probe_segment_z", probe_segment_z, &b_probe_segment_z);
   fChain->SetBranchAddress("probe_segment_px", probe_segment_px, &b_probe_segment_px);
   fChain->SetBranchAddress("probe_segment_py", probe_segment_py, &b_probe_segment_py);
   fChain->SetBranchAddress("probe_segment_pz", probe_segment_pz, &b_probe_segment_pz);
   fChain->SetBranchAddress("probe_segment_chiSquared", probe_segment_chiSquared, &b_probe_segment_chiSquared);
   fChain->SetBranchAddress("probe_segment_numberDoF", probe_segment_numberDoF, &b_probe_segment_numberDoF);
   fChain->SetBranchAddress("probe_segment_sector", probe_segment_sector, &b_probe_segment_sector);
   fChain->SetBranchAddress("probe_segment_chamberIndex", probe_segment_chamberIndex, &b_probe_segment_chamberIndex);
   fChain->SetBranchAddress("probe_segment_etaIndex", probe_segment_etaIndex, &b_probe_segment_etaIndex);
   fChain->SetBranchAddress("probe_segment_nPrecisionHits", probe_segment_nPrecisionHits, &b_probe_segment_nPrecisionHits);
   fChain->SetBranchAddress("probe_segment_nPhiLayers", probe_segment_nPhiLayers, &b_probe_segment_nPhiLayers);
   fChain->SetBranchAddress("probe_segment_nTrigEtaLayers", probe_segment_nTrigEtaLayers, &b_probe_segment_nTrigEtaLayers);
   fChain->SetBranchAddress("tp_dR", &tp_dR, &b_tp_dR);
   fChain->SetBranchAddress("tp_deta", &tp_deta, &b_tp_deta);
   fChain->SetBranchAddress("tp_dphi", &tp_dphi, &b_tp_dphi);
   fChain->SetBranchAddress("tp_extdR", &tp_extdR, &b_tp_extdR);
   fChain->SetBranchAddress("tp_extdeta", &tp_extdeta, &b_tp_extdeta);
   fChain->SetBranchAddress("tp_extdphi", &tp_extdphi, &b_tp_extdphi);
   fChain->SetBranchAddress("tp_extinndR", &tp_extinndR, &b_tp_extinndR);
   fChain->SetBranchAddress("tp_extinndeta", &tp_extinndeta, &b_tp_extinndeta);
   fChain->SetBranchAddress("tp_extinndphi", &tp_extinndphi, &b_tp_extinndphi);
   fChain->SetBranchAddress("tp_mass", &tp_mass, &b_tp_mass);
   fChain->SetBranchAddress("tp_vftlxy", &tp_vftlxy, &b_tp_vftlxy);
   fChain->SetBranchAddress("tp_vftchi2", &tp_vftchi2, &b_tp_vftchi2);
   fChain->SetBranchAddress("tp_vftndof", &tp_vftndof, &b_tp_vftndof);
   fChain->SetBranchAddress("probe_mesEFTAG_pass", &probe_mesEFTAG_pass, &b_probe_mesEFTAG_pass);
   fChain->SetBranchAddress("probe_mesEFTAG_dR", &probe_mesEFTAG_dR, &b_probe_mesEFTAG_dR);
   fChain->SetBranchAddress("probe_mesEFTAG_tpdR", &probe_mesEFTAG_tpdR, &b_probe_mesEFTAG_tpdR);
   fChain->SetBranchAddress("probe_mesEFTAG_pt", &probe_mesEFTAG_pt, &b_probe_mesEFTAG_pt);
   fChain->SetBranchAddress("probe_mesEFTAG_eta", &probe_mesEFTAG_eta, &b_probe_mesEFTAG_eta);
   fChain->SetBranchAddress("probe_mesEFTAG_phi", &probe_mesEFTAG_phi, &b_probe_mesEFTAG_phi);
   fChain->SetBranchAddress("probe_mesL1_pass", &probe_mesL1_pass, &b_probe_mesL1_pass);
   fChain->SetBranchAddress("probe_mesL1_dR", &probe_mesL1_dR, &b_probe_mesL1_dR);
   fChain->SetBranchAddress("probe_mesL1_tpdR", &probe_mesL1_tpdR, &b_probe_mesL1_tpdR);
   fChain->SetBranchAddress("probe_mesL1_pt", &probe_mesL1_pt, &b_probe_mesL1_pt);
   fChain->SetBranchAddress("probe_mesL1_eta", &probe_mesL1_eta, &b_probe_mesL1_eta);
   fChain->SetBranchAddress("probe_mesL1_phi", &probe_mesL1_phi, &b_probe_mesL1_phi);
   fChain->SetBranchAddress("probe_mesSA_pass", &probe_mesSA_pass, &b_probe_mesSA_pass);
   fChain->SetBranchAddress("probe_mesSA_dR", &probe_mesSA_dR, &b_probe_mesSA_dR);
   fChain->SetBranchAddress("probe_mesSA_tpdR", &probe_mesSA_tpdR, &b_probe_mesSA_tpdR);
   fChain->SetBranchAddress("probe_mesSA_pt", &probe_mesSA_pt, &b_probe_mesSA_pt);
   fChain->SetBranchAddress("probe_mesSA_eta", &probe_mesSA_eta, &b_probe_mesSA_eta);
   fChain->SetBranchAddress("probe_mesSA_phi", &probe_mesSA_phi, &b_probe_mesSA_phi);
   fChain->SetBranchAddress("probe_mesSA_etams", &probe_mesSA_etams, &b_probe_mesSA_etams);
   fChain->SetBranchAddress("probe_mesSA_phims", &probe_mesSA_phims, &b_probe_mesSA_phims);
   fChain->SetBranchAddress("probe_mesSA_etabe", &probe_mesSA_etabe, &b_probe_mesSA_etabe);
   fChain->SetBranchAddress("probe_mesSA_phibe", &probe_mesSA_phibe, &b_probe_mesSA_phibe);
   fChain->SetBranchAddress("probe_mesSA_tgcpt", &probe_mesSA_tgcpt, &b_probe_mesSA_tgcpt);
   fChain->SetBranchAddress("probe_mesSA_ptBarrelRadius", &probe_mesSA_ptBarrelRadius, &b_probe_mesSA_ptBarrelRadius);
   fChain->SetBranchAddress("probe_mesSA_ptBarrelSagitta", &probe_mesSA_ptBarrelSagitta, &b_probe_mesSA_ptBarrelSagitta);
   fChain->SetBranchAddress("probe_mesSA_ptEndcapAlpha", &probe_mesSA_ptEndcapAlpha, &b_probe_mesSA_ptEndcapAlpha);
   fChain->SetBranchAddress("probe_mesSA_ptEndcapBeta", &probe_mesSA_ptEndcapBeta, &b_probe_mesSA_ptEndcapBeta);
   fChain->SetBranchAddress("probe_mesSA_ptEndcapRadius", &probe_mesSA_ptEndcapRadius, &b_probe_mesSA_ptEndcapRadius);
   fChain->SetBranchAddress("probe_mesSA_ptCSC", &probe_mesSA_ptCSC, &b_probe_mesSA_ptCSC);
   fChain->SetBranchAddress("probe_mesSA_sAddress", &probe_mesSA_sAddress, &b_probe_mesSA_sAddress);
   fChain->SetBranchAddress("probe_mesSA_superPointR_BI", &probe_mesSA_superPointR_BI, &b_probe_mesSA_superPointR_BI);
   fChain->SetBranchAddress("probe_mesSA_superPointR_BM", &probe_mesSA_superPointR_BM, &b_probe_mesSA_superPointR_BM);
   fChain->SetBranchAddress("probe_mesSA_superPointR_BO", &probe_mesSA_superPointR_BO, &b_probe_mesSA_superPointR_BO);
   fChain->SetBranchAddress("probe_mesSA_superPointR_EI", &probe_mesSA_superPointR_EI, &b_probe_mesSA_superPointR_EI);
   fChain->SetBranchAddress("probe_mesSA_superPointR_EM", &probe_mesSA_superPointR_EM, &b_probe_mesSA_superPointR_EM);
   fChain->SetBranchAddress("probe_mesSA_superPointR_EO", &probe_mesSA_superPointR_EO, &b_probe_mesSA_superPointR_EO);
   fChain->SetBranchAddress("probe_mesSA_superPointR_EE", &probe_mesSA_superPointR_EE, &b_probe_mesSA_superPointR_EE);
   fChain->SetBranchAddress("probe_mesSA_superPointR_CSC", &probe_mesSA_superPointR_CSC, &b_probe_mesSA_superPointR_CSC);
   fChain->SetBranchAddress("probe_mesSA_superPointR_BEE", &probe_mesSA_superPointR_BEE, &b_probe_mesSA_superPointR_BEE);
   fChain->SetBranchAddress("probe_mesSA_superPointR_BME", &probe_mesSA_superPointR_BME, &b_probe_mesSA_superPointR_BME);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_BI", &probe_mesSA_superPointZ_BI, &b_probe_mesSA_superPointZ_BI);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_BM", &probe_mesSA_superPointZ_BM, &b_probe_mesSA_superPointZ_BM);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_BO", &probe_mesSA_superPointZ_BO, &b_probe_mesSA_superPointZ_BO);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_EI", &probe_mesSA_superPointZ_EI, &b_probe_mesSA_superPointZ_EI);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_EM", &probe_mesSA_superPointZ_EM, &b_probe_mesSA_superPointZ_EM);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_EO", &probe_mesSA_superPointZ_EO, &b_probe_mesSA_superPointZ_EO);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_EE", &probe_mesSA_superPointZ_EE, &b_probe_mesSA_superPointZ_EE);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_CSC", &probe_mesSA_superPointZ_CSC, &b_probe_mesSA_superPointZ_CSC);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_BEE", &probe_mesSA_superPointZ_BEE, &b_probe_mesSA_superPointZ_BEE);
   fChain->SetBranchAddress("probe_mesSA_superPointZ_BME", &probe_mesSA_superPointZ_BME, &b_probe_mesSA_superPointZ_BME);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_BI", &probe_mesSA_superPointSlope_BI, &b_probe_mesSA_superPointSlope_BI);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_BM", &probe_mesSA_superPointSlope_BM, &b_probe_mesSA_superPointSlope_BM);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_BO", &probe_mesSA_superPointSlope_BO, &b_probe_mesSA_superPointSlope_BO);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_EI", &probe_mesSA_superPointSlope_EI, &b_probe_mesSA_superPointSlope_EI);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_EM", &probe_mesSA_superPointSlope_EM, &b_probe_mesSA_superPointSlope_EM);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_EO", &probe_mesSA_superPointSlope_EO, &b_probe_mesSA_superPointSlope_EO);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_EE", &probe_mesSA_superPointSlope_EE, &b_probe_mesSA_superPointSlope_EE);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_CSC", &probe_mesSA_superPointSlope_CSC, &b_probe_mesSA_superPointSlope_CSC);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_BEE", &probe_mesSA_superPointSlope_BEE, &b_probe_mesSA_superPointSlope_BEE);
   fChain->SetBranchAddress("probe_mesSA_superPointSlope_BME", &probe_mesSA_superPointSlope_BME, &b_probe_mesSA_superPointSlope_BME);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_BI", &probe_mesSA_superPointIntercept_BI, &b_probe_mesSA_superPointIntercept_BI);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_BM", &probe_mesSA_superPointIntercept_BM, &b_probe_mesSA_superPointIntercept_BM);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_BO", &probe_mesSA_superPointIntercept_BO, &b_probe_mesSA_superPointIntercept_BO);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_EI", &probe_mesSA_superPointIntercept_EI, &b_probe_mesSA_superPointIntercept_EI);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_EM", &probe_mesSA_superPointIntercept_EM, &b_probe_mesSA_superPointIntercept_EM);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_EO", &probe_mesSA_superPointIntercept_EO, &b_probe_mesSA_superPointIntercept_EO);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_EE", &probe_mesSA_superPointIntercept_EE, &b_probe_mesSA_superPointIntercept_EE);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_CSC", &probe_mesSA_superPointIntercept_CSC, &b_probe_mesSA_superPointIntercept_CSC);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_BEE", &probe_mesSA_superPointIntercept_BEE, &b_probe_mesSA_superPointIntercept_BEE);
   fChain->SetBranchAddress("probe_mesSA_superPointIntercept_BME", &probe_mesSA_superPointIntercept_BME, &b_probe_mesSA_superPointIntercept_BME);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_BI", &probe_mesSA_superPointChi2_BI, &b_probe_mesSA_superPointChi2_BI);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_BM", &probe_mesSA_superPointChi2_BM, &b_probe_mesSA_superPointChi2_BM);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_BO", &probe_mesSA_superPointChi2_BO, &b_probe_mesSA_superPointChi2_BO);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_EI", &probe_mesSA_superPointChi2_EI, &b_probe_mesSA_superPointChi2_EI);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_EM", &probe_mesSA_superPointChi2_EM, &b_probe_mesSA_superPointChi2_EM);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_EO", &probe_mesSA_superPointChi2_EO, &b_probe_mesSA_superPointChi2_EO);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_EE", &probe_mesSA_superPointChi2_EE, &b_probe_mesSA_superPointChi2_EE);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_CSC", &probe_mesSA_superPointChi2_CSC, &b_probe_mesSA_superPointChi2_CSC);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_BEE", &probe_mesSA_superPointChi2_BEE, &b_probe_mesSA_superPointChi2_BEE);
   fChain->SetBranchAddress("probe_mesSA_superPointChi2_BME", &probe_mesSA_superPointChi2_BME, &b_probe_mesSA_superPointChi2_BME);
   fChain->SetBranchAddress("probe_mesSA_isRpcFailure", &probe_mesSA_isRpcFailure, &b_probe_mesSA_isRpcFailure);
   fChain->SetBranchAddress("probe_mesSA_isTgcFailure", &probe_mesSA_isTgcFailure, &b_probe_mesSA_isTgcFailure);
   fChain->SetBranchAddress("probe_mesSA_rpcHitX", &probe_mesSA_rpcHitX, &b_probe_mesSA_rpcHitX);
   fChain->SetBranchAddress("probe_mesSA_rpcHitY", &probe_mesSA_rpcHitY, &b_probe_mesSA_rpcHitY);
   fChain->SetBranchAddress("probe_mesSA_rpcHitZ", &probe_mesSA_rpcHitZ, &b_probe_mesSA_rpcHitZ);
   fChain->SetBranchAddress("probe_mesSA_rpcHitR", &probe_mesSA_rpcHitR, &b_probe_mesSA_rpcHitR);
   fChain->SetBranchAddress("probe_mesSA_rpcHitEta", &probe_mesSA_rpcHitEta, &b_probe_mesSA_rpcHitEta);
   fChain->SetBranchAddress("probe_mesSA_rpcHitPhi", &probe_mesSA_rpcHitPhi, &b_probe_mesSA_rpcHitPhi);
   fChain->SetBranchAddress("probe_mesSA_rpcHitStationNumber", &probe_mesSA_rpcHitStationNumber, &b_probe_mesSA_rpcHitStationNumber);
   fChain->SetBranchAddress("probe_mesSA_mdtHitIsOutlier",     &probe_mesSA_mdtHitIsOutlier,     &b_probe_mesSA_mdtHitIsOutlier);
   fChain->SetBranchAddress("probe_mesSA_mdtHitChamber",       &probe_mesSA_mdtHitChamber,       &b_probe_mesSA_mdtHitChamber);
   fChain->SetBranchAddress("probe_mesSA_mdtHitR",             &probe_mesSA_mdtHitR,             &b_probe_mesSA_mdtHitR);
   fChain->SetBranchAddress("probe_mesSA_mdtHitZ",             &probe_mesSA_mdtHitZ,             &b_probe_mesSA_mdtHitZ);
   fChain->SetBranchAddress("probe_mesSA_mdtHitPhi",           &probe_mesSA_mdtHitPhi,           &b_probe_mesSA_mdtHitPhi);
   fChain->SetBranchAddress("probe_mesSA_mdtHitResidual",      &probe_mesSA_mdtHitResidual,      &b_probe_mesSA_mdtHitResidual);

   fChain->SetBranchAddress("probe_mesSA_roadAw", &probe_mesSA_roadAw, &b_probe_mesSA_roadAw );
   fChain->SetBranchAddress("probe_mesSA_roadBw", &probe_mesSA_roadBw, &b_probe_mesSA_roadBw );
   fChain->SetBranchAddress("probe_mesSA_zMin",   &probe_mesSA_zMin,   &b_probe_mesSA_zMin   );
   fChain->SetBranchAddress("probe_mesSA_zMax",   &probe_mesSA_zMax,   &b_probe_mesSA_zMax   );
   fChain->SetBranchAddress("probe_mesSA_rMin",   &probe_mesSA_rMin,   &b_probe_mesSA_rMin   );
   fChain->SetBranchAddress("probe_mesSA_rMax",   &probe_mesSA_rMax,   &b_probe_mesSA_rMax   );
   fChain->SetBranchAddress("probe_mesSA_etaMin", &probe_mesSA_etaMin, &b_probe_mesSA_etaMin );
   fChain->SetBranchAddress("probe_mesSA_etaMax", &probe_mesSA_etaMax, &b_probe_mesSA_etaMax );

   fChain->SetBranchAddress("probe_mesSA_roadFtkAw", &probe_mesSA_roadFtkAw, &b_probe_mesSA_roadFtkAw );
   fChain->SetBranchAddress("probe_mesSA_roadFtkBw", &probe_mesSA_roadFtkBw, &b_probe_mesSA_roadFtkBw );
   fChain->SetBranchAddress("probe_mesSA_extFtkEta", &probe_mesSA_extFtkEta, &b_probe_mesSA_extFtkEta );
   fChain->SetBranchAddress("probe_mesSA_extFtkPhi", &probe_mesSA_extFtkPhi, &b_probe_mesSA_extFtkPhi );
   fChain->SetBranchAddress("probe_mesSA_extFtkR",   &probe_mesSA_extFtkR, &b_probe_mesSA_extFtkR );
   fChain->SetBranchAddress("probe_mesSA_extFtkZ",   &probe_mesSA_extFtkZ, &b_probe_mesSA_extFtkZ );
   fChain->SetBranchAddress("probe_mesSA_roadRpcAw", &probe_mesSA_roadRpcAw, &b_probe_mesSA_roadRpcAw );
   fChain->SetBranchAddress("probe_mesSA_roadRpcBw", &probe_mesSA_roadRpcBw, &b_probe_mesSA_roadRpcBw );
   fChain->SetBranchAddress("probe_mesSA_roadAlgo", &probe_mesSA_roadAlgo, &b_probe_mesSA_roadAlgo );
   fChain->SetBranchAddress("probe_mesSA_ptFtk", &probe_mesSA_ptFtk, &b_probe_mesSA_ptFtk );
   fChain->SetBranchAddress("probe_mesSA_etaFtk", &probe_mesSA_etaFtk, &b_probe_mesSA_etaFtk );
   fChain->SetBranchAddress("probe_mesSA_phiFtk", &probe_mesSA_phiFtk, &b_probe_mesSA_phiFtk );

   fChain->SetBranchAddress("probe_mesL2_pass",       &probe_mesL2_pass, &b_probe_mesL2_pass);
   fChain->SetBranchAddress("probe_mesL2_pt",         &probe_mesL2_pt, &b_probe_mesL2_pt);
   fChain->SetBranchAddress("probe_mesL2_eta",        &probe_mesL2_eta, &b_probe_mesL2_eta);
   fChain->SetBranchAddress("probe_mesL2_phi",        &probe_mesL2_phi, &b_probe_mesL2_phi);

   fChain->SetBranchAddress("probe_mesL2_ptFtk",      &probe_mesL2_ptFtk, &b_probe_mesL2_ptFtk);
   fChain->SetBranchAddress("probe_mesL2_etaFtk",     &probe_mesL2_etaFtk, &b_probe_mesL2_etaFtk);
   fChain->SetBranchAddress("probe_mesL2_phiFtk",     &probe_mesL2_phiFtk, &b_probe_mesL2_phiFtk);
   fChain->SetBranchAddress("probe_mesL2_roadAlgo",   &probe_mesL2_roadAlgo, &b_probe_mesL2_roadAlgo);

   fChain->SetBranchAddress("probe_mesL2_passSA",     &probe_mesL2_passSA, &b_probe_mesL2_passSA);
   fChain->SetBranchAddress("probe_mesL2_numSpSA",    &probe_mesL2_numSpSA, &b_probe_mesL2_numSpSA);
   fChain->SetBranchAddress("probe_mesL2_ptSA",       &probe_mesL2_ptSA, &b_probe_mesL2_ptSA);
   fChain->SetBranchAddress("probe_mesL2_etaSA",      &probe_mesL2_etaSA, &b_probe_mesL2_etaSA);
   fChain->SetBranchAddress("probe_mesL2_phiSA",      &probe_mesL2_phiSA, &b_probe_mesL2_phiSA);
   fChain->SetBranchAddress("probe_mesL2_sAddressSA", &probe_mesL2_sAddressSA, &b_probe_mesL2_sAddressSA);

   fChain->SetBranchAddress("probe_mesSA_roiEta", &probe_mesSA_roiEta, &b_probe_mesSA_roiEta );
   fChain->SetBranchAddress("probe_mesSA_roiPhi", &probe_mesSA_roiPhi, &b_probe_mesSA_roiPhi );

   fChain->SetBranchAddress("probe_mesCB_pass", &probe_mesCB_pass, &b_probe_mesCB_pass);
   fChain->SetBranchAddress("probe_mesCB_dR", &probe_mesCB_dR, &b_probe_mesCB_dR);
   fChain->SetBranchAddress("probe_mesCB_tpdR", &probe_mesCB_tpdR, &b_probe_mesCB_tpdR);
   fChain->SetBranchAddress("probe_mesCB_pt", &probe_mesCB_pt, &b_probe_mesCB_pt);
   fChain->SetBranchAddress("probe_mesCB_eta", &probe_mesCB_eta, &b_probe_mesCB_eta);
   fChain->SetBranchAddress("probe_mesCB_phi", &probe_mesCB_phi, &b_probe_mesCB_phi);
   fChain->SetBranchAddress("probe_mesFTF_pass", &probe_mesFTF_pass, &b_probe_mesFTF_pass);
   fChain->SetBranchAddress("probe_mesFTF_dR", &probe_mesFTF_dR, &b_probe_mesFTF_dR);
   fChain->SetBranchAddress("probe_mesFTF_pt", &probe_mesFTF_pt, &b_probe_mesFTF_pt);
   fChain->SetBranchAddress("probe_mesFTF_eta", &probe_mesFTF_eta, &b_probe_mesFTF_eta);
   fChain->SetBranchAddress("probe_mesFTF_phi", &probe_mesFTF_phi, &b_probe_mesFTF_phi);
   fChain->SetBranchAddress("probe_mesEF_pass", &probe_mesEF_pass, &b_probe_mesEF_pass);
   fChain->SetBranchAddress("probe_mesEF_dR", &probe_mesEF_dR, &b_probe_mesEF_dR);
   fChain->SetBranchAddress("probe_mesEF_tpdR", &probe_mesEF_tpdR, &b_probe_mesEF_tpdR);
   fChain->SetBranchAddress("probe_mesEF_pt", &probe_mesEF_pt, &b_probe_mesEF_pt);
   fChain->SetBranchAddress("probe_mesEF_eta", &probe_mesEF_eta, &b_probe_mesEF_eta);
   fChain->SetBranchAddress("probe_mesEF_phi", &probe_mesEF_phi, &b_probe_mesEF_phi);

   InitHist();
   Notify();
}

void RPC::InitEffHist(){
  m_chain.push_back("PROBE");
  m_chain.push_back("L1");
  m_chain.push_back("SA");
  m_chain.push_back("CB");
  m_chain.push_back("L2");
  m_chain.push_back("EF");
  for ( int i = 0; i< m_chain.size(); i++){
    m_probe_pt_mu4.push_back(new TH1F(Form("probe_pt_mu4_%s",     m_chain[i].data()), ";p_{T}[GeV];Counts", 50, 0,    20));
    m_probe_mu_mu4.push_back(new TH1F(Form("probe_mu_mu4_%s",     m_chain[i].data()), ";<#mu>;Counts",      50, 0,    60));
    m_probe_eta_mu4.push_back(new TH1F(Form("probe_eta_mu4_%s",   m_chain[i].data()), ";#eta;Counts",       50, -3,   3));
    m_probe_phi_mu4.push_back(new TH1F(Form("probe_phi_mu4_%s",   m_chain[i].data()), ";phi;Counts",        50, -3.5, 3.5));
    m_probe_qetapt_mu4.push_back(new TH2F(Form("probe_qetapt_mu4_%s", m_chain[i].data()), ";Q#eta;p_{T}[GeV];Counts", 30, -2.5, 2.5, 30, 0, 20));
    m_probe_etaphi_mu4.push_back(new TH2F(Form("probe_etaphi_mu4_%s", m_chain[i].data()), ";#eta;#phi;Counts", 30, -2.5, 2.5, 30, -3.2, 3.2));
  }

  m_ftk_size_pt_mu4  = new TH2F("m_ftk_size_pt_mu4",  ";Offline p_{T}[GeV];# of FTKs", 50, 0, 50, 20, 0 ,10);
  m_ftk_size_eta_mu4  = new TH2F("m_ftk_size_eta_mu4",  ";Offline p_{T}[GeV];# of FTKs", 50, 0, 50, 20, -2.5, 2.5);

  m_ftk_size_pt_mu50  = new TH2F("m_ftk_size_pt_mu50",  ";Offline p_{T}[GeV];# of FTKs", 100, 0, 100, 20, 0 ,10);
  m_ftk_size_eta_mu50  = new TH2F("m_ftk_size_eta_mu50",  ";Offline p_{T}[GeV];# of FTKs", 100, 0, 100, 20, -2.5, 2.5);

  m_eff_pt_mu4_L1  = new TH1F("m_eff_pt_mu4_L1",  ";p_{T}[GeV];L1 Efficiency", 50, 0,    20);
  m_eff_mu_mu4_L1  = new TH1F("m_eff_mu_mu4_L1",  ";<#mu>;L1 Efficiency",      50, 0,    60);
  m_eff_eta_mu4_L1 = new TH1F("m_eff_eta_mu4_L1", ";#eta;L1 Efficiency",       50, -3,   3);
  m_eff_phi_mu4_L1 = new TH1F("m_eff_phi_mu4_L1", ";#phi;L1 Efficiency",       50, -3.5, 3.5);

  m_eff_pt_mu4_L1SA  = new TH1F("m_eff_pt_mu4_L1SA",  ";p_{T}[GeV];SA/L1 Efficiency", 50, 0,    20);
  m_eff_mu_mu4_L1SA  = new TH1F("m_eff_mu_mu4_L1SA",  ";<#mu>;SA/L1 Efficiency",      50, 0,    60);
  m_eff_eta_mu4_L1SA = new TH1F("m_eff_eta_mu4_L1SA", ";#eta;SA/L1 Efficiency",       50, -3,   3);
  m_eff_phi_mu4_L1SA = new TH1F("m_eff_phi_mu4_L1SA", ";#phi;SA/L1 Efficiency",       50, -3.5, 3.5);

  m_eff_pt_mu4_SACB  = new TH1F("m_eff_pt_mu4_SACB",  ";p_{T}[GeV];CB/SA Efficiency", 50, 0,    20);
  m_eff_mu_mu4_SACB  = new TH1F("m_eff_mu_mu4_SACB",  ";<#mu>;CB/SA Efficiency",      50, 0,    60);
  m_eff_eta_mu4_SACB = new TH1F("m_eff_eta_mu4_SACB", ";#eta;CB/SA Efficiency",       50, -3,   3);
  m_eff_phi_mu4_SACB = new TH1F("m_eff_phi_mu4_SACB", ";#phi;CB/SA Efficiency",       50, -3.5, 3.5);

  m_eff_pt_mu4_L1L2  = new TH1F("m_eff_pt_mu4_L1L2",  ";p_{T}[GeV];L2/L1 Efficiency", 50, 0,    20);
  m_eff_mu_mu4_L1L2  = new TH1F("m_eff_mu_mu4_L1L2",  ";<#mu>;L2/L1 Efficiency",      50, 0,    60);
  m_eff_eta_mu4_L1L2 = new TH1F("m_eff_eta_mu4_L1L2", ";#eta;L2/L1 Efficiency",       50, -3,   3);
  m_eff_phi_mu4_L1L2 = new TH1F("m_eff_phi_mu4_L1L2", ";#phi;L2/L1 Efficiency",       50, -3.5, 3.5);

  m_eff_qetapt_mu4_L1L2 = new TH2F("m_eff_qetapt_mu4_L1L2", ";Q#eta;p_{T}[GeV];L2/L1 Efficiency", 30, -2.5, 2.5, 30, 0, 20);
  m_eff_etaphi_mu4_L1L2 = new TH2F("m_eff_etaphi_mu4_L1L2", ";#eta;#phi;L2/L1 Efficiency",        30, -2.5, 2.5, 30, 0, 20);

  m_eff_pt_mu4_L1CB  = new TH1F("m_eff_pt_mu4_L1CB",  ";p_{T}[GeV];CB/L1 Efficiency", 50, 0,    20);
  m_eff_mu_mu4_L1CB  = new TH1F("m_eff_mu_mu4_L1CB",  ";<#mu>;CB/L1 Efficiency",      50, 0,    60);
  m_eff_eta_mu4_L1CB = new TH1F("m_eff_eta_mu4_L1CB", ";#eta;CB/L1 Efficiency",       50, -3,   3);
  m_eff_phi_mu4_L1CB = new TH1F("m_eff_phi_mu4_L1CB", ";#phi;CB/L1 Efficiency",       50, -3.5, 3.5);

  h_probe_mu_mu4_L1  = new TH1F("h_probe_mu_mu4_L1", "probe;<#mu>;Entries",  30, 0, 60);
  h_probe_mu_mu4_SA  = new TH1F("h_probe_mu_mu4_SA", "probe;<#mu>;Entries",  30, 0, 60);
  h_eff_mu_mu4_L1    = new TH1F("h_eff_mu_mu4_L1",   "eff;<#mu>;Efficiency", 30, 0, 60);
  h_eff_mu_mu4_L1SA  = new TH1F("h_eff_mu_mu4_L1SA", "eff;<#mu>;Efficiency", 30, 0, 60);

  h_probe_pt_mu4_PROBE  = new TH1F("h_probe_pt_mu4_PROBE", "probe;Probe p_{T}[GeV];Entries",  50, 0, 20);
  h_probe_pt_mu4_L1  = new TH1F("h_probe_pt_mu4_L1",       "probe;Probe p_{T}[GeV];Entries",  50, 0, 20);
  h_probe_pt_mu4_SA  = new TH1F("h_probe_pt_mu4_SA",       "probe;Probe p_{T}[GeV];Entries",  50, 0, 20);
  h_probe_pt_mu4_CB  = new TH1F("h_probe_pt_mu4_CB",       "probe;Probe p_{T}[GeV];Entries",  50, 0, 20);
  h_probe_pt_mu4_EF  = new TH1F("h_probe_pt_mu4_EF",       "probe;Probe p_{T}[GeV];Entries",  50, 0, 20);
  h_eff_pt_mu4_L1    = new TH1F("h_eff_pt_mu4_L1",         "eff;Probe p_{T}[GeV];L1 Efficiency", 50, 0, 20);
  h_eff_pt_mu4_L1SA  = new TH1F("h_eff_pt_mu4_L1SA",       "eff;Probe p_{T}[GeV];SA/L1 Efficiency", 50, 0, 20);
  h_eff_pt_mu4_SACB  = new TH1F("h_eff_pt_mu4_SACB",       "eff;Probe p_{T}[GeV];CB/SA Efficiency", 50, 0, 20);
  h_eff_pt_mu4_CBEF  = new TH1F("h_eff_pt_mu4_CBEF",       "eff;Probe p_{T}[GeV];EF/CBEfficiency", 50, 0, 20);

  h_probe_phi_mu4_L1  = new TH1F("h_probe_phi_mu4_L1", "probe;Probe #phi;Entries", 50, -3.5, 3.5);
  h_probe_phi_mu4_SA  = new TH1F("h_probe_phi_mu4_SA", "probe;Probe #phi;Entries", 50, -3.5, 3.5);
  h_eff_phi_mu4_L1    = new TH1F("h_eff_phi_mu4_L1",   "eff;Probe #phi;Efficiency", 50, -3.5, 3.5);
  h_eff_phi_mu4_L1SA  = new TH1F("h_eff_phi_mu4_L1SA", "eff;Probe #phi;Efficiency", 50, -3.5, 3.5);

  h_probe_eta_mu4_L1  = new TH1F("h_probe_eta_mu4_L1", "probe;Probe #eta;Entries", 50, -2.5, 2.5);
  h_probe_eta_mu4_SA  = new TH1F("h_probe_eta_mu4_SA", "probe;Probe #eta;Entries", 50, -2.5, 2.5);
  h_eff_eta_mu4_L1    = new TH1F("h_eff_eta_mu4_L1",   "eff;Probe #eta;Efficiency", 50, -2.5, 2.5);
  h_eff_eta_mu4_L1SA  = new TH1F("h_eff_eta_mu4_L1SA", "eff;Probe #eta;Efficiency", 50, -2.5, 2.5);

  hh_probe_etaphi_mu4_L1 = new TH2F("hh_probe_etaphi_mu4_L1", "probe;Probe #eta;Probe #phi;Entries",  30, -2.5, 2.5, 30, -3.5, 3.5);
  hh_probe_etaphi_mu4_SA = new TH2F("hh_probe_etaphi_mu4_SA", "probe;Probe #eta;Probe #phi;Entries",  30, -2.5, 2.5, 30, -3.5, 3.5);
  hh_eff_etaphi_mu4_L1   = new TH2F("hh_eff_etaphi_mu4_L1",   "eff;Probe #eta;Probe #phi;Efficiency", 30, -2.5, 2.5, 30, -3.5, 3.5);
  hh_eff_etaphi_mu4_L1SA = new TH2F("hh_eff_etaphi_mu4_L1SA", "eff;Probe #eta;Probe #phi;Efficiency", 30, -2.5, 2.5, 30, -3.5, 3.5);

  hh_probe_qetapt_mu4_L1 = new TH2F("hh_probe_qetapt_mu4_L1", "probe;Probe Q#eta;Probe p_{T}[GeV];Entries",  30, -2.5, 2.5, 30, 0, 14);
  hh_probe_qetapt_mu4_SA = new TH2F("hh_probe_qetapt_mu4_SA", "probe;Probe Q#eta;Probe p_{T}[GeV];Entries",  30, -2.5, 2.5, 30, 0, 14);
  hh_eff_qetapt_mu4_L1   = new TH2F("hh_eff_qetapt_mu4_L1",   "eff;Probe Q#eta;Probe p_{T}[GeV];Efficiency", 30, -2.5, 2.5, 30, 0, 14);
  hh_eff_qetapt_mu4_L1SA = new TH2F("hh_eff_qetapt_mu4_L1SA", "eff;Probe Q#eta;Probe p_{T}[GeV];Efficiency", 30, -2.5, 2.5, 30, 0, 14);

  h_probe_mu_mu50_L1  = new TH1F("h_probe_mu_mu50_L1", "probe;<#mu>;Entries",  30, 0, 60);
  h_probe_mu_mu50_SA  = new TH1F("h_probe_mu_mu50_SA", "probe;<#mu>;Entries",  30, 0, 60);
  h_probe_mu_mu50_CB  = new TH1F("h_probe_mu_mu50_CB", "probe;<#mu>;Entries",  30, 0, 60);
  h_eff_mu_mu50_L1    = new TH1F("h_eff_mu_mu50_L1",   "eff;<#mu>;Efficiency", 30, 0, 60);
  h_eff_mu_mu50_L1SA  = new TH1F("h_eff_mu_mu50_L1SA", "eff;<#mu>;Efficiency", 30, 0, 60);
  h_eff_mu_mu50_SACB  = new TH1F("h_eff_mu_mu50_SACB", "eff;<#mu>;Efficiency", 50, 0, 60);

  h_probe_pt_mu50_L1  = new TH1F("h_probe_pt_mu50_L1", "probe;Probe p_{T}[GeV];Entries",  50, 0, 100);
  h_probe_pt_mu50_SA  = new TH1F("h_probe_pt_mu50_SA", "probe;Probe p_{T}[GeV];Entries",  50, 0, 100);
  h_probe_pt_mu50_CB  = new TH1F("h_probe_pt_mu50_CB", "probe;Probe p_{T}[GeV];Entries",  50, 0, 100);
  h_eff_pt_mu50_L1    = new TH1F("h_eff_pt_mu50_L1",   "eff;Probe p_{T}[GeV];Efficiency", 50, 0, 100);
  h_eff_pt_mu50_L1SA  = new TH1F("h_eff_pt_mu50_L1SA", "eff;Probe p_{T}[GeV];Efficiency", 50, 0, 100);
  h_eff_pt_mu50_SACB  = new TH1F("h_eff_pt_mu50_SACB", "eff;Probe p_{T}[GeV];Efficiency", 50, 0, 100);

  h_probe_phi_mu50_L1  = new TH1F("h_probe_phi_mu50_L1", "probe;Probe #phi;Entries", 50, -3.5, 3.5);
  h_probe_phi_mu50_SA  = new TH1F("h_probe_phi_mu50_SA", "probe;Probe #phi;Entries", 50, -3.5, 3.5);
  h_probe_phi_mu50_CB  = new TH1F("h_probe_phi_mu50_CB", "probe;Probe #phi;Entries", 50, -3.5, 3.5);
  h_eff_phi_mu50_L1    = new TH1F("h_eff_phi_mu50_L1",   "eff;Probe #phi;Efficiency", 50, -3.5, 3.5);
  h_eff_phi_mu50_L1SA  = new TH1F("h_eff_phi_mu50_L1SA", "eff;Probe #phi;Efficiency", 50, -3.5, 3.5);
  h_eff_phi_mu50_SACB  = new TH1F("h_eff_phi_mu50_SACB", "eff;Probe #phi;Efficiency", 50, -3.5, 3.5);

  h_probe_eta_mu50_L1  = new TH1F("h_probe_eta_mu50_L1", "probe;Probe #eta;Entries", 50, -2.5, 2.5);
  h_probe_eta_mu50_SA  = new TH1F("h_probe_eta_mu50_SA", "probe;Probe #eta;Entries", 50, -2.5, 2.5);
  h_probe_eta_mu50_CB  = new TH1F("h_probe_eta_mu50_CB", "probe;Probe #eta;Entries", 50, -2.5, 2.5);
  h_eff_eta_mu50_L1    = new TH1F("h_eff_eta_mu50_L1",   "eff;Probe #eta;Efficiency", 50, -2.5, 2.5);
  h_eff_eta_mu50_L1SA    = new TH1F("h_eff_eta_mu50_L1SA",   "eff;Probe #eta;Efficiency", 50, -2.5, 2.5);
  h_eff_eta_mu50_SACB  = new TH1F("h_eff_eta_mu50_SACB", "eff;Probe #eta;Efficiency", 50, -2.5, 2.5);

  hh_probe_etaphi_mu50_L1 = new TH2F("hh_probe_etaphi_mu50_L1", "probe;Probe #eta;Probe #phi;Entries", 30, -2.5, 2.5, 30, -3.5, 3.5);
  hh_probe_etaphi_mu50_SA = new TH2F("hh_probe_etaphi_mu50_SA", "probe;Probe #eta;Probe #phi;Entries", 30, -2.5, 2.5, 30, -3.5, 3.5);
  hh_eff_etaphi_mu50_L1   = new TH2F("hh_eff_etaphi_mu50_L1",   "eff;Probe #eta;Probe #phi;Efficiency",   30, -2.5, 2.5, 30, -3.5, 3.5);
  hh_eff_etaphi_mu50_L1SA = new TH2F("hh_eff_etaphi_mu50_L1SA", "eff;Probe #eta;Probe #phi;Efficiency",   30, -2.5, 2.5, 30, -3.5, 3.5);

  hh_probe_qetapt_mu50_L1 = new TH2F("hh_probe_qetapt_mu50_L1", "probe;Probe Q#eta;Probe p_{T}[GeV];Entries", 30, -2.5, 2.5, 30, 0, 100);
  hh_probe_qetapt_mu50_SA = new TH2F("hh_probe_qetapt_mu50_SA", "probe;Probe Q#eta;Probe p_{T}[GeV];Entries", 30, -2.5, 2.5, 30, 0, 100);
  hh_eff_qetapt_mu50_L1   = new TH2F("hh_eff_qetapt_mu50_L1",   "eff;Probe Q#eta;Probe p_{T}[GeV];Efficiency",   30, -2.5, 2.5, 30, 0, 100);
  hh_eff_qetapt_mu50_L1SA = new TH2F("hh_eff_qetapt_mu50_L1SA", "eff;Probe Q#eta;Probe p_{T}[GeV];Efficiency",   30, -2.5, 2.5, 30, 0, 100);

}

void RPC::InitSPHist(){
  h_NumberOfSP_LumiBlock   = new TH2F("h_NumberOfSP_LumiBlock",   "h_NumberOfSP_LumiBlock;LumiBlock;Number;Counts", 100, 0,    500, 100, 0, 5);
  h_NumberOfSP_eta         = new TH2F("h_NumberOfSP_eta",         "h_NumberOfSP_eta;eta;Number;Counts",             50,  -2.5, 2.5, 100, 0, 5);
  h_NumberOfSP_qeta        = new TH2F("h_NumberOfSP_qeta",        "h_NumberOfSP_qeta;qeta;Number;Counts",           50,  -2.5, 2.5, 100, 0, 5);
  h_NumberOfSP_pt_barrel   = new TH2F("h_NumberOfSP_pt_barrel",   "h_NumberOfSP_pt_barrel;pT[GeV];Number;Counts",   50,  0,    20, 100, 0, 5);
  h_NumberOfSP_pass_barrel = new TH2F("h_NumberOfSP_pass_barrel", "h_NumberOfSP_pass_barrel;pass;Number;Counts",    50,  -3,   2,   100, 0, 5);
}

void RPC::InitMdtHist(){
  h_NumberOfMdt_LumiBlock    = new TH2F("h_NumberOfMdt_LumiBlock",    "h_NumberOfMdt_LumiBlock;LumiBlock;Number;Counts",    100, 0, 500, 100, 0, 50);
  h_NumberOfMdt_LumiBlock_BI = new TH2F("h_NumberOfMdt_LumiBlock_BI", "h_NumberOfMdt_LumiBlock_BI;LumiBlock;Number (BI);Counts", 100, 0, 500, 100, 0, 50);
  h_NumberOfMdt_LumiBlock_BM = new TH2F("h_NumberOfMdt_LumiBlock_BM", "h_NumberOfMdt_LumiBlock_BM;LumiBlock;Number (BM);Counts", 100, 0, 500, 100, 0, 50);
  h_NumberOfMdt_LumiBlock_BO = new TH2F("h_NumberOfMdt_LumiBlock_BO", "h_NumberOfMdt_LumiBlock_BO;LumiBlock;Number (BO);Counts", 100, 0, 500, 100, 0, 50);

  h_NumberOfMdt_eta    = new TH2F("h_NumberOfMdt_eta",                "h_NumberOfMdt_eta;eta;Number;Counts",              100, -2.5, 2.5, 100, 0, 50);
  h_NumberOfMdt_eta_BI = new TH2F("h_NumberOfMdt_eta_BI",             "h_NumberOfMdt_eta_BI;eta;Number (BI);Counts",      100, -1.05, 1.05, 100, 0, 50);
  h_NumberOfMdt_eta_BM = new TH2F("h_NumberOfMdt_eta_BM",             "h_NumberOfMdt_eta_BM;eta;Number (BM);Counts",      100, -1.05, 1.05, 100, 0, 50);
  h_NumberOfMdt_eta_BO = new TH2F("h_NumberOfMdt_eta_BO",             "h_NumberOfMdt_eta_BO;eta;Number (BO);Counts",      100, -1.05, 1.05, 100, 0, 50);

  h_NumberOfMdt_qeta    = new TH2F("h_NumberOfMdt_qeta",              "h_NumberOfMdt_qeta;qeta;Number;Counts",            100, -1.05, 1.05, 100, 0, 50);
  h_NumberOfMdt_qeta_BI = new TH2F("h_NumberOfMdt_qeta_BI",           "h_NumberOfMdt_qeta_BI;qeta;Number (BI);Counts",    100, -1.05, 1.05, 100, 0, 50);
  h_NumberOfMdt_qeta_BM = new TH2F("h_NumberOfMdt_qeta_BM",           "h_NumberOfMdt_qeta_BM;qeta;Number (BM);Counts",    100, -1.05, 1.05, 100, 0, 50);
  h_NumberOfMdt_qeta_BO = new TH2F("h_NumberOfMdt_qeta_BO",           "h_NumberOfMdt_qeta_BO;qeta;Number (BO);Counts",    100, -1.05, 1.05, 100, 0, 50);

  h_NumberOfMdt_pt_barrel    = new TH2F("h_NumberOfMdt_pt_barrel",    "h_NumberOfMdt_pt_barrel;pT;Number;Counts",         100, 0,    100, 100, 0, 50);
  h_NumberOfMdt_pt_barrel_BI = new TH2F("h_NumberOfMdt_pt_barrel_BI", "h_NumberOfMdt_pt_barrel_BI;pT (BI);Number;Counts", 100, 0,    100, 100, 0, 50);
  h_NumberOfMdt_pt_barrel_BM = new TH2F("h_NumberOfMdt_pt_barrel_BM", "h_NumberOfMdt_pt_barrel_BM;pT (BM);Number;Counts", 100, 0,    100, 100, 0, 50);
  h_NumberOfMdt_pt_barrel_BO = new TH2F("h_NumberOfMdt_pt_barrel_BO", "h_NumberOfMdt_pt_barrel_BO;pT (BO);Number;Counts", 100, 0,    100, 100, 0, 50);

  h_ResidualMdt_Outlier_eta    = new TH2F("h_ResidualMdt_Outlier_eta",    "h_ResidualMdt_Outlier_eta;#eta;Residual;Counts",                              50, -1.05, 1.05, 50, -1000, 1000);
  h_ResidualMdt_Outlier_eta_BI = new TH2F("h_ResidualMdt_Outlier_eta_BI", "h_ResidualMdt_Outlier_eta_BI;#eta;(BI) Outlier MDT hit residual [mm];Counts", 50, -1.05, 1.05, 50, -1000, 1000);
  h_ResidualMdt_Outlier_eta_BM = new TH2F("h_ResidualMdt_Outlier_eta_BM", "h_ResidualMdt_Outlier_eta_BM;#eta;(BM) Outlier MDT hit residual [mm];Counts", 50, -1.05, 1.05, 50, -1000, 1000);
  h_ResidualMdt_Outlier_eta_BO = new TH2F("h_ResidualMdt_Outlier_eta_BO", "h_ResidualMdt_Outlier_eta_BO;#eta;(BO) Outlier MDT hit residual [mm];Counts", 50, -1.05, 1.05, 50, -1000, 1000);

  h_ResidualMdt_Inlier_eta    = new TH2F("h_ResidualMdt_Inlier_eta",    "h_ResidualMdt_Inlier_eta;#eta;Residual;Counts",                             50, -1.05, 1.05, 50, -1000, 1000);
  h_ResidualMdt_Inlier_eta_BI = new TH2F("h_ResidualMdt_Inlier_eta_BI", "h_ResidualMdt_Inlier_eta_BI;#eta;(BI) Inlier MDT hit residual [mm];Counts", 50, -1.05, 1.05, 50, -1000, 1000);
  h_ResidualMdt_Inlier_eta_BM = new TH2F("h_ResidualMdt_Inlier_eta_BM", "h_ResidualMdt_Inlier_eta_BM;#eta;(BM) Inlier MDT hit residual [mm];Counts", 50, -1.05, 1.05, 50, -1000, 1000);
  h_ResidualMdt_Inlier_eta_BO = new TH2F("h_ResidualMdt_Inlier_eta_BO", "h_ResidualMdt_Inlier_eta_BO;#eta;(BO) Inlier MDT hit residual [mm];Counts", 50, -1.05, 1.05, 50, -1000, 1000);

  h_ResidualMdt_Outlier_pt_barrel    = new TH2F("h_ResidualMdt_Outlier_pt_barrel",    "h_ResidualMdt_Outlier_pt_barrel;p_{T}[GeV];Outlier MDT hit residual [mm];Counts",         50, 0, 100, 50, -1000, 1000);
  h_ResidualMdt_Outlier_pt_barrel_BI = new TH2F("h_ResidualMdt_Outlier_pt_barrel_BI", "h_ResidualMdt_Outlier_pt_barrel_BI;p_{T}[GeV];(BI) Outlier MDT hit residual [mm];Counts", 50, 0, 100, 50, -1000, 1000);
  h_ResidualMdt_Outlier_pt_barrel_BM = new TH2F("h_ResidualMdt_Outlier_pt_barrel_BM", "h_ResidualMdt_Outlier_pt_barrel_BM;p_{T}[GeV];(BM) Outlier MDT hit residual [mm];Counts", 50, 0, 100, 50, -1000, 1000);
  h_ResidualMdt_Outlier_pt_barrel_BO = new TH2F("h_ResidualMdt_Outlier_pt_barrel_BO", "h_ResidualMdt_Outlier_pt_barrel_BO;p_{T}[GeV];(BO) Outlier MDT hit residual [mm];Counts", 50, 0, 100, 50, -1000, 1000);

  h_ResidualMdt_Inlier_pt_barrel    = new TH2F("h_ResidualMdt_Inlier_pt_barrel",    "h_ResidualMdt_Inlier_pt_barrel;p_{T}[GeV];Inlier MDT hit residual [mm];Counts",         50, 0, 100, 50, -1000, 1000);
  h_ResidualMdt_Inlier_pt_barrel_BI = new TH2F("h_ResidualMdt_Inlier_pt_barrel_BI", "h_ResidualMdt_Inlier_pt_barrel_BI;p_{T}[GeV];(BI) Inlier MDT hit residual [mm];Counts", 50, 0, 100, 50, -1000, 1000);
  h_ResidualMdt_Inlier_pt_barrel_BM = new TH2F("h_ResidualMdt_Inlier_pt_barrel_BM", "h_ResidualMdt_Inlier_pt_barrel_BM;p_{T}[GeV];(BM) Inlier MDT hit residual [mm];Counts", 50, 0, 100, 50, -1000, 1000);
  h_ResidualMdt_Inlier_pt_barrel_BO = new TH2F("h_ResidualMdt_Inlier_pt_barrel_BO", "h_ResidualMdt_Inlier_pt_barrel_BO;pT [GeV];(BO) Inlier MDT hit residual [mm];Counts",   50, 0, 100, 50, -1000, 1000);

}

void RPC::InitPtResidualHist(){
  h_PtResidual_pt1  = new TH2F("h_PtResidual_pt1",       "h_PtResidual_pt1;Probe p_{T} [GeV];p_{T} residual;Counts",                 100, 0,    20,  50,  -1.1, 1.1);
  h_PtResidual_pt2  = new TH2F("h_PtResidual_pt2",       "h_PtResidual_pt2;Probe p_{T} [GeV];p_{T} residual;Counts",                 100, 0,    20,  50,  -1.1, 1.1);
  h_PtResidual_pt3  = new TH2F("h_PtResidual_pt3",       "h_PtResidual_pt3;Probe p_{T} [GeV];p_{T} residual;Counts",                 100, 0,    20,  50,  -1.1, 1.1);
  h_PtResidual_ptFtk  = new TH2F("h_PtResidual_ptFtk", "h_PtResidual_ptFtk;Probe p_{T} [GeV];p_{T} residual;Counts",              100, 0,    20,  50,  -1.1, 1.1);
  h_PtResidual_ptLut  = new TH2F("h_PtResidual_ptLut", "h_PtResidual_ptLut;Probe p_{T} [GeV];p_{T} residual;Counts",              100, 0,    20,  50,  -1.1, 1.1);
  h_PtResidual_eta = new TH2F("h_PtResidual_eta",      "h_PtResidual_eta;Probe #eta;p_{T} residual;Counts",                       100, -2.5, 2.5, 50,  -1.1, 1.1);
  h_pt_vs_pt       = new TH2F("h_pt_vs_pt",            "h_pt_vs_pt;Probe Offline p_{T} [GeV];Probe L2MuonSA p_{T} [GeV];Counts",  100, 0,    20,  100, 0,  20);
  h_ptl2_vs_ptFtk       = new TH2F("h_ptl2_vs_ptFtk",  "h_ptl2_vs_ptFtk;Probe L2MuonSA p_{T} [GeV];Probe FTK p_{T} [GeV];Counts", 100, 0,    20,  100, 0,  20);
  h_pt_vs_ptFtk       = new TH2F("h_pt_vs_ptFtk",      "h_pt_vs_ptFtk;Probe Offline p_{T} [GeV];Probe FTK p_{T} [GeV];Counts",    100, 0,    20,  100, 0,  20);
}

void RPC::InitInEffHist(){
  h_InEff_pt     = new TH2F("h_InEff_pt",     "h_InEff_pt;Probe p_{T} [GeV];;Counts",                40, 0,    14,  50,  -20,  20);
  h_InEff_eta    = new TH2F("h_InEff_eta",    "h_InEff_eta;Probe #eta;;Counts",                      40, -2.5, 2.5, 50,  -20,  20);
  h_InEff_qeta   = new TH2F("h_InEff_qeta",   "h_InEff_qeta;Probe Q#eta;;Counts",                    40, -2.5, 2.5, 100, -20,  20);
  h_InEff_qetapt = new TH3F("h_InEff_qetapt", "h_InEff_qetapt;Probe p_{T} [GeV];Probe #eta;;Counts", 40, 0,    14,  100, -2.5, 2.5, 50, -20, 20);
}


void RPC::InitHist(){
  InitEffHist();
  InitInEffHist();
  InitSPHist();
  InitMdtHist();
  InitPtResidualHist();

  // Histgrams
  h_superPointRZ_BIS = new TH2F("h_superPointRZ_BIS", "h_superPointRZ_BIS;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_superPointRZ_BIL = new TH2F("h_superPointRZ_BIL", "h_superPointRZ_BIL;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_superPointRZ_BMS = new TH2F("h_superPointRZ_BMS", "h_superPointRZ_BMS;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_superPointRZ_BML = new TH2F("h_superPointRZ_BML", "h_superPointRZ_BML;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_superPointRZ_BOS = new TH2F("h_superPointRZ_BOS", "h_superPointRZ_BOS;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_superPointRZ_BOL = new TH2F("h_superPointRZ_BOL", "h_superPointRZ_BOL;Z;R;Counts", 100, -15, 15, 100, 0, 20);

  h_segmentRZ_BIS = new TH2F("h_segmentRZ_BIS", "h_segmentRZ_BIS;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_segmentRZ_BIL = new TH2F("h_segmentRZ_BIL", "h_segmentRZ_BIL;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_segmentRZ_BMS = new TH2F("h_segmentRZ_BMS", "h_segmentRZ_BMS;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_segmentRZ_BML = new TH2F("h_segmentRZ_BML", "h_segmentRZ_BML;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_segmentRZ_BOS = new TH2F("h_segmentRZ_BOS", "h_segmentRZ_BOS;Z;R;Counts", 100, -15, 15, 100, 0, 20);
  h_segmentRZ_BOL = new TH2F("h_segmentRZ_BOL", "h_segmentRZ_BOL;Z;R;Counts", 100, -15, 15, 100, 0, 20);

  h_residualRZ_BIS = new TH2F("h_residualRZ_BIS", "h_residualRZ_BIS;Z;R;Counts", 100, -0.1, 0.1, 100, -0.1, 0.1);
  h_residualRZ_BIL = new TH2F("h_residualRZ_BIL", "h_residualRZ_BIL;Z;R;Counts", 100, -0.1, 0.1, 100, -0.1, 0.1);
  h_residualRZ_BMS = new TH2F("h_residualRZ_BMS", "h_residualRZ_BMS;Z;R;Counts", 100, -0.1, 0.1, 100, -0.1, 0.1);
  h_residualRZ_BML = new TH2F("h_residualRZ_BML", "h_residualRZ_BML;Z;R;Counts", 100, -0.1, 0.1, 100, -0.1, 0.1);
  h_residualRZ_BOS = new TH2F("h_residualRZ_BOS", "h_residualRZ_BOS;Z;R;Counts", 100, -0.1, 0.1, 100, -0.1, 0.1);
  h_residualRZ_BOL = new TH2F("h_residualRZ_BOL", "h_residualRZ_BOL;Z;R;Counts", 100, -0.1, 0.1, 100, -0.1, 0.1);

  h_ResidualSegment_eta        = new TH2F("h_ResidualSegment_eta",        "h_ResidualSegment_eta;#eta;Residual;Counts",                          100, -1.05, 1.05, 100, -200, 200);
  h_ResidualSegment_eta_BI     = new TH2F("h_ResidualSegment_eta_BI",     "h_ResidualSegment_eta_BI;#eta;(BI) Segment residual [mm];Counts",     100, -1.05, 1.05, 100, -200, 200);
  h_ResidualSegment_eta_BI_2d     = new TH2F("h_ResidualSegment_eta_BI_2d",     "h_ResidualSegment_eta_BI_2d;(BI) RPC road residual;(BI) FTK road residual;Counts",     100, -200, 200, 100, -200, 200);
  h_ResidualSegment_eta_BM     = new TH2F("h_ResidualSegment_eta_BM",     "h_ResidualSegment_eta_BM;#eta;(BM) Segment residual [mm];Counts",     100, -1.05, 1.05, 100, -200, 200);
  h_ResidualSegment_eta_BO     = new TH2F("h_ResidualSegment_eta_BO",     "h_ResidualSegment_eta_BO;#eta;(BO) Segment residual [mm];Counts",     100, -1.05, 1.05, 100, -200, 200);

  h_ResidualSegment_eta_rpc_BI = new TH1F("h_ResidualSegment_eta_rpc_BI", "h_ResidualSegment_eta_rpc_BI;(BI) Segment residual [mm];Counts", 200, -1000, 1000);
  h_ResidualSegment_eta_rpc_BM = new TH1F("h_ResidualSegment_eta_rpc_BM", "h_ResidualSegment_eta_rpc_BM;(BM) Segment residual [mm];Counts", 200, -1000, 1000);
  h_ResidualSegment_eta_rpc_BO = new TH1F("h_ResidualSegment_eta_rpc_BO", "h_ResidualSegment_eta_rpc_BO;(BO) Segment residual [mm];Counts", 200, -1000, 1000);

  h_ResidualSegment_eta_ftk_BI = new TH1F("h_ResidualSegment_eta_ftk_BI", "h_ResidualSegment_eta_ftk_BI;(BI) Segment residual [mm];Counts", 200, -1000, 1000);
  h_ResidualSegment_eta_ftk_BM = new TH1F("h_ResidualSegment_eta_ftk_BM", "h_ResidualSegment_eta_ftk_BM;(BM) Segment residual [mm];Counts", 200, -1000, 1000);
  h_ResidualSegment_eta_ftk_BO = new TH1F("h_ResidualSegment_eta_ftk_BO", "h_ResidualSegment_eta_ftk_BO;(BO) Segment residual [mm];Counts", 200, -1000, 1000);
}

void RPC::End(){
  if(h_superPointRZ_BIL != 0) {
    delete h_superPointRZ_BIL; h_superPointRZ_BIL = 0;
  }
  if(h_superPointRZ_BIS != 0) {
    delete h_superPointRZ_BIS; h_superPointRZ_BIS = 0;
  }
  if(h_superPointRZ_BML != 0) {
    delete h_superPointRZ_BML; h_superPointRZ_BML = 0;
  }
  if(h_superPointRZ_BMS != 0) {
    delete h_superPointRZ_BMS; h_superPointRZ_BMS = 0;
  }
  if(h_superPointRZ_BOL != 0) {
    delete h_superPointRZ_BOL; h_superPointRZ_BOL = 0;
  }
  if(h_superPointRZ_BOS != 0) {
    delete h_superPointRZ_BOS; h_superPointRZ_BOS = 0;
  }

  if(h_segmentRZ_BIL != 0) {
    delete h_segmentRZ_BIL; h_segmentRZ_BIL = 0;
  }
  if(h_segmentRZ_BIS != 0) {
    delete h_segmentRZ_BIS; h_segmentRZ_BIS = 0;
  }
  if(h_segmentRZ_BML != 0) {
    delete h_segmentRZ_BML; h_segmentRZ_BML = 0;
  }
  if(h_segmentRZ_BMS != 0) {
    delete h_segmentRZ_BMS; h_segmentRZ_BMS = 0;
  }
  if(h_segmentRZ_BOL != 0) {
    delete h_segmentRZ_BOL; h_segmentRZ_BOL = 0;
  }
  if(h_segmentRZ_BOS != 0) {
    delete h_segmentRZ_BOS; h_segmentRZ_BOS = 0;
  }

  if(h_residualRZ_BIL != 0) {
    delete h_residualRZ_BIL; h_residualRZ_BIL = 0;
  }
  if(h_residualRZ_BIS != 0) {
    delete h_residualRZ_BIS; h_residualRZ_BIS = 0;
  }
  if(h_residualRZ_BML != 0) {
    delete h_residualRZ_BML; h_residualRZ_BML = 0;
  }
  if(h_residualRZ_BMS != 0) {
    delete h_residualRZ_BMS; h_residualRZ_BMS = 0;
  }
  if(h_residualRZ_BOL != 0) {
    delete h_residualRZ_BOL; h_residualRZ_BOL = 0;
  }
  if(h_residualRZ_BOS != 0) {
    delete h_residualRZ_BOS; h_residualRZ_BOS = 0;
  }

  if(h_NumberOfSP_eta != 0) {
    delete h_NumberOfSP_eta; h_NumberOfSP_eta = 0;
  }
  if(h_NumberOfSP_qeta != 0) {
    delete h_NumberOfSP_qeta; h_NumberOfSP_qeta = 0;
  }
  if(h_NumberOfSP_pt_barrel != 0) {
    delete h_NumberOfSP_pt_barrel; h_NumberOfSP_pt_barrel = 0;
  }
  if(h_NumberOfSP_pass_barrel != 0) {
    delete h_NumberOfSP_pass_barrel; h_NumberOfSP_pass_barrel = 0;
  }

  if(h_NumberOfMdt_eta != 0) {
    delete h_NumberOfMdt_eta; h_NumberOfMdt_eta = 0;
  }
  if(h_NumberOfMdt_qeta != 0) {
    delete h_NumberOfMdt_qeta; h_NumberOfMdt_qeta = 0;
  }
  if(h_NumberOfMdt_pt_barrel != 0) {
    delete h_NumberOfMdt_pt_barrel; h_NumberOfMdt_pt_barrel = 0;
  }
}

Bool_t RPC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RPC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}


Int_t RPC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}



#endif // #ifdef RPC_cxx
