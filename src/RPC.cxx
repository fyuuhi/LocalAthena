#define RPC_cxx
#include "RPC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TFile.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include <vector>
#include <iostream>

using namespace std;
//declare functions


#include "/gpfs/home/yfukuhar/RootUtils/rootlogon.C"

double Res(double param1, double param2);
bool GRLlist(int LumiBlock);


//==================================================================
//main function
//==================================================================
int main(int argc, char **argv){
  rootlogon();
  cout << "---start---" << endl;
  TColor::InvertPalette();
  // tree1
  TChain *tree1 = new TChain("t_tap", "t_tap");
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/data/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_jpzYFV4_GRL_F_tree_v1_349533_EXT0/hadd/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_jpzYFV4_GRL_F_tree_v1_349533_EXT0.root");
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/data/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v1_EXT0/hadd/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v1_EXT.root");
  tree1 -> Add("/gpfs/fs2001/yfukuhar/data/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v2_EXT0/hadd/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v2_EXT.root");
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/data/data18_0621/data18_0621.root");
  //tree1 -> Add("/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/");
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16c_Jpsimu6_default/mc16c_Jpsimu6_default.root");

  RPC t_349014(tree1); 

  t_349014.Loop(-1, 10000);
  cout << "[INFO]: Loop SUCCESS" << endl;

  //t_349014.DrawHist("../plot/data18_0621.pdf");
  t_349014.DrawHist("../plot/t_349014.pdf");
  cout << "[INFO]: DrawHist SUCCESS" << endl;

  t_349014.CalcEff();
  cout << "[INFO]: CalcEff SUCCESS" << endl;

  //t_349014.DrawEffHist("../plot/data18_0621_eff.pdf");
  t_349014.DrawEffHist("../plot/test_eff_all.pdf");
  cout << "[INFO]: DrawEffHist SUCCESS" << endl;

  t_349014.End();
  cout << "[INFO]: End SUCCESS" << endl;

  //for ( int i=0; i < 60;i++){
  //  t_349014.Display(i);
  //}
  //t_349014.Display(5);

  // tree2
  //TChain *tree2 = new TChain("t_tap", "t_tap");
  ////tree2 -> Add("/gpfs/fs2001/yfukuhar/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_jpzYFV3GRL_EXT0.root");
  ////tree2 -> Add("/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/_gpfs_home_yfukuhar_public_dataset_aod_mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.recon.AOD.e5441_s3126.21.0.32.noRpcHitWide.pool.root/mc16c_Jpsimu6_NoTag.root");
  ////tree2 -> Add("/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/_gpfs_home_yfukuhar_public_dataset_aod_mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.recon.AOD.e5441_s3126.21.0.32.noRpcHit.pool.root/mc16c_Jpsimu6_NoTag.root");
  //tree2 -> Add("/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/_gpfs_home_yfukuhar_public_dataset_aod_mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.recon.AOD.e5441_s3126.21.0.32.noRpcHitWide.pool.root/mc16c_Jpsimu6_NoTag.root");

  //RPC t_349533(tree2); 

  //t_349533.Loop(-1, 10000);
  //cout << "[INFO]: Loop SUCCESS" << endl;

  //t_349533.DrawHist("../plot/t_mc16c_noRpcHitWide.pdf");
  //cout << "[INFO]: DrawHist SUCCESS" << endl;

  //t_349014.End();
  //cout << "[INFO]: End SUCCESS" << endl;


  //delete tree1;
  //delete tree2;

  return 0;
}


void RPC::Loop( int Nevents, int DisplayNumber )
{
//   In a ROOT session, you can do:
//      root> .L RPC.C
//      root> RPC t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   int nLoop;
   if (Nevents == -1) {
     nLoop = fChain -> GetEntries();
   } else {
     nLoop = Nevents;
   }

   //Long64_t nentries = fChain->GetEntriesFast();
   double entries = fChain->GetEntries();
   cout << "[INFO]: Nentries:" << entries << endl;

   //const int N50 = 14;


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nLoop;jentry++) {
      int  ientry = LoadTree(jentry);
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (ientry < 0) break;
      if( ientry%DisplayNumber == 0){
        cout << "now event number -->> " << ientry << "/" << nLoop << " : " << ((double) ientry)/nLoop * 100. << "%" << endl;
      }


      //==================================================================
      //Analysis code for a entry
      //==================================================================
      // Check GRL
      //if (GRLlist(LumiBlock)){
      //  continue;
      //}
      FillProbeHist();

      tag_proc = NTagProc;
      switch (tag_proc) {
        case 1: //Jpsi until L2
          break;
        case 2: //Jpsi from L2:
          break;
        case 3: //Z

          // Check TAP
          if (!(probe_mesEFTAG_pass -> at(N50) > -1 && probe_mesL1_pass -> at(N50) > 0 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
            continue;
          }
          // Set superpoint and segment for each station
          double probe_segmentR_BIS = 0;
          double probe_segmentZ_BIS = 0;
          double probe_segmentR_BIL = 0;
          double probe_segmentZ_BIL = 0;
          double probe_segmentR_BMS = 0;
          double probe_segmentZ_BMS = 0;
          double probe_segmentR_BML = 0;
          double probe_segmentZ_BML = 0;
          double probe_segmentR_BOS = 0;
          double probe_segmentZ_BOS = 0;
          double probe_segmentR_BOL = 0;
          double probe_segmentZ_BOL = 0;
          for ( int i = 0; i < probe_segment_n; i++){
            if (probe_segment_chamberIndex[i] == 0) {
              probe_segmentR_BIS = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i])/1000.;
              probe_segmentZ_BIS = probe_segment_z[i]/1000.;
            } else if(probe_segment_chamberIndex[i] == 1) {
              probe_segmentR_BIL = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i])/1000.;
              probe_segmentZ_BIL = probe_segment_z[i]/1000.;
            } else if(probe_segment_chamberIndex[i] == 2) {
              probe_segmentR_BMS = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i])/1000.;
              probe_segmentZ_BMS = probe_segment_z[i]/1000.;
            } else if(probe_segment_chamberIndex[i] == 3) {
              probe_segmentR_BML = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i])/1000.;
              probe_segmentZ_BML = probe_segment_z[i]/1000.;
            } else if(probe_segment_chamberIndex[i] == 4) {
              probe_segmentR_BOS = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i])/1000.;
              probe_segmentZ_BOS = probe_segment_z[i]/1000.;
            } else if(probe_segment_chamberIndex[i] == 5) {
              probe_segmentR_BOL = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i])/1000.;
              probe_segmentZ_BOL = probe_segment_z[i]/1000.;
            }
          }

          // Check GRL
          if (!GRLlist(LumiBlock)){
            continue;
          }

          //// Check isRpcFailure
          //if (probe_mesSA_isRpcFailure -> at(N50) == 1){
          //  continue;
          //}


          int nMdtBI = 0;
          int nMdtBM = 0;
          int nMdtBO = 0;
          //cout << "NmdtHits1: " << probe_mesSA_mdtHitIsOutlier->size() << endl;
          //cout << "NmdtHits2: " << (probe_mesSA_mdtHitIsOutlier -> at(14)).size() << endl;
          for ( uint32_t i = 0; i < (probe_mesSA_mdtHitIsOutlier -> at(N50)).size();i++){
            if (probe_mesSA_mdtHitIsOutlier -> at(N50)[i] == 1) {
              continue;
            }
            //cout << "chamber: " << probe_mesSA_mdtHitChamber -> at(14)[i] << endl;
            switch (probe_mesSA_mdtHitChamber -> at(N50)[i]) {
              case 0:
                nMdtBI += 1;
                //cout << "BI" << endl;
                h_ResidualMdt_eta_BI    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
                if (abs(probe_eta) < 1.05) h_ResidualMdt_pt_barrel_BI -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
                break;
              case 1:
                nMdtBM += 1;
                //cout << "BM" << endl;
                h_ResidualMdt_eta_BM    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
                if (abs(probe_eta) < 1.05) h_ResidualMdt_pt_barrel_BM -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
                break;
              case 2:
                nMdtBO += 1;
                //cout << "BO" << endl;
                h_ResidualMdt_eta_BO    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
                if (abs(probe_eta) < 1.05) h_ResidualMdt_pt_barrel_BO -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
                break;
            }
            //cout << "mdtHits: "<< i << ": " << (probe_mesSA_mdtHitChamber -> at(14))[i] << endl;
          }

          h_NumberOfMdt_eta    -> Fill(probe_eta, nMdtBI + nMdtBM + nMdtBO);
          h_NumberOfMdt_eta_BI -> Fill(probe_eta, nMdtBI);
          h_NumberOfMdt_eta_BM -> Fill(probe_eta, nMdtBM);
          h_NumberOfMdt_eta_BO -> Fill(probe_eta, nMdtBO);

          h_NumberOfSP_eta->Fill(probe_eta, NumberOfSP());
          h_NumberOfSP_qeta->Fill(probe_charge*probe_eta, NumberOfSP());
          if (abs(probe_eta) < 1.05){
            //cout << NumberOfSP() << endl;
            h_NumberOfMdt_LumiBlock->Fill(LumiBlock, (probe_mesSA_mdtHitIsOutlier -> at(14)).size());
            h_NumberOfMdt_LumiBlock_BI->Fill(LumiBlock, nMdtBI);
            h_NumberOfMdt_LumiBlock_BM->Fill(LumiBlock, nMdtBM);
            h_NumberOfMdt_LumiBlock_BO->Fill(LumiBlock, nMdtBO);
            h_NumberOfSP_LumiBlock->Fill(LumiBlock, NumberOfSP());
            h_NumberOfSP_pt_barrel->Fill(probe_pt/1000., NumberOfSP());
            h_NumberOfSP_pass_barrel->Fill(probe_mesSA_pass->at(N50), NumberOfSP());
          }

          //BIS
          if (probe_mesSA_superPointR_BI->at(N50) != 0 &&
              probe_mesSA_superPointR_BI -> at(N50) / 1000. > -90 &&
              probe_mesSA_sAddress -> at(N50) == 2 ) {
            h_superPointRZ_BIS -> Fill( probe_mesSA_superPointZ_BI->at(N50)/1000, probe_mesSA_superPointR_BI->at(N50)/1000); 
            h_segmentRZ_BIS -> Fill( probe_segmentZ_BIS, probe_segmentR_BIS); 
            h_residualRZ_BIS -> Fill( Res(probe_mesSA_superPointZ_BI->at(N50)/1000,probe_segmentZ_BIS), Res(probe_mesSA_superPointR_BI->at(N50)/1000,probe_segmentR_BIS)); 
          }
          //BIL
          if (probe_mesSA_superPointR_BI->at(N50) != 0 &&
              probe_mesSA_superPointR_BI -> at(N50) / 1000. > -90 &&
              probe_mesSA_sAddress -> at(N50) == 0 ) {
            h_superPointRZ_BIL -> Fill( probe_mesSA_superPointZ_BI->at(N50)/1000, probe_mesSA_superPointR_BI->at(N50)/1000); 
            h_segmentRZ_BIL -> Fill( probe_segmentZ_BIL, probe_segmentR_BIL); 
            h_residualRZ_BIL -> Fill( Res(probe_mesSA_superPointZ_BI->at(N50)/1000,probe_segmentZ_BIL), Res(probe_mesSA_superPointR_BI->at(N50)/1000,probe_segmentR_BIL)); 
          }
          //BMS
          if (probe_mesSA_superPointR_BM->at(N50) != 0 &&
              probe_mesSA_superPointR_BM -> at(N50) / 1000. > -90 &&
              probe_mesSA_sAddress -> at(N50) == 2 ) {
            h_superPointRZ_BMS -> Fill( probe_mesSA_superPointZ_BM->at(N50)/1000, probe_mesSA_superPointR_BM->at(N50)/1000); 
            h_segmentRZ_BMS -> Fill( probe_segmentZ_BMS, probe_segmentR_BMS); 
            h_residualRZ_BMS -> Fill( Res(probe_mesSA_superPointZ_BM->at(N50)/1000,probe_segmentZ_BMS), Res(probe_mesSA_superPointR_BM->at(N50)/1000,probe_segmentR_BMS)); 
          }
          //BML
          if (probe_mesSA_superPointR_BM->at(N50) != 0 &&
              probe_mesSA_superPointR_BM -> at(N50) / 1000. > -90 &&
              probe_mesSA_sAddress -> at(N50) == 0 ) {
            h_superPointRZ_BML -> Fill( probe_mesSA_superPointZ_BM->at(N50)/1000, probe_mesSA_superPointR_BM->at(N50)/1000); 
            h_segmentRZ_BML -> Fill( probe_segmentZ_BML, probe_segmentR_BML); 
            h_residualRZ_BML -> Fill( Res(probe_mesSA_superPointZ_BM->at(N50)/1000,probe_segmentZ_BML), Res(probe_mesSA_superPointR_BM->at(N50)/1000,probe_segmentR_BML)); 
          }
          //BOS
          if (probe_mesSA_superPointR_BO->at(N50) != 0 &&
              probe_mesSA_superPointR_BO -> at(N50) / 1000. > -90 &&
              probe_mesSA_sAddress -> at(N50) == 2 ) {
            h_superPointRZ_BOS -> Fill( probe_mesSA_superPointZ_BO->at(N50)/1000, probe_mesSA_superPointR_BO->at(N50)/1000); 
            h_segmentRZ_BOS -> Fill( probe_segmentZ_BOS, probe_segmentR_BOS); 
            h_residualRZ_BOS -> Fill( Res(probe_mesSA_superPointZ_BO->at(N50)/1000,probe_segmentZ_BOS), Res(probe_mesSA_superPointR_BO->at(N50)/1000,probe_segmentR_BOS)); 
          }
          //BOL
          if (probe_mesSA_superPointR_BO->at(N50) != 0 &&
              probe_mesSA_superPointR_BO -> at(N50) / 1000. > -90 &&
              probe_mesSA_sAddress -> at(N50) == 0 ) {
            h_superPointRZ_BOL -> Fill( probe_mesSA_superPointZ_BO->at(N50)/1000, probe_mesSA_superPointR_BO->at(N50)/1000); 
            h_segmentRZ_BOL -> Fill( probe_segmentZ_BOL, probe_segmentR_BOL); 
            h_residualRZ_BOL -> Fill( Res(probe_mesSA_superPointZ_BO->at(N50)/1000,probe_segmentZ_BOL), Res(probe_mesSA_superPointR_BO->at(N50)/1000,probe_segmentR_BOL)); 
          }
          break;
      }

   } // end of each entry

} // end of RPC::Loop()



void RPC::DrawHist(TString pdf){
  //==================================================================
  //Set Canvas
  //==================================================================
  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1020, 700);
  c1->SetGrid();
  c1->SetRightMargin(0.20);
  c1->SetLeftMargin(0.23);
  c1->SetBottomMargin(0.20);

  c1 -> Print( pdf + "[", "pdf" );

  h_ResidualMdt_pt_barrel_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_pt_barrel_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_pt_barrel_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_eta_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_eta_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_eta_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );


  h_NumberOfMdt_eta->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  // Begin of a Fraction plot
  h_NumberOfMdt_LumiBlock->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  FractionOfnMDTs(h_NumberOfMdt_LumiBlock, c1, pdf);

  FractionOfnMDTs(h_NumberOfMdt_eta_BI, c1, pdf);
  FractionOfnMDTs(h_NumberOfMdt_eta_BM, c1, pdf);
  FractionOfnMDTs(h_NumberOfMdt_eta_BO, c1, pdf);

  h_NumberOfMdt_eta_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_NumberOfMdt_eta_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_NumberOfMdt_eta_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_NumberOfMdt_LumiBlock->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_NumberOfMdt_LumiBlock_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_NumberOfMdt_LumiBlock_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_NumberOfMdt_LumiBlock_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  // Begin of a Fraction plot
  h_NumberOfSP_LumiBlock->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  TProfile *pf_LumiBlock = h_NumberOfSP_LumiBlock->ProfileX();
  pf_LumiBlock->Draw();
  pf_LumiBlock->GetYaxis()->SetTitle("Average number of SP");
  c1 -> Print(pdf, "pdf" );

  TH1D* h_LumiBlock_all = h_NumberOfSP_LumiBlock->ProjectionX("_all");
  TH1D* h_LumiBlock_5 = h_NumberOfSP_LumiBlock->ProjectionX("5",h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(4.9),h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(5.1));
  TH1D* h_LumiBlock_4 = h_NumberOfSP_LumiBlock->ProjectionX("4",h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(3.9),h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(4.1));
  TH1D* h_LumiBlock_3 = h_NumberOfSP_LumiBlock->ProjectionX("3",h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(2.9),h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(3.1));
  TH1D* h_LumiBlock_2 = h_NumberOfSP_LumiBlock->ProjectionX("2",h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(1.9),h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(2.1));
  TH1D* h_LumiBlock_1 = h_NumberOfSP_LumiBlock->ProjectionX("1",h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(0.9),h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(1.1));
  TH1D* h_LumiBlock_0 = h_NumberOfSP_LumiBlock->ProjectionX("0",h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(-0.9),h_NumberOfSP_LumiBlock->GetYaxis()->FindBin(0.1));

  h_LumiBlock_5->Divide(h_LumiBlock_all);
  h_LumiBlock_4->Divide(h_LumiBlock_all);
  h_LumiBlock_3->Divide(h_LumiBlock_all);
  h_LumiBlock_2->Divide(h_LumiBlock_all);
  h_LumiBlock_1->Divide(h_LumiBlock_all);
  h_LumiBlock_0->Divide(h_LumiBlock_all);


  THStack *hs_LumiBlock = new THStack("hs_LumiBlock",";LumiBlock;Fraction of number of SPs");
  h_LumiBlock_5->SetFillColor(kCyan);//あらかじめFillColorをSetしておく
  h_LumiBlock_4->SetFillColor(kMagenta);//あらかじめFillColorをSetしておく
  h_LumiBlock_3->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_LumiBlock_2->SetFillColor(kBlue);
  h_LumiBlock_1->SetFillColor(kGreen);
  h_LumiBlock_0->SetFillColor(kYellow);
  hs_LumiBlock->Add(h_LumiBlock_5);
  hs_LumiBlock->Add(h_LumiBlock_4);
  hs_LumiBlock->Add(h_LumiBlock_3);
  hs_LumiBlock->Add(h_LumiBlock_2);
  hs_LumiBlock->Add(h_LumiBlock_1);
  hs_LumiBlock->Add(h_LumiBlock_0);

  hs_LumiBlock->Draw();

  TLegend *leg_LumiBlock = new TLegend(0.82,0.62,0.9,0.92);
  leg_LumiBlock->AddEntry(h_LumiBlock_0," n=0","f");
  leg_LumiBlock->AddEntry(h_LumiBlock_1," n=1","f");
  leg_LumiBlock->AddEntry(h_LumiBlock_2," n=2","f");
  leg_LumiBlock->AddEntry(h_LumiBlock_3," n=3","f");
  leg_LumiBlock->AddEntry(h_LumiBlock_4," n=4","f");
  leg_LumiBlock->AddEntry(h_LumiBlock_5," n=5","f");
  leg_LumiBlock->Draw();
  c1 -> Print(pdf, "pdf" );
  delete leg_LumiBlock;
  // End of a Fraction plot


  h_NumberOfSP_eta->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  TProfile *pf_eta = h_NumberOfSP_eta->ProfileX();
  pf_eta->Draw();
  pf_eta->GetYaxis()->SetTitle("Average number of SP");
  c1 -> Print(pdf, "pdf" );

  TH1D* h_eta_all = h_NumberOfSP_eta->ProjectionX("_all");
  TH1D* h_eta_5 = h_NumberOfSP_eta->ProjectionX("5",h_NumberOfSP_eta->GetYaxis()->FindBin(4.9),h_NumberOfSP_eta->GetYaxis()->FindBin(5.1));
  TH1D* h_eta_4 = h_NumberOfSP_eta->ProjectionX("4",h_NumberOfSP_eta->GetYaxis()->FindBin(3.9),h_NumberOfSP_eta->GetYaxis()->FindBin(4.1));
  TH1D* h_eta_3 = h_NumberOfSP_eta->ProjectionX("3",h_NumberOfSP_eta->GetYaxis()->FindBin(2.9),h_NumberOfSP_eta->GetYaxis()->FindBin(3.1));
  TH1D* h_eta_2 = h_NumberOfSP_eta->ProjectionX("2",h_NumberOfSP_eta->GetYaxis()->FindBin(1.9),h_NumberOfSP_eta->GetYaxis()->FindBin(2.1));
  TH1D* h_eta_1 = h_NumberOfSP_eta->ProjectionX("1",h_NumberOfSP_eta->GetYaxis()->FindBin(0.9),h_NumberOfSP_eta->GetYaxis()->FindBin(1.1));
  TH1D* h_eta_0 = h_NumberOfSP_eta->ProjectionX("0",h_NumberOfSP_eta->GetYaxis()->FindBin(-0.9),h_NumberOfSP_eta->GetYaxis()->FindBin(0.1));

  h_eta_5->Divide(h_eta_all);
  h_eta_4->Divide(h_eta_all);
  h_eta_3->Divide(h_eta_all);
  h_eta_2->Divide(h_eta_all);
  h_eta_1->Divide(h_eta_all);
  h_eta_0->Divide(h_eta_all);


  THStack *hs_eta = new THStack("hs_eta",";eta;Fraction of number of SPs");
  h_eta_5->SetFillColor(kCyan);//あらかじめFillColorをSetしておく
  h_eta_4->SetFillColor(kMagenta);//あらかじめFillColorをSetしておく
  h_eta_3->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_eta_2->SetFillColor(kBlue);
  h_eta_1->SetFillColor(kGreen);
  h_eta_0->SetFillColor(kYellow);
  hs_eta->Add(h_eta_5);
  hs_eta->Add(h_eta_4);
  hs_eta->Add(h_eta_3);
  hs_eta->Add(h_eta_2);
  hs_eta->Add(h_eta_1);
  hs_eta->Add(h_eta_0);

  hs_eta->Draw();

  TLegend *leg_eta = new TLegend(0.82,0.62,0.9,0.92);
  leg_eta->AddEntry(h_eta_0," n=0","f");
  leg_eta->AddEntry(h_eta_1," n=1","f");
  leg_eta->AddEntry(h_eta_2," n=2","f");
  leg_eta->AddEntry(h_eta_3," n=3","f");
  leg_eta->AddEntry(h_eta_4," n=4","f");
  leg_eta->AddEntry(h_eta_5," n=5","f");
  leg_eta->Draw();
  c1 -> Print(pdf, "pdf" );
  delete leg_eta;

  h_NumberOfSP_qeta->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  TProfile *pf_qeta = h_NumberOfSP_qeta->ProfileX();
  pf_qeta->Draw();
  pf_qeta->GetYaxis()->SetTitle("Average number of SP");
  c1 -> Print(pdf, "pdf" );

  TH1D* h_qeta_all = h_NumberOfSP_qeta->ProjectionX("_all");
  TH1D* h_qeta_5 = h_NumberOfSP_qeta->ProjectionX("5",h_NumberOfSP_qeta->GetYaxis()->FindBin(4.9),h_NumberOfSP_qeta->GetYaxis()->FindBin(5.1));
  TH1D* h_qeta_4 = h_NumberOfSP_qeta->ProjectionX("4",h_NumberOfSP_qeta->GetYaxis()->FindBin(3.9),h_NumberOfSP_qeta->GetYaxis()->FindBin(4.1));
  TH1D* h_qeta_3 = h_NumberOfSP_qeta->ProjectionX("3",h_NumberOfSP_qeta->GetYaxis()->FindBin(2.9),h_NumberOfSP_qeta->GetYaxis()->FindBin(3.1));
  TH1D* h_qeta_2 = h_NumberOfSP_qeta->ProjectionX("2",h_NumberOfSP_qeta->GetYaxis()->FindBin(1.9),h_NumberOfSP_qeta->GetYaxis()->FindBin(2.1));
  TH1D* h_qeta_1 = h_NumberOfSP_qeta->ProjectionX("1",h_NumberOfSP_qeta->GetYaxis()->FindBin(0.9),h_NumberOfSP_qeta->GetYaxis()->FindBin(1.1));
  TH1D* h_qeta_0 = h_NumberOfSP_qeta->ProjectionX("0",h_NumberOfSP_qeta->GetYaxis()->FindBin(-0.9),h_NumberOfSP_qeta->GetYaxis()->FindBin(0.1));

  h_qeta_5->Divide(h_qeta_all);
  h_qeta_4->Divide(h_qeta_all);
  h_qeta_3->Divide(h_qeta_all);
  h_qeta_2->Divide(h_qeta_all);
  h_qeta_1->Divide(h_qeta_all);
  h_qeta_0->Divide(h_qeta_all);


  THStack *hs_qeta = new THStack("hs_qeta",";qeta;Fraction of number of SPs");
  h_qeta_5->SetFillColor(kCyan);//あらかじめFillColorをSetしておく
  h_qeta_4->SetFillColor(kMagenta);//あらかじめFillColorをSetしておく
  h_qeta_3->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_qeta_2->SetFillColor(kBlue);
  h_qeta_1->SetFillColor(kGreen);
  h_qeta_0->SetFillColor(kYellow);
  hs_qeta->Add(h_qeta_5);
  hs_qeta->Add(h_qeta_4);
  hs_qeta->Add(h_qeta_3);
  hs_qeta->Add(h_qeta_2);
  hs_qeta->Add(h_qeta_1);
  hs_qeta->Add(h_qeta_0);

  hs_qeta->Draw();

  TLegend *leg_qeta = new TLegend(0.82,0.62,0.9,0.92);
  leg_qeta->AddEntry(h_qeta_0," n=0","f");
  leg_qeta->AddEntry(h_qeta_1," n=1","f");
  leg_qeta->AddEntry(h_qeta_2," n=2","f");
  leg_qeta->AddEntry(h_qeta_3," n=3","f");
  leg_qeta->AddEntry(h_qeta_4," n=4","f");
  leg_qeta->AddEntry(h_qeta_5," n=5","f");
  leg_qeta->Draw();
  c1 -> Print(pdf, "pdf" );
  delete leg_qeta;

  h_NumberOfSP_pt_barrel->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  TProfile *pf_pt = h_NumberOfSP_pt_barrel->ProfileX();
  pf_pt->Draw();
  pf_pt->GetYaxis()->SetTitle("Average number of SP");
  c1 -> Print(pdf, "pdf" );

  TH1D* h_pt_all = h_NumberOfSP_pt_barrel->ProjectionX("_all");
  TH1D* h_pt_5 = h_NumberOfSP_pt_barrel->ProjectionX("5",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(4.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(5.1));
  TH1D* h_pt_4 = h_NumberOfSP_pt_barrel->ProjectionX("4",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(3.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(4.1));
  TH1D* h_pt_3 = h_NumberOfSP_pt_barrel->ProjectionX("3",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(2.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(3.1));
  TH1D* h_pt_2 = h_NumberOfSP_pt_barrel->ProjectionX("2",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(1.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(2.1));
  TH1D* h_pt_1 = h_NumberOfSP_pt_barrel->ProjectionX("1",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(0.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(1.1));
  TH1D* h_pt_0 = h_NumberOfSP_pt_barrel->ProjectionX("0",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(-0.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(0.1));

  h_pt_5->Divide(h_pt_all);
  h_pt_4->Divide(h_pt_all);
  h_pt_3->Divide(h_pt_all);
  h_pt_2->Divide(h_pt_all);
  h_pt_1->Divide(h_pt_all);
  h_pt_0->Divide(h_pt_all);


  THStack *hs_pt = new THStack("hs_pt",";pT [GeV];Fraction of number of SPs");
  h_pt_5->SetFillColor(kCyan);//あらかじめFillColorをSetしておく
  h_pt_4->SetFillColor(kMagenta);//あらかじめFillColorをSetしておく
  h_pt_3->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_pt_2->SetFillColor(kBlue);
  h_pt_1->SetFillColor(kGreen);
  h_pt_0->SetFillColor(kYellow);
  hs_pt->Add(h_pt_5);
  hs_pt->Add(h_pt_4);
  hs_pt->Add(h_pt_3);
  hs_pt->Add(h_pt_2);
  hs_pt->Add(h_pt_1);
  hs_pt->Add(h_pt_0);

  hs_pt->Draw();

  TLegend *leg_pt = new TLegend(0.82,0.62,0.9,0.92);
  leg_pt->AddEntry(h_pt_0," n=0","f");
  leg_pt->AddEntry(h_pt_1," n=1","f");
  leg_pt->AddEntry(h_pt_2," n=2","f");
  leg_pt->AddEntry(h_pt_3," n=3","f");
  leg_pt->AddEntry(h_pt_4," n=4","f");
  leg_pt->AddEntry(h_pt_5," n=5","f");
  leg_pt->Draw();
  c1 -> Print(pdf, "pdf" );
  delete leg_pt;

  h_NumberOfSP_pass_barrel->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  TProfile *pf_pass = h_NumberOfSP_pass_barrel->ProfileX();
  pf_pass->GetYaxis()->SetTitle("Average number of SP");
  pf_pass->Draw();
  c1 -> Print(pdf, "pdf" );


  TH1D* h_pass_all = h_NumberOfSP_pass_barrel->ProjectionX("_all");
  TH1D* h_pass_3 = h_NumberOfSP_pass_barrel->ProjectionX("pass_3",h_NumberOfSP_pass_barrel->GetYaxis()->FindBin(2.9),h_NumberOfSP_pass_barrel->GetYaxis()->FindBin(3.1));
  TH1D* h_pass_2 = h_NumberOfSP_pass_barrel->ProjectionX("pass_2",h_NumberOfSP_pass_barrel->GetYaxis()->FindBin(1.9),h_NumberOfSP_pass_barrel->GetYaxis()->FindBin(2.1));
  TH1D* h_pass_1 = h_NumberOfSP_pass_barrel->ProjectionX("pass_1",h_NumberOfSP_pass_barrel->GetYaxis()->FindBin(0.9),h_NumberOfSP_pass_barrel->GetYaxis()->FindBin(1.1));
  TH1D* h_pass_0 = h_NumberOfSP_pass_barrel->ProjectionX("pass_0",h_NumberOfSP_pass_barrel->GetYaxis()->FindBin(-0.9),h_NumberOfSP_pass_barrel->GetYaxis()->FindBin(0.1));

  h_pass_3->Divide(h_pass_all);
  h_pass_2->Divide(h_pass_all);
  h_pass_1->Divide(h_pass_all);
  h_pass_0->Divide(h_pass_all);

  THStack *hs_pass = new THStack("hs_pass",";pass;Fraction of number of SPs");
  h_pass_3->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_pass_2->SetFillColor(kBlue);
  h_pass_1->SetFillColor(kGreen);
  h_pass_0->SetFillColor(kYellow);
  hs_pass->Add(h_pass_3);
  hs_pass->Add(h_pass_2);
  hs_pass->Add(h_pass_1);
  hs_pass->Add(h_pass_0);

  hs_pass->Draw();
  c1 -> Print(pdf, "pdf" );

  h_superPointRZ_BIL -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_segmentRZ_BIL -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BIL  -> SetTitle("R;Z;BIL Residual");
  h_residualRZ_BIL -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BIL->ProjectionX("name",h_residualRZ_BIL->GetYaxis()->FindBin(-0.1),h_residualRZ_BIL->GetYaxis()->FindBin(0.1)) -> Draw();
  //h_residualRZ_BIL->ProjectionX("name",h_residualRZ_BIL->GetYaxis()->FindBin(-0.1),h_residualRZ_BIL->GetYaxis()->FindBin(0.1)) -> DrawNormalized("",1./(h_residualRZ_BIL->ProjectionX("name",h_residualRZ_BIL->GetYaxis()->FindBin(-0.1),h_residualRZ_BIL->GetYaxis()->FindBin(0.1))->GetBinWidth(1)));
  c1 -> Print(pdf, "pdf" );

  h_superPointRZ_BIS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_segmentRZ_BIS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BIS  -> SetTitle("R;Z;BIS Residual");
  h_residualRZ_BIS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BIS->ProjectionX("name",h_residualRZ_BIS->GetYaxis()->FindBin(-0.1),h_residualRZ_BIS->GetYaxis()->FindBin(0.1)) -> Draw();
  //h_residualRZ_BIS->ProjectionX("name",h_residualRZ_BIS->GetYaxis()->FindBin(-0.1),h_residualRZ_BIS->GetYaxis()->FindBin(0.1)) -> DrawNormalized("",1./(h_residualRZ_BIS->ProjectionX("name",h_residualRZ_BIS->GetYaxis()->FindBin(-0.1),h_residualRZ_BIS->GetYaxis()->FindBin(0.1))->GetBinWidth(1)));
  c1 -> Print(pdf, "pdf" );

  h_superPointRZ_BML -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_segmentRZ_BML -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BML  -> SetTitle("R;Z;BML Residual");
  h_residualRZ_BML -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BML->ProjectionX("name",h_residualRZ_BML->GetYaxis()->FindBin(-0.1),h_residualRZ_BML->GetYaxis()->FindBin(0.1)) -> Draw();
  //h_residualRZ_BML->ProjectionX("name",h_residualRZ_BML->GetYaxis()->FindBin(-0.1),h_residualRZ_BML->GetYaxis()->FindBin(0.1)) -> DrawNormalized("",1./(h_residualRZ_BML->ProjectionX("name",h_residualRZ_BML->GetYaxis()->FindBin(-0.1),h_residualRZ_BML->GetYaxis()->FindBin(0.1))->GetBinWidth(1)));
  c1 -> Print(pdf, "pdf" );

  h_superPointRZ_BMS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_segmentRZ_BMS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BMS  -> SetTitle("R;Z;BML Residual");
  h_residualRZ_BMS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BMS->ProjectionX("name",h_residualRZ_BMS->GetYaxis()->FindBin(-0.1),h_residualRZ_BMS->GetYaxis()->FindBin(0.1)) -> Draw();
  //h_residualRZ_BMS->ProjectionX("name",h_residualRZ_BMS->GetYaxis()->FindBin(-0.1),h_residualRZ_BMS->GetYaxis()->FindBin(0.1)) -> DrawNormalized("",1./(h_residualRZ_BMS->ProjectionX("name",h_residualRZ_BMS->GetYaxis()->FindBin(-0.1),h_residualRZ_BMS->GetYaxis()->FindBin(0.1))->GetBinWidth(1)));
  c1 -> Print(pdf, "pdf" );

  h_superPointRZ_BOL -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_segmentRZ_BOL -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BOL  -> SetTitle("R;Z;BOL Residual");
  h_residualRZ_BOL -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BOL->ProjectionX("name",h_residualRZ_BOL->GetYaxis()->FindBin(-0.1),h_residualRZ_BOL->GetYaxis()->FindBin(0.1)) -> Draw();
  //h_residualRZ_BOL->ProjectionX("name",h_residualRZ_BOL->GetYaxis()->FindBin(-0.1),h_residualRZ_BOL->GetYaxis()->FindBin(0.1)) -> DrawNormalized("",1./(h_residualRZ_BOL->ProjectionX("name",h_residualRZ_BOL->GetYaxis()->FindBin(-0.1),h_residualRZ_BOL->GetYaxis()->FindBin(0.1))->GetBinWidth(1)));
  c1 -> Print(pdf, "pdf" );

  h_superPointRZ_BOS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_segmentRZ_BOS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BOS  -> SetTitle("R;Z;BOS Residual");
  h_residualRZ_BOS -> Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_residualRZ_BOS->ProjectionX("name",h_residualRZ_BOS->GetYaxis()->FindBin(-0.1),h_residualRZ_BOS->GetYaxis()->FindBin(0.1)) -> Draw();
  //h_residualRZ_BOS->ProjectionX("name",h_residualRZ_BOS->GetYaxis()->FindBin(-0.1),h_residualRZ_BOS->GetYaxis()->FindBin(0.1)) -> DrawNormalized("",1./(h_residualRZ_BOS->ProjectionX("name",h_residualRZ_BOS->GetYaxis()->FindBin(-0.1),h_residualRZ_BOS->GetYaxis()->FindBin(0.1))->GetBinWidth(1)));
  c1 -> Print(pdf, "pdf" );

  c1 -> Print( pdf + "]", "pdf" );

  delete c1;
}


double Res(double param1, double param2){
  if (param2 == 0.) {
    return -99999.;
  } else if (param2 != 0.){
    return (param1-param2)/param2;
  }
  return -99999.;
}

//int NumberOfSP(vector<double> SPs){
//  int number = 0;
//  for (int i = 0; i! = SPs.size(); ++i){
//    if (SPs[i] > 0.) {
//      number += 1;
//    }
//  }
//
//}

int RPC::NumberOfSP(){
  int number = 0;
  if (abs(probe_mesSA_superPointZ_BI -> at(N50)) > 0. && probe_mesSA_superPointZ_BI -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BM -> at(N50)) > 0. && probe_mesSA_superPointZ_BM -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BO -> at(N50)) > 0. && probe_mesSA_superPointZ_BO -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EI -> at(N50)) > 0. && probe_mesSA_superPointZ_EI -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EM -> at(N50)) > 0. && probe_mesSA_superPointZ_EM -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EO -> at(N50)) > 0. && probe_mesSA_superPointZ_EO -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EE -> at(N50)) > 0. && probe_mesSA_superPointZ_EE -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_CSC -> at(N50)) > 0. && probe_mesSA_superPointZ_CSC -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BEE -> at(N50)) > 0. && probe_mesSA_superPointZ_BEE -> at(N50) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BME -> at(N50)) > 0. && probe_mesSA_superPointZ_BME -> at(N50) > -99990.) {
    number += 1;
  }
  return number;
}

void RPC::Display(Long64_t entry)
{
  Long64_t tmp_nb;
  Long64_t ientry = LoadTree(entry);
  tmp_nb = fChain->GetEntry(entry);
  
  //cout << ientry << endl;
  //cout << tmp_nb << endl;

  TCanvas *c2 = new TCanvas("c2", "c2", 10, 10, 1020, 700);
  c2->SetGrid();
  c2->SetRightMargin(0.20);
  c2->SetLeftMargin(0.23);
  c2->SetBottomMargin(0.20);

  cout << "size: " << (probe_mesSA_rpcHitEta -> at(N50)).size() << endl;
  cout << "eta: " << (probe_mesSA_eta -> at(N50)) << endl;
  //cout << probe_mesSA_pt -> at(N50) << endl;

  //fChain->Show(entry);

  const Int_t n = (probe_mesSA_rpcHitEta -> at(N50)).size();
  TGraph *gr = new TGraph(n); //各点が(0,0)で初期化される

  for (Int_t i=0;i<n;++i) {
    gr->SetPoint(i , (probe_mesSA_rpcHitZ -> at(N50))[i] / 1000. , (probe_mesSA_rpcHitR -> at(N50))[i] / 1000.); //SetPoint(点番号,x座標,y座標)
  }

  gr -> GetXaxis()->SetLimits(0,20);
  gr -> GetYaxis()->SetRangeUser(0,20);
  gr->Draw("AP");
  
  c2->Print("../plot/test0531.pdf");

  delete gr;
  delete c2;
}

bool GRLlist(int LumiBlock){
 bool lb_1 = (LumiBlock > 110 && LumiBlock < 124);
 bool lb_2= (LumiBlock > 126 && LumiBlock < 130);
 bool lb_3= (LumiBlock > 169 && LumiBlock < 239);

 return (lb_1 || lb_2 || lb_3);
}

void RPC::FractionOfnMDTs(TH2F* h_NumberOfMdt, TCanvas* c1, TString pdf){

  // Begin of a Fraction plot
  TProfile *pf_Mdt = h_NumberOfMdt->ProfileX();
  pf_Mdt->Draw();
  pf_Mdt->GetYaxis()->SetTitle("Average number of Mdt");
  c1 -> Print(pdf, "pdf" );

  TH1D* h_all_Mdt = h_NumberOfMdt->ProjectionX("_all");
  TH1D* h_6_Mdt = h_NumberOfMdt->ProjectionX("6",h_NumberOfMdt->GetYaxis()->FindBin(24.9),h_NumberOfMdt->GetYaxis()->FindBin(30.1));
  TH1D* h_5_Mdt = h_NumberOfMdt->ProjectionX("5",h_NumberOfMdt->GetYaxis()->FindBin(19.9),h_NumberOfMdt->GetYaxis()->FindBin(25.1));
  TH1D* h_4_Mdt = h_NumberOfMdt->ProjectionX("4",h_NumberOfMdt->GetYaxis()->FindBin(14.9),h_NumberOfMdt->GetYaxis()->FindBin(19.9));
  TH1D* h_3_Mdt = h_NumberOfMdt->ProjectionX("3",h_NumberOfMdt->GetYaxis()->FindBin(9.9),h_NumberOfMdt->GetYaxis()->FindBin(14.9));
  TH1D* h_2_Mdt = h_NumberOfMdt->ProjectionX("2",h_NumberOfMdt->GetYaxis()->FindBin(4.9),h_NumberOfMdt->GetYaxis()->FindBin(9.9));
  TH1D* h_1_Mdt = h_NumberOfMdt->ProjectionX("1",h_NumberOfMdt->GetYaxis()->FindBin(-0.9),h_NumberOfMdt->GetYaxis()->FindBin(4.9));

  h_6_Mdt->Divide(h_all_Mdt);
  h_5_Mdt->Divide(h_all_Mdt);
  h_4_Mdt->Divide(h_all_Mdt);
  h_3_Mdt->Divide(h_all_Mdt);
  h_2_Mdt->Divide(h_all_Mdt);
  h_1_Mdt->Divide(h_all_Mdt);

  //TString ytitle = h_NumberOfMdt -> GetYaxis()-> GetTitle();

  THStack *hs_Mdt = new THStack("hs_Mdt",";LumiBlock;Fraction of number of Mdts");
  h_6_Mdt->SetFillColor(kYellow);//あらかじめFillColorをSetしておく
  h_5_Mdt->SetFillColor(kCyan);//あらかじめFillColorをSetしておく
  h_4_Mdt->SetFillColor(kMagenta);//あらかじめFillColorをSetしておく
  h_3_Mdt->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_2_Mdt->SetFillColor(kBlue);
  h_1_Mdt->SetFillColor(kGreen);
  hs_Mdt->Add(h_1_Mdt);
  hs_Mdt->Add(h_2_Mdt);
  hs_Mdt->Add(h_3_Mdt);
  hs_Mdt->Add(h_4_Mdt);
  hs_Mdt->Add(h_5_Mdt);
  hs_Mdt->Add(h_6_Mdt);

  hs_Mdt->Draw();

  TLegend *leg_Mdt = new TLegend(0.82,0.62,0.95,0.92);
  leg_Mdt->AddEntry(h_1_Mdt," n=0~5","f");
  leg_Mdt->AddEntry(h_2_Mdt," n=5~10","f");
  leg_Mdt->AddEntry(h_3_Mdt," n=10~15","f");
  leg_Mdt->AddEntry(h_4_Mdt," n=15~20","f");
  leg_Mdt->AddEntry(h_5_Mdt," n=20~25","f");
  leg_Mdt->AddEntry(h_6_Mdt," n=25~30","f");
  leg_Mdt->Draw();
  c1 -> Print(pdf, "pdf" );
  delete leg_Mdt;
  delete hs_Mdt;
  // End of a Fraction plot
}


void RPC::FillProbeHist(){
  switch (tag_proc) {
    case 1: //Jpsi until L2
      // Check TAP
      if (!(probe_mesEFTAG_pass -> at(N4) > -1 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
        return;
      }
      // Fill L1 probe hists
      if (probe_mesL1_pass -> at(N4) > -1){
        h_probe_pt_mu4_L1 -> Fill(probe_pt/1000.);
        hh_probe_qetapt_mu4_L1 -> Fill(probe_charge*probe_eta, probe_pt/1000.);
        // Plateau cut
        if (probe_pt/1000. > 8.0) {
          h_probe_eta_mu4_L1 -> Fill(probe_eta);
          h_probe_phi_mu4_L1 -> Fill(probe_phi);
          hh_probe_etaphi_mu4_L1 -> Fill(probe_eta, probe_phi);
        }
      }
      // Fill SA probe hists
      if (probe_mesL1_pass -> at(N4) > -1 && probe_mesSA_pass -> at(N4) > -1){
        h_probe_pt_mu4_SA -> Fill(probe_pt/1000.);
        hh_probe_qetapt_mu4_SA -> Fill(probe_charge*probe_eta, probe_pt/1000.);
        // Plateau cut
        if (probe_pt/1000. > 8.0) {
          h_probe_eta_mu4_SA -> Fill(probe_eta);
          h_probe_phi_mu4_SA -> Fill(probe_phi);
          hh_probe_etaphi_mu4_SA -> Fill(probe_eta, probe_phi);
        }
      }

      break;
    case 2: //Jpsi from L2:
      break;
    case 3: //Z
      // Check TAP
      if (!(probe_mesEFTAG_pass -> at(N50) > -1 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
        return;
      }
      // Fill L1 probe hists
      if (probe_mesL1_pass -> at(N50) > -1){
        h_probe_pt_mu50_L1 -> Fill(probe_pt/1000.);
        hh_probe_qetapt_mu50_L1 -> Fill(probe_charge*probe_eta, probe_pt/1000.);
        // Plateau cut
        if (probe_pt/1000. > 8.0) {
          h_probe_eta_mu50_L1 -> Fill(probe_eta);
          h_probe_phi_mu50_L1 -> Fill(probe_phi);
          hh_probe_etaphi_mu50_L1 -> Fill(probe_eta, probe_phi);
        }
      }
      // Fill SA probe hists
      if (probe_mesL1_pass -> at(N50) > -1 && probe_mesSA_pass -> at(N50) > -1){
        h_probe_pt_mu50_SA -> Fill(probe_pt/1000.);
        hh_probe_qetapt_mu50_SA -> Fill(probe_charge*probe_eta, probe_pt/1000.);
        // Plateau cut
        if (probe_pt/1000. > 8.0) {
          h_probe_eta_mu50_SA -> Fill(probe_eta);
          h_probe_phi_mu50_SA -> Fill(probe_phi);
          hh_probe_etaphi_mu50_SA -> Fill(probe_eta, probe_phi);
        }
      }

      break;
  }
  return;
}

void RPC::CalcEff(){
  // mu4
  CalcHistToHist( h_probe_pt_mu4_SA,       h_probe_pt_mu4_L1,       h_eff_pt_mu4_L1SA);
  CalcHistToHist( h_probe_eta_mu4_SA,      h_probe_eta_mu4_L1,      h_eff_eta_mu4_L1SA);
  CalcHistToHist( h_probe_phi_mu4_SA,      h_probe_phi_mu4_L1,      h_eff_phi_mu4_L1SA);
  CalcHistToHist( hh_probe_qetapt_mu4_SA,  hh_probe_qetapt_mu4_L1,  hh_eff_qetapt_mu4_L1SA);
  CalcHistToHist( hh_probe_etaphi_mu4_SA,  hh_probe_etaphi_mu4_L1,  hh_eff_etaphi_mu4_L1SA);

  // mu50
  CalcHistToHist( h_probe_pt_mu50_SA,      h_probe_pt_mu50_L1,      h_eff_pt_mu50_L1SA);
  CalcHistToHist( h_probe_eta_mu50_SA,     h_probe_eta_mu50_L1,     h_eff_eta_mu50_L1SA);
  CalcHistToHist( h_probe_phi_mu50_SA,     h_probe_phi_mu50_L1,     h_eff_phi_mu50_L1SA);
  CalcHistToHist( hh_probe_qetapt_mu50_SA, hh_probe_qetapt_mu50_L1, hh_eff_qetapt_mu50_L1SA);
  CalcHistToHist( hh_probe_etaphi_mu50_SA, hh_probe_etaphi_mu50_L1, hh_eff_etaphi_mu50_L1SA);
}

void RPC::DrawEffHist(TString pdf){
  //==================================================================
  //Set Canvas
  //==================================================================
  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1020, 700);
  c1->SetGrid();
  c1->SetRightMargin(0.20);
  c1->SetLeftMargin(0.23);
  c1->SetBottomMargin(0.20);

  c1 -> Print( pdf + "[", "pdf" );

  h_eff_pt_mu4_L1SA->Draw();
  h_eff_pt_mu4_L1SA->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_eta_mu4_L1SA->Draw();
  h_eff_eta_mu4_L1SA->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_phi_mu4_L1SA->Draw();
  h_eff_phi_mu4_L1SA->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  hh_eff_etaphi_mu4_L1SA->Draw("colz");
  hh_eff_etaphi_mu4_L1SA->GetZaxis()->SetRangeUser(0,1.0);
  c1 -> Print( pdf, "pdf" );

  hh_eff_qetapt_mu4_L1SA->Draw("colz");
  hh_eff_qetapt_mu4_L1SA->GetZaxis()->SetRangeUser(0,1.0);
  c1 -> Print( pdf, "pdf" );


  h_eff_pt_mu50_L1SA->Draw();
  h_eff_phi_mu50_L1SA->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_eta_mu50_L1SA->Draw();
  h_eff_eta_mu50_L1SA->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_phi_mu50_L1SA->Draw();
  h_eff_phi_mu50_L1SA->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  hh_eff_etaphi_mu50_L1SA->Draw("colz");
  hh_eff_etaphi_mu50_L1SA->GetZaxis()->SetRangeUser(0,1.0);
  c1 -> Print( pdf, "pdf" );

  hh_eff_qetapt_mu50_L1SA->Draw("colz");
  hh_eff_qetapt_mu50_L1SA->GetZaxis()->SetRangeUser(0,1.0);
  c1 -> Print( pdf, "pdf" );

  c1 -> Print( pdf + "]", "pdf" );
}

void RPC::CalcHistToHist( TH1F* h1, TH1F* h2, TH1F* hout ) {
 
  double entry[3], error[3];
  const int nbinX = hout->GetXaxis()->GetNbins();

  if ( h1->GetXaxis()->GetNbins() != nbinX || h2->GetXaxis()->GetNbins() != nbinX ) {
    cerr << "Error::CalcHistToHist() : Number of bin mismatched" << endl;
    return;
  }

  for (int ibinX=0; ibinX<nbinX; ibinX++) {
    entry[0] = h1->GetBinContent(ibinX+1);
    error[0] = h1->GetBinError(ibinX+1);
    entry[1] = h2->GetBinContent(ibinX+1);
    error[1] = h2->GetBinError(ibinX+1);
    entry[2] = ( entry[1]!=0 )? entry[0]/entry[1]:0;
    error[2] = ( entry[1]!=0 )? sqrt( ( entry[0]*( 1. - 2*entry[2] ) + ( entry[1]*entry[2]*entry[2] ) )/(entry[1]*entry[1]) ):0;
    if ( entry[2] >= 1 ) {
      error[2] = 0;
      }
    hout->SetBinContent(ibinX+1, entry[2]);
    hout->SetBinError(ibinX+1, error[2]);
  }

  return;
}

void RPC::CalcHistToHist( TH2F* h1, TH2F* h2, TH2F* hout ) {

  double entry[3], error[3];
  const int nbinX = hout->GetXaxis()->GetNbins();
  const int nbinY = hout->GetYaxis()->GetNbins();

  if ( h1->GetXaxis()->GetNbins() != nbinX || h2->GetXaxis()->GetNbins() != nbinX ) {
    cerr << "Error::CalcHistToHist() : Number of Xbin mismatched" << endl;
    return;
  }
  if ( h1->GetYaxis()->GetNbins() != nbinY || h2->GetYaxis()->GetNbins() != nbinY ) {
    cerr << "Error::CalcHistToHist() : Number of Ybin mismatched" << endl;
    return;
  }

  for (int ibinX=0; ibinX<nbinX; ibinX++) {
    for (int ibinY=0; ibinY<nbinY; ibinY++) {
      entry[0] = h1->GetBinContent(ibinX+1, ibinY+1);
      error[0] = h1->GetBinError(ibinX+1, ibinY+1);
      entry[1] = h2->GetBinContent(ibinX+1, ibinY+1);
      error[1] = h2->GetBinError(ibinX+1, ibinY+1);
      error[0] = sqrt( entry[0] + 0.25 ) + 0.5;
      error[1] = sqrt( entry[0] + 0.25 ) + 0.5;
      entry[2] = ( entry[1]!=0 )? entry[0]/entry[1]:0;
      error[2] = ( entry[1]!=0 )? sqrt( entry[0]*entry[0]*error[1]*error[1] + entry[1]*entry[1]*error[0]*error[0] )/( entry[1]*entry[1] ):0;
      hout->SetBinContent(ibinX+1, ibinY+1, entry[2]);
      hout->SetBinError(ibinX+1, ibinY+1, error[2]);
    }
  }

  return;
}
