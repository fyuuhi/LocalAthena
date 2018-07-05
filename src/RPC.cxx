#define RPC_cxx
#include "RPC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TH1F.h"
#include "TF1.h"
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
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/data/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v2_EXT0/hadd/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v2_EXT.root");
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/data/data18_0621/data18_0621.root");
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16_youhei_Zmumu_default/mc16_youhei_Zmumu_default.root");
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16_youhei_Zmumu_noRpcHit/mc16_youhei_Zmumu_noRpcHit.root");
  tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16_youhei_Zmumu_noRpcHitWide/mc16_youhei_Zmumu_noRpcHitWide.root");
  //tree1 -> Add("/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/");
  //tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16c_Jpsimu6_default/mc16c_Jpsimu6_default.root");

  RPC t_349014(tree1); 

  t_349014.Loop(-1, 10000);
  cout << "[INFO]: Loop SUCCESS" << endl;

  //t_349014.DrawHist("../plot/data18_0621.pdf");
  t_349014.DrawHist("../plot/mc16_youhei_Zmumu_noRpcHitWide.pdf");
  //t_349014.DrawHist("../plot/mc16_youhei_Zmumu_noRpcHit.pdf");
  //t_349014.DrawHist("../plot/mc16_youhei_Zmumu_default.pdf");
  cout << "[INFO]: DrawHist SUCCESS" << endl;

  //t_349014.CalcEff();
  //cout << "[INFO]: CalcEff SUCCESS" << endl;

  ////t_349014.DrawEffHist("../plot/data18_0621_eff.pdf");
  //t_349014.DrawEffHist("../plot/test_eff_all.pdf");
  //cout << "[INFO]: DrawEffHist SUCCESS" << endl;
  t_349014.Display(8, "../plot/test0628.pdf");

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
          // Check TAP
          if (!(probe_mesEFTAG_pass -> at(N4) > -1 && probe_mesL1_pass -> at(N4) > 0 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
            continue;
          }
          break;
        case 2: //Jpsi from L2:
          break;
        case 3: //Z

          // Check TAP
          if (!(probe_mesEFTAG_pass -> at(N50) > -1 && probe_mesL1_pass -> at(N50) > 0 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
            continue;
          }
          // Set superpoint and segment for each station
          for ( int i = 0; i < probe_segment_n; i++){
            double R = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i]);
            double Z = probe_segment_z[i];
            if (probe_segment_chamberIndex[i] == 0 || probe_segment_chamberIndex[i] == 1) {
              h_ResidualSegment_eta_BI -> Fill(probe_eta, calc_residual(probe_mesSA_roadAw -> at(N50)[0], probe_mesSA_roadBw -> at(N50)[0], Z, R));
              cout << "BI: " << calc_residual(probe_mesSA_roadAw -> at(N50)[0], probe_mesSA_roadBw -> at(N50)[0], Z, R) << endl;
            } else if(probe_segment_chamberIndex[i] == 2 || probe_segment_chamberIndex[i] == 3) {
              h_ResidualSegment_eta_BM -> Fill(probe_eta, calc_residual(probe_mesSA_roadAw -> at(N50)[1], probe_mesSA_roadBw -> at(N50)[1], Z, R));
              cout <<"BM: " <<  calc_residual(probe_mesSA_roadAw -> at(N50)[0], probe_mesSA_roadBw -> at(N50)[0], Z, R) << endl;
            } else if(probe_segment_chamberIndex[i] == 4 || probe_segment_chamberIndex[i] == 5) {
              h_ResidualSegment_eta_BO -> Fill(probe_eta, calc_residual(probe_mesSA_roadAw -> at(N50)[2], probe_mesSA_roadBw -> at(N50)[2], Z, R));
              cout <<"BO: " <<  calc_residual(probe_mesSA_roadAw -> at(N50)[0], probe_mesSA_roadBw -> at(N50)[0], Z, R) << endl;
            }
          }


          //// Check isRpcFailure
          //if (probe_mesSA_isRpcFailure -> at(N50) == 1){
          //  continue;
          //}

          FillMdtHist();
          FillSPHist();

          //BIS
          //if (probe_mesSA_superPointR_BI->at(N50) != 0 &&
          //    probe_mesSA_superPointR_BI -> at(N50) / 1000. > -90 &&
          //    probe_mesSA_sAddress -> at(N50) == 2 ) {
          //  h_superPointRZ_BIS -> Fill( probe_mesSA_superPointZ_BI->at(N50)/1000, probe_mesSA_superPointR_BI->at(N50)/1000); 
          //  h_segmentRZ_BIS -> Fill( probe_segmentZ_BIS, probe_segmentR_BIS); 
          //  h_residualRZ_BIS -> Fill( Res(probe_mesSA_superPointZ_BI->at(N50)/1000,probe_segmentZ_BIS), Res(probe_mesSA_superPointR_BI->at(N50)/1000,probe_segmentR_BIS));
          //}
          ////BIL
          //if (probe_mesSA_superPointR_BI->at(N50) != 0 &&
          //    probe_mesSA_superPointR_BI -> at(N50) / 1000. > -90 &&
          //    probe_mesSA_sAddress -> at(N50) == 0 ) {
          //  h_superPointRZ_BIL -> Fill( probe_mesSA_superPointZ_BI->at(N50)/1000, probe_mesSA_superPointR_BI->at(N50)/1000); 
          //  h_segmentRZ_BIL -> Fill( probe_segmentZ_BIL, probe_segmentR_BIL); 
          //  h_residualRZ_BIL -> Fill( Res(probe_mesSA_superPointZ_BI->at(N50)/1000,probe_segmentZ_BIL), Res(probe_mesSA_superPointR_BI->at(N50)/1000,probe_segmentR_BIL));
          //}
          ////BMS
          //if (probe_mesSA_superPointR_BM->at(N50) != 0 &&
          //    probe_mesSA_superPointR_BM -> at(N50) / 1000. > -90 &&
          //    probe_mesSA_sAddress -> at(N50) == 2 ) {
          //  h_superPointRZ_BMS -> Fill( probe_mesSA_superPointZ_BM->at(N50)/1000, probe_mesSA_superPointR_BM->at(N50)/1000); 
          //  h_segmentRZ_BMS -> Fill( probe_segmentZ_BMS, probe_segmentR_BMS); 
          //  h_residualRZ_BMS -> Fill( Res(probe_mesSA_superPointZ_BM->at(N50)/1000,probe_segmentZ_BMS), Res(probe_mesSA_superPointR_BM->at(N50)/1000,probe_segmentR_BMS));
          //}
          ////BML
          //if (probe_mesSA_superPointR_BM->at(N50) != 0 &&
          //    probe_mesSA_superPointR_BM -> at(N50) / 1000. > -90 &&
          //    probe_mesSA_sAddress -> at(N50) == 0 ) {
          //  h_superPointRZ_BML -> Fill( probe_mesSA_superPointZ_BM->at(N50)/1000, probe_mesSA_superPointR_BM->at(N50)/1000); 
          //  h_segmentRZ_BML -> Fill( probe_segmentZ_BML, probe_segmentR_BML); 
          //  h_residualRZ_BML -> Fill( Res(probe_mesSA_superPointZ_BM->at(N50)/1000,probe_segmentZ_BML), Res(probe_mesSA_superPointR_BM->at(N50)/1000,probe_segmentR_BML));
          //}
          ////BOS
          //if (probe_mesSA_superPointR_BO->at(N50) != 0 &&
          //    probe_mesSA_superPointR_BO -> at(N50) / 1000. > -90 &&
          //    probe_mesSA_sAddress -> at(N50) == 2 ) {
          //  h_superPointRZ_BOS -> Fill( probe_mesSA_superPointZ_BO->at(N50)/1000, probe_mesSA_superPointR_BO->at(N50)/1000); 
          //  h_segmentRZ_BOS -> Fill( probe_segmentZ_BOS, probe_segmentR_BOS); 
          //  h_residualRZ_BOS -> Fill( Res(probe_mesSA_superPointZ_BO->at(N50)/1000,probe_segmentZ_BOS), Res(probe_mesSA_superPointR_BO->at(N50)/1000,probe_segmentR_BOS));
          //}
          ////BOL
          //if (probe_mesSA_superPointR_BO->at(N50) != 0 &&
          //    probe_mesSA_superPointR_BO -> at(N50) / 1000. > -90 &&
          //    probe_mesSA_sAddress -> at(N50) == 0 ) {
          //  h_superPointRZ_BOL -> Fill( probe_mesSA_superPointZ_BO->at(N50)/1000, probe_mesSA_superPointR_BO->at(N50)/1000); 
          //  h_segmentRZ_BOL -> Fill( probe_segmentZ_BOL, probe_segmentR_BOL); 
          //  h_residualRZ_BOL -> Fill( Res(probe_mesSA_superPointZ_BO->at(N50)/1000,probe_segmentZ_BOL), Res(probe_mesSA_superPointR_BO->at(N50)/1000,probe_segmentR_BOL));
          //}
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

  h_ResidualSegment_eta_BI -> Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualSegment_eta_BM -> Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualSegment_eta_BO -> Draw("colz");
  c1 -> Print(pdf, "pdf" );

  // Inlier
  h_ResidualMdt_Inlier_pt_barrel_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Inlier_pt_barrel_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Inlier_pt_barrel_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Inlier_eta_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Inlier_eta_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Inlier_eta_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );


  // Outlier
  h_ResidualMdt_Outlier_pt_barrel_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Outlier_pt_barrel_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Outlier_pt_barrel_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Outlier_eta_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Outlier_eta_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_ResidualMdt_Outlier_eta_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );


  h_NumberOfMdt_eta->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  DrawFractionOfnMDTs(h_NumberOfMdt_pt_barrel_BI, c1, pdf);
  DrawFractionOfnMDTs(h_NumberOfMdt_pt_barrel_BM, c1, pdf);
  DrawFractionOfnMDTs(h_NumberOfMdt_pt_barrel_BO, c1, pdf);
  DrawFractionOfnMDTs(h_NumberOfMdt_eta, c1, pdf);
  DrawFractionOfnMDTs(h_NumberOfMdt_eta_BI, c1, pdf);
  DrawFractionOfnMDTs(h_NumberOfMdt_eta_BM, c1, pdf);
  DrawFractionOfnMDTs(h_NumberOfMdt_eta_BO, c1, pdf);

  h_NumberOfMdt_eta_BI->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_NumberOfMdt_eta_BM->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  h_NumberOfMdt_eta_BO->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  //h_NumberOfMdt_LumiBlock->Draw("colz");
  //c1 -> Print(pdf, "pdf" );
  //h_NumberOfMdt_LumiBlock_BI->Draw("colz");
  //c1 -> Print(pdf, "pdf" );
  //h_NumberOfMdt_LumiBlock_BM->Draw("colz");
  //c1 -> Print(pdf, "pdf" );
  //h_NumberOfMdt_LumiBlock_BO->Draw("colz");
  //c1 -> Print(pdf, "pdf" );



  h_NumberOfSP_eta->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_NumberOfSP_qeta->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_NumberOfSP_pt_barrel->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  h_NumberOfSP_pass_barrel->Draw("colz");
  c1 -> Print(pdf, "pdf" );

  DrawFractionOfnSPs(h_NumberOfSP_LumiBlock, c1, pdf);


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

void RPC::DrawFractionOfnSPs(TH2F* h_NumberOfSP, TCanvas* c1, TString pdf){
  // Begin of a Fraction plot
  TProfile *pf_SP = h_NumberOfSP->ProfileX();
  pf_SP->Draw();
  pf_SP->GetYaxis()->SetTitle("Average number of SP");
  c1 -> Print(pdf, "pdf" );

  TH1D* h_SP_all = h_NumberOfSP->ProjectionX("_all");
  TH1D* h_SP_5 = h_NumberOfSP->ProjectionX("5",h_NumberOfSP->GetYaxis()->FindBin(4.9),h_NumberOfSP->GetYaxis()->FindBin(5.1));
  TH1D* h_SP_4 = h_NumberOfSP->ProjectionX("4",h_NumberOfSP->GetYaxis()->FindBin(3.9),h_NumberOfSP->GetYaxis()->FindBin(4.1));
  TH1D* h_SP_3 = h_NumberOfSP->ProjectionX("3",h_NumberOfSP->GetYaxis()->FindBin(2.9),h_NumberOfSP->GetYaxis()->FindBin(3.1));
  TH1D* h_SP_2 = h_NumberOfSP->ProjectionX("2",h_NumberOfSP->GetYaxis()->FindBin(1.9),h_NumberOfSP->GetYaxis()->FindBin(2.1));
  TH1D* h_SP_1 = h_NumberOfSP->ProjectionX("1",h_NumberOfSP->GetYaxis()->FindBin(0.9),h_NumberOfSP->GetYaxis()->FindBin(1.1));
  TH1D* h_SP_0 = h_NumberOfSP->ProjectionX("0",h_NumberOfSP->GetYaxis()->FindBin(-0.9),h_NumberOfSP->GetYaxis()->FindBin(0.1));

  h_SP_5->Divide(h_SP_all);
  h_SP_4->Divide(h_SP_all);
  h_SP_3->Divide(h_SP_all);
  h_SP_2->Divide(h_SP_all);
  h_SP_1->Divide(h_SP_all);
  h_SP_0->Divide(h_SP_all);


  THStack *hs_SP = new THStack("hs_SP",Form(";%s;Fraction of nSPs",h_NumberOfSP->GetXaxis()->GetTitle()));
  h_SP_5->SetFillColor(kCyan);//あらかじめFillColorをSetしておく
  h_SP_4->SetFillColor(kMagenta);//あらかじめFillColorをSetしておく
  h_SP_3->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_SP_2->SetFillColor(kBlue);
  h_SP_1->SetFillColor(kGreen);
  h_SP_0->SetFillColor(kYellow);
  hs_SP->Add(h_SP_5);
  hs_SP->Add(h_SP_4);
  hs_SP->Add(h_SP_3);
  hs_SP->Add(h_SP_2);
  hs_SP->Add(h_SP_1);
  hs_SP->Add(h_SP_0);


  hs_SP->Draw();

  TLegend *leg_SP = new TLegend(0.82,0.62,0.9,0.92);
  leg_SP->AddEntry(h_SP_0," n=0","f");
  leg_SP->AddEntry(h_SP_1," n=1","f");
  leg_SP->AddEntry(h_SP_2," n=2","f");
  leg_SP->AddEntry(h_SP_3," n=3","f");
  leg_SP->AddEntry(h_SP_4," n=4","f");
  leg_SP->AddEntry(h_SP_5," n=5","f");
  leg_SP->Draw();
  c1 -> Print(pdf, "pdf" );
  delete hs_SP;
  delete leg_SP;
  // End of a Fraction plot

}
void RPC::FillSPHist(){
  h_NumberOfSP_eta->Fill(probe_eta, NumberOfSP());
  h_NumberOfSP_qeta->Fill(probe_charge*probe_eta, NumberOfSP());
  if (abs(probe_eta) < 1.05){
    //cout << NumberOfSP() << endl;
    h_NumberOfSP_LumiBlock->Fill(LumiBlock, NumberOfSP());
    h_NumberOfSP_pt_barrel->Fill(probe_pt/1000., NumberOfSP());
    h_NumberOfSP_pass_barrel->Fill(probe_mesSA_pass->at(N50), NumberOfSP());
  }
}


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

void RPC::Display(Long64_t begin_entry, Long64_t limit_entry, TString pdf)
{
   // Prepare Canvas
   TCanvas *c2 = new TCanvas("c2", "c2", 10, 10, 1020, 700);
   c2->SetGrid();
   c2->SetRightMargin(0.20);
   c2->SetLeftMargin(0.23);
   c2->SetBottomMargin(0.20);
   c2->Print(pdf + "[", "pdf");

   c2->Print(pdf, "pdf");

   // Prepare Loop
   if (fChain == 0) return;
   int nLoop = fChain -> GetEntries();
   //Long64_t nentries = fChain->GetEntriesFast();
   double entries = fChain->GetEntries();
   cout << "[INFO]: Nentries:" << entries << endl;
   Long64_t nbytes = 0, nb = 0;

   int current_entry = 0;
   // Start Loop
   for (Long64_t jentry = begin_entry; jentry<nLoop;jentry++) {
     int  ientry = LoadTree(jentry);
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if (ientry < 0) break;
     // Check GRL
      //if (GRLlist(LumiBlock)){
      //  continue;
      //}

      tag_proc = NTagProc;
      switch (tag_proc) {
        case 1: //Jpsi until L2
          // Check TAP
          if (!(probe_mesEFTAG_pass -> at(N4) > -1 && probe_mesL1_pass -> at(N4) > 0 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
            continue;
          }
          break;
        case 2: //Jpsi from L2:
          break;
        case 3: //Z
          // Check TAP
          if (!(probe_mesEFTAG_pass -> at(N50) > -1 && probe_mesL1_pass -> at(N50) > 0 && probe_mesSA_pass -> at(N50) == -2 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
            continue;
          }
          // Check Barrel
          if (abs(probe_eta) > 1.05){
            continue;
          }
          // Check NumberOfSP
          if (NumberOfSP() != 3){
            continue;
          }
          // Check size of mdtHits
          if ((probe_mesSA_mdtHitChamber -> at(N50)).size() == 0){
            continue;
          }
          // Check sAddress to avoid Spercial sector
          if (probe_mesSA_sAddress -> at(N50) == 1 || probe_mesSA_sAddress -> at(N50) == 3){
            continue;
          }

          cout << "size: " << (probe_mesSA_mdtHitChamber -> at(N50)).size() << endl;
          cout << "eta: " << (probe_mesSA_eta -> at(N50)) << endl;

          // Size of Mdt hits and RPC hits
          const Int_t nMDT = (probe_mesSA_mdtHitChamber -> at(N50)).size();
          const Int_t nRPC = (probe_mesSA_rpcHitR -> at(N50)).size();

          // Offline segment
          // Set superpoint and segment for each station
          TGraph gr_segment = TGraph(probe_segment_n);
          for ( int i = 0; i < probe_segment_n; i++){
            if (probe_segment_chamberIndex[i] >= 0 && probe_segment_chamberIndex[i] <= 5) {
              double R = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i]);
              double Z = probe_segment_z[i];
              gr_segment.SetPoint(i, Z / 1000., R / 1000.);
            }
          }
          // RPC hits
          TGraph gr_RPC = TGraph((probe_mesSA_rpcHitX ->at(N50)).size());
          for (int i = 0; i< nRPC; i++){
            if (probe_mesSA_rpcHitStationNumber->at(N50)[i] == 1. ||
                probe_mesSA_rpcHitStationNumber->at(N50)[i] == 2. ||
                probe_mesSA_rpcHitStationNumber->at(N50)[i] == 5. ||
                probe_mesSA_rpcHitStationNumber->at(N50)[i] == 6.) {
              gr_RPC.SetPoint(i, (probe_mesSA_rpcHitZ -> at(N50)[i]) / 1000., (probe_mesSA_rpcHitR -> at(N50)[i]) / 1000.);
            }
          }

          TGraph gr_SP = TGraph(NumberOfSP());

          // SetPoint each SuperPoint
          if (abs(probe_mesSA_superPointZ_BI -> at(N50)) > 0. && probe_mesSA_superPointZ_BI -> at(N50) > -99990.){
            cout << "SP_BI: " << probe_mesSA_superPointR_BI -> at(N50) << endl;
            //cout << "sAdress: " << probe_mesSA_sAddress -> at(N50) << endl;
            gr_SP.SetPoint(0 , (probe_mesSA_superPointZ_BI -> at(N50)) / 1000. , (probe_mesSA_superPointR_BI -> at(N50)) / 1000.);
          }
          if (abs(probe_mesSA_superPointZ_BM -> at(N50)) > 0. && probe_mesSA_superPointZ_BM -> at(N50) > -99990.){
            cout << "SP_BM: " << probe_mesSA_superPointR_BM -> at(N50) << endl;
            gr_SP.SetPoint(1 , (probe_mesSA_superPointZ_BM -> at(N50)) / 1000. , (probe_mesSA_superPointR_BM -> at(N50)) / 1000.);
          }
          if (abs(probe_mesSA_superPointZ_BO -> at(N50)) > 0. && probe_mesSA_superPointZ_BO -> at(N50) > -99990.){
            cout << "SP_BO: " << probe_mesSA_superPointR_BO -> at(N50) << endl;
            gr_SP.SetPoint(2 , (probe_mesSA_superPointZ_BO -> at(N50)) / 1000. , (probe_mesSA_superPointR_BO -> at(N50)) / 1000.);
          }

          TGraph gr = TGraph(nMDT); //各点が(0,0)で初期化される
          // Inlier Mdt Hit
          TGraph gr_MdtHit_Inlier_BI = TGraph(nMDT); //各点が(0,0)で初期化される
          TGraph gr_MdtHit_Inlier_BM = TGraph(nMDT); //各点が(0,0)で初期化される
          TGraph gr_MdtHit_Inlier_BO = TGraph(nMDT); //各点が(0,0)で初期化される
          // Outlier Mdt Hit
          TGraph gr_MdtHit_Outlier_BI = TGraph(nMDT); //各点が(0,0)で初期化される
          TGraph gr_MdtHit_Outlier_BM = TGraph(nMDT); //各点が(0,0)で初期化される
          TGraph gr_MdtHit_Outlier_BO = TGraph(nMDT); //各点が(0,0)で初期化される

          double Z_BI=0, R_BI=0;
          double Z_BM=0, R_BM=0;
          double Z_BO=0, R_BO=0;
          double Zmin_BI=0, Zmax_BI=0, Rmin_BI=0, Rmax_BI=0;
          double Zmin_BM=0, Zmax_BM=0, Rmin_BM=0, Rmax_BM=0;
          double Zmin_BO=0, Zmax_BO=0, Rmin_BO=0, Rmax_BO=0;
          int nBI=0;
          int nBM=0;
          int nBO=0;
          int nBI_Inlier=0;
          int nBM_Inlier=0;
          int nBO_Inlier=0;
          int nBI_Outlier=0;
          int nBM_Outlier=0;
          int nBO_Outlier=0;
          double deltaZ = 1.;
          double deltaR = 0.4;

          // Setpoint each Mdt Hit
          for (Int_t i=0;i<nMDT;++i) {
            cout << probe_mesSA_mdtHitChamber -> at(N50)[i] << endl;
            cout << "Z: " << probe_mesSA_mdtHitZ -> at(N50)[i]/1000. << ", R: " << probe_mesSA_mdtHitR -> at(N50)[i]/1000. << endl;
            cout << nBI << ":" << nBM << ":" << nBO << endl;
            gr.SetPoint(i , (probe_mesSA_mdtHitZ -> at(N50))[i] / 1000. , (probe_mesSA_mdtHitR -> at(N50))[i] / 1000.);

            //BI
            if (probe_mesSA_mdtHitChamber -> at(N50)[i] == 0){
              if (probe_mesSA_mdtHitIsOutlier -> at(N50)[i] == 1) {
                gr_MdtHit_Outlier_BI.SetPoint(i , (probe_mesSA_mdtHitZ -> at(N50))[i] / 1000. , (probe_mesSA_mdtHitR -> at(N50))[i] / 1000.);
                nBI_Outlier += 1;
              } else{
                gr_MdtHit_Inlier_BI.SetPoint(i , (probe_mesSA_mdtHitZ -> at(N50))[i] / 1000. , (probe_mesSA_mdtHitR -> at(N50))[i] / 1000.);
                nBI_Inlier += 1;
              }
              Z_BI += probe_mesSA_mdtHitZ -> at(N50)[i]/1000.;
              R_BI += probe_mesSA_mdtHitR -> at(N50)[i]/1000.;
              nBI += 1;
            }
            //BM
            if (probe_mesSA_mdtHitChamber -> at(N50)[i] == 1){
              if (probe_mesSA_mdtHitIsOutlier -> at(N50)[i] == 1) {
                gr_MdtHit_Outlier_BM.SetPoint(i , (probe_mesSA_mdtHitZ -> at(N50))[i] / 1000. , (probe_mesSA_mdtHitR -> at(N50))[i] / 1000.);
                nBM_Outlier += 1;
              } else{
                gr_MdtHit_Inlier_BM.SetPoint(i , (probe_mesSA_mdtHitZ -> at(N50))[i] / 1000. , (probe_mesSA_mdtHitR -> at(N50))[i] / 1000.);
                nBM_Inlier += 1;
              }
              Z_BM += probe_mesSA_mdtHitZ -> at(N50)[i]/1000.;
              R_BM += probe_mesSA_mdtHitR -> at(N50)[i]/1000.;
              nBM += 1;
            }
            //BO
            if (probe_mesSA_mdtHitChamber -> at(N50)[i] == 2){
              if (probe_mesSA_mdtHitIsOutlier -> at(N50)[i] == 1) {
                gr_MdtHit_Outlier_BO.SetPoint(i , (probe_mesSA_mdtHitZ -> at(N50))[i] / 1000. , (probe_mesSA_mdtHitR -> at(N50))[i] / 1000.);
                nBO_Outlier += 1;
              } else{
                gr_MdtHit_Inlier_BO.SetPoint(i , (probe_mesSA_mdtHitZ -> at(N50))[i] / 1000. , (probe_mesSA_mdtHitR -> at(N50))[i] / 1000.);
                nBO_Inlier += 1;
              }
              Z_BO += probe_mesSA_mdtHitZ -> at(N50)[i]/1000.;
              R_BO += probe_mesSA_mdtHitR -> at(N50)[i]/1000.;
              nBO += 1;
            }
          }

          // Avoid nBI==0
          if (nBI != 0) {
            Zmax_BI = Z_BI/nBI + deltaZ;
            Zmin_BI = Z_BI/nBI - deltaZ;
            Rmax_BI = R_BI/nBI + deltaR;
            Rmin_BI = R_BI/nBI - deltaR;
          }else{
            R_BI = 5;
            deltaR = 1;
            Z_BI = (R_BI - probe_mesSA_roadBw->at(N50)[0]/1000.) / (probe_mesSA_roadAw->at(N50)[0]);
            Zmax_BI = Z_BI + deltaZ;
            Zmin_BI = Z_BI - deltaZ;
            Rmax_BI = R_BI + deltaR;
            Rmin_BI = R_BI - deltaR;
          }
          // Avoid nBM==0
          if (nBM != 0) {
            Zmax_BM = Z_BM/nBM + deltaZ;
            Zmin_BM = Z_BM/nBM - deltaZ;
            Rmax_BM = R_BM/nBM + deltaR;
            Rmin_BM = R_BM/nBM - deltaR;
          }else{
            R_BM = 5;
            deltaR = 1;
            Z_BM = (R_BM - probe_mesSA_roadBw->at(N50)[0]/1000.) / (probe_mesSA_roadAw->at(N50)[0]);
            Zmax_BM = Z_BM + deltaZ;
            Zmin_BM = Z_BM - deltaZ;
            Rmax_BM = R_BM + deltaR;
            Rmin_BM = R_BM - deltaR;
          }
          // Avoid nBO==0
          if (nBO != 0) {
            Zmax_BO = Z_BO/nBO + deltaZ;
            Zmin_BO = Z_BO/nBO - deltaZ;
            Rmax_BO = R_BO/nBO + deltaR;
            Rmin_BO = R_BO/nBO - deltaR;
          }else{
            R_BO = 5;
            deltaR = 1;
            Z_BO = (R_BO - probe_mesSA_roadBw->at(N50)[0]/1000.) / (probe_mesSA_roadAw->at(N50)[0]);
            Zmax_BO = Z_BO + deltaZ;
            Zmin_BO = Z_BO - deltaZ;
            Rmax_BO = R_BO + deltaR;
            Rmin_BO = R_BO - deltaR;
          }
          // Road
          TF1 f_road_BI = TF1("f_road_BI", "[0]*x+[1]", -20, 20);
          cout << "road: " << probe_mesSA_roadAw -> at(N50)[0] << ": " << probe_mesSA_roadBw->at(N50)[0] << endl;
          cout << "road: " << probe_mesSA_roadAw -> at(N50)[1] << ": " << probe_mesSA_roadBw->at(N50)[1] << endl;
          cout << "road: " << probe_mesSA_roadAw -> at(N50)[2] << ": " << probe_mesSA_roadBw->at(N50)[2] << endl;
          f_road_BI.SetTitle(";Z [m];R [m]");
          f_road_BI.SetParameter(0,probe_mesSA_roadAw -> at(N50)[0]);
          f_road_BI.SetParameter(1,probe_mesSA_roadBw -> at(N50)[0]/1000.);
          f_road_BI.SetLineColor(15);
          f_road_BI.SetLineWidth(2);
          f_road_BI.SetLineStyle(2);

          TF1 f_road_BI_plus = TF1("f_road_BI_plus", "[0]*x+[1]", -20, 20);
          f_road_BI_plus.SetParameter(0,probe_mesSA_roadAw -> at(N50)[0]);
          f_road_BI_plus.SetParameter(1,( (probe_mesSA_roadBw -> at(N50)[0]/1000.) + rWidthToBw(probe_mesSA_roadAw->at(N50)[0],0.4)));
          f_road_BI_plus.SetLineColor(15);
          f_road_BI_plus.SetLineWidth(2);
          f_road_BI_plus.SetLineStyle(1);

          TF1 f_road_BI_minus = TF1("f_road_BI_minus", "[0]*x+[1]", -20, 20);
          f_road_BI_minus.SetParameter(0,probe_mesSA_roadAw -> at(N50)[0]);
          f_road_BI_minus.SetParameter(1,( (probe_mesSA_roadBw -> at(N50)[0]/1000.) - rWidthToBw(probe_mesSA_roadAw->at(N50)[0],0.4)));
          f_road_BI_minus.SetLineColor(15);
          f_road_BI_minus.SetLineWidth(2);
          f_road_BI_minus.SetLineStyle(1);


          TF1 f_road_BM = TF1("f_road_BM", "[0]*x+[1]", -20, 20);
          f_road_BM.SetTitle(";Z [m];R [m]");
          f_road_BM.SetParameter(0,probe_mesSA_roadAw -> at(N50)[1]);
          f_road_BM.SetParameter(1,probe_mesSA_roadBw -> at(N50)[1]/1000.);
          f_road_BM.SetLineColor(15);
          f_road_BM.SetLineWidth(2);
          f_road_BM.SetLineStyle(2);

          TF1 f_road_BM_plus = TF1("f_road_BM_plus", "[0]*x+[1]", -20, 20);
          f_road_BM_plus.SetParameter(0,probe_mesSA_roadAw -> at(N50)[1]);
          f_road_BM_plus.SetParameter(1,( (probe_mesSA_roadBw -> at(N50)[1]/1000.) + rWidthToBw(probe_mesSA_roadAw->at(N50)[1],0.2)));
          f_road_BM_plus.SetLineColor(15);
          f_road_BM_plus.SetLineWidth(2);
          f_road_BM_plus.SetLineStyle(1);

          TF1 f_road_BM_minus = TF1("f_road_BM_minus", "[0]*x+[1]", -20, 20);
          f_road_BM_minus.SetParameter(0,probe_mesSA_roadAw -> at(N50)[1]);
          f_road_BM_minus.SetParameter(1,( (probe_mesSA_roadBw -> at(N50)[1]/1000.) - rWidthToBw(probe_mesSA_roadAw->at(N50)[1],0.2)));
          f_road_BM_minus.SetLineColor(15);
          f_road_BM_minus.SetLineWidth(2);
          f_road_BM_minus.SetLineStyle(1);

          TF1 f_road_BO = TF1("f_road_BO", "[0]*x+[1]", -20, 20);
          f_road_BO.SetTitle(";Z [m];R [m]");
          f_road_BO.SetParameter(0,probe_mesSA_roadAw -> at(N50)[2]);
          f_road_BO.SetParameter(1,probe_mesSA_roadBw -> at(N50)[2]/1000.);
          f_road_BO.SetLineColor(15);
          f_road_BO.SetLineWidth(2);
          f_road_BO.SetLineStyle(2);

          TF1 f_road_BO_plus = TF1("f_road_BO_plus", "[0]*x+[1]", -20, 20);
          f_road_BO_plus.SetParameter(0,probe_mesSA_roadAw -> at(N50)[2]);
          f_road_BO_plus.SetParameter(1,( (probe_mesSA_roadBw -> at(N50)[2]/1000.) + rWidthToBw(probe_mesSA_roadAw->at(N50)[2],0.4)));
          f_road_BO_plus.SetLineColor(15);
          f_road_BO_plus.SetLineWidth(2);
          f_road_BO_plus.SetLineStyle(1);

          TF1 f_road_BO_minus = TF1("f_road_BO_minus", "[0]*x+[1]", -20, 20);
          f_road_BO_minus.SetParameter(0,probe_mesSA_roadAw -> at(N50)[2]);
          f_road_BO_minus.SetParameter(1,( (probe_mesSA_roadBw -> at(N50)[2]/1000.) - rWidthToBw(probe_mesSA_roadAw->at(N50)[2],0.4)));
          f_road_BO_minus.SetLineColor(15);
          f_road_BO_minus.SetLineWidth(2);
          f_road_BO_minus.SetLineStyle(1);


          // RoI
          TF1 f_roi = TF1("f_roi", "[0]*x", -20, 20);
          f_roi.SetTitle(";Z [m];R [m]");
          f_roi.SetParameter(0, tan((2*atan(exp(-probe_mesSA_roiEta->at(N50))))));
          f_roi.SetLineColor(kYellow+2);
          f_roi.SetLineWidth(3);
          f_roi.SetLineStyle(10);

          // Set Legend
          TLegend leg = TLegend(0.805,0.22,0.99,0.95);
          leg.SetTextSize(0.035);
          leg.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
          leg.AddEntry(&gr,Form("MDT hit (%d)",nMDT),"p");
          leg.AddEntry(&gr_RPC,Form("RPC hit (%d)",nRPC),"p");
          leg.AddEntry(&gr_SP,Form("SuperPoint (%d)",NumberOfSP()),"p");
          leg.AddEntry(&f_road_BI,"Road","l");
          leg.AddEntry(&f_roi,"RoI center","l");

          // Set Legend
          TLegend leg_left = TLegend(-0.025,0.12,0.10,0.8);
          leg_left.SetTextSize(0.03);
          leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe p_{T}}{: %4.3f [GeV]}",probe_pt/1000.),"");
          leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe p_{T}}{: %4.3f [GeV]}",abs(probe_mesSA_pt->at(N50))),"");
          leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe #eta}{: %4.3f}",probe_eta),"");
          leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe #eta}{: %4.3f}",probe_mesSA_eta->at(N50)),"");
          leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe #phi}{: %4.3f}",probe_phi),"");
          leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe #phi}{: %4.3f}",probe_mesSA_phi->at(N50)),"");
          leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe #eta}{: %d}",probe_mesSA_pass->at(N50)),"");


          // Draw Display
          double Zmax = 20;
          double Zmin = -20;
          TH1 *frame = c2->DrawFrame(Zmin,0,Zmax,12);
          frame->GetXaxis()->SetTitle("Z [m]");
          frame->GetYaxis()->SetTitle("R [m]");
          gr_segment.Draw("P");
          f_roi.Draw("same");
          gr_segment.SetMarkerStyle(21);
          gr_segment.SetMarkerSize(2);
          gr_segment.SetMarkerColor(6);
          gr_segment.GetXaxis()->SetLimits(Zmin,Zmax);
          gr_segment.SetTitle(";Z [m];R [m]");
          gr_segment.GetYaxis()->SetRangeUser(0,12);
          f_road_BI.Draw("same");
          f_road_BM.Draw("same");
          f_road_BO.Draw("same");
          gr.SetMarkerStyle(24);
          gr.SetMarkerSize(1);
          gr.SetMarkerColor(2);
          gr.Draw("P");
          leg.Draw();
          leg_left.Draw();
          gr_RPC.SetMarkerStyle(22);
          gr_RPC.SetMarkerSize(1);
          gr_RPC.SetMarkerColor(kCyan-3);
          gr_RPC.Draw("P, same");
          gr_SP.SetMarkerStyle(8);
          gr_SP.SetMarkerSize(2);
          gr_SP.SetMarkerColor(4);
          gr_SP.Draw("P, same");
          c2->Print(pdf, "pdf");
          c2->RedrawAxis();
          delete frame;

          // Set Legend for BI
          TLegend leg_BI = TLegend(0.805,0.22,0.99,0.95);
          leg_BI.SetHeader("Barrel Inner","C");
          leg_BI.SetTextSize(0.035);
          leg_BI.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
          leg_BI.AddEntry(&gr_MdtHit_Inlier_BI,  Form("#splitline{Inlier}{MDT hit (%d)}",nBI_Inlier),  "p");
          leg_BI.AddEntry(&gr_MdtHit_Outlier_BI, Form("#splitline{Outlier}{MDT hit (%d)}",nBI_Outlier), "p");
          leg_BI.AddEntry(&gr_RPC,               "RPC hit",         "p");
          leg_BI.AddEntry(&gr_SP,                "SuperPoint",      "p");
          leg_BI.AddEntry(&f_road_BI,            "Road",            "l");
          leg_BI.AddEntry(&f_roi,"RoI center","l");
          TH1* frame_BI = c2->DrawFrame(Zmin_BI,Rmin_BI,Zmax_BI,Rmax_BI);
          frame_BI->GetXaxis()->SetTitle("Z [m]");
          frame_BI->GetYaxis()->SetTitle("R [m]");
          f_road_BI.Draw("same");
          f_road_BI_plus.Draw("same");
          f_road_BI_minus.Draw("same");
          f_roi.Draw("same");
          cout << "Zmin_BI: " << Zmin_BI <<endl;
          cout << "Zmax_BI: " << Zmax_BI <<endl;
          f_road_BI.GetXaxis()->SetLimits(Zmin_BI,Zmax_BI);
          f_road_BI.GetYaxis()->SetRangeUser(Rmin_BI,Rmax_BI);
          gr_MdtHit_Inlier_BI.SetMarkerColor(kGreen + 2);
          gr_MdtHit_Inlier_BI.SetMarkerStyle(24);
          gr_MdtHit_Inlier_BI.SetMarkerSize(1);
          gr_MdtHit_Inlier_BI.Draw("P, same");
          gr_segment.Draw("P, same");

          leg_BI.Draw();

          gr_MdtHit_Outlier_BI.SetMarkerColor(kRed);
          gr_MdtHit_Outlier_BI.SetMarkerStyle(24);
          gr_MdtHit_Outlier_BI.SetMarkerSize(1);
          gr_MdtHit_Outlier_BI.Draw("P,same");

          gr_SP.Draw("P, same");
          c2->Print(pdf, "pdf");
          c2->RedrawAxis();
          delete frame_BI;

          // Set Legend for BM
          TLegend leg_BM = TLegend(0.805,0.22,0.99,0.95);
          leg_BM.SetHeader("Barrel Middle","C");
          leg_BM.SetTextSize(0.035);
          leg_BM.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
          leg_BM.AddEntry(&gr_MdtHit_Inlier_BM,  Form("#splitline{Inlier}{MDT hit (%d)}",nBM_Inlier),  "p");
          leg_BM.AddEntry(&gr_MdtHit_Outlier_BM, Form("#splitline{Outlier}{MDT hit (%d)}",nBM_Outlier), "p");
          leg_BM.AddEntry(&gr_RPC,"RPC hit","p");
          leg_BM.AddEntry(&gr_SP,"SuperPoint","p");
          leg_BM.AddEntry(&f_road_BM,"Road","l");
          leg_BM.AddEntry(&f_roi,"RoI center","l");
          TH1* frame_BM = c2->DrawFrame(Zmin_BM,Rmin_BM,Zmax_BM,Rmax_BM);
          frame_BM->GetXaxis()->SetTitle("Z [m]");
          frame_BM->GetYaxis()->SetTitle("R [m]");
          f_road_BM.Draw("same");
          f_road_BM_plus.Draw("same");
          f_road_BM_minus.Draw("same");
          f_roi.Draw("same");
          f_road_BM.GetXaxis()->SetLimits(Zmin_BM,Zmax_BM);
          f_road_BM.GetYaxis()->SetRangeUser(Rmin_BM,Rmax_BM);
          gr_MdtHit_Inlier_BM.SetMarkerColor(kGreen + 2);
          gr_MdtHit_Inlier_BM.SetMarkerStyle(24);
          gr_MdtHit_Inlier_BM.SetMarkerSize(1);
          gr_MdtHit_Inlier_BM.Draw("P, same");
          gr_segment.Draw("P, same");

          leg_BM.Draw();

          gr_RPC.Draw("P, same");

          gr_MdtHit_Outlier_BM.SetMarkerColor(kRed);
          gr_MdtHit_Outlier_BM.SetMarkerStyle(24);
          gr_MdtHit_Outlier_BM.SetMarkerSize(1);
          gr_MdtHit_Outlier_BM.Draw("P,same");
          gr_SP.Draw("P, same");
          c2->Print(pdf, "pdf");
          c2->RedrawAxis();
          delete frame_BM;

          // Set Legend for BO
          TLegend leg_BO = TLegend(0.81,0.22,0.99,0.95);
          leg_BO.SetHeader("Barrel Outer","C");
          leg_BO.SetTextSize(0.035);
          leg_BO.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
          leg_BO.AddEntry(&gr_MdtHit_Inlier_BO,  Form("#splitline{Inlier}{MDT hit (%d)}",nBO_Inlier),  "p");
          leg_BO.AddEntry(&gr_MdtHit_Outlier_BO, Form("#splitline{Outlier}{MDT hit (%d)}",nBO_Outlier), "p");
          leg_BO.AddEntry(&gr_RPC,"RPC hit","p");
          leg_BO.AddEntry(&gr_SP,"SuperPoint","p");
          leg_BO.AddEntry(&f_road_BO,"Road","l");
          leg_BO.AddEntry(&f_roi,"RoI center","l");
          TH1* frame_BO = c2->DrawFrame(Zmin_BO,Rmin_BO,Zmax_BO,Rmax_BO);
          frame_BO->GetXaxis()->SetTitle("Z [m]");
          frame_BO->GetYaxis()->SetTitle("R [m]");
          f_road_BO.Draw("same");
          f_road_BO_plus.Draw("same");
          f_road_BO_minus.Draw("same");
          f_roi.Draw("same");
          f_road_BO.GetXaxis()->SetLimits(Zmin_BO,Zmax_BO);
          f_road_BO.GetYaxis()->SetRangeUser(Rmin_BO,Rmax_BO);
          gr_MdtHit_Inlier_BO.SetMarkerColor(kGreen + 2);
          gr_MdtHit_Inlier_BO.SetMarkerStyle(24);
          gr_MdtHit_Inlier_BO.SetMarkerSize(1);
          gr_MdtHit_Inlier_BO.Draw("P, same");
          gr_segment.Draw("P, same");

          leg_BO.Draw();

          gr_RPC.Draw("P, same");

          gr_MdtHit_Outlier_BO.SetMarkerColor(kRed);
          gr_MdtHit_Outlier_BO.SetMarkerStyle(24);
          gr_MdtHit_Outlier_BO.SetMarkerSize(1);
          gr_MdtHit_Outlier_BO.Draw("P,same");
          gr_SP.Draw("P, same");
          c2->Print(pdf, "pdf");
          c2->RedrawAxis();
          delete frame_BO;

          current_entry += 1;
          break;
      }

      if (current_entry == limit_entry) {
        break;
      }
      c2->Clear();
   } // end of each entry
   c2->Print(pdf + "]", "pdf");
   delete c2;
}

bool GRLlist(int LumiBlock){
  bool lb_1 = (LumiBlock > 110 && LumiBlock < 124);
  bool lb_2= (LumiBlock > 126 && LumiBlock < 130);
  bool lb_3= (LumiBlock > 169 && LumiBlock < 239);

  return (lb_1 || lb_2 || lb_3);
}

void RPC::FillMdtHist(){
  int nMdtBI = 0;
  int nMdtBM = 0;
  int nMdtBO = 0;
  int nMdt = 0;
  for ( uint32_t i = 0; i < (probe_mesSA_mdtHitIsOutlier -> at(N50)).size();i++){
    // Count each MDT hit regardless of IsOutlier or not
    nMdt += 1;

    // Separate histograms between Inlier and Outlier
    switch (probe_mesSA_mdtHitIsOutlier -> at(N50)[i]) {
      case 0: // Inlier
        // Separate each MDT station
        switch (probe_mesSA_mdtHitChamber -> at(N50)[i]) {
          case 0:
            nMdtBI += 1;
            //cout << "BI" << endl;
            h_ResidualMdt_Inlier_eta_BI    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
            if (abs(probe_eta) < 1.05) h_ResidualMdt_Inlier_pt_barrel_BI -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
            break;
          case 1:
            nMdtBM += 1;
            //cout << "BM" << endl;
            h_ResidualMdt_Inlier_eta_BM    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
            if (abs(probe_eta) < 1.05) h_ResidualMdt_Inlier_pt_barrel_BM -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
            break;
          case 2:
            nMdtBO += 1;
            //cout << "BO" << endl;
            h_ResidualMdt_Inlier_eta_BO    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
            if (abs(probe_eta) < 1.05) h_ResidualMdt_Inlier_pt_barrel_BO -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
            break;
        }
        break;
      case 1: // Outlier
        // Separate each MDT station
        switch (probe_mesSA_mdtHitChamber -> at(N50)[i]) {
          case 0:
            nMdtBI += 1;
            //cout << "BI" << endl;
            h_ResidualMdt_Outlier_eta_BI    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
            if (abs(probe_eta) < 1.05) h_ResidualMdt_Outlier_pt_barrel_BI -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
            break;
          case 1:
            nMdtBM += 1;
            //cout << "BM" << endl;
            h_ResidualMdt_Outlier_eta_BM    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
            if (abs(probe_eta) < 1.05) h_ResidualMdt_Outlier_pt_barrel_BM -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
            break;
          case 2:
            nMdtBO += 1;
            //cout << "BO" << endl;
            h_ResidualMdt_Outlier_eta_BO    -> Fill(probe_eta, probe_mesSA_mdtHitResidual -> at(N50)[i]);
            if (abs(probe_eta) < 1.05) h_ResidualMdt_Outlier_pt_barrel_BO -> Fill(probe_pt/1000., probe_mesSA_mdtHitResidual -> at(N50)[i]);
            break;
        }
        break;
    }
  }

  h_NumberOfMdt_eta    -> Fill(probe_eta, nMdt);
  h_NumberOfMdt_eta_BI -> Fill(probe_eta, nMdtBI);
  h_NumberOfMdt_eta_BM -> Fill(probe_eta, nMdtBM);
  h_NumberOfMdt_eta_BO -> Fill(probe_eta, nMdtBO);

  // Select barrel region because of pT histograms
  if (abs(probe_eta) < 1.05){
    //h_NumberOfMdt_LumiBlock->Fill(LumiBlock, (probe_mesSA_mdtHitIsOutlier -> at(14)).size());
    h_NumberOfMdt_pt_barrel_BI->Fill(probe_pt/1000., nMdtBI);
    h_NumberOfMdt_pt_barrel_BM->Fill(probe_pt/1000., nMdtBM);
    h_NumberOfMdt_pt_barrel_BO->Fill(probe_pt/1000., nMdtBO);
  }
}

void RPC::DrawFractionOfnMDTs(TH2F* h_NumberOfMdt, TCanvas* c1, TString pdf){

  // Begin of a Fraction plot
  TProfile *pf_Mdt = h_NumberOfMdt->ProfileX();
  pf_Mdt->Draw();
  pf_Mdt->GetYaxis()->SetTitle("Average number of Mdt");
  c1 -> Print(pdf, "pdf" );

  const int delta = 0.1;

  TH1D* h_all_Mdt = h_NumberOfMdt->ProjectionX("_all");
  TH1D* h_1_Mdt = h_NumberOfMdt->ProjectionX("1", h_NumberOfMdt->GetYaxis()->FindBin(0. - delta),  h_NumberOfMdt->GetYaxis()->FindBin(5. + delta));
  TH1D* h_2_Mdt = h_NumberOfMdt->ProjectionX("2", h_NumberOfMdt->GetYaxis()->FindBin(6. - delta),  h_NumberOfMdt->GetYaxis()->FindBin(10. + delta));
  TH1D* h_3_Mdt = h_NumberOfMdt->ProjectionX("3", h_NumberOfMdt->GetYaxis()->FindBin(11. - delta), h_NumberOfMdt->GetYaxis()->FindBin(15. + delta));
  TH1D* h_4_Mdt = h_NumberOfMdt->ProjectionX("4", h_NumberOfMdt->GetYaxis()->FindBin(16. - delta), h_NumberOfMdt->GetYaxis()->FindBin(20. + delta));
  TH1D* h_5_Mdt = h_NumberOfMdt->ProjectionX("5", h_NumberOfMdt->GetYaxis()->FindBin(21. - delta), h_NumberOfMdt->GetYaxis()->FindBin(25. + delta));
  TH1D* h_6_Mdt = h_NumberOfMdt->ProjectionX("6", h_NumberOfMdt->GetYaxis()->FindBin(26. - delta), h_NumberOfMdt->GetYaxis()->FindBin(30. + delta));

  h_6_Mdt->Divide(h_all_Mdt);
  h_5_Mdt->Divide(h_all_Mdt);
  h_4_Mdt->Divide(h_all_Mdt);
  h_3_Mdt->Divide(h_all_Mdt);
  h_2_Mdt->Divide(h_all_Mdt);
  h_1_Mdt->Divide(h_all_Mdt);

  //TString ytitle = h_NumberOfMdt -> GetYaxis()-> GetTitle();

  THStack *hs_Mdt = new THStack("hs_Mdt",Form(";%s;Fraction of nMDTs",h_NumberOfMdt->GetXaxis()->GetTitle()));
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

  TLegend *leg_Mdt = new TLegend(0.81,0.22,0.99,0.95);
  leg_Mdt->AddEntry(h_1_Mdt," n = 0~5","f");
  leg_Mdt->AddEntry(h_2_Mdt," n = 5~10","f");
  leg_Mdt->AddEntry(h_3_Mdt," n = 10~15","f");
  leg_Mdt->AddEntry(h_4_Mdt," n = 15~20","f");
  leg_Mdt->AddEntry(h_5_Mdt," n = 20~25","f");
  leg_Mdt->AddEntry(h_6_Mdt," n = 25~30","f");
  leg_Mdt->Draw();
  c1 -> Print(pdf, "pdf" );
  delete leg_Mdt;
  delete hs_Mdt;
  // End of a Fraction plot
}


void RPC::FillProbeHist(){
  if (abs(probe_eta) > 1.05){ return;}
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
  h_eff_pt_mu50_L1SA->GetYaxis()->SetRangeUser(0,1.1);
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

double RPC::rWidthToBw(double aw, double rWidth)
{
  return rWidth*sqrt(aw*aw+1);
}

double RPC::calc_residual(double aw,double bw,double x,double y)
{
  const double ZERO_LIMIT = 1e-4;
  if( fabs(aw) < ZERO_LIMIT ) return y-bw;
  double ia  = 1/aw;
  double iaq = ia*ia;
  double dz  = x - (y-bw)*ia;
  return dz/sqrt(1.+iaq);
}
