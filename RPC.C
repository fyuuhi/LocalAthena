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
#include "TCanvas.h"
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


//==================================================================
//main function
//==================================================================
int main(int argc, char **argv){
  rootlogon();
  cout << "---start---" << endl;
  TColor::InvertPalette();
  //BarrelResidual( -1, 100000, "data18" );
  //BarrelResidual( -1, 100000, "data15" );
  //ForwardResidual(-1, 100000, "data16" );
  //ForwardResidual(-1, 100000, "data17" );
  //app.Run();

  // tree1
  TChain *tree1 = new TChain("t_tap", "t_tap");
  tree1 -> Add("/gpfs/fs2001/yfukuhar/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0.root");

  RPC t_349014(tree1); 

  t_349014.Loop(-1, 10000);
  cout << "[INFO]: Loop SUCCESS" << endl;

  t_349014.DrawHist("t_349014.pdf");
  cout << "[INFO]: DrawHist SUCCESS" << endl;

  t_349014.End();
  cout << "[INFO]: End SUCCESS" << endl;


  // tree2
  TChain *tree2 = new TChain("t_tap", "t_tap");
  tree2 -> Add("/gpfs/fs2001/yfukuhar/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_jpzYFV3GRL_EXT0.root");

  RPC t_349533(tree2); 

  t_349533.Loop(-1, 10000);
  cout << "[INFO]: Loop SUCCESS" << endl;

  t_349533.DrawHist("t_349533.pdf");
  cout << "[INFO]: DrawHist SUCCESS" << endl;

  t_349014.End();
  cout << "[INFO]: End SUCCESS" << endl;


  delete tree1;
  delete tree2;

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
          // Set suerpoint and segment for each station
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

          h_NumberOfSP_eta->Fill(probe_eta, NumberOfSP());
          if (abs(probe_eta) < 1.05){
            //cout << NumberOfSP() << endl;
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
  c1->SetRightMargin(0.20);
  c1->SetLeftMargin(0.23);
  c1->SetBottomMargin(0.20);

  c1 -> Print( pdf + "[", "pdf" );

  h_NumberOfSP_eta->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  TProfile *pf_eta = h_NumberOfSP_eta->ProfileX();
  pf_eta->Draw();
  pf_eta->GetYaxis()->SetTitle("Average number of SP");
  c1 -> Print(pdf, "pdf" );

  TH1D* h_eta_all = h_NumberOfSP_eta->ProjectionX("_all");
  TH1D* h_eta_3 = h_NumberOfSP_eta->ProjectionX("eta_3",h_NumberOfSP_eta->GetYaxis()->FindBin(2.9),h_NumberOfSP_eta->GetYaxis()->FindBin(3.1));
  TH1D* h_eta_2 = h_NumberOfSP_eta->ProjectionX("eta_2",h_NumberOfSP_eta->GetYaxis()->FindBin(1.9),h_NumberOfSP_eta->GetYaxis()->FindBin(2.1));
  TH1D* h_eta_1 = h_NumberOfSP_eta->ProjectionX("eta_1",h_NumberOfSP_eta->GetYaxis()->FindBin(0.9),h_NumberOfSP_eta->GetYaxis()->FindBin(1.1));
  TH1D* h_eta_0 = h_NumberOfSP_eta->ProjectionX("eta_0",h_NumberOfSP_eta->GetYaxis()->FindBin(-0.9),h_NumberOfSP_eta->GetYaxis()->FindBin(0.1));

  h_eta_3->Divide(h_eta_all);
  h_eta_2->Divide(h_eta_all);
  h_eta_1->Divide(h_eta_all);
  h_eta_0->Divide(h_eta_all);


  THStack *hs_eta = new THStack("hs_eta",";eta;Fraction of number of SP");
  h_eta_3->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_eta_2->SetFillColor(kBlue);
  h_eta_1->SetFillColor(kGreen);
  h_eta_0->SetFillColor(kYellow);
  hs_eta->Add(h_eta_3);
  hs_eta->Add(h_eta_2);
  hs_eta->Add(h_eta_1);
  hs_eta->Add(h_eta_0);

  hs_eta->Draw();
  c1 -> Print(pdf, "pdf" );


  h_NumberOfSP_pt_barrel->Draw("colz");
  c1 -> Print(pdf, "pdf" );
  TProfile *pf_pt = h_NumberOfSP_pt_barrel->ProfileX();
  pf_pt->Draw();
  pf_pt->GetYaxis()->SetTitle("Average number of SP");
  c1 -> Print(pdf, "pdf" );

  TH1D* h_pt_all = h_NumberOfSP_pt_barrel->ProjectionX("_all");
  TH1D* h_pt_3 = h_NumberOfSP_pt_barrel->ProjectionX("pt_3",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(2.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(3.1));
  TH1D* h_pt_2 = h_NumberOfSP_pt_barrel->ProjectionX("pt_2",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(1.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(2.1));
  TH1D* h_pt_1 = h_NumberOfSP_pt_barrel->ProjectionX("pt_1",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(0.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(1.1));
  TH1D* h_pt_0 = h_NumberOfSP_pt_barrel->ProjectionX("pt_0",h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(-0.9),h_NumberOfSP_pt_barrel->GetYaxis()->FindBin(0.1));

  h_pt_3->Divide(h_pt_all);
  h_pt_2->Divide(h_pt_all);
  h_pt_1->Divide(h_pt_all);
  h_pt_0->Divide(h_pt_all);


  THStack *hs_pt = new THStack("hs_pt",";pT [GeV];Fraction of number of SP");
  h_pt_3->SetFillColor(kRed);//あらかじめFillColorをSetしておく
  h_pt_2->SetFillColor(kBlue);
  h_pt_1->SetFillColor(kGreen);
  h_pt_0->SetFillColor(kYellow);
  hs_pt->Add(h_pt_3);
  hs_pt->Add(h_pt_2);
  hs_pt->Add(h_pt_1);
  hs_pt->Add(h_pt_0);

  hs_pt->Draw();
  c1 -> Print(pdf, "pdf" );

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

  THStack *hs_pass = new THStack("hs_pass",";pass;Fraction of number of SP");
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
  return number;
}
