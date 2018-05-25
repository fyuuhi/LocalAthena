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




using namespace std;
//declare functions

#include "/gpfs/home/yfukuhar/RootUtils/rootlogon.C"


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

  TChain *tree1 = new TChain("t_tap", "t_tap");
  tree1 -> Add("/gpfs/fs2001/yfukuhar/data/hadd_data18_v3_mu26ivm_ok/user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0/hadd_data18_v3_mu26ivm_ok_user.yfukuhar.00349014.physics_Main.YFTAP.f926_m1955_jpzYFV3GRL_EXT0.root");

  RPC t_349014(tree1); 

  t_349014.Loop();
  cout << "[INFO]: Loop SUCCESS" << endl;

  t_349014.DrawHist();
  cout << "[INFO]: DrawHist SUCCESS" << endl;

  return 0;
}


void RPC::Loop()
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

   Long64_t nentries = fChain->GetEntriesFast();
   double entries = fChain->GetEntries();
   cout << "[INFO]: Nentries:" << entries << endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      int  ientry = LoadTree(jentry);
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (ientry < 0) break;

      //==================================================================
      //Analysis code for a entry
      //==================================================================
      switch (tag_proc) {
        case 1: //Jpsi until L2
          break;
        case 2: //Jpsi from L2:
          break;
        case 3: //Z
          if (probe_mesSA_superPointR_BI->at(0) != 0 && probe_mesSA_superPointR_BI -> at(0) / 1000. > -90) {
            h_superPointRZ_BI -> Fill( probe_mesSA_superPointZ_BI->at(0)/1000, probe_mesSA_superPointR_BI->at(0)/1000); 
          }
          if (probe_mesSA_superPointR_BM->at(0) != 0 && probe_mesSA_superPointR_BM -> at(0) / 1000. > -90) {
            h_superPointRZ_BI -> Fill( probe_mesSA_superPointZ_BM->at(0)/1000, probe_mesSA_superPointR_BM->at(0)/1000); 
          }
          if (probe_mesSA_superPointR_BO->at(0) != 0 && probe_mesSA_superPointR_BO -> at(0) / 1000. > -90) {
            h_superPointRZ_BI -> Fill( probe_mesSA_superPointZ_BO->at(0)/1000, probe_mesSA_superPointR_BO->at(0)/1000); 
          }
          break;
      }

   } // end of each entry

} // end of RPC::Loop()






void RPC::DrawHist(){
  //==================================================================
  //Set Canvas
  //==================================================================
  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1020, 700);
  c1->SetRightMargin(0.20);
  c1->SetLeftMargin(0.23);
  c1->SetBottomMargin(0.20);

  TString pdf = "test.pdf";
  c1 -> Print( pdf + "[", "pdf" );

  h_superPointRZ_BI -> Draw("colz");

  h_superPointRZ_BM -> Draw("colz same");

  h_superPointRZ_BO -> Draw("colz same");
  c1 -> Print(pdf, "pdf" );

  c1 -> Print( pdf + "]", "pdf" );
}

