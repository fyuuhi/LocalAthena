#define RPC_cxx
#include "RPC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH3F.h"
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
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace std;

#include "/gpfs/home/yfukuhar/RootUtils/rootlogon.C"
#include <plog/Log.h> // Step1: include the header.
#include <plog/Appenders/ConsoleAppender.h>
#include <plog/Appenders/ColorConsoleAppender.h>

//#include "/gpfs/fs7001/yfukuhar/workspace/TestCalcEffTool/CalcEffTool/CalcEfficiency/CalcEfficiency/easy_logger.h"
//#include "/gpfs/home/yfukuhar/work/RpcL2MuonSA/RpcL2MuonSA/src/easy_logger.h"
//#include "/gpfs/home/yfukuhar/work/RpcL2MuonSA/RpcL2MuonSA/src/tmp.h"

double Res(double param1, double param2);
bool GRLlist(int LumiBlock);

//==================================================================
//main function (argv[1]: PDF_LABEL, argv[2]: INPUT_NTUPLE, argv[3]: IS_DRAW, arg[4]: IS_EVENTDISPLAY, arg[5]: BEGIN_ENTRY, arg[6]: LIMIT_ENTRY, arg[7]: TAP_TYPE, arg[8]: trig_chain)
//==================================================================
int main(int argc, char *argv[]){
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  //plog::init(plog::info, &consoleAppender);
  plog::init(plog::verbose, &consoleAppender);
  //plog::init(plog::debug, &consoleAppender);

  rootlogon();
  LOGI << "---start---";
  TColor::InvertPalette();
  gStyle->SetLegendBorderSize(0);

  TChain *tree1 = new TChain("t_tap", "t_tap");

  TString PdfLabel = argv[1];
  tree1 -> Add(argv[2]);
  TString isDraw = argv[3];
  TString isEventDisplay = argv[4];
  Int_t begin_entry = atoi(argv[5]);
  Int_t limit_entry = atoi(argv[6]);
  Int_t tap_type = atoi(argv[7]);
  Int_t trig_chain = atoi(argv[8]);
  LOGI << begin_entry << "," << limit_entry;

  TFile *fout = new TFile(Form("outroot/%s.root", PdfLabel.Data()), "recreate");

  // == Core Begin ==
  RPC t_349014(tree1); 

  if (isDraw == "true"){
    //t_349014.Loop(limit_entry, 10000);
    t_349014.Loop(limit_entry, 1);
    LOGI << "Loop SUCCESS";

    t_349014.DrawHist("../plot/DrawHist_" + PdfLabel + ".pdf");
    LOGI << "DrawHist SUCCESS";

    t_349014.CalcEff();
    LOGI << "CalcEff SUCCESS";

    t_349014.DrawEffHist("../plot/DrawEffHist_" + PdfLabel + ".pdf");
    LOGI << "DrawEffHist SUCCESS";

    t_349014.DrawInEffHist("../plot/DrawInEffHist_" + PdfLabel + ".pdf");
    LOGI << "DrawInEffHist SUCCESS";
  }

  if (isEventDisplay == "true"){
    t_349014.Display(tap_type, trig_chain, begin_entry, limit_entry, "../plot/Display_" + PdfLabel + ".pdf");
    LOGI << "Display SUCCESS";
  }
  // == Core End ==

  // Output root file
  fout -> Write();

  t_349014.End();
  LOGI << "End SUCCESS";

  return 0;
}

void RPC::Loop( int Nevents, int DisplayNumber )
{
  const double ZERO_LIMIT = 1e-5;
  if (fChain == 0) return;

  int nLoop;
  if (Nevents == -1) {
    nLoop = fChain -> GetEntries();
  } else {
    nLoop = Nevents;
  }

  //Long64_t nentries = fChain->GetEntriesFast();
  double entries = fChain->GetEntries();
  LOGI << "Nentries:" << entries;

  //const int N50 = 14;


  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nLoop;jentry++) {
    int  ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (ientry < 0) break;
    if( ientry%DisplayNumber == 0){
      LOGI << "now event number -->> " << ientry << "/" << nLoop << " : " << ((double) ientry)/nLoop * 100. << "%";
    }
    LOGD << "EventNumber = " << EventNumber;

    //==================================================================
    //Analysis code for a entry
    //==================================================================

    // Check GRL
    if (RunNumber == 349533){
      if (GRLlist(LumiBlock)){
        continue;
      }
    }

    FillProbeHist();
    FillInEffHist(3, 0);

    if ( tag_proc == NTagProc && probe_mesL1_pass->at(N4) > 0){
      FillMdtHist();
      FillSPHist();
      FillPtResidualHist();
    }


    ///*
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
        LOGD << "l2sa : pass/pt/eta/phi/sAddress/roadAlgo: " << probe_mesSA_pass->at(N50) << "/" << probe_mesSA_pt->at(N50) << "/" << probe_mesSA_eta->at(N50) << "/" << probe_mesSA_phi->at(N50) << "/" << probe_mesSA_sAddress->at(N50) << "/" << probe_mesSA_roadAlgo->at(N50);

        // Check TAP
        if (!(probe_mesEFTAG_pass -> at(N50) > -1 && probe_mesL1_pass -> at(N50) > 0 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
          continue;
        }
        // Set superpoint and segment for each station
        for ( int i = 0; i < probe_segment_n; i++){
          double R = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i]);
          double Z = probe_segment_z[i];
          LOGD << "segment R/Z: " << R << "/" << Z;

          //BI
          if (probe_segment_chamberIndex[i] == 0 || probe_segment_chamberIndex[i] == 1) {
            // normal
            LOGD << "size: " << (probe_mesSA_roadAw->at(N50)).size() << "/" << (probe_mesSA_roadFtkAw->at(N50)).size() << "/" << (probe_mesSA_roadRpcAw->at(N50)).size();
            h_ResidualSegment_eta_BI -> Fill(probe_eta, calc_residual(probe_mesSA_roadAw -> at(N50)[0], probe_mesSA_roadBw -> at(N50)[0], Z, R));
           if ((probe_mesSA_roadFtkAw->at(N50)).size()>0)  h_ResidualSegment_eta_BI_2d -> Fill(calc_residual(probe_mesSA_roadFtkAw -> at(N50)[0], probe_mesSA_roadFtkBw -> at(N50)[0], Z, R), calc_residual(probe_mesSA_roadRpcAw -> at(N50)[0], probe_mesSA_roadRpcBw -> at(N50)[0], Z, R));

            // RPC
            if (probe_mesSA_roadAlgo->at(N50) == 2){
              LOGD << "BI RPC road Aw/Bw: " << probe_mesSA_roadRpcAw->at(N50)[0] << "/" << probe_mesSA_roadRpcBw->at(N50)[0];
              LOGD << "BI RPC road residual: " << calc_residual(probe_mesSA_roadRpcAw -> at(N50)[0], probe_mesSA_roadRpcBw -> at(N50)[0], Z, R);
              h_ResidualSegment_eta_rpc_BI -> Fill(calc_residual(probe_mesSA_roadRpcAw -> at(N50)[0], probe_mesSA_roadRpcBw -> at(N50)[0], Z, R));
            }

            // FTK
            if (probe_mesSA_roadAlgo->at(N50) == 1){
              LOGD << "BI FTK road Aw/Bw: " << probe_mesSA_roadFtkAw->at(N50)[0] << "/" << probe_mesSA_roadFtkBw->at(N50)[0];
              LOGD << "BI FTK road residual: " << calc_residual(probe_mesSA_roadFtkAw -> at(N50)[0], probe_mesSA_roadFtkBw -> at(N50)[0], Z, R);
              h_ResidualSegment_eta_ftk_BI -> Fill(calc_residual(probe_mesSA_roadFtkAw -> at(N50)[0], probe_mesSA_roadFtkBw -> at(N50)[0], Z, R));
            }
            //cout << "BI: " << calc_residual(probe_mesSA_roadAw -> at(N50)[0], probe_mesSA_roadBw -> at(N50)[0], Z, R) << endl;

          } else if(probe_segment_chamberIndex[i] == 2 || probe_segment_chamberIndex[i] == 3) { //BM
            // normal
            h_ResidualSegment_eta_BM -> Fill(probe_eta, calc_residual(probe_mesSA_roadAw -> at(N50)[1], probe_mesSA_roadBw -> at(N50)[1], Z, R));

            // RPC
            if (probe_mesSA_roadAlgo->at(N50) == 2){
              LOGD << "BM RPC road Aw/Bw: " << probe_mesSA_roadRpcAw->at(N50)[1] << "/" << probe_mesSA_roadRpcBw->at(N50)[1];
              h_ResidualSegment_eta_rpc_BM -> Fill(calc_residual(probe_mesSA_roadRpcAw -> at(N50)[1], probe_mesSA_roadRpcBw -> at(N50)[1], Z, R));
            }

            // FTK
            if (probe_mesSA_roadAlgo->at(N50) == 1){
              LOGD << "BM FTK road Aw/Bw: " << probe_mesSA_roadFtkAw->at(N50)[1] << "/" << probe_mesSA_roadFtkBw->at(N50)[1];
              h_ResidualSegment_eta_ftk_BM -> Fill(calc_residual(probe_mesSA_roadFtkAw -> at(N50)[1], probe_mesSA_roadFtkBw -> at(N50)[1], Z, R));
            }
            //cout <<"BM: " <<  calc_residual(probe_mesSA_roadAw -> at(N50)[0], probe_mesSA_roadBw -> at(N50)[0], Z, R) << endl;

          } else if(probe_segment_chamberIndex[i] == 4 || probe_segment_chamberIndex[i] == 5) { //BO
            h_ResidualSegment_eta_BO -> Fill(probe_eta, calc_residual(probe_mesSA_roadAw -> at(N50)[2], probe_mesSA_roadBw -> at(N50)[2], Z, R));

            // RPC
            if (probe_mesSA_roadAlgo->at(N50) == 2){
              LOGD << "BO RPC road Aw/Bw: " << probe_mesSA_roadRpcAw->at(N50)[2] << "/" << probe_mesSA_roadRpcBw->at(N50)[2];
              h_ResidualSegment_eta_rpc_BO -> Fill(calc_residual(probe_mesSA_roadRpcAw -> at(N50)[2], probe_mesSA_roadRpcBw -> at(N50)[2], Z, R));
            }

            // FTK
            if (probe_mesSA_roadAlgo->at(N50) == 1){
              LOGD << "BO FTK road Aw/Bw: " << probe_mesSA_roadFtkAw->at(N50)[2] << "/" << probe_mesSA_roadFtkBw->at(N50)[2];
              h_ResidualSegment_eta_ftk_BO -> Fill(calc_residual(probe_mesSA_roadFtkAw -> at(N50)[2], probe_mesSA_roadFtkBw -> at(N50)[2], Z, R));
            }
            //cout <<"BO: " <<  calc_residual(probe_mesSA_roadAw -> at(N50)[0], probe_mesSA_roadBw -> at(N50)[0], Z, R) << endl;
          }
        }

        //// Check isRpcFailure
        //if (probe_mesSA_isRpcFailure -> at(N50) == 1){
        //  continue;
        //}

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
    //*/

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

  DrawPtResidualHist(c1, pdf);

  h_ResidualSegment_eta_BI -> Draw("colz");
  c1 -> Print(pdf, "pdf" );

  TH1D* tBI = h_ResidualSegment_eta_BI->ProjectionY("eta_BI",h_ResidualSegment_eta_BI->GetXaxis()->FindBin(-2),h_ResidualSegment_eta_BI->GetXaxis()->FindBin(2));
  tBI->Draw();
  c1 -> Print(pdf, "pdf" );

  h_ResidualSegment_eta_BM -> Draw("colz");
  c1 -> Print(pdf, "pdf" );

  TH1D* tBM = h_ResidualSegment_eta_BM->ProjectionY("eta_BM",h_ResidualSegment_eta_BM->GetXaxis()->FindBin(-2),h_ResidualSegment_eta_BM->GetXaxis()->FindBin(2));
  tBM->Draw();
  c1 -> Print(pdf, "pdf" );

  h_ResidualSegment_eta_BO -> Draw("colz");
  c1 -> Print(pdf, "pdf" );

  TH1D* tBO = h_ResidualSegment_eta_BO->ProjectionY("eta_BI",h_ResidualSegment_eta_BO->GetXaxis()->FindBin(-2),h_ResidualSegment_eta_BO->GetXaxis()->FindBin(2));
  tBO->Draw();
  c1 -> Print(pdf, "pdf" );

  h_ResidualSegment_eta_BI_2d -> Draw("colz");
  c1 -> Print(pdf, "pdf" );


  LOGD << "GetNbinsX(BI): " << h_ResidualSegment_eta_rpc_BI->GetNbinsX();
  LOGD << "0(BI): " << h_ResidualSegment_eta_rpc_BI->GetBinContent(0);
  LOGD << "1(BI): " << h_ResidualSegment_eta_rpc_BI->GetBinContent(1);
  LOGD << "GetNbinsX(BI): " << h_ResidualSegment_eta_ftk_BI->GetNbinsX();
  LOGD << "0(BI): " << h_ResidualSegment_eta_ftk_BI->GetBinContent(0);
  LOGD << "1(BI): " << h_ResidualSegment_eta_ftk_BI->GetBinContent(1);
  h_ResidualSegment_eta_rpc_BI->SetLineColor(kRed);
  //h_ResidualSegment_eta_rpc_BI->SetBinContent(1, h_ResidualSegment_eta_rpc_BI->GetBinContent(0)+h_ResidualSegment_eta_rpc_BI->GetBinContent(1));
  //h_ResidualSegment_eta_rpc_BI->SetBinContent(h_ResidualSegment_eta_rpc_BI->GetNbinsX(), h_ResidualSegment_eta_rpc_BI->GetBinContent(h_ResidualSegment_eta_rpc_BI->GetNbinsX())+h_ResidualSegment_eta_rpc_BI->GetBinContent(h_ResidualSegment_eta_rpc_BI->GetNbinsX()+1));
  h_ResidualSegment_eta_rpc_BI->Draw();
  h_ResidualSegment_eta_ftk_BI->SetLineColor(kBlue);
  //h_ResidualSegment_eta_ftk_BI->SetBinContent(1, h_ResidualSegment_eta_ftk_BI->GetBinContent(0)+h_ResidualSegment_eta_ftk_BI->GetBinContent(1));
  //h_ResidualSegment_eta_ftk_BI->SetBinContent(h_ResidualSegment_eta_ftk_BI->GetNbinsX(), h_ResidualSegment_eta_ftk_BI->GetBinContent(h_ResidualSegment_eta_ftk_BI->GetNbinsX())+h_ResidualSegment_eta_ftk_BI->GetBinContent(h_ResidualSegment_eta_ftk_BI->GetNbinsX()+1));
  h_ResidualSegment_eta_ftk_BI->Draw("same");
  TLegend *leg2 = new TLegend(0.00,0.00,0.3,0.15);
  leg2->AddEntry(h_ResidualSegment_eta_rpc_BI,"RPC road","f");
  leg2->AddEntry(h_ResidualSegment_eta_ftk_BI,"FTK road","f");
  leg2->Draw();
  c1 -> Print(pdf, "pdf" );
  LOGD << "0(BI): " << h_ResidualSegment_eta_rpc_BI->GetBinContent(0);
  LOGD << "1(BI): " << h_ResidualSegment_eta_rpc_BI->GetBinContent(1);
  LOGD << "0(BI): " << h_ResidualSegment_eta_ftk_BI->GetBinContent(0);
  LOGD << "1(BI): " << h_ResidualSegment_eta_ftk_BI->GetBinContent(1);

  h_ResidualSegment_eta_rpc_BM->SetLineColor(kRed);
  //h_ResidualSegment_eta_rpc_BM->SetBinContent(1, h_ResidualSegment_eta_rpc_BM->GetBinContent(0)+h_ResidualSegment_eta_rpc_BM->GetBinContent(1));
  h_ResidualSegment_eta_rpc_BM->Draw();
  h_ResidualSegment_eta_ftk_BM->SetLineColor(kBlue);
  //h_ResidualSegment_eta_ftk_BM->SetBinContent(1, h_ResidualSegment_eta_ftk_BM->GetBinContent(0)+h_ResidualSegment_eta_ftk_BM->GetBinContent(1));
  h_ResidualSegment_eta_ftk_BM->Draw("same");
  leg2->Draw();
  c1 -> Print(pdf, "pdf" );

  h_ResidualSegment_eta_rpc_BO->SetLineColor(kRed);
  //h_ResidualSegment_eta_rpc_BO->SetBinContent(1, h_ResidualSegment_eta_rpc_BO->GetBinContent(0)+h_ResidualSegment_eta_rpc_BO->GetBinContent(1));
  h_ResidualSegment_eta_rpc_BO->Draw();
  h_ResidualSegment_eta_ftk_BO->SetLineColor(kBlue);
  //h_ResidualSegment_eta_ftk_BO->SetBinContent(1, h_ResidualSegment_eta_ftk_BO->GetBinContent(0)+h_ResidualSegment_eta_ftk_BO->GetBinContent(1));
  h_ResidualSegment_eta_ftk_BO->Draw("same");
  leg2->Draw();
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



  //h_NumberOfSP_eta->Draw("colz");
  //c1 -> Print(pdf, "pdf" );

  //h_NumberOfSP_qeta->Draw("colz");
  //c1 -> Print(pdf, "pdf" );

  //h_NumberOfSP_pt_barrel->Draw("colz");
  //c1 -> Print(pdf, "pdf" );

  //h_NumberOfSP_pass_barrel->Draw("colz");
  //c1 -> Print(pdf, "pdf" );

  DrawFractionOfnSPs(h_NumberOfSP_LumiBlock, c1, pdf);
  DrawFractionOfnSPs(h_NumberOfSP_pt_barrel, c1, pdf);
  DrawFractionOfnSPs(h_NumberOfSP_eta, c1, pdf);
  DrawFractionOfnSPs(h_NumberOfSP_qeta, c1, pdf);


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

void RPC::FillInEffHist(int tap_type, int NTrigChain ){
  const double ZERO_LIMIT = 1e-5;

  LOGI << "FillInEffHist";


  // Check TAP
  LOGI << "tag_proc/tap_type: " << tag_proc << "/" << tap_type;
  if (tag_proc != tap_type){ return;}
  LOGI << "FillInEffHist0";
  if (!(probe_mesEFTAG_pass -> at(NTrigChain) > -1 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
    return;
  }

  LOGI << "FillInEffHist1";

  bool isBarrel = (abs(probe_eta) < 1.05);

  // Check L1 pass and SA not pass
  if (probe_mesL1_pass -> at(NTrigChain) <= -1) {return;}
  //if (probe_mesSA_pass -> at(NTrigChain) > 0) {return;}

  LOGI << "FillInEffHist2";

  // pT
  if(isBarrel) {
    if (probe_mesSA_pass -> at(NTrigChain) > 0) { // pass
      LOGI << "pT pass";
      h_InEff_pt -> Fill(probe_pt/1000., 1.);
    } else { // not pass
      if(abs(probe_mesSA_pt -> at(NTrigChain)) > ZERO_LIMIT){ // pT was calculated
        h_InEff_pt -> Fill(probe_pt/1000., -1.);
      } else { // p was not calculated
        h_InEff_pt -> Fill(probe_pt/1000., -2.);
      }
    }
  }

  // eta and qeta
  if (probe_pt/1000. > 4.){
    if (probe_mesSA_pass -> at(NTrigChain) > 0) { // pass
      h_InEff_eta -> Fill(probe_eta, 1.);
      h_InEff_qeta -> Fill(probe_eta*probe_charge, 1.);
    } else { // not pass
      if(abs(probe_mesSA_pt -> at(NTrigChain)) > ZERO_LIMIT){ // pT was calculated
        h_InEff_eta -> Fill(probe_eta, -1.);
        h_InEff_qeta -> Fill(probe_eta*probe_charge, -1.);
      } else { // pT was not calculated
        if (NumberOfSP() > 1.5){
          h_InEff_eta -> Fill(probe_eta, -2.);
          h_InEff_qeta -> Fill(probe_eta*probe_charge, -2.);
        } else if (NumberOfSP() < 1.5){
          h_InEff_eta -> Fill(probe_eta, -3.);
          h_InEff_qeta -> Fill(probe_eta*probe_charge, -3.);
        }
      }
    }
  }


  // Plateau cut
  //if (probe_pt/1000. > 8.0) {
  //  if(isBarrel) h_probe_mu_mu4_SA -> Fill(AverageInteractionsPerCrossing);
  //  h_probe_eta_mu4_SA -> Fill(probe_eta);
  //  if(isBarrel) h_probe_phi_mu4_SA -> Fill(probe_phi);
  //  hh_probe_etaphi_mu4_SA -> Fill(probe_eta, probe_phi);
  //}

  return;
}

void RPC::DrawInEffHist(TString pdf){
  //==================================================================
  //Set Canvas
  //==================================================================
  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1020, 700);
  c1->SetGrid();
  c1->SetRightMargin(0.20);
  c1->SetLeftMargin(0.23);
  c1->SetBottomMargin(0.20);

  c1 -> Print( pdf + "[", "pdf" );


  // pT
  h_InEff_pt->GetXaxis()->SetRangeUser(0,14);
  h_InEff_pt->GetYaxis()->SetRangeUser(-2,2);
  h_InEff_pt->Draw("colz");
  c1 -> Print( pdf, "pdf" );

  // eta
  h_InEff_eta->Draw("colz");
  h_InEff_eta->GetYaxis()->SetRangeUser(-4,4);
  c1 -> Print( pdf, "pdf" );

  // qeta
  h_InEff_qeta->Draw("colz");
  h_InEff_qeta->GetYaxis()->SetRangeUser(-4,4);
  c1 -> Print( pdf, "pdf" );

  // fraction of pt
  TH1D* h_pt_all            = h_InEff_pt->ProjectionX("_all");
  TH1D* h_pt_pass           = h_InEff_pt->ProjectionX("pass", h_InEff_pt->GetYaxis()->FindBin(0.9),  h_InEff_pt->GetYaxis()->FindBin(1.1));
  TH1D* h_pt_notpass        = h_InEff_pt->ProjectionX("notpass", h_InEff_pt->GetYaxis()->FindBin(-2.1),  h_InEff_pt->GetYaxis()->FindBin(-0.9));
  TH1D* h_pt_notpass_cal    = h_InEff_pt->ProjectionX("notpass_cal", h_InEff_pt->GetYaxis()->FindBin(-1.1),  h_InEff_pt->GetYaxis()->FindBin(-0.9));
  TH1D* h_pt_notpass_notcal = h_InEff_pt->ProjectionX("notpass_notcal", h_InEff_pt->GetYaxis()->FindBin(-2.1),  h_InEff_pt->GetYaxis()->FindBin(-1.9));

  h_pt_pass -> Draw();
  c1 -> Print( pdf, "pdf" );

  h_pt_notpass -> Draw();
  c1 -> Print( pdf, "pdf" );

  h_pt_notpass_cal -> Draw();
  c1 -> Print( pdf, "pdf" );

  h_pt_notpass_notcal -> Draw();
  c1 -> Print( pdf, "pdf" );

  TH1D *tmp_h_pt_pass           = (TH1D*)h_pt_pass -> Clone();
  TH1D *tmp_h_pt_notpass        = (TH1D*)h_pt_notpass -> Clone();
  TH1D *tmp_h_pt_notpass_cal    = (TH1D*)h_pt_notpass_cal -> Clone();
  TH1D *tmp_h_pt_notpass_notcal = (TH1D*)h_pt_notpass_notcal -> Clone();

  // divided by pass-hist and not-pass-hist
  h_pt_pass->Divide(h_pt_all);
  h_pt_notpass_cal->Divide(h_pt_all);
  h_pt_notpass_notcal->Divide(h_pt_all);

  THStack *hs_pt = new THStack("hs_pt",Form(";%s;Fraction",h_InEff_pt->GetXaxis()->GetTitle()));
  h_pt_pass->SetFillColor(kGreen);
  h_pt_notpass_cal->SetFillColor(kBlue);
  h_pt_notpass_notcal->SetFillColor(kRed);
  hs_pt->Add(h_pt_pass);
  hs_pt->Add(h_pt_notpass_cal);
  hs_pt->Add(h_pt_notpass_notcal);

  gStyle->SetHistTopMargin(0);
  hs_pt->Draw();
  hs_pt->SetMaximum(1.1);
  hs_pt->Draw();

  TLegend *leg_pt = new TLegend(0.00,0.00,0.57,0.15);
  leg_pt->AddEntry(h_pt_pass,"p_{T} was calculated and passed in L2MuonSA","f");
  leg_pt->AddEntry(h_pt_notpass_cal,"p_{T} was calculated and not passed in L2MuonSA","f");
  leg_pt->AddEntry(h_pt_notpass_notcal,"p_{T} was not calculated and not passed in L2MuonSA","f");
  leg_pt->Draw();

  c1 -> Print(pdf, "pdf" );

  // divided by not-pass-hist
  tmp_h_pt_notpass_cal->Divide(tmp_h_pt_notpass);
  tmp_h_pt_notpass_notcal->Divide(tmp_h_pt_notpass);

  THStack *tmp_hs_pt = new THStack("tmp_hs_pt",Form(";%s;Fraction",h_InEff_pt->GetXaxis()->GetTitle()));
  tmp_h_pt_notpass_cal->SetFillColor(kBlue);
  tmp_h_pt_notpass_notcal->SetFillColor(kRed);
  tmp_hs_pt->Add(tmp_h_pt_notpass_cal);
  tmp_hs_pt->Add(tmp_h_pt_notpass_notcal);

  gStyle->SetHistTopMargin(0);
  tmp_hs_pt->Draw();
  tmp_hs_pt->SetMaximum(1.1);
  tmp_hs_pt->Draw();

  TLegend *tmp_leg_pt = new TLegend(0.00,0.00,0.57,0.15);
  tmp_leg_pt->AddEntry(h_pt_notpass_cal,"p_{T} was calculated and not passed in L2MuonSA","f");
  tmp_leg_pt->AddEntry(h_pt_notpass_notcal,"p_{T} was not calculated and not passed in L2MuonSA","f");
  tmp_leg_pt->Draw();


  c1 -> Print(pdf, "pdf" );
  delete hs_pt;

  // fraction of eta
  TH1D* h_eta_all            = h_InEff_eta->ProjectionX("_all");
  TH1D* h_eta_pass           = h_InEff_eta->ProjectionX("pass",           h_InEff_eta->GetYaxis()->FindBin(0.9),  h_InEff_eta->GetYaxis()->FindBin(1.1));
  TH1D* h_eta_notpass        = h_InEff_eta->ProjectionX("notpass",        h_InEff_eta->GetYaxis()->FindBin(-3.1), h_InEff_eta->GetYaxis()->FindBin(-0.9));
  TH1D* h_eta_notpass_cal    = h_InEff_eta->ProjectionX("notpass_cal",    h_InEff_eta->GetYaxis()->FindBin(-1.1), h_InEff_eta->GetYaxis()->FindBin(-0.9));
  TH1D* h_eta_notpass_notcal_SP2 = h_InEff_eta->ProjectionX("notpass_notcal_SP2", h_InEff_eta->GetYaxis()->FindBin(-2.1), h_InEff_eta->GetYaxis()->FindBin(-1.9));
  TH1D* h_eta_notpass_notcal_SP1 = h_InEff_eta->ProjectionX("notpass_notcal_SP1", h_InEff_eta->GetYaxis()->FindBin(-3.1), h_InEff_eta->GetYaxis()->FindBin(-2.9));

  TH1D *tmp_h_eta_pass           = (TH1D*)h_eta_pass -> Clone();
  TH1D *tmp_h_eta_notpass        = (TH1D*)h_eta_notpass -> Clone();
  TH1D *tmp_h_eta_notpass_cal    = (TH1D*)h_eta_notpass_cal -> Clone();
  TH1D *tmp_h_eta_notpass_notcal_SP1 = (TH1D*)h_eta_notpass_notcal_SP1 -> Clone();
  TH1D *tmp_h_eta_notpass_notcal_SP2 = (TH1D*)h_eta_notpass_notcal_SP2 -> Clone();

  h_eta_pass->Divide(h_eta_all);
  h_eta_notpass_cal->Divide(h_eta_all);
  h_eta_notpass_notcal_SP1->Divide(h_eta_all);
  h_eta_notpass_notcal_SP2->Divide(h_eta_all);

  THStack *hs_eta = new THStack("hs_eta",Form(";%s;Fraction",h_InEff_eta->GetXaxis()->GetTitle()));
  h_eta_pass->SetFillColor(kGreen);
  h_eta_pass->SetFillStyle(3244);
  h_eta_notpass_cal->SetFillColor(kBlue);
  h_eta_notpass_notcal_SP1->SetFillColor(kRed);
  h_eta_notpass_notcal_SP2->SetFillColor(kGreen+2);
  hs_eta->Add(h_eta_pass);
  hs_eta->Add(h_eta_notpass_cal);
  hs_eta->Add(h_eta_notpass_notcal_SP1);
  hs_eta->Add(h_eta_notpass_notcal_SP2);

  gStyle->SetHistTopMargin(0);
  hs_eta->Draw();
  hs_eta->SetMaximum(1.1);
  hs_eta->Draw();

  leg_pt->Draw();

  TLegend *leg_eta = new TLegend(0.00,0.00,0.57,0.15);
  leg_eta->AddEntry(h_eta_pass,"Passed: pT > 0","f");
  leg_eta->AddEntry(h_eta_notpass_cal,"Not Passed: pT > 0","f");
  leg_eta->AddEntry(h_eta_notpass_notcal_SP1,"Not Passed: pT = 0, nSP = 0 or 1","f");
  leg_eta->AddEntry(h_eta_notpass_notcal_SP2,"Not Passed: pT = 0, nSP > 1","f");
  leg_eta->Draw();

  c1 -> Print(pdf, "pdf" );

  // divided by not-pass-hist
  tmp_h_eta_notpass_cal->Divide(tmp_h_eta_notpass);
  tmp_h_eta_notpass_notcal_SP1->Divide(tmp_h_eta_notpass);
  tmp_h_eta_notpass_notcal_SP2->Divide(tmp_h_eta_notpass);

  THStack *tmp_hs_eta = new THStack("tmp_hs_eta",Form(";%s;Fraction",h_InEff_eta->GetXaxis()->GetTitle()));
  tmp_h_eta_notpass_cal->SetFillColor(kBlue);
  tmp_h_eta_notpass_notcal_SP1->SetFillColor(kRed);
  tmp_h_eta_notpass_notcal_SP2->SetFillColor(kGreen+2);
  tmp_hs_eta->Add(tmp_h_eta_notpass_cal);
  tmp_hs_eta->Add(tmp_h_eta_notpass_notcal_SP1);
  tmp_hs_eta->Add(tmp_h_eta_notpass_notcal_SP2);

  gStyle->SetHistTopMargin(0);
  tmp_hs_eta->Draw();
  tmp_hs_eta->SetMaximum(1.1);
  tmp_hs_eta->Draw();


  TLegend *tmp_leg_eta = new TLegend(0.00,0.00,0.57,0.15);
  tmp_leg_eta->AddEntry(h_eta_notpass_cal,"Not Passed: pT > 0","f");
  tmp_leg_eta->AddEntry(h_eta_notpass_notcal_SP1,"Not Passed: pT = 0, nSP = 0 or 1","f");
  tmp_leg_eta->AddEntry(h_eta_notpass_notcal_SP2,"Not Passed: pT = 0, nSP > 1","f");
  tmp_leg_eta->Draw();

  tmp_leg_eta->Draw();

  c1 -> Print(pdf, "pdf" );
  delete hs_eta;

  // fraction of qeta
  TH1D* h_qeta_all            = h_InEff_qeta->ProjectionX("_all");
  TH1D* h_qeta_pass           = h_InEff_qeta->ProjectionX("pass",           h_InEff_qeta->GetYaxis()->FindBin(0.9),  h_InEff_qeta->GetYaxis()->FindBin(1.1));
  TH1D* h_qeta_notpass        = h_InEff_qeta->ProjectionX("notpass",        h_InEff_qeta->GetYaxis()->FindBin(-3.1), h_InEff_qeta->GetYaxis()->FindBin(-0.9));
  TH1D* h_qeta_notpass_cal    = h_InEff_qeta->ProjectionX("notpass_cal",    h_InEff_qeta->GetYaxis()->FindBin(-1.1), h_InEff_qeta->GetYaxis()->FindBin(-0.9));
  TH1D* h_qeta_notpass_notcal_SP2 = h_InEff_qeta->ProjectionX("notpass_notcal_SP2", h_InEff_qeta->GetYaxis()->FindBin(-2.1), h_InEff_qeta->GetYaxis()->FindBin(-1.9));
  TH1D* h_qeta_notpass_notcal_SP1 = h_InEff_qeta->ProjectionX("notpass_notcal_SP1", h_InEff_qeta->GetYaxis()->FindBin(-3.1), h_InEff_qeta->GetYaxis()->FindBin(-2.9));

  TH1D *tmp_h_qeta_pass           = (TH1D*)h_qeta_pass -> Clone();
  TH1D *tmp_h_qeta_notpass        = (TH1D*)h_qeta_notpass -> Clone();
  TH1D *tmp_h_qeta_notpass_cal    = (TH1D*)h_qeta_notpass_cal -> Clone();
  TH1D *tmp_h_qeta_notpass_notcal_SP1 = (TH1D*)h_qeta_notpass_notcal_SP1 -> Clone();
  TH1D *tmp_h_qeta_notpass_notcal_SP2 = (TH1D*)h_qeta_notpass_notcal_SP2 -> Clone();

  h_qeta_pass->Divide(h_qeta_all);
  h_qeta_notpass_cal->Divide(h_qeta_all);
  h_qeta_notpass_notcal_SP1->Divide(h_qeta_all);
  h_qeta_notpass_notcal_SP2->Divide(h_qeta_all);

  THStack *hs_qeta = new THStack("hs_qeta",Form(";%s;Fraction",h_InEff_qeta->GetXaxis()->GetTitle()));
  h_qeta_pass->SetFillStyle(3244);
  h_qeta_pass->SetFillColor(kGreen);
  h_qeta_notpass_cal->SetFillColor(kBlue);
  h_qeta_notpass_notcal_SP1->SetFillColor(kRed);
  h_qeta_notpass_notcal_SP2->SetFillColor(kGreen+2);
  hs_qeta->Add(h_qeta_pass);
  hs_qeta->Add(h_qeta_notpass_cal);
  hs_qeta->Add(h_qeta_notpass_notcal_SP1);
  hs_qeta->Add(h_qeta_notpass_notcal_SP2);

  gStyle->SetHistTopMargin(0);
  hs_qeta->Draw();
  hs_qeta->SetMaximum(1.1);
  hs_qeta->Draw();


  TLegend *leg_qeta = new TLegend(0.00,0.00,0.57,0.15);
  leg_qeta->AddEntry(h_qeta_pass,"Passed: pT > 0","f");
  leg_qeta->AddEntry(h_qeta_notpass_cal,"Not Passed: pT > 0","f");
  leg_qeta->AddEntry(h_qeta_notpass_notcal_SP1,"Not Passed: pT = 0, nSP = 0 or 1","f");
  leg_qeta->AddEntry(h_qeta_notpass_notcal_SP2,"Not Passed: pT = 0, nSP > 1","f");
  leg_qeta->Draw();


  c1 -> Print(pdf, "pdf" );

  // divided by not-pass-hist
  tmp_h_qeta_notpass_cal->Divide(tmp_h_qeta_notpass);
  tmp_h_qeta_notpass_notcal_SP1->Divide(tmp_h_qeta_notpass);
  tmp_h_qeta_notpass_notcal_SP2->Divide(tmp_h_qeta_notpass);

  THStack *tmp_hs_qeta = new THStack("tmp_hs_qeta",Form(";%s;Fraction",h_InEff_qeta->GetXaxis()->GetTitle()));
  tmp_h_qeta_notpass_cal->SetFillColor(kBlue);
  tmp_h_qeta_notpass_notcal_SP1->SetFillColor(kRed);
  tmp_h_qeta_notpass_notcal_SP2->SetFillColor(kGreen+2);
  tmp_hs_qeta->Add(tmp_h_qeta_notpass_cal);
  tmp_hs_qeta->Add(tmp_h_qeta_notpass_notcal_SP1);
  tmp_hs_qeta->Add(tmp_h_qeta_notpass_notcal_SP2);

  gStyle->SetHistTopMargin(0);
  tmp_hs_qeta->Draw();
  tmp_hs_qeta->SetMaximum(1.1);
  tmp_hs_qeta->Draw();


  TLegend *tmp_leg_qeta = new TLegend(0.00,0.00,0.57,0.15);
  tmp_leg_qeta->AddEntry(h_qeta_notpass_cal,"Not Passed: pT > 0","f");
  tmp_leg_qeta->AddEntry(h_qeta_notpass_notcal_SP1,"Not Passed: pT = 0, nSP = 0 or 1","f");
  tmp_leg_qeta->AddEntry(h_qeta_notpass_notcal_SP2,"Not Passed: pT = 0, nSP > 1","f");
  tmp_leg_qeta->Draw();

  tmp_leg_qeta->Draw();

  c1 -> Print(pdf, "pdf" );
  delete hs_qeta;


  /*
  // fraction of qeta
  TH1D* h_qeta_all = h_InEff_qeta->ProjectionX("_all");
  TH1D* h_qeta_pass = h_InEff_qeta->ProjectionX("pass", h_InEff_qeta->GetYaxis()->FindBin(0.9),  h_InEff_qeta->GetYaxis()->FindBin(1.1));
  TH1D* h_qeta_notpass = h_InEff_qeta->ProjectionX("notpass", h_InEff_qeta->GetYaxis()->FindBin(-2.1),  h_InEff_qeta->GetYaxis()->FindBin(-0.9));
  TH1D* h_qeta_notpass_cal = h_InEff_qeta->ProjectionX("notpass_cal", h_InEff_qeta->GetYaxis()->FindBin(-1.1),  h_InEff_qeta->GetYaxis()->FindBin(-0.9));
  TH1D* h_qeta_notpass_notcal = h_InEff_qeta->ProjectionX("notpass_notcal", h_InEff_qeta->GetYaxis()->FindBin(-2.1),  h_InEff_qeta->GetYaxis()->FindBin(-1.9));

  //h_qeta_pass -> Draw();
  //c1 -> Print( pdf, "pdf" );

  //h_qeta_notpass -> Draw();
  //c1 -> Print( pdf, "pdf" );
  //leg_pt->Draw();

  //h_qeta_notpass_cal -> Draw();
  //c1 -> Print( pdf, "pdf" );

  //h_qeta_notpass_notcal -> Draw();
  //c1 -> Print( pdf, "pdf" );

  TH1D *tmp_h_qeta_pass           = (TH1D*)h_qeta_pass -> Clone();
  TH1D *tmp_h_qeta_notpass        = (TH1D*)h_qeta_notpass -> Clone();
  TH1D *tmp_h_qeta_notpass_cal    = (TH1D*)h_qeta_notpass_cal -> Clone();
  TH1D *tmp_h_qeta_notpass_notcal = (TH1D*)h_qeta_notpass_notcal -> Clone();

  h_qeta_pass->Divide(h_qeta_all);
  h_qeta_notpass_cal->Divide(h_qeta_all);
  h_qeta_notpass_notcal->Divide(h_qeta_all);

  THStack *hs_qeta = new THStack("hs_qeta",Form(";%s;Fraction",h_InEff_qeta->GetXaxis()->GetTitle()));
  h_qeta_pass->SetFillColor(kGreen);
  h_qeta_notpass_cal->SetFillColor(kBlue);
  h_qeta_notpass_notcal->SetFillColor(kRed);
  hs_qeta->Add(h_qeta_pass);
  hs_qeta->Add(h_qeta_notpass_cal);
  hs_qeta->Add(h_qeta_notpass_notcal);

  gStyle->SetHistTopMargin(0);
  hs_qeta->Draw();
  hs_qeta->SetMaximum(1.1);
  hs_qeta->Draw();

  leg_pt->Draw();

  c1 -> Print(pdf, "pdf" );

  // divided by not-pass-hist
  tmp_h_qeta_notpass_cal->Divide(tmp_h_qeta_notpass);
  tmp_h_qeta_notpass_notcal->Divide(tmp_h_qeta_notpass);

  THStack *tmp_hs_qeta = new THStack("tmp_hs_qeta",Form(";%s;Fraction",h_InEff_qeta->GetXaxis()->GetTitle()));
  tmp_h_qeta_notpass_cal->SetFillColor(kBlue);
  tmp_h_qeta_notpass_notcal->SetFillColor(kRed);
  tmp_hs_qeta->Add(tmp_h_qeta_notpass_cal);
  tmp_hs_qeta->Add(tmp_h_qeta_notpass_notcal);

  gStyle->SetHistTopMargin(0);
  tmp_hs_qeta->Draw();
  tmp_hs_qeta->SetMaximum(1.1);
  tmp_hs_qeta->Draw();

  tmp_leg_pt->Draw();

  c1 -> Print(pdf, "pdf" );
  delete hs_qeta;
  */

  delete leg_pt;
  c1 -> Print( pdf + "]", "pdf" );
}

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
  int trig = N50;
  if (abs(probe_mesSA_superPointZ_BI -> at(trig)) > 0. && probe_mesSA_superPointZ_BI -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BM -> at(trig)) > 0. && probe_mesSA_superPointZ_BM -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BO -> at(trig)) > 0. && probe_mesSA_superPointZ_BO -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EI -> at(trig)) > 0. && probe_mesSA_superPointZ_EI -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EM -> at(trig)) > 0. && probe_mesSA_superPointZ_EM -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EO -> at(trig)) > 0. && probe_mesSA_superPointZ_EO -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EE -> at(trig)) > 0. && probe_mesSA_superPointZ_EE -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_CSC -> at(trig)) > 0. && probe_mesSA_superPointZ_CSC -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BEE -> at(trig)) > 0. && probe_mesSA_superPointZ_BEE -> at(trig) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BME -> at(trig)) > 0. && probe_mesSA_superPointZ_BME -> at(trig) > -99990.) {
    number += 1;
  }
  return number;
}

int RPC::NumberOfSP(int NTrigChain){
  int number = 0;
  if (abs(probe_mesSA_superPointZ_BI -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_BI -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BM -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_BM -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BO -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_BO -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EI -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_EI -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EM -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_EM -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EO -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_EO -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_EE -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_EE -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_CSC -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_CSC -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BEE -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_BEE -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  if (abs(probe_mesSA_superPointZ_BME -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_BME -> at(NTrigChain) > -99990.) {
    number += 1;
  }
  return number;
}

void RPC::Display(int tap_type, int trig_chain, Long64_t begin_entry, Long64_t limit_entry, TString pdf)
{
  // Prepare Canvas
  TCanvas *c2 = new TCanvas("c2", "c2", 10, 10, 1020, 700);
  c2->SetGrid();
  c2->SetRightMargin(0.20);
  c2->SetLeftMargin(0.23);
  c2->SetBottomMargin(0.20);
  c2->Print(pdf + "[", "pdf");

  c2->Print(pdf, "pdf");

  const double ZERO_LIMIT = 1e-5;

  if (limit_entry == -1) {
    begin_entry = 0;
    limit_entry = fChain -> GetEntries();
    LOGI << "-------";
  }
  // Prepare Loop
  if (fChain == 0) return;
  int nLoop = fChain -> GetEntries();
  //Long64_t nentries = fChain->GetEntriesFast();
  double entries = fChain->GetEntries();
  LOGI << "Nentries:" << entries;
  Long64_t nbytes = 0, nb = 0;

  int current_entry = 0;
  // Start Loop
  for (Long64_t jentry = begin_entry; jentry<nLoop;jentry++) {
    int  ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    LOGI << "now event number -->> " << ientry << "/" << nLoop << " : " << ((double) ientry)/nLoop * 100. << "%";
    if (ientry < 0) break;
    // Check GRL
    if (RunNumber == 349533){
      if (GRLlist(LumiBlock)){
        continue;
      }
    }

    if (tag_proc != tap_type) {
      continue;
    }

    int NTrigChain = trig_chain;

    if (!(probe_mesEFTAG_pass -> at(NTrigChain) > -1 && probe_mesL1_pass -> at(NTrigChain) > 0 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
      continue;
    }

    //if ( probe_mesSA_pass -> at(NTrigChain) > 0 ) {
    //  continue;
    //}

    // Check sAddress to avoid Special sector
    //if (probe_mesSA_sAddress -> at(NTrigChain) == 2 || probe_mesSA_sAddress -> at(NTrigChain) == 3 || probe_mesSA_sAddress -> at(NTrigChain) == -1){
    //  continue;
    //}

    //if ( EventNumber != 162940171){
    //  continue;
    //}

    //offline pT cut
    //if ( !((probe_pt/1000. > 4.0) && (probe_pt/1000. < 6.0)) ){
    //  continue;
    //}

    //offline pT cut
    //if ( !((probe_pt/1000. > 10.0)) ){
    //  continue;
    //}

    // L2MuonSA pT cut
    //if ( abs(probe_mesSA_pt->at(NTrigChain)) > ZERO_LIMIT ){
    //  continue;
    //}

    // Check Barrel
    if (abs(probe_eta) > 1.05){
      continue;
    }

    //if (!(abs(probe_mesSA_superPointR_BM -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_EI -> at(NTrigChain) > 0.)){
    // continue;
    //}

    //double qeta = probe_eta*probe_charge;
    //if ( !(qeta > -1.3 && qeta < -0.9)){
    //  continue;
    //}


    // Check NumberOfSP
    //if (NumberOfSP() != 3){
    //  continue;
    //}
    // Check size of mdtHits
    //if ((probe_mesSA_mdtHitChamber -> at(NTrigChain)).size() == 0){
    //  continue;
    //}
    // Check sAddress to avoid Spercial sector
    //if (probe_mesSA_sAddress -> at(N50) == 1 || probe_mesSA_sAddress -> at(N50) == 3 || probe_mesSA_sAddress -> at(N50) == 0){
    //  continue;
    //}


    enum Station{ Inner, Middle, Outer};
    enum Sector{ Large, Small, LargeSpecial, SmallSpecial};
    double LargeR[3]        = {4.95, 7.,   9.5};
    double SmallR[3]        = {5.4,  7.,   9.5};
    double LargeSpecialR[3] = {4.55, 7.95, 10.45};
    double SmallSpecialR[3] = {4.55, 7.95, 10.45};
    double deltaR = 0.3;
    double deltaZ = 0.3;

    LOGD << "EventNumber = " << EventNumber;

    LOGD << "NumberOfSP(): " << NumberOfSP(NTrigChain);

    TString label_for_sector;
    if (probe_mesSA_sAddress->at(NTrigChain) == 0){
      label_for_sector = "Large";
    } else if (probe_mesSA_sAddress->at(NTrigChain) == 1){
      label_for_sector = "Small";
    } else if (probe_mesSA_sAddress->at(NTrigChain) == 2){
      label_for_sector = "Large Special";
    } else if (probe_mesSA_sAddress->at(NTrigChain) == 3){
      label_for_sector = "Small Special";
    } else if (probe_mesSA_sAddress->at(NTrigChain) == -1){
      label_for_sector = "Endcap";
    }


    LOGD << "Sector: " << label_for_sector;
    LOGD << probe_mesSA_sAddress->at(NTrigChain);
    //cout << "Rmin_BI: " << Rmin_BI << ", Rmin_BM: " << Rmin_BM << ", Rmin_BO: " << Rmin_BO << endl;
    //cout << "Rmax_BI: " << Rmax_BI << ", Rmax_BM: " << Rmax_BM << ", Rmax_BO: " << Rmax_BO << endl;
    //cout << "Zmin_BI: " << Zmin_BI << ", Zmin_BM: " << Zmin_BM << ", Zmin_BO: " << Zmin_BO << endl;
    //cout << "Zmax_BI: " << Zmax_BI << ", Zmax_BM: " << Zmax_BM << ", Zmax_BO: " << Zmax_BO << endl;

    //Offline segment for all station
    const Int_t nOffSeg = probe_segment_n;
    TGraph gr_segment = TGraph(0);
    TGraph gr_segment_EI = TGraph(0);
    for ( int i = 0; i < nOffSeg; i++){
      double R = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i]) / 1000.;
      double Z = probe_segment_z[i] / 1000.;
      gr_segment.SetPoint(i, Z, R);
      if ( probe_segment_chamberIndex[i] == 7 || probe_segment_chamberIndex[i] == 8){
        gr_segment_EI.SetPoint(i, Z, R);
      }
    }

    // RPC hits
    const Int_t nRPC = (probe_mesSA_rpcHitR -> at(NTrigChain)).size();
    TGraph gr_RPC = TGraph(nRPC);
    for (int i = 0; i< nRPC; i++){
      if (probe_mesSA_rpcHitStationNumber->at(NTrigChain)[i] > 0.) {
        gr_RPC.SetPoint(i, (probe_mesSA_rpcHitZ -> at(NTrigChain)[i]) / 1000., (probe_mesSA_rpcHitR -> at(NTrigChain)[i]) / 1000.);
      }
    }

    // SuperPoint
    TGraph gr_SP = TGraph(NumberOfSP(NTrigChain));
    TGraph gr_SP_EI = TGraph(NumberOfSP(NTrigChain));
    TGraph gr_SP_BI = TGraph(NumberOfSP(NTrigChain));
    int iSP = 0;
    if (abs(probe_mesSA_superPointR_BI -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_BI -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_BI -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BI -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_BM -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_BM -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_BM -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BM -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_BO -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_BO -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_BO -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BO -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_EI -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_EI -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_EI -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_EI -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_EM -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_EM -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_EM -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_EM -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_EO -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_EO -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_EO -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_EO -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_EE -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_EE -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_EE -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_EE -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_CSC -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_CSC -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_CSC -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_CSC -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_BEE -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_BEE -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_BEE -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BEE -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    if (abs(probe_mesSA_superPointR_BME -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_BME -> at(NTrigChain) > -99990.){
      gr_SP.SetPoint(iSP , (probe_mesSA_superPointZ_BME -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BME -> at(NTrigChain)) / 1000.);
      iSP += 1;
    }
    // special EI
    if (abs(probe_mesSA_superPointR_EI -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_EI -> at(NTrigChain) > -99990.){
      gr_SP_EI.SetPoint(0 , (probe_mesSA_superPointZ_EI -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_EI -> at(NTrigChain)) / 1000.);
    }
    // special BI
    if (abs(probe_mesSA_superPointR_BI -> at(NTrigChain)) > 0. && probe_mesSA_superPointR_BI -> at(NTrigChain) > -99990.){
      gr_SP_BI.SetPoint(0 , (probe_mesSA_superPointZ_BI -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BI -> at(NTrigChain)) / 1000.);
    }

    // MDT HITs
    int nMDT_Inlier=0;
    int nMDT_Outlier=0;
    int nMDT_BI_Inlier=0;
    int nMDT_BM_Inlier=0;
    int nMDT_BO_Inlier=0;
    int nMDT_BI_Outlier=0;
    int nMDT_BM_Outlier=0;
    int nMDT_BO_Outlier=0;

    int nMDT = (probe_mesSA_mdtHitChamber -> at(NTrigChain)).size();
    LOGD << "nMDT: " << nMDT;

    for (int iMDT = 0; iMDT < nMDT; iMDT++){
      // All station
      if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 0) {
        nMDT_Inlier += 1;
      } else if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 1) {
        nMDT_Outlier += 1;
      }
      // BI
      if (probe_mesSA_mdtHitChamber -> at(NTrigChain)[iMDT] == 0){
        if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 0) {
          nMDT_BI_Inlier += 1;
        } else if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 1) {
          nMDT_BI_Outlier += 1;
        }
      }
      // BM
      if (probe_mesSA_mdtHitChamber -> at(NTrigChain)[iMDT] == 1){
        if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 0) {
          nMDT_BM_Inlier += 1;
        } else if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 1) {
          nMDT_BM_Outlier += 1;
        }
      }
      // BO
      if (probe_mesSA_mdtHitChamber -> at(NTrigChain)[iMDT] == 2){
        if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 0) {
          nMDT_BO_Inlier += 1;
        } else if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 1) {
          nMDT_BO_Outlier += 1;
        }
      }
    }
    TGraph gr_MdtHit_Inlier = TGraph(nMDT_Inlier);
    TGraph gr_MdtHit_Outlier_1 = TGraph(nMDT_Outlier);
    TGraph gr_MdtHit_Outlier_2 = TGraph(nMDT_Outlier);
    TGraph gr_MdtHit_Outlier_3 = TGraph(nMDT_Outlier);

    for (int iMDT = 0; iMDT < nMDT; iMDT++){
      if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 0) {
        gr_MdtHit_Inlier.SetPoint(iMDT, (probe_mesSA_mdtHitZ -> at(NTrigChain)[iMDT]) / 1000., (probe_mesSA_mdtHitR -> at(NTrigChain)[iMDT] / 1000.));
      } else if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 1){
        gr_MdtHit_Outlier_1.SetPoint(iMDT, (probe_mesSA_mdtHitZ -> at(NTrigChain)[iMDT] / 1000.), (probe_mesSA_mdtHitR -> at(NTrigChain)[iMDT] / 1000.));
      } else if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 2){
        gr_MdtHit_Outlier_2.SetPoint(iMDT, (probe_mesSA_mdtHitZ -> at(NTrigChain)[iMDT] / 1000.), (probe_mesSA_mdtHitR -> at(NTrigChain)[iMDT] / 1000.));
      } else if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[iMDT] == 3){
        gr_MdtHit_Outlier_3.SetPoint(iMDT, (probe_mesSA_mdtHitZ -> at(NTrigChain)[iMDT] / 1000.), (probe_mesSA_mdtHitR -> at(NTrigChain)[iMDT] / 1000.));
      }
    }

    LOGD << "DEBUG2: " << nOffSeg << ": " << nRPC << ": " << nMDT;

    // road
    LOGD << "road: " << probe_mesSA_roadAlgo -> at(NTrigChain) << ": " << probe_mesSA_roadAlgo->at(NTrigChain);
    LOGD << "road: " << probe_mesSA_roadAw -> at(NTrigChain)[0] << ": " << probe_mesSA_roadBw->at(NTrigChain)[0];
    LOGD << "road: " << probe_mesSA_roadAw -> at(NTrigChain)[1] << ": " << probe_mesSA_roadBw->at(NTrigChain)[1];
    LOGD << "road: " << probe_mesSA_roadAw -> at(NTrigChain)[2] << ": " << probe_mesSA_roadBw->at(NTrigChain)[2];
    LOGD << "roadFtk: " << probe_mesSA_roadFtkAw -> at(NTrigChain)[0] << ": " << probe_mesSA_roadFtkBw->at(NTrigChain)[0];
    LOGD << "roadFtk: " << probe_mesSA_roadFtkAw -> at(NTrigChain)[1] << ": " << probe_mesSA_roadFtkBw->at(NTrigChain)[1];
    LOGD << "roadFtk: " << probe_mesSA_roadFtkAw -> at(NTrigChain)[2] << ": " << probe_mesSA_roadFtkBw->at(NTrigChain)[2];
    LOGD << "roadRpc: " << probe_mesSA_roadRpcAw -> at(NTrigChain)[0] << ": " << probe_mesSA_roadRpcBw->at(NTrigChain)[0];
    LOGD << "roadRpc: " << probe_mesSA_roadRpcAw -> at(NTrigChain)[1] << ": " << probe_mesSA_roadRpcBw->at(NTrigChain)[1];
    LOGD << "roadRpc: " << probe_mesSA_roadRpcAw -> at(NTrigChain)[2] << ": " << probe_mesSA_roadRpcBw->at(NTrigChain)[2];

    TF1 f_roadRpc_BI = TF1("f_roadRpc_BI", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BI.SetTitle(";Z [m];R [m]");
    f_roadRpc_BI.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[0]);
    f_roadRpc_BI.SetParameter(1,probe_mesSA_roadRpcBw -> at(NTrigChain)[0]/1000.);
    f_roadRpc_BI.SetLineColor(kCyan);
    f_roadRpc_BI.SetLineWidth(2);
    f_roadRpc_BI.SetLineStyle(2);

    //double aw_rpc_BI = probe_mesSA_roadRpcAw->at(NTrigChain)[0];
    //double bw_rpc_BI = probe_mesSA_roadRpcBw->at(NTrigChain)[0]/1000.;
    //double x1_rpc_BI = (4 - bw_rpc_BI)/aw_rpc_BI;
    //double x2_rpc_BI = (6 - bw_rpc_BI)/aw_rpc_BI;
    //double xMin_rpc_BI = (x1_rpc_BI<x2_rpc_BI)? x1_rpc_BI:x2_rpc_BI;
    //double xMax_rpc_BI = (x1_rpc_BI<x2_rpc_BI)? x2_rpc_BI:x1_rpc_BI;
    //cout << "xMin,Max: " << xMin_rpc_BI << ": " << xMax_rpc_BI << endl;
    //f_roadRpc_BI.SetRange(xMin_rpc_BI,xMax_rpc_BI,4,6);

    TF1 f_roadRpc_BI_plus = TF1("f_roadRpc_BI_plus", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BI_plus.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[0]);
    f_roadRpc_BI_plus.SetParameter(1,( (probe_mesSA_roadRpcBw -> at(NTrigChain)[0]/1000.) + rWidthToBw(probe_mesSA_roadRpcAw->at(NTrigChain)[0],0.4)));
    f_roadRpc_BI_plus.SetLineColor(kCyan+0);
    f_roadRpc_BI_plus.SetLineWidth(2);
    f_roadRpc_BI_plus.SetLineStyle(1);

    TF1 f_roadRpc_BI_minus = TF1("f_roadRpc_BI_minus", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BI_minus.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[0]);
    f_roadRpc_BI_minus.SetParameter(1,( (probe_mesSA_roadRpcBw -> at(NTrigChain)[0]/1000.) - rWidthToBw(probe_mesSA_roadRpcAw->at(NTrigChain)[0],0.4)));
    f_roadRpc_BI_minus.SetLineColor(kCyan+0);
    f_roadRpc_BI_minus.SetLineWidth(2);
    f_roadRpc_BI_minus.SetLineStyle(1);

    TF1 f_roadFtk_BI = TF1("f_roadFtk_BI", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BI.SetTitle(";Z [m];R [m]");
    f_roadFtk_BI.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[0]);
    f_roadFtk_BI.SetParameter(1,probe_mesSA_roadFtkBw -> at(NTrigChain)[0]/1000.);
    f_roadFtk_BI.SetLineColor(kMagenta);
    f_roadFtk_BI.SetLineWidth(2);
    f_roadFtk_BI.SetLineStyle(2);

    //double aw_ftk_BI = probe_mesSA_roadFtkAw->at(NTrigChain)[0];
    //double bw_ftk_BI = probe_mesSA_roadFtkBw->at(NTrigChain)[0]/1000.;
    //double x1_ftk_BI = (4 - bw_ftk_BI)/aw_ftk_BI;
    //double x2_ftk_BI = (6 - bw_ftk_BI)/aw_ftk_BI;
    //double xMin_ftk_BI = (x1_ftk_BI<x2_ftk_BI)? x1_ftk_BI:x2_ftk_BI;
    //double xMax_ftk_BI = (x1_ftk_BI<x2_ftk_BI)? x2_ftk_BI:x1_ftk_BI;
    //cout << "xMin,Max: " << xMin_ftk_BI << ": " << xMax_ftk_BI << endl;
    //f_roadFtk_BI.SetRange(xMin_ftk_BI,xMax_ftk_BI,4,6);

    TF1 f_roadFtk_BI_plus = TF1("f_roadFtk_BI_plus", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BI_plus.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[0]);
    f_roadFtk_BI_plus.SetParameter(1,( (probe_mesSA_roadFtkBw -> at(NTrigChain)[0]/1000.) + rWidthToBw(probe_mesSA_roadFtkAw->at(NTrigChain)[0],0.4)));
    f_roadFtk_BI_plus.SetLineColor(kMagenta+0);
    f_roadFtk_BI_plus.SetLineWidth(2);
    f_roadFtk_BI_plus.SetLineStyle(1);

    TF1 f_roadFtk_BI_minus = TF1("f_roadFtk_BI_minus", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BI_minus.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[0]);
    f_roadFtk_BI_minus.SetParameter(1,( (probe_mesSA_roadFtkBw -> at(NTrigChain)[0]/1000.) - rWidthToBw(probe_mesSA_roadFtkAw->at(NTrigChain)[0],0.4)));
    f_roadFtk_BI_minus.SetLineColor(kMagenta+0);
    f_roadFtk_BI_minus.SetLineWidth(2);
    f_roadFtk_BI_minus.SetLineStyle(1);


    TF1 f_roadRpc_BM = TF1("f_roadRpc_BM", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BM.SetTitle(";Z [m];R [m]");
    f_roadRpc_BM.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[1]);
    f_roadRpc_BM.SetParameter(1,probe_mesSA_roadRpcBw -> at(NTrigChain)[1]/1000.);
    f_roadRpc_BM.SetLineColor(kCyan+1);
    f_roadRpc_BM.SetLineWidth(2);
    f_roadRpc_BM.SetLineStyle(2);

    //double aw_rpc_BM = probe_mesSA_roadRpcAw->at(NTrigChain)[1];
    //double bw_rpc_BM = probe_mesSA_roadRpcBw->at(NTrigChain)[1]/1000.;
    //double x1_rpc_BM = (6 - bw_rpc_BM)/aw_rpc_BM;
    //double x2_rpc_BM = (8 - bw_rpc_BM)/aw_rpc_BM;
    //double xMin_rpc_BM = (x1_rpc_BM<x2_rpc_BM)? x1_rpc_BM:x2_rpc_BM;
    //double xMax_rpc_BM = (x1_rpc_BM<x2_rpc_BM)? x2_rpc_BM:x1_rpc_BM;
    //f_roadRpc_BM.SetRange(xMin_rpc_BM,xMax_rpc_BM,6,8);

    TF1 f_roadRpc_BM_plus = TF1("f_roadRpc_BM_plus", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BM_plus.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[1]);
    f_roadRpc_BM_plus.SetParameter(1,( (probe_mesSA_roadRpcBw -> at(NTrigChain)[1]/1000.) + rWidthToBw(probe_mesSA_roadRpcAw->at(NTrigChain)[1],0.2)));
    f_roadRpc_BM_plus.SetLineColor(kCyan+1);
    f_roadRpc_BM_plus.SetLineWidth(2);
    f_roadRpc_BM_plus.SetLineStyle(1);

    TF1 f_roadRpc_BM_minus = TF1("f_roadRpc_BM_minus", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BM_minus.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[1]);
    f_roadRpc_BM_minus.SetParameter(1,( (probe_mesSA_roadRpcBw -> at(NTrigChain)[1]/1000.) - rWidthToBw(probe_mesSA_roadRpcAw->at(NTrigChain)[1],0.2)));
    f_roadRpc_BM_minus.SetLineColor(kCyan+1);
    f_roadRpc_BM_minus.SetLineWidth(2);
    f_roadRpc_BM_minus.SetLineStyle(1);

    TF1 f_roadFtk_BM = TF1("f_roadFtk_BM", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BM.SetTitle(";Z [m];R [m]");
    f_roadFtk_BM.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[1]);
    f_roadFtk_BM.SetParameter(1,probe_mesSA_roadFtkBw -> at(NTrigChain)[1]/1000.);
    f_roadFtk_BM.SetLineColor(kMagenta+1);
    f_roadFtk_BM.SetLineWidth(2);
    f_roadFtk_BM.SetLineStyle(2);

    //double aw_ftk_BM = probe_mesSA_roadFtkAw->at(NTrigChain)[1];
    //double bw_ftk_BM = probe_mesSA_roadFtkBw->at(NTrigChain)[1]/1000.;
    //double x1_ftk_BM = (6 - bw_ftk_BM)/aw_ftk_BM;
    //double x2_ftk_BM = (8 - bw_ftk_BM)/aw_ftk_BM;
    //double xMin_ftk_BM = (x1_ftk_BM<x2_ftk_BM)? x1_ftk_BM:x2_ftk_BM;
    //double xMax_ftk_BM = (x1_ftk_BM<x2_ftk_BM)? x2_ftk_BM:x1_ftk_BM;
    //f_roadFtk_BM.SetRange(xMin_ftk_BM,xMax_ftk_BM,6,8);

    TF1 f_roadFtk_BM_plus = TF1("f_roadFtk_BM_plus", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BM_plus.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[1]);
    f_roadFtk_BM_plus.SetParameter(1,( (probe_mesSA_roadFtkBw -> at(NTrigChain)[1]/1000.) + rWidthToBw(probe_mesSA_roadFtkAw->at(NTrigChain)[1],0.2)));
    f_roadFtk_BM_plus.SetLineColor(kMagenta+1);
    f_roadFtk_BM_plus.SetLineWidth(2);
    f_roadFtk_BM_plus.SetLineStyle(1);

    TF1 f_roadFtk_BM_minus = TF1("f_roadFtk_BM_minus", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BM_minus.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[1]);
    f_roadFtk_BM_minus.SetParameter(1,( (probe_mesSA_roadFtkBw -> at(NTrigChain)[1]/1000.) - rWidthToBw(probe_mesSA_roadFtkAw->at(NTrigChain)[1],0.2)));
    f_roadFtk_BM_minus.SetLineColor(kMagenta+1);
    f_roadFtk_BM_minus.SetLineWidth(2);
    f_roadFtk_BM_minus.SetLineStyle(1);

    TF1 f_roadRpc_BO = TF1("f_roadRpc_BO", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BO.SetTitle(";Z [m];R [m]");
    f_roadRpc_BO.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[2]);
    f_roadRpc_BO.SetParameter(1,probe_mesSA_roadRpcBw -> at(NTrigChain)[2]/1000.);
    f_roadRpc_BO.SetLineColor(kCyan+2);
    f_roadRpc_BO.SetLineWidth(2);
    f_roadRpc_BO.SetLineStyle(2);

    //double aw_rpc_BO = probe_mesSA_roadRpcAw->at(NTrigChain)[2];
    //double bw_rpc_BO = probe_mesSA_roadRpcBw->at(NTrigChain)[2]/1000.;
    //double x1_rpc_BO = (8 - bw_rpc_BO)/aw_rpc_BO;
    //double x2_rpc_BO = (11 - bw_rpc_BO)/aw_rpc_BO;
    //double xMin_rpc_BO = (x1_rpc_BO<x2_rpc_BO)? x1_rpc_BO:x2_rpc_BO;
    //double xMax_rpc_BO = (x1_rpc_BO<x2_rpc_BO)? x2_rpc_BO:x1_rpc_BO;
    //f_roadRpc_BO.SetRange(xMin_rpc_BO,xMax_rpc_BO,8,11);

    TF1 f_roadRpc_BO_plus = TF1("f_roadRpc_BO_plus", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BO_plus.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[2]);
    f_roadRpc_BO_plus.SetParameter(1,( (probe_mesSA_roadRpcBw -> at(NTrigChain)[2]/1000.) + rWidthToBw(probe_mesSA_roadRpcAw->at(NTrigChain)[2],0.4)));
    f_roadRpc_BO_plus.SetLineColor(kCyan+2);
    f_roadRpc_BO_plus.SetLineWidth(2);
    f_roadRpc_BO_plus.SetLineStyle(1);

    TF1 f_roadRpc_BO_minus = TF1("f_roadRpc_BO_minus", "[0]*x+[1]", -9999, 9999);
    f_roadRpc_BO_minus.SetParameter(0,probe_mesSA_roadRpcAw -> at(NTrigChain)[2]);
    f_roadRpc_BO_minus.SetParameter(1,( (probe_mesSA_roadRpcBw -> at(NTrigChain)[2]/1000.) - rWidthToBw(probe_mesSA_roadRpcAw->at(NTrigChain)[2],0.4)));
    f_roadRpc_BO_minus.SetLineColor(kCyan+2);
    f_roadRpc_BO_minus.SetLineWidth(2);
    f_roadRpc_BO_minus.SetLineStyle(1);

    TF1 f_roadFtk_BO = TF1("f_roadFtk_BO", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BO.SetTitle(";Z [m];R [m]");
    f_roadFtk_BO.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[2]);
    f_roadFtk_BO.SetParameter(1,probe_mesSA_roadFtkBw -> at(NTrigChain)[2]/1000.);
    f_roadFtk_BO.SetLineColor(kMagenta+2);
    f_roadFtk_BO.SetLineWidth(2);
    f_roadFtk_BO.SetLineStyle(2);

    //double aw_ftk_BO = probe_mesSA_roadFtkAw->at(NTrigChain)[2];
    //double bw_ftk_BO = probe_mesSA_roadFtkBw->at(NTrigChain)[2]/1000.;
    //double x1_ftk_BO = (8 - bw_ftk_BO)/aw_ftk_BO;
    //double x2_ftk_BO = (11 - bw_ftk_BO)/aw_ftk_BO;
    //double xMin_ftk_BO = (x1_ftk_BO<x2_ftk_BO)? x1_ftk_BO:x2_ftk_BO;
    //double xMax_ftk_BO = (x1_ftk_BO<x2_ftk_BO)? x2_ftk_BO:x1_ftk_BO;
    //f_roadFtk_BO.SetRange(xMin_ftk_BO,xMax_ftk_BO,8,11);

    TF1 f_roadFtk_BO_plus = TF1("f_roadFtk_BO_plus", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BO_plus.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[2]);
    f_roadFtk_BO_plus.SetParameter(1,( (probe_mesSA_roadFtkBw -> at(NTrigChain)[2]/1000.) + rWidthToBw(probe_mesSA_roadFtkAw->at(NTrigChain)[2],0.4)));
    f_roadFtk_BO_plus.SetLineColor(kMagenta+2);
    f_roadFtk_BO_plus.SetLineWidth(2);
    f_roadFtk_BO_plus.SetLineStyle(1);

    TF1 f_roadFtk_BO_minus = TF1("f_roadFtk_BO_minus", "[0]*x+[1]", -9999, 9999);
    f_roadFtk_BO_minus.SetParameter(0,probe_mesSA_roadFtkAw -> at(NTrigChain)[2]);
    f_roadFtk_BO_minus.SetParameter(1,( (probe_mesSA_roadFtkBw -> at(NTrigChain)[2]/1000.) - rWidthToBw(probe_mesSA_roadFtkAw->at(NTrigChain)[2],0.4)));
    f_roadFtk_BO_minus.SetLineColor(kMagenta+2);
    f_roadFtk_BO_minus.SetLineWidth(2);
    f_roadFtk_BO_minus.SetLineStyle(1);

    TF1 f_road_EI = TF1("f_road_EI", "[0]*x+[1]", -9999, 9999);
    f_road_EI.SetTitle(";Z [m];R [m]");
    f_road_EI.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[3]);
    f_road_EI.SetParameter(1,probe_mesSA_roadBw -> at(NTrigChain)[3]/1000.);
    f_road_EI.SetLineColor(30);
    f_road_EI.SetLineWidth(2);
    f_road_EI.SetLineStyle(2);

    TF1 f_road_EI_plus = TF1("f_road_EI_plus", "[0]*x+[1]", -9999, 9999);
    f_road_EI_plus.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[3]);
    f_road_EI_plus.SetParameter(1,( (probe_mesSA_roadBw -> at(NTrigChain)[3]/1000.) + rWidthToBw(probe_mesSA_roadAw->at(NTrigChain)[3],0.4)));
    f_road_EI_plus.SetLineColor(30);
    f_road_EI_plus.SetLineWidth(2);
    f_road_EI_plus.SetLineStyle(1);

    TF1 f_road_EI_minus = TF1("f_road_EI_minus", "[0]*x+[1]", -9999, 9999);
    f_road_EI_minus.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[3]);
    f_road_EI_minus.SetParameter(1,( (probe_mesSA_roadBw -> at(NTrigChain)[3]/1000.) - rWidthToBw(probe_mesSA_roadAw->at(NTrigChain)[3],0.4)));
    f_road_EI_minus.SetLineColor(30);
    f_road_EI_minus.SetLineWidth(2);
    f_road_EI_minus.SetLineStyle(1);

    // RoI
    TF1 f_roi = TF1("f_roi", "[0]*x", -20, 20);
    f_roi.SetTitle(";Z [m];R [m]");
    f_roi.SetParameter(0, tan((2*atan(exp(-probe_mesSA_roiEta->at(NTrigChain))))));
    f_roi.SetLineColor(kYellow+2);
    f_roi.SetLineWidth(1);
    f_roi.SetLineStyle(9);

    // MdtRegion
    double Z_BI[4];
    double Z_BM[4];
    double Z_BO[4];
    double R_BI[4];
    double R_BM[4];
    double R_BO[4];
    double SlopeMin_BI = tan(2*atan(exp(- probe_mesSA_etaMin -> at(NTrigChain) [0])));
    double SlopeMax_BI = tan(2*atan(exp(- probe_mesSA_etaMax -> at(NTrigChain) [0])));
    double SlopeMin_BM = tan(2*atan(exp(- probe_mesSA_etaMin -> at(NTrigChain) [1])));
    double SlopeMax_BM = tan(2*atan(exp(- probe_mesSA_etaMax -> at(NTrigChain) [1])));
    double SlopeMin_BO = tan(2*atan(exp(- probe_mesSA_etaMin -> at(NTrigChain) [2])));
    double SlopeMax_BO = tan(2*atan(exp(- probe_mesSA_etaMax -> at(NTrigChain) [2])));
    Z_BI[0] = (probe_mesSA_rMax -> at(NTrigChain)[0] / 1000.) / SlopeMin_BI;
    Z_BI[1] = (probe_mesSA_rMax -> at(NTrigChain)[0] / 1000.) / SlopeMax_BI;
    Z_BI[2] = (probe_mesSA_rMin -> at(NTrigChain)[0] / 1000.) / SlopeMax_BI;
    Z_BI[3] = (probe_mesSA_rMin -> at(NTrigChain)[0] / 1000.) / SlopeMin_BI;
    R_BI[0] = (probe_mesSA_rMax -> at(NTrigChain)[0] / 1000.);
    R_BI[1] = (probe_mesSA_rMax -> at(NTrigChain)[0] / 1000.);
    R_BI[2] = (probe_mesSA_rMin -> at(NTrigChain)[0] / 1000.);
    R_BI[3] = (probe_mesSA_rMin -> at(NTrigChain)[0] / 1000.);

    Z_BM[0] = (probe_mesSA_rMax -> at(NTrigChain)[1] / 1000.) / SlopeMin_BM;
    Z_BM[1] = (probe_mesSA_rMax -> at(NTrigChain)[1] / 1000.) / SlopeMax_BM;
    Z_BM[2] = (probe_mesSA_rMin -> at(NTrigChain)[1] / 1000.) / SlopeMax_BM;
    Z_BM[3] = (probe_mesSA_rMin -> at(NTrigChain)[1] / 1000.) / SlopeMin_BM;
    R_BM[0] = (probe_mesSA_rMax -> at(NTrigChain)[1] / 1000.);
    R_BM[1] = (probe_mesSA_rMax -> at(NTrigChain)[1] / 1000.);
    R_BM[2] = (probe_mesSA_rMin -> at(NTrigChain)[1] / 1000.);
    R_BM[3] = (probe_mesSA_rMin -> at(NTrigChain)[1] / 1000.);

    Z_BO[0] = (probe_mesSA_rMax -> at(NTrigChain)[2] / 1000.) / SlopeMin_BO;
    Z_BO[1] = (probe_mesSA_rMax -> at(NTrigChain)[2] / 1000.) / SlopeMax_BO;
    Z_BO[2] = (probe_mesSA_rMin -> at(NTrigChain)[2] / 1000.) / SlopeMax_BO;
    Z_BO[3] = (probe_mesSA_rMin -> at(NTrigChain)[2] / 1000.) / SlopeMin_BO;
    R_BO[0] = (probe_mesSA_rMax -> at(NTrigChain)[2] / 1000.);
    R_BO[1] = (probe_mesSA_rMax -> at(NTrigChain)[2] / 1000.);
    R_BO[2] = (probe_mesSA_rMin -> at(NTrigChain)[2] / 1000.);
    R_BO[3] = (probe_mesSA_rMin -> at(NTrigChain)[2] / 1000.);

    LOGD << "Z_BI: " << Z_BI[0] << ":, " <<  Z_BI[1] << ":, " << Z_BI[2] << ":, " << Z_BI[3]; 
    LOGD << "R_BI: " << R_BI[0] << ":, " <<  R_BI[1] << ":, " << R_BI[2] << ":, " << R_BI[3]; 
    LOGD << "Z_BM: " << Z_BM[0] << ":, " <<  Z_BM[1] << ":, " << Z_BM[2] << ":, " << Z_BM[3]; 
    LOGD << "R_BM: " << R_BM[0] << ":, " <<  R_BM[1] << ":, " << R_BM[2] << ":, " << R_BM[3]; 
    LOGD << "Z_BO: " << Z_BO[0] << ":, " <<  Z_BO[1] << ":, " << Z_BO[2] << ":, " << Z_BO[3]; 
    LOGD << "R_BO: " << R_BO[0] << ":, " <<  R_BO[1] << ":, " << R_BO[2] << ":, " << R_BO[3]; 


    // find min, max of Z and R
    std::vector<double> vecZ;
    std::vector<double> vecR;
    std::vector<double> vecZ_BI;
    std::vector<double> vecR_BI;
    std::vector<double> vecZ_BM;
    std::vector<double> vecR_BM;
    std::vector<double> vecZ_BO;
    std::vector<double> vecR_BO;
    for ( int i = 0; i < 4; i++){
      vecZ.push_back(Z_BI[i]);
      vecZ.push_back(Z_BM[i]);
      vecZ.push_back(Z_BO[i]);
      vecR.push_back(R_BI[i]);
      vecR.push_back(R_BM[i]);
      vecR.push_back(R_BO[i]);

      vecZ_BI.push_back(Z_BI[i]);
      vecR_BI.push_back(R_BI[i]);

      vecZ_BM.push_back(Z_BM[i]);
      vecR_BM.push_back(R_BM[i]);

      vecZ_BO.push_back(Z_BO[i]);
      vecR_BO.push_back(R_BO[i]);
    }
    double Zmin = *std::min_element(vecZ.begin(), vecZ.end());
    double Zmax = *std::max_element(vecZ.begin(), vecZ.end());
    double Rmin = *std::min_element(vecR.begin(), vecR.end());
    double Rmax = *std::max_element(vecR.begin(), vecR.end());

    double Zmin_BI = *std::min_element(vecZ_BI.begin(), vecZ_BI.end());
    double Zmax_BI = *std::max_element(vecZ_BI.begin(), vecZ_BI.end());
    double Rmin_BI = *std::min_element(vecR_BI.begin(), vecR_BI.end());
    double Rmax_BI = *std::max_element(vecR_BI.begin(), vecR_BI.end());

    double Zmin_BM = *std::min_element(vecZ_BM.begin(), vecZ_BM.end());
    double Zmax_BM = *std::max_element(vecZ_BM.begin(), vecZ_BM.end());
    double Rmin_BM = *std::min_element(vecR_BM.begin(), vecR_BM.end());
    double Rmax_BM = *std::max_element(vecR_BM.begin(), vecR_BM.end());

    double Zmin_BO = *std::min_element(vecZ_BO.begin(), vecZ_BO.end());
    double Zmax_BO = *std::max_element(vecZ_BO.begin(), vecZ_BO.end());
    double Rmin_BO = *std::min_element(vecR_BO.begin(), vecR_BO.end());
    double Rmax_BO = *std::max_element(vecR_BO.begin(), vecR_BO.end());

    LOGD << "Zmin: " << Zmin << ", Zmax: " <<  Zmax << ", Rmin: " << Rmin << ", Rmax: " << Rmax;
    LOGD << "Zmin_BI: " << Zmin_BI << ", Zmax_BI: " <<  Zmax_BI << ", Rmin_BI: " << Rmin_BI << ", Rmax_BI: " << Rmax_BI;
    LOGD << "Zmin_BM: " << Zmin_BM << ", Zmax_BM: " <<  Zmax_BM << ", Rmin_BM: " << Rmin_BM << ", Rmax_BM: " << Rmax_BM;
    LOGD << "Zmin_BO: " << Zmin_BO << ", Zmax_BO: " <<  Zmax_BO << ", Rmin_BO: " << Rmin_BO << ", Rmax_BO: " << Rmax_BO;

    Zmin -= (abs(Zmin) >= ZERO_LIMIT) ? abs(Zmax - Zmin) / 2. : deltaZ;
    Zmax += (abs(Zmax) >= ZERO_LIMIT) ? abs(Zmax - Zmin) / 2. : deltaZ;
    Rmin -= (abs(Rmin) >= ZERO_LIMIT) ? abs(Rmax - Rmin) / 2. : deltaR;
    Rmax += (abs(Rmax) >= ZERO_LIMIT) ? abs(Rmax - Rmin) / 2. : deltaR;

    Zmin_BI -= (abs(Zmin_BI) >= ZERO_LIMIT) ? abs(Zmax_BI - Zmin_BI) / 2. : deltaZ;
    Zmax_BI += (abs(Zmax_BI) >= ZERO_LIMIT) ? abs(Zmax_BI - Zmin_BI) / 2. : deltaZ;
    Rmin_BI -= (abs(Rmin_BI) >= ZERO_LIMIT) ? abs(Rmax_BI - Rmin_BI) / 2. : deltaR;
    Rmax_BI += (abs(Rmax_BI) >= ZERO_LIMIT) ? abs(Rmax_BI - Rmin_BI) / 2. : deltaR;

    Zmin_BM -= (abs(Zmin_BM) >= ZERO_LIMIT) ? abs(Zmax_BM - Zmin_BM) / 2. : deltaZ;
    Zmax_BM += (abs(Zmax_BM) >= ZERO_LIMIT) ? abs(Zmax_BM - Zmin_BM) / 2. : deltaZ;
    Rmin_BM -= (abs(Rmin_BM) >= ZERO_LIMIT) ? abs(Rmax_BM - Rmin_BM) / 2. : deltaR;
    Rmax_BM += (abs(Rmax_BM) >= ZERO_LIMIT) ? abs(Rmax_BM - Rmin_BM) / 2. : deltaR;

    Zmin_BO -= (abs(Zmin_BO) >= ZERO_LIMIT) ? abs(Zmax_BO - Zmin_BO) / 2. : deltaZ;
    Zmax_BO += (abs(Zmax_BO) >= ZERO_LIMIT) ? abs(Zmax_BO - Zmin_BO) / 2. : deltaZ;
    Rmin_BO -= (abs(Rmin_BO) >= ZERO_LIMIT) ? abs(Rmax_BO - Rmin_BO) / 2. : deltaR;
    Rmax_BO += (abs(Rmax_BO) >= ZERO_LIMIT) ? abs(Rmax_BO - Rmin_BO) / 2. : deltaR;

    LOGD << "Zmin: " << Zmin << ", Zmax: " <<  Zmax << ", Rmin: " << Rmin << ", Rmax: " << Rmax;
    LOGD << "Zmin_BI: " << Zmin_BI << ", Zmax_BI: " <<  Zmax_BI << ", Rmin_BI: " << Rmin_BI << ", Rmax_BI: " << Rmax_BI;
    LOGD << "Zmin_BM: " << Zmin_BM << ", Zmax_BM: " <<  Zmax_BM << ", Rmin_BM: " << Rmin_BM << ", Rmax_BM: " << Rmax_BM;
    LOGD << "Zmin_BO: " << Zmin_BO << ", Zmax_BO: " <<  Zmax_BO << ", Rmin_BO: " << Rmin_BO << ", Rmax_BO: " << Rmax_BO;

    TGraph MdtRegion_BI = TGraph(4,Z_BI,R_BI);
    TGraph MdtRegion_BM = TGraph(4,Z_BM,R_BM);
    TGraph MdtRegion_BO = TGraph(4,Z_BO,R_BO);

    MdtRegion_BI.SetFillStyle(3002);
    MdtRegion_BI.SetFillColorAlpha(13,0.40);
    MdtRegion_BM.SetFillStyle(3002);
    MdtRegion_BM.SetFillColorAlpha(13,0.40);
    MdtRegion_BO.SetFillStyle(3002);
    MdtRegion_BO.SetFillColorAlpha(13,0.40);

    // Set Legend
    TLegend leg = TLegend(0.805,0.22,0.99,0.95);
    leg.SetTextSize(0.035);
    leg.SetHeader(Form("#splitline{(1) Barrel All}{%s}", label_for_sector.Data()),"C");
    leg.AddEntry(&gr_segment,Form("#splitline{Offline}{segment (%d)}",nOffSeg),"p");
    leg.AddEntry(&gr_SP,Form("SuperPoint (%d)",NumberOfSP(NTrigChain)),"p");
    leg.AddEntry(&gr_MdtHit_Inlier,Form("#splitline{Inlier}{MDT hit (%d)}",nMDT_Inlier),"p");
    leg.AddEntry(&gr_MdtHit_Outlier_1,Form("#splitline{Outlier}{MDT hit (%d)}",nMDT_Outlier),"p");
    leg.AddEntry(&gr_RPC,Form("RPC hit (%d)",nRPC),"p");
    leg.AddEntry(&MdtRegion_BI,"MDT Region","f");
    leg.AddEntry(&f_roi,"RoI direction","l");

    // Set Legend
    TLegend leg_left = TLegend(-0.027,0.08,0.10,0.9);
    leg_left.SetTextSize(0.03);
    // Display passe or not
    TLegendEntry* l_pass;
    if (probe_mesSA_pass->at(NTrigChain) == 1){
      l_pass = leg_left.AddEntry((TObject*)0,"#splitline{Passed}{in L2MuonSA}","");
      l_pass -> SetTextColor(kGreen+2);
    } else{
      l_pass = leg_left.AddEntry((TObject*)0,"#splitline{NOT Passed}{in L2MuonSA}","");
      l_pass -> SetTextColor(kRed);
    }
    // Display passe or not
    TLegendEntry* l_rpc;
    if (probe_mesSA_isRpcFailure->at(NTrigChain) == 0){
      l_rpc = leg_left.AddEntry((TObject*)0,"RPC Fit Success","");
      l_rpc -> SetTextColor(kGreen+2);
    } else{
      l_rpc = leg_left.AddEntry((TObject*)0,"RPC Fit Failure","");
      l_rpc -> SetTextColor(kRed);
    }
    // Display road
    TLegendEntry* l_road;
    if (probe_mesSA_roadAlgo->at(NTrigChain) == 1){
      l_road = leg_left.AddEntry((TObject*)0,"FTK road used","");
      l_road -> SetTextColor(kMagenta+2);
    } else{
      l_road = leg_left.AddEntry((TObject*)0,"RPC road used","");
      l_road -> SetTextColor(kCyan+2);
    }


    leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe pT}{ : %4.3f [GeV]}",probe_pt/1000.),"");
    leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe pT}{ : %4.3f [GeV]}",abs(probe_mesSA_pt->at(NTrigChain))),"");
    leg_left.AddEntry((TObject*)0,Form("#splitline{FTK probe pT}{ : %4.3f [GeV]}",abs(probe_mesSA_ptFtk->at(NTrigChain))/1000.),"");
    leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe #eta}{ : %4.3f}",probe_eta),"");
    leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe #eta}{ : %4.3f}",probe_mesSA_eta->at(NTrigChain)),"");
    leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe #phi}{ : %4.3f}",probe_phi),"");
    leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe #phi}{ : %4.3f}",probe_mesSA_phi->at(NTrigChain)),"");
    leg_left.AddEntry((TObject*)0,Form("#splitline{RoI #eta}{ : %4.3f}",probe_mesSA_roiEta->at(NTrigChain)),"");
    leg_left.AddEntry((TObject*)0,Form("#splitline{RoI #phi}{ : %4.3f}",probe_mesSA_roiPhi->at(NTrigChain)),"");


    // Draw Display
    TH1 *frame = c2->DrawFrame(Zmin,Rmin,Zmax,Rmax);
    frame->GetXaxis()->SetTitle("Z [m]");
    frame->GetYaxis()->SetTitle("R [m]");
    gr_segment.Draw("P");
    f_roi.Draw("same");
    MdtRegion_BI.Draw("f, same");
    MdtRegion_BM.Draw("f, same");
    MdtRegion_BO.Draw("f, same");
    f_roadRpc_BI.Draw(" same");
    //f_roadRpc_BI_plus.Draw(" same");
    //f_roadRpc_BI_minus.Draw(" same");
    f_roadFtk_BI.Draw(" same");
    //f_roadFtk_BI_plus.Draw(" same");
    //f_roadFtk_BI_minus.Draw(" same");
    f_roadRpc_BM.Draw(" same");
    f_roadFtk_BM.Draw(" same");
    f_roadRpc_BO.Draw(" same");
    f_roadFtk_BO.Draw(" same");
    gr_segment.SetMarkerStyle(21);
    gr_segment.SetMarkerSize(2);
    gr_segment.SetMarkerColor(6);
    gr_segment_EI.SetMarkerStyle(21);
    gr_segment_EI.SetMarkerSize(2);
    gr_segment_EI.SetMarkerColor(kBlack);
    gr_segment.GetXaxis()->SetLimits(Zmin,Zmax);
    gr_segment.SetTitle(";Z [m];R [m]");
    gr_MdtHit_Inlier.SetMarkerStyle(24);
    gr_MdtHit_Inlier.SetMarkerSize(1);
    gr_MdtHit_Inlier.SetMarkerColor(kGreen+2);
    gr_MdtHit_Inlier.Draw("P");
    gr_MdtHit_Outlier_1.SetMarkerStyle(24);
    gr_MdtHit_Outlier_1.SetMarkerSize(1);
    gr_MdtHit_Outlier_1.SetMarkerColor(kRed);
    gr_MdtHit_Outlier_1.Draw("P");
    gr_MdtHit_Outlier_2.SetMarkerStyle(24);
    gr_MdtHit_Outlier_2.SetMarkerSize(1);
    gr_MdtHit_Outlier_2.SetMarkerColor(kBlue);
    gr_MdtHit_Outlier_2.Draw("P");
    gr_MdtHit_Outlier_3.SetMarkerStyle(24);
    gr_MdtHit_Outlier_3.SetMarkerSize(1);
    gr_MdtHit_Outlier_3.SetMarkerColor(kBlack);
    gr_MdtHit_Outlier_3.Draw("P");
    leg.Draw();
    leg_left.Draw();
    gr_RPC.SetMarkerStyle(22);
    gr_RPC.SetMarkerSize(1);
    gr_RPC.SetMarkerColor(kCyan-3);
    gr_RPC.Draw("P, same");
    gr_SP.SetMarkerStyle(8);
    gr_SP.SetMarkerSize(4);
    gr_SP.SetMarkerColor(4);
    gr_SP.Draw("P, same");
    gr_SP_EI.SetMarkerStyle(8);
    gr_SP_EI.SetMarkerSize(4);
    gr_SP_EI.SetMarkerColor(kGreen+2);
    gr_SP_EI.Draw("P, same");
    gr_SP_BI.SetMarkerStyle(8);
    gr_SP_BI.SetMarkerSize(4);
    gr_SP_BI.SetMarkerColor(kRed+2);
    gr_SP_BI.Draw("P, same");
    gr_segment.Draw("P, same");
    gr_segment_EI.Draw("P, same");

    TText eventInfo = TText(0.05,0.02,Form("EventNumber = %llu, RunNumber = %d, LumiBlock = %d, AverageInteractionsPerCrossing = %5.3f",EventNumber, RunNumber, LumiBlock, AverageInteractionsPerCrossing));
    eventInfo.SetNDC();
    eventInfo.SetTextSize(0.03);
    eventInfo.Draw();
    c2->Print(pdf, "pdf");
    c2->RedrawAxis();
    delete frame;

    ///*
    // Set Legend for BI
    TLegend leg_BI = TLegend(0.805,0.22,0.99,0.95);
    leg_BI.SetHeader(Form("#splitline{(2) Barrel Inner}{%s}", label_for_sector.Data()),"C");
    leg_BI.SetTextSize(0.035);
    leg_BI.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
    leg_BI.AddEntry(&gr_SP,Form("#splitline{SuperPoint}{(%4.3f, %4.3f)}",(probe_mesSA_superPointZ_BI -> at(NTrigChain))/1000., (probe_mesSA_superPointR_BI -> at(NTrigChain))/1000.),      "p");
    leg_BI.AddEntry(&gr_MdtHit_Inlier,  Form("#splitline{Inlier}{MDT hit (%d)}",nMDT_BI_Inlier),  "p");
    leg_BI.AddEntry(&gr_MdtHit_Outlier_1, Form("#splitline{Outlier}{MDT hit (%d)}",nMDT_BI_Outlier), "p");
    leg_BI.AddEntry(&f_roadRpc_BI,            "Rpc Road",            "l");
    leg_BI.AddEntry(&f_roadFtk_BI,            "Ftk Road",            "l");
    leg_BI.AddEntry(&MdtRegion_BI,"MDT Region","f");
    leg_BI.AddEntry(&f_roi,"RoI direction","l");
    TH1* frame_BI = c2->DrawFrame(Zmin_BI,Rmin_BI,Zmax_BI,Rmax_BI);
    frame_BI->GetXaxis()->SetTitle("Z [m]");
    frame_BI->GetYaxis()->SetTitle("R [m]");
    f_roadRpc_BI.Draw("same");
    f_roadRpc_BI_plus.Draw("same");
    f_roadRpc_BI_minus.Draw("same");
    f_roadFtk_BI.Draw("same");
    f_roadFtk_BI_plus.Draw("same");
    f_roadFtk_BI_minus.Draw("same");
    //MdtRegion_BI.Draw("f, same");
    f_roi.Draw("same");
    gr_MdtHit_Outlier_3.Draw("P, same");
    gr_MdtHit_Outlier_2.Draw("P, same");
    gr_MdtHit_Outlier_1.Draw("P, same");
    gr_MdtHit_Inlier.Draw("P, same");
    gr_segment.Draw("P, same");

    leg_BI.Draw();
    // Set Legend
    TLegend leg_left_BI = TLegend(-0.027,0.12,0.10,0.8);
    leg_left_BI.SetTextSize(0.03);
    // Display passe or not
    TLegendEntry* l_pass_BI;
    if (probe_mesSA_pass->at(NTrigChain) == 1){
      l_pass_BI = leg_left_BI.AddEntry((TObject*)0,"#splitline{Passed}{in L2MuonSA}","");
      l_pass_BI -> SetTextColor(kGreen+2);
    } else{
      l_pass_BI = leg_left_BI.AddEntry((TObject*)0,"#splitline{NOT Passed}{in L2MuonSA}","");
      l_pass_BI -> SetTextColor(kRed);
    }
    // Display passe or not
    TLegendEntry* l_rpc_BI;
    if (probe_mesSA_isRpcFailure->at(NTrigChain) == 0){
      l_rpc_BI = leg_left_BI.AddEntry((TObject*)0,"RPC Fit Success","");
      l_rpc_BI -> SetTextColor(kGreen+2);
    } else{
      l_rpc_BI = leg_left_BI.AddEntry((TObject*)0,"RPC Fit Failure","");
      l_rpc_BI -> SetTextColor(kRed);
    }
    // Display road
    TLegendEntry* l_road_BI;
    if (probe_mesSA_roadAlgo->at(NTrigChain) == 1){
      l_road_BI = leg_left_BI.AddEntry((TObject*)0,"FTK road used","");
      l_road_BI -> SetTextColor(kMagenta+2);
    } else{
      l_road_BI = leg_left_BI.AddEntry((TObject*)0,"RPC road used","");
      l_road_BI -> SetTextColor(kCyan+2);
    }

    leg_left_BI.AddEntry((TObject*)0,Form("#splitline{Offline probe p_{T}}{ : %4.3f [GeV]}",probe_pt/1000.),"");
    leg_left_BI.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe p_{T}}{ : %4.3f [GeV]}",abs(probe_mesSA_pt->at(NTrigChain))),"");
    leg_left_BI.AddEntry((TObject*)0,Form("#splitline{FTK probe p_{T}}{ : %4.3f [GeV]}",abs(probe_mesSA_ptFtk->at(NTrigChain))/1000.),"");
    leg_left_BI.AddEntry((TObject*)0,Form("#splitline{FTK Road}{%4.3f, %4.3f}", probe_mesSA_roadFtkAw->at(NTrigChain)[0], probe_mesSA_roadFtkBw->at(NTrigChain)[0]/1000.),"");
    leg_left_BI.AddEntry((TObject*)0,Form("#splitline{RPC Road}{%4.3f, %4.3f}", probe_mesSA_roadRpcAw->at(NTrigChain)[0], probe_mesSA_roadRpcBw->at(NTrigChain)[0]/1000.),"");
    leg_left_BI.Draw();

    gr_SP.Draw("P, same");
    eventInfo.Draw();
    c2->Print(pdf, "pdf");
    c2->RedrawAxis();
    delete frame_BI;

    // Set Legend for BM
    TLegend leg_BM = TLegend(0.805,0.22,0.99,0.95);
    leg_BM.SetHeader(Form("#splitline{(3) Barrel Middle}{%s}", label_for_sector.Data()),"C");
    leg_BM.SetTextSize(0.035);
    leg_BM.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
    leg_BM.AddEntry(&gr_SP,Form("#splitline{SuperPoint}{(%4.3f, %4.3f)}",(probe_mesSA_superPointZ_BM -> at(NTrigChain))/1000., (probe_mesSA_superPointR_BM -> at(NTrigChain))/1000.),      "p");
    leg_BM.AddEntry(&gr_MdtHit_Inlier,  Form("#splitline{Inlier}{MDT hit (%d)}",nMDT_BM_Inlier),  "p");
    leg_BM.AddEntry(&gr_MdtHit_Outlier_1, Form("#splitline{Outlier}{MDT hit (%d)}",nMDT_BM_Outlier), "p");
    leg_BM.AddEntry(&gr_RPC,"RPC hit","p");
    leg_BM.AddEntry(&f_roadRpc_BM,"Rpc Road","l");
    leg_BM.AddEntry(&f_roadFtk_BM,"Ftk Road","l");
    leg_BM.AddEntry(&MdtRegion_BI,"MDT Region","f");
    leg_BM.AddEntry(&f_roi,"RoI direction","l");
    TH1* frame_BM = c2->DrawFrame(Zmin_BM,Rmin_BM,Zmax_BM,Rmax_BM);
    frame_BM->GetXaxis()->SetTitle("Z [m]");
    frame_BM->GetYaxis()->SetTitle("R [m]");
    f_roadRpc_BM.Draw("same");
    f_roadRpc_BM_plus.Draw("same");
    f_roadRpc_BM_minus.Draw("same");
    f_roadFtk_BM.Draw("same");
    f_roadFtk_BM_plus.Draw("same");
    f_roadFtk_BM_minus.Draw("same");
    //MdtRegion_BM.Draw("f, same");
    f_roi.Draw("same");
    gr_MdtHit_Outlier_3.Draw("P,same");
    gr_MdtHit_Outlier_2.Draw("P,same");
    gr_MdtHit_Outlier_1.Draw("P,same");
    gr_MdtHit_Inlier.Draw("P, same");
    gr_segment.Draw("P, same");

    leg_BM.Draw();
    // Set Legend
    TLegend leg_left_BM = TLegend(-0.027,0.12,0.10,0.8);
    leg_left_BM.SetTextSize(0.03);
    // Display passe or not
    TLegendEntry* l_pass_BM;
    if (probe_mesSA_pass->at(NTrigChain) == 1){
      l_pass_BM = leg_left_BM.AddEntry((TObject*)0,"#splitline{Passed}{in L2MuonSA}","");
      l_pass_BM -> SetTextColor(kGreen+2);
    } else{
      l_pass_BM = leg_left_BM.AddEntry((TObject*)0,"#splitline{NOT Passed}{in L2MuonSA}","");
      l_pass_BM -> SetTextColor(kRed);
    }
    // Display passe or not
    TLegendEntry* l_rpc_BM;
    if (probe_mesSA_isRpcFailure->at(NTrigChain) == 0){
      l_rpc_BM = leg_left_BM.AddEntry((TObject*)0,"RPC Fit Success","");
      l_rpc_BM -> SetTextColor(kGreen+2);
    } else{
      l_rpc_BM = leg_left_BM.AddEntry((TObject*)0,"RPC Fit Failure","");
      l_rpc_BM -> SetTextColor(kRed);
    }
    // Display road
    TLegendEntry* l_road_BM;
    if (probe_mesSA_roadAlgo->at(NTrigChain) == 1){
      l_road_BM = leg_left_BM.AddEntry((TObject*)0,"FTK road used","");
      l_road_BM -> SetTextColor(kMagenta+2);
    } else{
      l_road_BM = leg_left_BM.AddEntry((TObject*)0,"RPC road used","");
      l_road_BM -> SetTextColor(kCyan+2);
    }

    leg_left_BM.AddEntry((TObject*)0,Form("#splitline{Offline probe p_{T}}{ : %4.3f [GeV]}",probe_pt/1000.),"");
    leg_left_BM.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe p_{T}}{ : %4.3f [GeV]}",abs(probe_mesSA_pt->at(NTrigChain))),"");
    leg_left_BM.AddEntry((TObject*)0,Form("#splitline{FTK probe p_{T}}{ : %4.3f [GeV]}",abs(probe_mesSA_ptFtk->at(NTrigChain))/1000.),"");
    leg_left_BM.AddEntry((TObject*)0,Form("#splitline{FTK Road}{%4.3f, %4.3f}", probe_mesSA_roadFtkAw->at(NTrigChain)[1], probe_mesSA_roadFtkBw->at(NTrigChain)[1]/1000.),"");
    leg_left_BM.AddEntry((TObject*)0,Form("#splitline{RPC Road}{%4.3f, %4.3f}", probe_mesSA_roadRpcAw->at(NTrigChain)[1], probe_mesSA_roadRpcBw->at(NTrigChain)[1]/1000.),"");
    leg_left_BM.Draw();

    gr_RPC.Draw("P, same");

    gr_SP.Draw("P, same");
    eventInfo.Draw();
    c2->Print(pdf, "pdf");
    c2->RedrawAxis();
    delete frame_BM;

    // Set Legend for BO
    TLegend leg_BO = TLegend(0.81,0.22,0.99,0.95);
    leg_BO.SetHeader(Form("#splitline{(4) Barrel Outer}{%s}", label_for_sector.Data()),"C");
    leg_BO.SetTextSize(0.035);
    leg_BO.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
    leg_BO.AddEntry(&gr_SP,Form("#splitline{SuperPoint}{(%4.3f, %4.3f)}",(probe_mesSA_superPointZ_BO -> at(NTrigChain))/1000., (probe_mesSA_superPointR_BO -> at(NTrigChain))/1000.),      "p");
    leg_BO.AddEntry(&gr_MdtHit_Inlier,  Form("#splitline{Inlier}{MDT hit (%d)}",nMDT_BO_Inlier),  "p");
    leg_BO.AddEntry(&gr_MdtHit_Outlier_1, Form("#splitline{Outlier}{MDT hit (%d)}",nMDT_BO_Outlier), "p");
    leg_BO.AddEntry(&gr_RPC,"RPC hit","p");
    leg_BO.AddEntry(&f_roadRpc_BO,"Rpc Road","l");
    leg_BO.AddEntry(&f_roadFtk_BO,"Ftk Road","l");
    leg_BO.AddEntry(&MdtRegion_BI,"MDT Region","f");
    leg_BO.AddEntry(&f_roi,"RoI direction","l");
    TH1* frame_BO = c2->DrawFrame(Zmin_BO,Rmin_BO,Zmax_BO,Rmax_BO);
    frame_BO->GetXaxis()->SetTitle("Z [m]");
    frame_BO->GetYaxis()->SetTitle("R [m]");
    f_roadRpc_BO.Draw("same");
    f_roadRpc_BO_plus.Draw("same");
    f_roadRpc_BO_minus.Draw("same");
    f_roadFtk_BO.Draw("same");
    f_roadFtk_BO_plus.Draw("same");
    f_roadFtk_BO_minus.Draw("same");
    //MdtRegion_BO.Draw("f, same");
    f_roi.Draw("same");
    gr_MdtHit_Outlier_3.Draw("P,same");
    gr_MdtHit_Outlier_2.Draw("P,same");
    gr_MdtHit_Outlier_1.Draw("P,same");
    gr_MdtHit_Inlier.Draw("P, same");
    gr_segment.Draw("P, same");

    leg_BO.Draw();
    // Set Legend
    TLegend leg_left_BO = TLegend(-0.027,0.12,0.10,0.8);
    leg_left_BO.SetTextSize(0.03);
    // Display passe or not
    TLegendEntry* l_pass_BO;
    if (probe_mesSA_pass->at(NTrigChain) == 1){
      l_pass_BO = leg_left_BO.AddEntry((TObject*)0,"#splitline{Passed}{in L2MuonSA}","");
      l_pass_BO -> SetTextColor(kGreen+2);
    } else{
      l_pass_BO = leg_left_BO.AddEntry((TObject*)0,"#splitline{NOT Passed}{in L2MuonSA}","");
      l_pass_BO -> SetTextColor(kRed);
    }
    // Display passe or not
    TLegendEntry* l_rpc_BO;
    if (probe_mesSA_isRpcFailure->at(NTrigChain) == 0){
      l_rpc_BO = leg_left_BO.AddEntry((TObject*)0,"RPC Fit Success","");
      l_rpc_BO -> SetTextColor(kGreen+2);
    } else{
      l_rpc_BO = leg_left_BO.AddEntry((TObject*)0,"RPC Fit Failure","");
      l_rpc_BO -> SetTextColor(kRed);
    }
    // Display road
    TLegendEntry* l_road_BO;
    if (probe_mesSA_roadAlgo->at(NTrigChain) == 1){
      l_road_BO = leg_left_BO.AddEntry((TObject*)0,"FTK road used","");
      l_road_BO -> SetTextColor(kMagenta+2);
    } else{
      l_road_BO = leg_left_BO.AddEntry((TObject*)0,"RPC road used","");
      l_road_BO -> SetTextColor(kCyan+2);
    }

    leg_left_BO.AddEntry((TObject*)0,Form("#splitline{Offline probe p_{T}}{ : %4.3f [GeV]}",probe_pt/1000.),"");
    leg_left_BO.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe p_{T}}{ : %4.3f [GeV]}",abs(probe_mesSA_pt->at(NTrigChain))),"");
    leg_left_BO.AddEntry((TObject*)0,Form("#splitline{FTK probe p_{T}}{ : %4.3f [GeV]}",abs(probe_mesSA_ptFtk->at(NTrigChain))/1000.),"");
    leg_left_BO.AddEntry((TObject*)0,Form("#splitline{FTK Road}{%4.3f, %4.3f}", probe_mesSA_roadFtkAw->at(NTrigChain)[2], probe_mesSA_roadFtkBw->at(NTrigChain)[2]/1000.),"");
    leg_left_BO.AddEntry((TObject*)0,Form("#splitline{RPC Road}{%4.3f, %4.3f}", probe_mesSA_roadRpcAw->at(NTrigChain)[2], probe_mesSA_roadRpcBw->at(NTrigChain)[2]/1000.),"");
    leg_left_BO.Draw();

    gr_RPC.Draw("P, same");

    gr_SP.Draw("P, same");
    eventInfo.Draw();
    c2->Print(pdf, "pdf");
    c2->RedrawAxis();

    delete frame_BO;
    //*/

    current_entry += 1;
    LOGI << "===" << begin_entry << ": " << current_entry << ": " << limit_entry;
    if (current_entry > limit_entry) {
      LOGI << "END!!";
      break;
    }
    c2->Clear();
  } // end of each entry
  c2->Print(pdf + "]", "pdf");
  delete c2;
}


// Size of Mdt hits and RPC hits
//const Int_t nMDT = (probe_mesSA_mdtHitChamber -> at(NTrigChain)).size();
//const Int_t nRPC = (probe_mesSA_rpcHitR -> at(NTrigChain)).size();
//
//// Offline segment
//// Set superpoint and segment for each station
//TGraph gr_segment = TGraph(probe_segment_n);
//for ( int i = 0; i < probe_segment_n; i++){
//  if (probe_segment_chamberIndex[i] >= 0 && probe_segment_chamberIndex[i] <= 5) {
//    double R = TMath::Sqrt(probe_segment_x[i]*probe_segment_x[i] + probe_segment_y[i]*probe_segment_y[i]);
//    double Z = probe_segment_z[i];
//    gr_segment.SetPoint(i, Z / 1000., R / 1000.);
//  }
//}
//
//// RPC hits
//TGraph gr_RPC = TGraph((probe_mesSA_rpcHitX ->at(NTrigChain)).size());
//for (int i = 0; i< nRPC; i++){
//  if (probe_mesSA_rpcHitStationNumber->at(NTrigChain)[i] == 1. ||
//      probe_mesSA_rpcHitStationNumber->at(NTrigChain)[i] == 2. ||
//      probe_mesSA_rpcHitStationNumber->at(NTrigChain)[i] == 5. ||
//      probe_mesSA_rpcHitStationNumber->at(NTrigChain)[i] == 6.) {
//    gr_RPC.SetPoint(i, (probe_mesSA_rpcHitZ -> at(NTrigChain)[i]) / 1000., (probe_mesSA_rpcHitR -> at(NTrigChain)[i]) / 1000.);
//  }
//}
//
//// SetPoint each SuperPoint
//TGraph gr_SP = TGraph(NumberOfSP());
//if (abs(probe_mesSA_superPointZ_BI -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_BI -> at(NTrigChain) > -99990.){
//  cout << "SP_BI: " << probe_mesSA_superPointR_BI -> at(NTrigChain) << endl;
//  //cout << "sAdress: " << probe_mesSA_sAddress -> at(NTrigChain) << endl;
//  gr_SP.SetPoint(0 , (probe_mesSA_superPointZ_BI -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BI -> at(NTrigChain)) / 1000.);
//}
//if (abs(probe_mesSA_superPointZ_BM -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_BM -> at(NTrigChain) > -99990.){
//  cout << "SP_BM: " << probe_mesSA_superPointR_BM -> at(NTrigChain) << endl;
//  gr_SP.SetPoint(1 , (probe_mesSA_superPointZ_BM -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BM -> at(NTrigChain)) / 1000.);
//}
//if (abs(probe_mesSA_superPointZ_BO -> at(NTrigChain)) > 0. && probe_mesSA_superPointZ_BO -> at(NTrigChain) > -99990.){
//  cout << "SP_BO: " << probe_mesSA_superPointR_BO -> at(NTrigChain) << endl;
//  gr_SP.SetPoint(2 , (probe_mesSA_superPointZ_BO -> at(NTrigChain)) / 1000. , (probe_mesSA_superPointR_BO -> at(NTrigChain)) / 1000.);
//}
//
//TGraph gr = TGraph(nMDT); //各点が(0,0)で初期化される
//// Inlier Mdt Hit
//TGraph gr_MdtHit_Inlier_BI = TGraph(nMDT); //各点が(0,0)で初期化される
//TGraph gr_MdtHit_Inlier_BM = TGraph(nMDT); //各点が(0,0)で初期化される
//TGraph gr_MdtHit_Inlier_BO = TGraph(nMDT); //各点が(0,0)で初期化される
//// Outlier Mdt Hit
//TGraph gr_MdtHit_Outlier_BI = TGraph(nMDT); //各点が(0,0)で初期化される
//TGraph gr_MdtHit_Outlier_BM = TGraph(nMDT); //各点が(0,0)で初期化される
//TGraph gr_MdtHit_Outlier_BO = TGraph(nMDT); //各点が(0,0)で初期化される
//// Outlier2 Mdt Hit
//TGraph gr_MdtHit_Outlier2_BI = TGraph(nMDT); //各点が(0,0)で初期化される
//TGraph gr_MdtHit_Outlier2_BM = TGraph(nMDT); //各点が(0,0)で初期化される
//TGraph gr_MdtHit_Outlier2_BO = TGraph(nMDT); //各点が(0,0)で初期化される
//
//
//double Z_BI=0, R_BI=0;
//double Z_BM=0, R_BM=0;
//double Z_BO=0, R_BO=0;
//double Zmin_BI=0, Zmax_BI=0, Rmin_BI=0, Rmax_BI=0;
//double Zmin_BM=0, Zmax_BM=0, Rmin_BM=0, Rmax_BM=0;
//double Zmin_BO=0, Zmax_BO=0, Rmin_BO=0, Rmax_BO=0;
//int nBI=0;
//int nBM=0;
//int nBO=0;
//int nBI_Inlier=0;
//int nBM_Inlier=0;
//int nBO_Inlier=0;
//int nBI_Outlier=0;
//int nBM_Outlier=0;
//int nBO_Outlier=0;
//double deltaZ = 1.;
//double deltaR = 0.4;
//
//// Setpoint each Mdt Hit
//for (Int_t i=0;i<nMDT;++i) {
//  cout << probe_mesSA_mdtHitChamber -> at(NTrigChain)[i] << endl;
//  cout << "Z: " << probe_mesSA_mdtHitZ -> at(NTrigChain)[i]/1000. << ", R: " << probe_mesSA_mdtHitR -> at(NTrigChain)[i]/1000. << endl;
//  cout << nBI << ":" << nBM << ":" << nBO << endl;
//  gr.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//
//  //BI
//  if (probe_mesSA_mdtHitChamber -> at(NTrigChain)[i] == 0){
//    if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[i] == 1) {
//      gr_MdtHit_Outlier_BI.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//      nBI_Outlier += 1;
//    } else if (probe_mesSA_mdtHitIsOutlier->at(NTrigChain)[i] == 0){
//      gr_MdtHit_Inlier_BI.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//      nBI_Inlier += 1;
//    } else if (probe_mesSA_mdtHitIsOutlier->at(NTrigChain)[i] == 2){
//      gr_MdtHit_Outlier2_BI.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//    }
//    Z_BI += probe_mesSA_mdtHitZ -> at(NTrigChain)[i]/1000.;
//    R_BI += probe_mesSA_mdtHitR -> at(NTrigChain)[i]/1000.;
//    nBI += 1;
//  }
//  //BM
//  if (probe_mesSA_mdtHitChamber -> at(NTrigChain)[i] == 1){
//    if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[i] == 1) {
//      gr_MdtHit_Outlier_BM.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//      nBM_Outlier += 1;
//    } else if (probe_mesSA_mdtHitIsOutlier->at(NTrigChain)[i] == 0){
//      gr_MdtHit_Inlier_BM.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//      nBM_Inlier += 1;
//    } else if (probe_mesSA_mdtHitIsOutlier->at(NTrigChain)[i] == 2){
//      gr_MdtHit_Outlier2_BM.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//    }
//    Z_BM += probe_mesSA_mdtHitZ -> at(NTrigChain)[i]/1000.;
//    R_BM += probe_mesSA_mdtHitR -> at(NTrigChain)[i]/1000.;
//    nBM += 1;
//  }
//  //BO
//  if (probe_mesSA_mdtHitChamber -> at(NTrigChain)[i] == 2){
//    if (probe_mesSA_mdtHitIsOutlier -> at(NTrigChain)[i] == 1) {
//      gr_MdtHit_Outlier_BO.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//      nBO_Outlier += 1;
//    } else if (probe_mesSA_mdtHitIsOutlier->at(NTrigChain)[i] == 0){
//      gr_MdtHit_Inlier_BO.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//      nBO_Inlier += 1;
//    } else if (probe_mesSA_mdtHitIsOutlier->at(NTrigChain)[i] == 2){
//      gr_MdtHit_Outlier2_BO.SetPoint(i , (probe_mesSA_mdtHitZ -> at(NTrigChain))[i] / 1000. , (probe_mesSA_mdtHitR -> at(NTrigChain))[i] / 1000.);
//    }
//    Z_BO += probe_mesSA_mdtHitZ -> at(NTrigChain)[i]/1000.;
//    R_BO += probe_mesSA_mdtHitR -> at(NTrigChain)[i]/1000.;
//    nBO += 1;
//  }
//}
//
//// Avoid nBI==0
//if (nBI != 0) {
//  Zmax_BI = Z_BI/nBI + deltaZ;
//  Zmin_BI = Z_BI/nBI - deltaZ;
//  Rmax_BI = R_BI/nBI + deltaR;
//  Rmin_BI = R_BI/nBI - deltaR;
//}else{
//  R_BI = 5;
//  deltaR = 1;
//  Z_BI = (R_BI - probe_mesSA_roadBw->at(NTrigChain)[0]/1000.) / (probe_mesSA_roadAw->at(NTrigChain)[0]);
//  Zmax_BI = Z_BI + deltaZ;
//  Zmin_BI = Z_BI - deltaZ;
//  Rmax_BI = R_BI + deltaR;
//  Rmin_BI = R_BI - deltaR;
//}
//// Avoid nBM==0
//if (nBM != 0) {
//  Zmax_BM = Z_BM/nBM + deltaZ;
//  Zmin_BM = Z_BM/nBM - deltaZ;
//  Rmax_BM = R_BM/nBM + deltaR;
//  Rmin_BM = R_BM/nBM - deltaR;
//}else{
//  R_BM = 5;
//  deltaR = 1;
//  Z_BM = (R_BM - probe_mesSA_roadBw->at(NTrigChain)[0]/1000.) / (probe_mesSA_roadAw->at(NTrigChain)[0]);
//  Zmax_BM = Z_BM + deltaZ;
//  Zmin_BM = Z_BM - deltaZ;
//  Rmax_BM = R_BM + deltaR;
//  Rmin_BM = R_BM - deltaR;
//}
//// Avoid nBO==0
//if (nBO != 0) {
//  Zmax_BO = Z_BO/nBO + deltaZ;
//  Zmin_BO = Z_BO/nBO - deltaZ;
//  Rmax_BO = R_BO/nBO + deltaR;
//  Rmin_BO = R_BO/nBO - deltaR;
//}else{
//  R_BO = 5;
//  deltaR = 1;
//  Z_BO = (R_BO - probe_mesSA_roadBw->at(NTrigChain)[0]/1000.) / (probe_mesSA_roadAw->at(NTrigChain)[0]);
//  Zmax_BO = Z_BO + deltaZ;
//  Zmin_BO = Z_BO - deltaZ;
//  Rmax_BO = R_BO + deltaR;
//  Rmin_BO = R_BO - deltaR;
//}
//// Road
//TF1 f_road_BI = TF1("f_road_BI", "[0]*x+[1]", -20, 20);
//cout << "road: " << probe_mesSA_roadAw -> at(NTrigChain)[0] << ": " << probe_mesSA_roadBw->at(NTrigChain)[0] << endl;
//cout << "road: " << probe_mesSA_roadAw -> at(NTrigChain)[1] << ": " << probe_mesSA_roadBw->at(NTrigChain)[1] << endl;
//cout << "road: " << probe_mesSA_roadAw -> at(NTrigChain)[2] << ": " << probe_mesSA_roadBw->at(NTrigChain)[2] << endl;
//f_road_BI.SetTitle(";Z [m];R [m]");
//f_road_BI.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[0]);
//f_road_BI.SetParameter(1,probe_mesSA_roadBw -> at(NTrigChain)[0]/1000.);
//f_road_BI.SetLineColor(15);
//f_road_BI.SetLineWidth(2);
//f_road_BI.SetLineStyle(2);
//
//TF1 f_road_BI_plus = TF1("f_road_BI_plus", "[0]*x+[1]", -20, 20);
//f_road_BI_plus.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[0]);
//f_road_BI_plus.SetParameter(1,( (probe_mesSA_roadBw -> at(NTrigChain)[0]/1000.) + rWidthToBw(probe_mesSA_roadAw->at(NTrigChain)[0],0.4)));
//f_road_BI_plus.SetLineColor(15);
//f_road_BI_plus.SetLineWidth(2);
//f_road_BI_plus.SetLineStyle(1);
//
//TF1 f_road_BI_minus = TF1("f_road_BI_minus", "[0]*x+[1]", -20, 20);
//f_road_BI_minus.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[0]);
//f_road_BI_minus.SetParameter(1,( (probe_mesSA_roadBw -> at(NTrigChain)[0]/1000.) - rWidthToBw(probe_mesSA_roadAw->at(NTrigChain)[0],0.4)));
//f_road_BI_minus.SetLineColor(15);
//f_road_BI_minus.SetLineWidth(2);
//f_road_BI_minus.SetLineStyle(1);
//
//
//TF1 f_road_BM = TF1("f_road_BM", "[0]*x+[1]", -20, 20);
//f_road_BM.SetTitle(";Z [m];R [m]");
//f_road_BM.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[1]);
//f_road_BM.SetParameter(1,probe_mesSA_roadBw -> at(NTrigChain)[1]/1000.);
//f_road_BM.SetLineColor(15);
//f_road_BM.SetLineWidth(2);
//f_road_BM.SetLineStyle(2);
//
//TF1 f_road_BM_plus = TF1("f_road_BM_plus", "[0]*x+[1]", -20, 20);
//f_road_BM_plus.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[1]);
//f_road_BM_plus.SetParameter(1,( (probe_mesSA_roadBw -> at(NTrigChain)[1]/1000.) + rWidthToBw(probe_mesSA_roadAw->at(NTrigChain)[1],0.2)));
//f_road_BM_plus.SetLineColor(15);
//f_road_BM_plus.SetLineWidth(2);
//f_road_BM_plus.SetLineStyle(1);
//
//TF1 f_road_BM_minus = TF1("f_road_BM_minus", "[0]*x+[1]", -20, 20);
//f_road_BM_minus.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[1]);
//f_road_BM_minus.SetParameter(1,( (probe_mesSA_roadBw -> at(NTrigChain)[1]/1000.) - rWidthToBw(probe_mesSA_roadAw->at(NTrigChain)[1],0.2)));
//f_road_BM_minus.SetLineColor(15);
//f_road_BM_minus.SetLineWidth(2);
//f_road_BM_minus.SetLineStyle(1);
//
//TF1 f_road_BO = TF1("f_road_BO", "[0]*x+[1]", -20, 20);
//f_road_BO.SetTitle(";Z [m];R [m]");
//f_road_BO.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[2]);
//f_road_BO.SetParameter(1,probe_mesSA_roadBw -> at(NTrigChain)[2]/1000.);
//f_road_BO.SetLineColor(15);
//f_road_BO.SetLineWidth(2);
//f_road_BO.SetLineStyle(2);
//
//TF1 f_road_BO_plus = TF1("f_road_BO_plus", "[0]*x+[1]", -20, 20);
//f_road_BO_plus.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[2]);
//f_road_BO_plus.SetParameter(1,( (probe_mesSA_roadBw -> at(NTrigChain)[2]/1000.) + rWidthToBw(probe_mesSA_roadAw->at(NTrigChain)[2],0.4)));
//f_road_BO_plus.SetLineColor(15);
//f_road_BO_plus.SetLineWidth(2);
//f_road_BO_plus.SetLineStyle(1);
//
//TF1 f_road_BO_minus = TF1("f_road_BO_minus", "[0]*x+[1]", -20, 20);
//f_road_BO_minus.SetParameter(0,probe_mesSA_roadAw -> at(NTrigChain)[2]);
//f_road_BO_minus.SetParameter(1,( (probe_mesSA_roadBw -> at(NTrigChain)[2]/1000.) - rWidthToBw(probe_mesSA_roadAw->at(NTrigChain)[2],0.4)));
//f_road_BO_minus.SetLineColor(15);
//f_road_BO_minus.SetLineWidth(2);
//f_road_BO_minus.SetLineStyle(1);
//
//
//// RoI
//TF1 f_roi = TF1("f_roi", "[0]*x", -20, 20);
//f_roi.SetTitle(";Z [m];R [m]");
//f_roi.SetParameter(0, tan((2*atan(exp(-probe_mesSA_roiEta->at(NTrigChain))))));
//f_roi.SetLineColor(kYellow+2);
//f_roi.SetLineWidth(1);
//f_roi.SetLineStyle(9);
//
//TF1 f_roi_min = TF1("f_roi_min", "[0]*x", -20, 20);
//f_roi_min.SetParameter(0, tan((2*atan(exp(- (probe_mesSA_roiEta->at(NTrigChain) - 0.05) )))));
//f_roi_min.SetLineColor(kYellow+2);
//f_roi_min.SetLineWidth(2);
//f_roi_min.SetLineStyle(1);
//
//TF1 f_roi_max = TF1("f_roi_max", "[0]*x", -20, 20);
//f_roi_max.SetParameter(0, tan((2*atan(exp(- (probe_mesSA_roiEta->at(NTrigChain) + 0.05) )))));
//f_roi_max.SetLineColor(kYellow+2);
//f_roi_max.SetLineWidth(2);
//f_roi_max.SetLineStyle(1);
//
//// Set Legend
//TLegend leg = TLegend(0.805,0.22,0.99,0.95);
//leg.SetTextSize(0.035);
//leg.SetHeader(Form("#splitline{(1) Barrel All}{%s}", label_for_sector.Data()),"C");
//leg.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
//leg.AddEntry(&gr,Form("MDT hit (%d)",nMDT),"p");
//leg.AddEntry(&gr_RPC,Form("RPC hit (%d)",nRPC),"p");
//leg.AddEntry(&gr_SP,Form("SuperPoint (%d)",NumberOfSP()),"p");
////leg.AddEntry(&f_road_BI,"Road","l");
//leg.AddEntry(&f_roi,"RoI","l");
//
//// Set Legend
//TLegend leg_left = TLegend(-0.027,0.12,0.10,0.8);
//leg_left.SetTextSize(0.03);
//// Display passe or not
//TLegendEntry* l_pass;
//if (probe_mesSA_pass->at(NTrigChain) == 1){
//  l_pass = leg_left.AddEntry((TObject*)0,"#splitline{Passed}{in L2MuonSA}","");
//  l_pass -> SetTextColor(kGreen+2);
//} else{
//  l_pass = leg_left.AddEntry((TObject*)0,"#splitline{NOT Passed}{in L2MuonSA}","");
//  l_pass -> SetTextColor(kRed);
//}
//// Display passe or not
//TLegendEntry* l_rpc;
//if (probe_mesSA_isRpcFailure->at(NTrigChain) == 0){
//  l_rpc = leg_left.AddEntry((TObject*)0,"RPC Fit Success","");
//  l_rpc -> SetTextColor(kGreen+2);
//} else{
//  l_rpc = leg_left.AddEntry((TObject*)0,"RPC Fit Failure","");
//  l_rpc -> SetTextColor(kRed);
//}
//
//
//leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe p_{T}}{: %4.3f [GeV]}",probe_pt/1000.),"");
//leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe p_{T}}{: %4.3f [GeV]}",abs(probe_mesSA_pt->at(NTrigChain))),"");
//leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe #eta}{: %4.3f}",probe_eta),"");
//leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe #eta}{: %4.3f}",probe_mesSA_eta->at(NTrigChain)),"");
//leg_left.AddEntry((TObject*)0,Form("#splitline{Offline probe #phi}{: %4.3f}",probe_phi),"");
//leg_left.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe #phi}{: %4.3f}",probe_mesSA_phi->at(NTrigChain)),"");
//
//
//// Draw Display
//double Zmax = 20;
//double Zmin = -20;
//TH1 *frame = c2->DrawFrame(Zmin,0,Zmax,12);
//frame->GetXaxis()->SetTitle("Z [m]");
//frame->GetYaxis()->SetTitle("R [m]");
//gr_segment.Draw("P");
//f_roi.Draw("same");
//gr_segment.SetMarkerStyle(21);
//gr_segment.SetMarkerSize(2);
//gr_segment.SetMarkerColor(6);
//gr_segment.GetXaxis()->SetLimits(Zmin,Zmax);
//gr_segment.SetTitle(";Z [m];R [m]");
//gr_segment.GetYaxis()->SetRangeUser(0,12);
//f_road_BI.Draw("same");
//f_road_BM.Draw("same");
//f_road_BO.Draw("same");
//gr.SetMarkerStyle(24);
//gr.SetMarkerSize(1);
//gr.SetMarkerColor(2);
//gr.Draw("P");
//leg.Draw();
//leg_left.Draw();
//gr_RPC.SetMarkerStyle(22);
//gr_RPC.SetMarkerSize(1);
//gr_RPC.SetMarkerColor(kCyan-3);
//gr_RPC.Draw("P, same");
//gr_SP.SetMarkerStyle(8);
//gr_SP.SetMarkerSize(2);
//gr_SP.SetMarkerColor(4);
//gr_SP.Draw("P, same");
//
//TText eventInfo = TText(0.32,0.03,Form("EventNumber=%d, RunNumber=%d, LumiBlock=%d",EventNumber, RunNumber, LumiBlock));
//eventInfo.SetNDC();
//eventInfo.SetTextSize(0.03);
//eventInfo.Draw();
//c2->Print(pdf, "pdf");
//c2->RedrawAxis();
//delete frame;
//
//// Set Legend for BI
//TLegend leg_BI = TLegend(0.805,0.22,0.99,0.95);
//leg_BI.SetHeader(Form("#splitline{(2) Barrel Inner}{%s}", label_for_sector.Data()),"C");
//leg_BI.SetTextSize(0.035);
//leg_BI.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
//leg_BI.AddEntry(&gr_MdtHit_Inlier_BI,  Form("#splitline{Inlier}{MDT hit (%d)}",nBI_Inlier),  "p");
//leg_BI.AddEntry(&gr_MdtHit_Outlier_BI, Form("#splitline{Outlier}{MDT hit (%d)}",nBI_Outlier), "p");
//leg_BI.AddEntry(&gr_RPC,               "RPC hit",         "p");
//leg_BI.AddEntry(&gr_SP,                "SuperPoint",      "p");
//leg_BI.AddEntry(&f_road_BI,            "Road",            "l");
//leg_BI.AddEntry(&f_roi,"RoI","l");
//TH1* frame_BI = c2->DrawFrame(Zmin_BI,Rmin_BI,Zmax_BI,Rmax_BI);
//frame_BI->GetXaxis()->SetTitle("Z [m]");
//frame_BI->GetYaxis()->SetTitle("R [m]");
//f_road_BI.Draw("same");
//f_road_BI_plus.Draw("same");
//f_road_BI_minus.Draw("same");
//f_roi.Draw("same");
//f_roi_min.Draw("same");
//f_roi_max.Draw("same");
/////cout << "Zmin_BI: " << Zmin_BI <<endl;
/////cout << "Zmax_BI: " << Zmax_BI <<endl;
//f_road_BI.GetXaxis()->SetLimits(Zmin_BI,Zmax_BI);
//f_road_BI.GetYaxis()->SetRangeUser(Rmin_BI,Rmax_BI);
//gr_MdtHit_Inlier_BI.SetMarkerColor(kGreen + 2);
//gr_MdtHit_Inlier_BI.SetMarkerStyle(24);
//gr_MdtHit_Inlier_BI.SetMarkerSize(1);
//gr_MdtHit_Inlier_BI.Draw("P, same");
//gr_segment.Draw("P, same");
//
//leg_BI.Draw();
//// Set Legend
//TLegend leg_left_BI = TLegend(-0.027,0.12,0.10,0.8);
//leg_left_BI.SetTextSize(0.03);
//// Display passe or not
//TLegendEntry* l_pass_BI;
//if (probe_mesSA_pass->at(NTrigChain) == 1){
//  l_pass_BI = leg_left_BI.AddEntry((TObject*)0,"#splitline{Passed}{in L2MuonSA}","");
//  l_pass_BI -> SetTextColor(kGreen+2);
//} else{
//  l_pass_BI = leg_left_BI.AddEntry((TObject*)0,"#splitline{NOT Passed}{in L2MuonSA}","");
//  l_pass_BI -> SetTextColor(kRed);
//}
//// Display passe or not
//TLegendEntry* l_rpc_BI;
//if (probe_mesSA_isRpcFailure->at(NTrigChain) == 0){
//  l_rpc_BI = leg_left_BI.AddEntry((TObject*)0,"RPC Fit Success","");
//  l_rpc_BI -> SetTextColor(kGreen+2);
//} else{
//  l_rpc_BI = leg_left_BI.AddEntry((TObject*)0,"RPC Fit Failure","");
//  l_rpc_BI -> SetTextColor(kRed);
//}
//
//leg_left_BI.AddEntry((TObject*)0,Form("#splitline{Offline probe p_{T}}{: %4.3f [GeV]}",probe_pt/1000.),"");
//leg_left_BI.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe p_{T}}{: %4.3f [GeV]}",abs(probe_mesSA_pt->at(NTrigChain))),"");
//leg_left_BI.AddEntry((TObject*)0,Form("#splitline{RoI #eta}{: %4.3f}",probe_mesSA_roiEta->at(NTrigChain)),"");
//leg_left_BI.AddEntry((TObject*)0,Form("#splitline{Road Slope #eta}{: %4.3f}", -TMath::Log(tan((atan(probe_mesSA_roadAw->at(NTrigChain)[0]))/2.))),"");
//leg_left_BI.AddEntry((TObject*)0,Form("#splitline{Road Slope}{: %4.3f}", probe_mesSA_roadAw->at(NTrigChain)[0]),"");
//leg_left_BI.AddEntry((TObject*)0,Form("#splitline{Road Z-Intercept}{: %4.3f}",probe_mesSA_roadBw->at(NTrigChain)[0]),"");
//leg_left_BI.Draw();
//
//gr_MdtHit_Outlier_BI.SetMarkerColor(kRed);
//gr_MdtHit_Outlier_BI.SetMarkerStyle(24);
//gr_MdtHit_Outlier_BI.SetMarkerSize(1);
//gr_MdtHit_Outlier_BI.Draw("P,same");
//
//gr_MdtHit_Outlier2_BI.SetMarkerColor(12);
//gr_MdtHit_Outlier2_BI.SetMarkerStyle(24);
//gr_MdtHit_Outlier2_BI.SetMarkerSize(1);
//gr_MdtHit_Outlier2_BI.Draw("P,same");
//
//gr_SP.Draw("P, same");
//eventInfo.Draw();
//c2->Print(pdf, "pdf");
//c2->RedrawAxis();
//delete frame_BI;
//
//// Set Legend for BM
//TLegend leg_BM = TLegend(0.805,0.22,0.99,0.95);
//leg_BM.SetHeader(Form("#splitline{(3) Barrel Middle}{%s}", label_for_sector.Data()),"C");
//leg_BM.SetTextSize(0.035);
//leg_BM.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
//leg_BM.AddEntry(&gr_MdtHit_Inlier_BM,  Form("#splitline{Inlier}{MDT hit (%d)}",nBM_Inlier),  "p");
//leg_BM.AddEntry(&gr_MdtHit_Outlier_BM, Form("#splitline{Outlier}{MDT hit (%d)}",nBM_Outlier), "p");
//leg_BM.AddEntry(&gr_RPC,"RPC hit","p");
//leg_BM.AddEntry(&gr_SP,"SuperPoint","p");
//leg_BM.AddEntry(&f_road_BM,"Road","l");
//leg_BM.AddEntry(&f_roi,"RoI","l");
//TH1* frame_BM = c2->DrawFrame(Zmin_BM,Rmin_BM,Zmax_BM,Rmax_BM);
//frame_BM->GetXaxis()->SetTitle("Z [m]");
//frame_BM->GetYaxis()->SetTitle("R [m]");
//f_road_BM.Draw("same");
//f_road_BM_plus.Draw("same");
//f_road_BM_minus.Draw("same");
//f_roi.Draw("same");
//f_roi_min.Draw("same");
//f_roi_max.Draw("same");
//f_road_BM.GetXaxis()->SetLimits(Zmin_BM,Zmax_BM);
//f_road_BM.GetYaxis()->SetRangeUser(Rmin_BM,Rmax_BM);
//gr_MdtHit_Inlier_BM.SetMarkerColor(kGreen + 2);
//gr_MdtHit_Inlier_BM.SetMarkerStyle(24);
//gr_MdtHit_Inlier_BM.SetMarkerSize(1);
//gr_MdtHit_Inlier_BM.Draw("P, same");
//gr_segment.Draw("P, same");
//
//leg_BM.Draw();
//// Set Legend
//TLegend leg_left_BM = TLegend(-0.027,0.12,0.10,0.8);
//leg_left_BM.SetTextSize(0.03);
//// Display passe or not
//TLegendEntry* l_pass_BM;
//if (probe_mesSA_pass->at(NTrigChain) == 1){
//  l_pass_BM = leg_left_BM.AddEntry((TObject*)0,"#splitline{Passed}{in L2MuonSA}","");
//  l_pass_BM -> SetTextColor(kGreen+2);
//} else{
//  l_pass_BM = leg_left_BM.AddEntry((TObject*)0,"#splitline{NOT Passed}{in L2MuonSA}","");
//  l_pass_BM -> SetTextColor(kRed);
//}
//// Display passe or not
//TLegendEntry* l_rpc_BM;
//if (probe_mesSA_isRpcFailure->at(NTrigChain) == 0){
//  l_rpc_BM = leg_left_BM.AddEntry((TObject*)0,"RPC Fit Success","");
//  l_rpc_BM -> SetTextColor(kGreen+2);
//} else{
//  l_rpc_BM = leg_left_BM.AddEntry((TObject*)0,"RPC Fit Failure","");
//  l_rpc_BM -> SetTextColor(kRed);
//}
//
//leg_left_BM.AddEntry((TObject*)0,Form("#splitline{Offline probe p_{T}}{: %4.3f [GeV]}",probe_pt/1000.),"");
//leg_left_BM.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe p_{T}}{: %4.3f [GeV]}",abs(probe_mesSA_pt->at(NTrigChain))),"");
//leg_left_BM.AddEntry((TObject*)0,Form("#splitline{RoI #eta}{: %4.3f}",probe_mesSA_roiEta->at(NTrigChain)),"");
//leg_left_BM.AddEntry((TObject*)0,Form("#splitline{Road Slope}{: %4.3f}", -TMath::Log(tan((atan(probe_mesSA_roadAw->at(NTrigChain)[1]))/2.))),"");
//leg_left_BM.AddEntry((TObject*)0,Form("#splitline{Road Slope}{: %4.3f}", probe_mesSA_roadAw->at(NTrigChain)[1]),"");
//leg_left_BM.AddEntry((TObject*)0,Form("#splitline{Road Z-Intercept}{: %4.3f}",probe_mesSA_roadBw->at(NTrigChain)[1]),"");
//leg_left_BM.Draw();
//
//gr_RPC.Draw("P, same");
//
//gr_MdtHit_Outlier_BM.SetMarkerColor(kRed);
//gr_MdtHit_Outlier_BM.SetMarkerStyle(24);
//gr_MdtHit_Outlier_BM.SetMarkerSize(1);
//gr_MdtHit_Outlier_BM.Draw("P,same");
//
//gr_MdtHit_Outlier2_BM.SetMarkerColor(12);
//gr_MdtHit_Outlier2_BM.SetMarkerStyle(24);
//gr_MdtHit_Outlier2_BM.SetMarkerSize(1);
//gr_MdtHit_Outlier2_BM.Draw("P,same");
//gr_SP.Draw("P, same");
//eventInfo.Draw();
//c2->Print(pdf, "pdf");
//c2->RedrawAxis();
//delete frame_BM;
//
//// Set Legend for BO
//TLegend leg_BO = TLegend(0.81,0.22,0.99,0.95);
//leg_BO.SetHeader(Form("#splitline{(4) Barrel Outer}{%s}", label_for_sector.Data()),"C");
//leg_BO.SetTextSize(0.035);
//leg_BO.AddEntry(&gr_segment,"#splitline{Offline}{segment}","p");
//leg_BO.AddEntry(&gr_MdtHit_Inlier_BO,  Form("#splitline{Inlier}{MDT hit (%d)}",nBO_Inlier),  "p");
//leg_BO.AddEntry(&gr_MdtHit_Outlier_BO, Form("#splitline{Outlier}{MDT hit (%d)}",nBO_Outlier), "p");
//leg_BO.AddEntry(&gr_RPC,"RPC hit","p");
//leg_BO.AddEntry(&gr_SP,"SuperPoint","p");
//leg_BO.AddEntry(&f_road_BO,"Road","l");
//leg_BO.AddEntry(&f_roi,"RoI","l");
//TH1* frame_BO = c2->DrawFrame(Zmin_BO,Rmin_BO,Zmax_BO,Rmax_BO);
//frame_BO->GetXaxis()->SetTitle("Z [m]");
//frame_BO->GetYaxis()->SetTitle("R [m]");
//f_road_BO.Draw("same");
//f_road_BO_plus.Draw("same");
//f_road_BO_minus.Draw("same");
//f_roi.Draw("same");
//f_roi_min.Draw("same");
//f_roi_max.Draw("same");
//f_road_BO.GetXaxis()->SetLimits(Zmin_BO,Zmax_BO);
//f_road_BO.GetYaxis()->SetRangeUser(Rmin_BO,Rmax_BO);
//gr_MdtHit_Inlier_BO.SetMarkerColor(kGreen + 2);
//gr_MdtHit_Inlier_BO.SetMarkerStyle(24);
//gr_MdtHit_Inlier_BO.SetMarkerSize(1);
//gr_MdtHit_Inlier_BO.Draw("P, same");
//gr_segment.Draw("P, same");
//
//leg_BO.Draw();
//// Set Legend
//TLegend leg_left_BO = TLegend(-0.027,0.12,0.10,0.8);
//leg_left_BO.SetTextSize(0.03);
//// Display passe or not
//TLegendEntry* l_pass_BO;
//if (probe_mesSA_pass->at(NTrigChain) == 1){
//  l_pass_BO = leg_left_BO.AddEntry((TObject*)0,"#splitline{Passed}{in L2MuonSA}","");
//  l_pass_BO -> SetTextColor(kGreen+2);
//} else{
//  l_pass_BO = leg_left_BO.AddEntry((TObject*)0,"#splitline{NOT Passed}{in L2MuonSA}","");
//  l_pass_BO -> SetTextColor(kRed);
//}
//// Display passe or not
//TLegendEntry* l_rpc_BO;
//if (probe_mesSA_isRpcFailure->at(NTrigChain) == 0){
//  l_rpc_BO = leg_left_BO.AddEntry((TObject*)0,"RPC Fit Success","");
//  l_rpc_BO -> SetTextColor(kGreen+2);
//} else{
//  l_rpc_BO = leg_left_BO.AddEntry((TObject*)0,"RPC Fit Failure","");
//  l_rpc_BO -> SetTextColor(kRed);
//}
//
//leg_left_BO.AddEntry((TObject*)0,Form("#splitline{Offline probe p_{T}}{: %4.3f [GeV]}",probe_pt/1000.),"");
//leg_left_BO.AddEntry((TObject*)0,Form("#splitline{L2MuonSA probe p_{T}}{: %4.3f [GeV]}",abs(probe_mesSA_pt->at(NTrigChain))),"");
//leg_left_BO.AddEntry((TObject*)0,Form("#splitline{RoI #eta}{: %4.3f}",probe_mesSA_roiEta->at(NTrigChain)),"");
//leg_left_BO.AddEntry((TObject*)0,Form("#splitline{Road Slope}{: %4.3f}", -TMath::Log(tan((atan(probe_mesSA_roadAw->at(NTrigChain)[2]))/2.))),"");
//leg_left_BO.AddEntry((TObject*)0,Form("#splitline{Road Slope}{: %4.3f}", probe_mesSA_roadAw->at(NTrigChain)[2]),"");
//leg_left_BO.AddEntry((TObject*)0,Form("#splitline{Road Z-Intercept}{: %4.3f}",probe_mesSA_roadBw->at(NTrigChain)[2]),"");
//leg_left_BO.Draw();
//
//gr_RPC.Draw("P, same");
//
//gr_MdtHit_Outlier_BO.SetMarkerColor(kRed);
//gr_MdtHit_Outlier_BO.SetMarkerStyle(24);
//gr_MdtHit_Outlier_BO.SetMarkerSize(1);
//gr_MdtHit_Outlier_BO.Draw("P,same");
//
//gr_MdtHit_Outlier2_BO.SetMarkerColor(12);
//gr_MdtHit_Outlier2_BO.SetMarkerStyle(24);
//gr_MdtHit_Outlier2_BO.SetMarkerSize(1);
//gr_MdtHit_Outlier2_BO.Draw("P,same");
//
//gr_SP.Draw("P, same");
//eventInfo.Draw();
//c2->Print(pdf, "pdf");
//c2->RedrawAxis();
//delete frame_BO;
//
//current_entry += 1;
//cout << "===" << begin_entry << ": " << current_entry << ": " << limit_entry << endl;
//if (current_entry > limit_entry) {
//  cout << "=================" << endl;
//  break;
//}
//c2->Clear();
//} // end of each entry
//c2->Print(pdf + "]", "pdf");
//delete c2;
//}

bool GRLlist(int LumiBlock){
  bool lb_1 = (LumiBlock > 110 && LumiBlock < 124);
  bool lb_2= (LumiBlock > 126 && LumiBlock < 130);
  bool lb_3= (LumiBlock > 169 && LumiBlock < 239);

  return (lb_1 || lb_2 || lb_3);
}

void RPC::FillMdtHist(){ // Only eta region is applied in this function, not tag_proc, pass, ....
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
  const double ZERO_LIMIT = 1e-5;

  double qeta = probe_eta*probe_charge;
  //bool isBarrel = true;
  bool isQetaCut = ((qeta > -0.4) && (qeta < -0.0));
  //bool isBarrel = (abs(probe_eta) < 1.05) && (isQetaCut);
  //bool isBarrel = true;
  bool isBarrel = (abs(probe_eta) < 1.05);

  LOGD << "size: " << (probe_mesL2_pt->at(N4)).size();
  for ( int i = 0; i < (probe_mesL2_pass->at(N4)).size(); i++ ){
    if ( probe_mesL2_sAddressSA->at(N4)[i] > -1  && probe_mesL2_pass->at(N4)[i] > -3){
      m_ftk_size_pt_mu4 -> Fill( probe_pt/1000., (probe_mesL2_pt->at(N4)).size() );

      m_ftk_size_pt_mu50 -> Fill( probe_pt/1000., (probe_mesL2_pt->at(14)).size() );

    }
  }

  //if ( !(qeta > -1.3 && qeta < -0.9)){
  //  return;
  //}

  // Special For FTK
  if ( tag_proc == 3 ){
    //if (!(sumReqdRL1<tp_extdR && 0.2<tp_extdR)){
    //  ;
    //}
    LOGD << "size: " << (probe_mesL2_pass->at(N4)).size();
    LOGD << "l2sa : pass/pt/eta/phi/sAddress/roi: " << probe_mesSA_pass->at(N4) << "/" << probe_mesSA_pt->at(N4) << "/" << probe_mesSA_eta->at(N4) << "/" << probe_mesSA_phi->at(N4)<< "/" << probe_mesSA_sAddress->at(N4);
    for ( int i = 0; i < (probe_mesL2_pass->at(N4)).size(); i++ ){
      if ( probe_mesL2_sAddressSA->at(N4)[i] > -1){
        LOGD << "l2sa Barrel: pass/pt/eta/phi/sAddress: " << probe_mesL2_pass->at(N4)[i] << "/" << probe_mesL2_ptSA->at(N4)[i] << "/" << probe_mesL2_etaSA->at(N4)[i] << "/" << probe_mesL2_phiSA->at(N4)[i] << "/" << probe_mesL2_sAddressSA->at(N4)[i];
        LOGD << "l2sa Barrel: ptFtk/etaFtk/phiFtk/roadAlgo: " << probe_mesL2_ptFtk->at(N4)[i] << "/" << probe_mesL2_etaFtk->at(N4)[i] << "/" << probe_mesL2_phiFtk->at(N4)[i] << "/" << probe_mesL2_roadAlgo->at(N4)[i];
      } else {
        LOGD << "l2sa Endcap: pass/pt/eta/phi/sAddress: " << probe_mesL2_pass->at(N4)[i] << "/" << probe_mesL2_ptSA->at(N4)[i] << "/" << probe_mesL2_etaSA->at(N4)[i] << "/" << probe_mesL2_phiSA->at(N4)[i] << "/" << probe_mesL2_sAddressSA->at(N4)[i];
        LOGD << "l2sa Endcap: ptFtk/etaFtk/phiFtk/roadAlgo: " << probe_mesL2_ptFtk->at(N4)[i] << "/" << probe_mesL2_etaFtk->at(N4)[i] << "/" << probe_mesL2_phiFtk->at(N4)[i] << "/" << probe_mesL2_roadAlgo->at(N4)[i];
      }
    }
    //LOGD << "vec size" << pass.size();
    //LOGD << "pass" << pass;



    //CB threshold 
    //'4GeV_v15a'             : [ [0,1.05,1.5,2.0,9.9], [  3.86,  3.77,  3.69,  3.70] ],
    //'5GeV_v15a'             : [ [0,1.05,1.5,2.0,9.9], [  4.9,  4.8,  4.8,  4.8] ],
    //'6GeV_v15a'             : [ [0,1.05,1.5,2.0,9.9], [  5.87,  5.79,  5.70,  5.62] ],
    double CBthre = 0;
    if ( N4 == 0 ){
      CBthre = 3.86;
    } else if ( N4 == 1){
      CBthre = 5.87;
    }

    if ( abs(probe_eta) < 1.05 ){
      LOGD << "barrel";
      m_probe_pt_mu4[0]->Fill(probe_pt/1000.); // PROBE
      if ( probe_mesL1_pass->at(N4) > -1){
        LOGD << "L1Pass";
        m_probe_pt_mu4[1]->Fill(probe_pt/1000.); // L1
        m_probe_qetapt_mu4[1]->Fill(probe_charge*probe_eta, probe_pt/1000.); // L1
        if ( CBthre < probe_pt/1000. ){
          m_probe_eta_mu4[1]->Fill(probe_eta); // L1
          m_probe_phi_mu4[1]->Fill(probe_phi); // L1
          m_probe_etaphi_mu4[1]->Fill(probe_eta, probe_phi); // L1
        }
        // SA and CB
        LOGD << "pass: " << probe_mesSA_pass->at(N4);
        LOGD << "roadAlgo: " << probe_mesSA_roadAlgo->at(N4);
        LOGD << "ptOff: " << probe_pt/1000.;
        LOGD << "ptSA/ptFtk: " << abs(probe_mesSA_pt->at(N4)) << "/" << probe_mesSA_ptFtk->at(N4)/1000.;
        LOGD << "NumberOfSP(): " <<  NumberOfSP();
        bool isSA_cut = false;
        //bool isSA_cut = (probe_mesSA_roadAlgo->at(N4)==1 && (NumberOfSP()>0 && (probe_mesSA_ptFtk->at(N4))/1000. > 3.38));
        //bool isSA_cut = probe_mesSA_pass->at(N4)>-1;
        //bool isSA_cut = false;
        //bool isSA_cut = true;
        //if ( probe_mesSA_roadAlgo->at(N4)==1 && NumberOfSP()>0 && (probe_mesSA_ptFtk->at(N4))/1000. > 3.38 ) 


        // L2SA threshold mu4 or mu6 => 3.38 or 5.17
        if ( probe_mesSA_roadAlgo->at(N4)==1 && abs(probe_mesSA_pt->at(N4)) > 3.38 ) {
          LOGD << "pass SA as FTK";
          isSA_cut = true;
         } else if ( probe_mesSA_roadAlgo->at(N4)==2 && abs(probe_mesSA_pt->at(N4)) > 3.38 ) {
           LOGD << "pass SA as RPC";
          isSA_cut = true;
         } else if ( abs(probe_mesSA_pt->at(N4)) > 3.38 ) {
           LOGD << "pass SA in any way";
          isSA_cut = true;
         } else {
           LOGD << "nopass SA";
         }

        if (probe_mesL1_pass->at(N4) > -1 && probe_mesSA_pass->at(N4) > 0){
          LOGD << "SAPass";
          m_probe_pt_mu4[2]->Fill(probe_pt/1000.); // SA
          if (probe_pt/1000. < 5){
            LOGD << "mu4cb";
          }
          if ( probe_mesCB_pass->at(N4) > -1 ){
            LOGD << "CBPass, " << "ptCB: " << abs(probe_mesCB_pt->at(N4)/1000.);
            m_probe_pt_mu4[3]->Fill(probe_pt/1000.); // CB
          } else {
            LOGD << "noCBPass, " << "ptCB: " << abs(probe_mesCB_pt->at(N4)/1000.);
          }
        }

        // L2
        int SApass = -1;
        int CBThrepass = -1;
        int L2pass = -1;
        LOGD << "chain: " << N4 << "--> thre: " << CBthre;
        for ( int i = 0; i < (probe_mesL2_pass->at(N4)).size(); i++ ){
          //double L2pt = abs(probe_mesL2_ptFtk->at(N4)[i]);
          double SApt       = abs(probe_mesL2_ptSA->at(N4)[i]);
          double L2ptFtk    = abs(probe_mesL2_ptFtk->at(N4)[i])/1000.;
          double L2roadAlgo = probe_mesL2_roadAlgo->at(N4)[i];
          LOGD << "";
          LOGD << "roiNumber/SApt/SAeta/SAphi: " << probe_mesL2_passSA->at(N4)[i] << "/" << abs(probe_mesL2_ptSA->at(N4)[i]) << "/" << probe_mesL2_etaSA->at(N4)[i] << "/" << probe_mesL2_phiSA->at(N4)[i];
          LOGD << "L2ptFtk/L2etaFtk/L2phiFtk: " << abs(probe_mesL2_ptFtk->at(N4)[i])/1000. << "/" << probe_mesL2_etaFtk->at(N4)[i] << "/" << probe_mesL2_phiFtk->at(N4)[i];
          LOGD << "L2roadAlgo: " << L2roadAlgo;
          LOGD << "L2pass: " << probe_mesL2_pass->at(N4)[i];
          if ( probe_mesL2_roadAlgo->at(N4)[i] == 1 ){ // FTK road
            LOGD << "FTK L2Muon";
            if ( abs(probe_mesL2_ptSA->at(N4)[i]) > CBthre ){
              L2pass = 1;
              LOGD << "L2pass as FTK L2Muon";
            } else {
              LOGD << "not L2pass as FTK L2Muon";
            }
          } else if ( probe_mesL2_roadAlgo->at(N4)[i] == 2) { // RPC road
              LOGD << "RPC L2MuonSA --> l2muComb";
          } else if ( probe_mesL2_roadAlgo->at(N4)[i] == -2) { // muComb from RPC road
            LOGD << "L2muComb";
            if ( probe_mesL2_pass->at(N4)[i] > 0 ){
              //L2pass = 1;
              LOGD << "L2pass as L2muComb";
            } else {
              LOGD << "not L2pass as L2muComb";
            }
          } else if ( probe_mesL2_roadAlgo->at(N4)[i] == 0) { // Endcap
            LOGD << "Endcap L2MuonSA";
          } else if ( probe_mesL2_roadAlgo->at(N4)[i] == -3) { // No FTK sample i.e. No dynamic variable of roadAlgo
            LOGD << "No FTK sample L2MuonSA";
          }

          //if ( L2pt > CBthre || probe_mesCB_pass->at(N4) > -1 ){ // CB thre of HLT_mu4 in barrel is 3.86 GeV
          //  CBThrepass = 1;
          //}
        }

        LOGD << "";
        LOGD << "L2pass result: " << L2pass;
        //if ( SApass == 1 && CBThrepass == 1) 
        if ( L2pass > -1) {
          m_probe_pt_mu4[4]->Fill(probe_pt/1000.);
          m_probe_qetapt_mu4[4]->Fill(probe_charge*probe_eta, probe_pt/1000.);
          if ( CBthre < probe_pt/1000. ){
            m_probe_eta_mu4[4]->Fill(probe_eta); // L1
            m_probe_phi_mu4[4]->Fill(probe_phi); // L1
            m_probe_etaphi_mu4[4]->Fill(probe_eta, probe_phi);
          }
          LOGD << "m_probe_pt_mu4[4]->Fill()";
        } else {
          LOGD << "not m_probe_pt_mu4[4]->Fill()";
        }
        LOGD << "SApass: " << SApass << ", probe_mesSA_pass: " << probe_mesSA_pass->at(N4);
        //LOGD << "CBThrepass: " << CBThrepass;
      }

      //if(isBarrel) h_probe_pt_mu4_PROBE -> Fill(probe_pt/1000.);
    }
  }

  // mu cut
  //if (AverageInteractionsPerCrossing < 30){return;}

  switch (tag_proc) {
    case 3: //Jpsi until L2
      // Check TAP
      //if (!(probe_mesEFTAG_pass -> at(N4) > -1 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
      //  return;
      //}

      // offline distribution
      if(isBarrel) h_probe_pt_mu4_PROBE -> Fill(probe_pt/1000.);
      // Fill L1 probe hists
      if (probe_mesL1_pass -> at(N4) > -1){
        if(isBarrel) h_probe_pt_mu4_L1 -> Fill(probe_pt/1000.);
        hh_probe_qetapt_mu4_L1 -> Fill(probe_charge*probe_eta, probe_pt/1000.);
        // Plateau cut
        if (probe_pt/1000. > 3.38) {
          if(isBarrel) h_probe_mu_mu4_L1 -> Fill(AverageInteractionsPerCrossing);
          h_probe_eta_mu4_L1 -> Fill(probe_eta);
          if(isBarrel) h_probe_phi_mu4_L1 -> Fill(probe_phi);
          hh_probe_etaphi_mu4_L1 -> Fill(probe_eta, probe_phi);
        }
        // for LUT-SP check
        //if(abs(probe_mesSA_pt->at(N4)) < ZERO_LIMIT && NumberOfSP() > 1.5) hh_probe_etaphi_mu4_SA -> Fill(probe_mesSA_eta->at(N4), probe_mesSA_phi->at(N4));
      }
      // Fill SA probe hists

      if (probe_mesL1_pass -> at(N4) > -1 && probe_mesSA_pass -> at(N4) > -1){
        if(isBarrel) h_probe_pt_mu4_SA -> Fill(probe_pt/1000.);
        hh_probe_qetapt_mu4_SA -> Fill(probe_charge*probe_eta, probe_pt/1000.);
        // Plateau cut
        if (probe_pt/1000. > 3.38) {
          if(isBarrel) h_probe_mu_mu4_SA -> Fill(AverageInteractionsPerCrossing);
          h_probe_eta_mu4_SA -> Fill(probe_eta);
          if(isBarrel) h_probe_phi_mu4_SA -> Fill(probe_phi);
          hh_probe_etaphi_mu4_SA -> Fill(probe_eta, probe_phi);
        }
      }
      // Fill CB probe hists
      if (probe_mesL1_pass -> at(N4) > -1 && probe_mesSA_pass -> at(N4) > -1 && probe_mesCB_pass -> at(N4) > -1){
        if(isBarrel) h_probe_pt_mu4_CB -> Fill(probe_pt/1000.);
      }
      // Fill EF probe hists
      if (probe_mesL1_pass -> at(N4) > -1 && probe_mesSA_pass -> at(N4) > -1 && probe_mesCB_pass -> at(N4) > -1 && probe_mesEF_pass -> at(N4) > -1){
        if(isBarrel) h_probe_pt_mu4_EF -> Fill(probe_pt/1000.);
      }

      //FillPtResidualHist();

      break;
    case 2: //Jpsi from L2:
      break;
    case 1: //Z
      // Check TAP

      if (!(probe_mesEFTAG_pass -> at(N50) > -1 && ( sumReqdRL1<tp_extdR && 0.2<tp_extdR ) && ( sumReqdREF<tp_dR ))){
        return;
      }
      // Fill L1 probe hists
      if (probe_mesL1_pass -> at(N50) > -1){
        if(isBarrel) h_probe_pt_mu50_L1 -> Fill(probe_pt/1000.);
        hh_probe_qetapt_mu50_L1 -> Fill(probe_charge*probe_eta, probe_pt/1000.);
        // Plateau cut
        if (probe_pt/1000. > 8.0) {
          if(isBarrel) h_probe_mu_mu50_L1 -> Fill(AverageInteractionsPerCrossing);
          h_probe_eta_mu50_L1 -> Fill(probe_eta);
          if(isBarrel) h_probe_phi_mu50_L1 -> Fill(probe_phi);
          hh_probe_etaphi_mu50_L1 -> Fill(probe_eta, probe_phi);
        }
      }
      // Fill SA probe hists
      if (probe_mesL1_pass -> at(N50) > -1 && probe_mesSA_pass -> at(N50) > -1){
        if(isBarrel) h_probe_pt_mu50_SA -> Fill(probe_pt/1000.);
        hh_probe_qetapt_mu50_SA -> Fill(probe_charge*probe_eta, probe_pt/1000.);
        // Plateau cut
        if (probe_pt/1000. > 8.0) {
          if(isBarrel) h_probe_mu_mu50_SA -> Fill(AverageInteractionsPerCrossing);
          h_probe_eta_mu50_SA -> Fill(probe_eta);
          if(isBarrel) h_probe_phi_mu50_SA -> Fill(probe_phi);
          hh_probe_etaphi_mu50_SA -> Fill(probe_eta, probe_phi);
        }
      }
      //cout << probe_mesSA_pass -> at(N50) << ", " << probe_mesCB_pass -> at(N50) << endl;

      // Fill CB probe hists
      if (probe_mesL1_pass -> at(N50) > -1 && probe_mesSA_pass -> at(N50) > -1 && probe_mesCB_pass -> at(N50) > -1){
        h_probe_pt_mu50_CB -> Fill(probe_pt/1000.);
        if (probe_pt/1000. > 8.0) {
          if(isBarrel) h_probe_mu_mu50_CB -> Fill(AverageInteractionsPerCrossing);
          if(isBarrel) h_probe_phi_mu50_CB -> Fill(probe_phi);
          h_probe_eta_mu50_CB -> Fill(probe_eta);
        }
      }

      break;
  }
  return;
}

void RPC::CalcEff(){

  CalcHistToHist( m_probe_pt_mu4[1], m_probe_pt_mu4[0], m_eff_pt_mu4_L1);
  CalcHistToHist( m_probe_pt_mu4[2], m_probe_pt_mu4[1], m_eff_pt_mu4_L1SA);
  CalcHistToHist( m_probe_pt_mu4[3], m_probe_pt_mu4[2], m_eff_pt_mu4_SACB);
  CalcHistToHist( m_probe_pt_mu4[4], m_probe_pt_mu4[1], m_eff_pt_mu4_L1L2);
  CalcHistToHist( m_probe_pt_mu4[3], m_probe_pt_mu4[1], m_eff_pt_mu4_L1CB);

  CalcHistToHist( m_probe_eta_mu4[4], m_probe_eta_mu4[1], m_eff_eta_mu4_L1L2);
  CalcHistToHist( m_probe_phi_mu4[4], m_probe_phi_mu4[1], m_eff_phi_mu4_L1L2);
  CalcHistToHist( m_probe_qetapt_mu4[4], m_probe_qetapt_mu4[1], m_eff_qetapt_mu4_L1L2);
  CalcHistToHist( m_probe_etaphi_mu4[4], m_probe_etaphi_mu4[1], m_eff_etaphi_mu4_L1L2);

  // mu4
  CalcHistToHist( h_probe_mu_mu4_SA,       h_probe_mu_mu4_L1,       h_eff_mu_mu4_L1SA);
  CalcHistToHist( h_probe_pt_mu4_L1,       h_probe_pt_mu4_PROBE,    h_eff_pt_mu4_L1);
  CalcHistToHist( h_probe_pt_mu4_SA,       h_probe_pt_mu4_L1,       h_eff_pt_mu4_L1SA);
  CalcHistToHist( h_probe_pt_mu4_CB,       h_probe_pt_mu4_SA,       h_eff_pt_mu4_SACB);
  CalcHistToHist( h_probe_pt_mu4_EF,       h_probe_pt_mu4_CB,       h_eff_pt_mu4_CBEF);
  CalcHistToHist( h_probe_eta_mu4_SA,      h_probe_eta_mu4_L1,      h_eff_eta_mu4_L1SA);
  CalcHistToHist( h_probe_phi_mu4_SA,      h_probe_phi_mu4_L1,      h_eff_phi_mu4_L1SA);
  CalcHistToHist( hh_probe_qetapt_mu4_SA,  hh_probe_qetapt_mu4_L1,  hh_eff_qetapt_mu4_L1SA);
  CalcHistToHist( hh_probe_etaphi_mu4_SA,  hh_probe_etaphi_mu4_L1,  hh_eff_etaphi_mu4_L1SA);

  // mu50
  CalcHistToHist( h_probe_mu_mu50_SA,      h_probe_mu_mu50_L1,      h_eff_mu_mu50_L1SA);
  CalcHistToHist( h_probe_pt_mu50_SA,      h_probe_pt_mu50_L1,      h_eff_pt_mu50_L1SA);
  CalcHistToHist( h_probe_eta_mu50_SA,     h_probe_eta_mu50_L1,     h_eff_eta_mu50_L1SA);
  CalcHistToHist( h_probe_phi_mu50_SA,     h_probe_phi_mu50_L1,     h_eff_phi_mu50_L1SA);
  CalcHistToHist( h_probe_mu_mu50_CB,      h_probe_mu_mu50_SA,      h_eff_mu_mu50_SACB);
  CalcHistToHist( h_probe_pt_mu50_CB,      h_probe_pt_mu50_SA,      h_eff_pt_mu50_SACB);
  CalcHistToHist( h_probe_eta_mu50_CB,     h_probe_eta_mu50_SA,     h_eff_eta_mu50_SACB);
  CalcHistToHist( h_probe_phi_mu50_CB,     h_probe_phi_mu50_SA,     h_eff_phi_mu50_SACB);
  CalcHistToHist( hh_probe_qetapt_mu50_SA, hh_probe_qetapt_mu50_L1, hh_eff_qetapt_mu50_L1SA);
  CalcHistToHist( hh_probe_etaphi_mu50_SA, hh_probe_etaphi_mu50_L1, hh_eff_etaphi_mu50_L1SA);
}

void RPC::DrawEffHist(TString pdf){
  LOGI << "---start---";
  //==================================================================
  //Set Canvas
  //==================================================================
  TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 1020, 700);
  c1->SetGrid();
  c1->SetRightMargin(0.20);
  c1->SetLeftMargin(0.23);
  c1->SetBottomMargin(0.20);

  c1 -> Print( pdf + "[", "pdf" );

  m_ftk_size_pt_mu4->Draw("colz");
  c1 -> Print( pdf, "pdf" );

  TProfile *pf4 = m_ftk_size_pt_mu4->ProfileX();
  pf4->Draw();
  pf4->GetYaxis()->SetTitle("Average number of FTKs");
  c1 -> Print(pdf, "pdf" );

  //m_ftk_size_eta_mu4->Draw("colz");
  //c1 -> Print( pdf, "pdf" );

  m_ftk_size_pt_mu50->Draw("colz");
  c1 -> Print( pdf, "pdf" );

  TProfile *pf50 = m_ftk_size_pt_mu50->ProfileX();
  pf50->Draw();
  pf50->GetYaxis()->SetTitle("Average number of FTKs");
  c1 -> Print(pdf, "pdf" );

  //m_ftk_size_eta_mu50->Draw("colz");
  //c1 -> Print( pdf, "pdf" );

  m_eff_pt_mu4_L1->Draw();
  m_eff_pt_mu4_L1->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  m_eff_pt_mu4_L1SA->Draw();
  m_eff_pt_mu4_L1SA->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  m_eff_pt_mu4_SACB->Draw();
  m_eff_pt_mu4_SACB->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  m_eff_pt_mu4_L1L2->Draw();
  m_eff_pt_mu4_L1L2->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  m_eff_pt_mu4_L1CB->Draw();
  m_eff_pt_mu4_L1CB->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  m_eff_pt_mu4_L1->Draw();
  m_eff_pt_mu4_L1->SetMarkerColor(kRed);
  m_eff_pt_mu4_L1->SetLineColor(kRed);
  m_eff_pt_mu4_L1SA->Draw("same");
  m_eff_pt_mu4_L1SA->SetMarkerColor(kYellow+2);
  m_eff_pt_mu4_L1SA->SetLineColor(kYellow+2);
  m_eff_pt_mu4_SACB->Draw("same");
  m_eff_pt_mu4_SACB->SetMarkerColor(kPink+2);
  m_eff_pt_mu4_SACB->SetLineColor(kPink+2);
  m_eff_pt_mu4_L1L2->Draw("same");
  m_eff_pt_mu4_L1L2->SetMarkerColor(kBlue);
  m_eff_pt_mu4_L1L2->SetLineColor(kBlue);
  m_eff_pt_mu4_L1CB->Draw("same");
  m_eff_pt_mu4_L1CB->SetMarkerColor(kCyan+2);
  m_eff_pt_mu4_L1CB->SetLineColor(kCyan+2);
  c1 -> Print( pdf, "pdf" );

  m_eff_eta_mu4_L1L2->Draw();
  m_eff_eta_mu4_L1L2->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  m_eff_phi_mu4_L1L2->Draw();
  m_eff_phi_mu4_L1L2->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  m_eff_qetapt_mu4_L1L2->Draw("colz");
  m_eff_qetapt_mu4_L1L2->GetZaxis()->SetRangeUser(0,1.0);
  c1 -> Print( pdf, "pdf" );

  m_eff_etaphi_mu4_L1L2->Draw("colz");
  m_eff_etaphi_mu4_L1L2->GetZaxis()->SetRangeUser(0,1.0);
  c1 -> Print( pdf, "pdf" );


  h_eff_mu_mu4_L1SA->Draw();
  h_eff_mu_mu4_L1SA->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

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

  h_eff_pt_mu4_SACB->Draw();
  h_eff_pt_mu4_SACB->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_pt_mu4_CBEF->Draw();
  h_eff_pt_mu4_CBEF->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_mu_mu50_L1SA->Draw();
  h_eff_mu_mu50_L1SA->GetYaxis()->SetRangeUser(0,1.1);
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

  h_eff_mu_mu50_SACB->Draw();
  h_eff_mu_mu50_SACB->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_pt_mu50_SACB->Draw();
  h_eff_pt_mu50_SACB->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_eta_mu50_SACB->Draw();
  h_eff_eta_mu50_SACB->GetYaxis()->SetRangeUser(0,1.1);
  c1 -> Print( pdf, "pdf" );

  h_eff_phi_mu50_SACB->Draw();
  h_eff_phi_mu50_SACB->GetYaxis()->SetRangeUser(0,1.1);
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

double RPC::calc_pTresidual( double offline_pt, double trig_pt){
  return ( 1. - (offline_pt) / (trig_pt));
}

void RPC::FillPtResidualHist(){ // Only eta region cut here
  //if ( ((probe_mesSA_pt -> at(N4)) == 0) && ((probe_mesSA_pt -> at(N4)) > -9999) ){ return;}
  //if ( ((probe_mesSA_pt -> at(N4)) == 0) && ((probe_mesSA_pt -> at(N4)) > -9999) ){ return;}

  double pTresidual[3];
  pTresidual[0] = 9999;
  pTresidual[1] = 9999;
  pTresidual[2] = 9999;
  double pTresidualFtk = 9999;
  double pTresidualLut = 9999;
  pTresidual[0] = calc_pTresidual(probe_pt/1000., abs(probe_mesSA_pt -> at(N4)));
  //cout << probe_mesSA_pass -> at(N4) << ": " << (probe_mesSA_ptFtk->at(N4))/1000. << ":" << probe_pt/1000. << ": " << abs(probe_mesSA_pt->at(N4)) << ": "  << pTresidual << endl;
  //if (abs(probe_mesSA_ptFtk->at(N4))>0.)
  if (probe_mesSA_ptFtk->size()>0. && probe_mesSA_ptFtk->at(N4) && probe_mesSA_roadAlgo->at(N4) == 1){
    pTresidualFtk = calc_pTresidual(probe_pt/1000., abs(probe_mesSA_ptFtk -> at(N4))/1000.);
    pTresidual[1] = calc_pTresidual(probe_pt/1000., abs(probe_mesSA_pt -> at(N4)));
    LOGD << "ptFtk is used as L2SA pt: residual = " << pTresidual;
  } else if (probe_mesSA_ptFtk->size()>0. && probe_mesSA_ptFtk->at(N4) && probe_mesSA_roadAlgo->at(N4) == 2) {
    pTresidualLut = calc_pTresidual(probe_pt/1000., abs(probe_mesSA_pt -> at(N4)));
    pTresidual[2] = calc_pTresidual(probe_pt/1000., abs(probe_mesSA_pt -> at(N4)));
    LOGD << "ptFtk is not available";
  }

  LOGD << "pTresidual: " << pTresidual[0] << "/" << pTresidual[1] << "/" << pTresidual[2] << endl;

  if (probe_mesSA_ptFtk->size()>0.){
    LOGD << "probe_mesSA_ptFtk.size() > 0";
    LOGD << "pass/ptFtk/probe_pt/mesSA_pt/residual: " << probe_mesSA_pass -> at(N4) << "/" << (probe_mesSA_ptFtk->at(N4))/1000. << "/" << probe_pt/1000. << "/" << abs(probe_mesSA_pt->at(N4)) << ":"  << pTresidual[0];
    LOGD << "roadAlgo: " << probe_mesSA_roadAlgo -> at(N4) << ": " << probe_mesSA_roadAlgo->at(N4);
  } else {
    LOGD << "probe_mesSA_ptFtk.size() = 0";
  }


  if (abs(probe_eta) < 1.05){
    h_PtResidual_pt1 -> Fill(probe_pt/1000., pTresidual[0]);
    h_PtResidual_pt2 -> Fill(probe_pt/1000., pTresidual[1]);
    h_PtResidual_pt3 -> Fill(probe_pt/1000., pTresidual[2]);
    h_PtResidual_ptFtk -> Fill(probe_pt/1000., pTresidualFtk);
    h_PtResidual_ptLut -> Fill(probe_pt/1000., pTresidualLut);
    h_pt_vs_pt -> Fill(probe_pt/1000., abs(probe_mesSA_pt->at(N4)));

    if (probe_mesSA_ptFtk->size()>0. && probe_mesSA_ptFtk->at(N4) && probe_mesSA_roadAlgo->at(N4) == 1){
      h_ptl2_vs_ptFtk -> Fill(abs(probe_mesSA_pt->at(N4)), abs(probe_mesSA_ptFtk->at(N4))/1000.);
      h_pt_vs_ptFtk -> Fill(probe_pt/1000., abs(probe_mesSA_ptFtk->at(N4))/1000.);
    }
  }

  //h_PtResidual_eta -> Fill(probe_eta, pTresidual);
}

void RPC::DrawPtResidualHist(TCanvas* c1, TString pdf){
  // default
  h_PtResidual_pt1 -> Draw("colz");
  c1 -> Print(pdf, "pdf");

  TH1D* tmp_h_pt1 = h_PtResidual_pt1->ProjectionY("pt1", h_PtResidual_pt1->GetXaxis()->FindBin(0.),  h_PtResidual_pt1->GetXaxis()->FindBin(20));
  tmp_h_pt1->Draw();
  c1 -> Print(pdf, "pdf");


  TH1D* tmp_h_pt2 = h_PtResidual_pt2->ProjectionY("pt2", h_PtResidual_pt2->GetXaxis()->FindBin(0.),  h_PtResidual_pt2->GetXaxis()->FindBin(20));
  tmp_h_pt2->SetLineColor(kBlue);
  TH1D* tmp_h_pt3 = h_PtResidual_pt3->ProjectionY("pt3", h_PtResidual_pt3->GetXaxis()->FindBin(0.),  h_PtResidual_pt3->GetXaxis()->FindBin(20));
  tmp_h_pt3->SetLineColor(kRed);
  tmp_h_pt3->Draw();
  tmp_h_pt2->Draw("same");
  TLegend *leg2 = new TLegend(0.2,0.00,0.4,0.15);
  leg2->AddEntry(tmp_h_pt2,"FTK road","l");
  leg2->AddEntry(tmp_h_pt3,"RPC road","l");
  leg2->Draw();
  c1 -> Print(pdf, "pdf");

  // ftk
  h_PtResidual_ptFtk -> Draw("colz");
  c1 -> Print(pdf, "pdf");

  TH1D* tmp_h_ptFtk = h_PtResidual_ptFtk->ProjectionY("pass1", h_PtResidual_ptFtk->GetXaxis()->FindBin(0.),  h_PtResidual_ptFtk->GetXaxis()->FindBin(20));
  tmp_h_ptFtk->Draw();
  c1 -> Print(pdf, "pdf");

  // rpc
  h_PtResidual_ptLut -> Draw("colz");
  c1 -> Print(pdf, "pdf");

  TH1D* tmp_h_ptLut = h_PtResidual_ptLut->ProjectionY("pass2", h_PtResidual_ptLut->GetXaxis()->FindBin(0.),  h_PtResidual_ptLut->GetXaxis()->FindBin(20));
  tmp_h_ptLut->Draw();
  c1 -> Print(pdf, "pdf");

  tmp_h_ptFtk->SetLineColor(kBlue);
  tmp_h_ptFtk->SetMarkerColor(kBlue);
  tmp_h_ptLut->SetLineColor(kRed);
  tmp_h_ptFtk->Draw("");
  tmp_h_ptLut->Draw("same");
  leg2->Draw();
  c1 -> Print(pdf, "pdf");

  h_pt_vs_pt -> Draw();
  h_pt_vs_ptFtk -> Draw("same");
  h_pt_vs_pt->SetMarkerColor(kRed);
  h_pt_vs_ptFtk->SetMarkerColor(kBlue);
  c1 -> Print(pdf, "pdf");

  h_ptl2_vs_ptFtk -> Draw("colz");
  c1 -> Print(pdf, "pdf");

  h_pt_vs_ptFtk -> Draw("colz");
  c1 -> Print(pdf, "pdf");

  h_PtResidual_eta -> Draw("colz");
  c1 -> Print(pdf, "pdf");
}
