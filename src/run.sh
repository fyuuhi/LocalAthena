#!/bin/sh

PDF_LABEL="youhei_Zmumu_AOD_default"
INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/hadd_youhei_Zmumu_AOD_default.root"

#PDF_LABEL="youhei_Zmumu_AOD_noRpcHit"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/hadd_youhei_Zmumu_AOD_noRpcHit.root"

#PDF_LABEL="youhei_Zmumu_AOD_noRpcHitWide"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/hadd_youhei_Zmumu_AOD_noRpcHitWide.root"

#PDF_LABEL="youhei_Zmumu_AOD_outlier2"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/_home_yfukuhar_gpfs_data_youhei_Zmumu_AOD_Zmumu_outlier2_AOD.pool.root/test0605_01.root"

#PDF_LABEL="youhei_Zmumu_AOD_noRpcHitWideRoI"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/_home_youhei_L2MuonSA_dataset_aod_test.noRpcHitWideRoI.pool.root/test0605_01.root"

#PDF_LABEL="youhei_Zmumu_AOD_noRpcHitWideRoI2"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/_home_youhei_L2MuonSA_dataset_aod_test.noRpcHitWideRoI.pool.root.2/test0605_01.root"

#PDF_LABEL="yfukuhar_Zmumu_AOD_noRpcHitWide"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/hadd_yfukuhar_Zmumu_AOD_noRpcHitWide.root"

#PDF_LABEL="yfukuhar_Zmumu_AOD_noRpcHitWide1"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/hadd_yfukuhar_Zmumu_AOD_noRpcHitWide1.root"

#PDF_LABEL="yfukuhar_Zmumu_AOD_noRpcHitWide2"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/hadd_yfukuhar_Zmumu_AOD_noRpcHitWide2.root"

#PDF_LABEL="youhei_Zmumu_AOD_MdtRegion"
#INPUT_NTUPLE="/gpfs/home/yfukuhar/work/CalcEffTool/run/hadd_yfukuhar_Zmumu_AOD_MdtRegion.root"

IS_DRAW="true"
IS_EVENTDISPLAY="false"
BEGIN_ENTRY=1
LIMIT_ENTRY=-1
#LIMIT_ENTRY=10

COMMAND="./RPC ${PDF_LABEL} ${INPUT_NTUPLE} ${IS_DRAW} ${IS_EVENTDISPLAY} ${BEGIN_ENTRY} ${LIMIT_ENTRY}"
LOG="log/log_"${PDF_LABEL}

eval ${COMMAND} 2>&1 | tee ${LOG}

#eof


  #tree1 -> Add("/gpfs/fs2001/yfukuhar/data/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_jpzYFV4_GRL_F_tree_v1_349533_EXT0/hadd/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_jpzYFV4_GRL_F_tree_v1_349533_EXT0.root");
  #tree1 -> Add("/gpfs/fs2001/yfukuhar/data/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v1_EXT0/hadd/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v1_EXT.root");
  #tree1 -> Add("/gpfs/fs2001/yfukuhar/data/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v2_EXT0/hadd/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v2_EXT.root");
  #tree1 -> Add("/gpfs/fs2001/yfukuhar/data/data18_0621/data18_0621.root");
  #tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16_youhei_Zmumu_default/mc16_youhei_Zmumu_default.root");
  #tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16_youhei_Zmumu_noRpcHit/mc16_youhei_Zmumu_noRpcHit.root");
  #tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16_youhei_Zmumu_noRpcHitWide/mc16_youhei_Zmumu_noRpcHitWide.root");
  #tree1 -> Add("/gpfs/home/yfukuhar/work/CalcEffTool/run/Output/");
  #tree1 -> Add("/gpfs/fs2001/yfukuhar/CalcEffPlotMakerOrigin/data/mc16c_Jpsimu6_default/mc16c_Jpsimu6_default.root");
  #tree1 -> Add("/gpfs/fs2001/yfukuhar/data/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v3_EXT0/hadd/user.yfukuhar.00349533.physics_Main.YFTAP.f929_m1955_GRL_False_349533_mdtHit_v3_EXT0.root");
