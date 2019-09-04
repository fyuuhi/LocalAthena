
### triger flags config ###
#from RecExConfig.RecFlags import rec
#rec.doWriteAOD=True
#rec.doAOD=True
#rec.doFTK=True
#
#TriggerFlags.doFTK=True
#
#TriggerFlags.readHLTconfigFromXML=False
#TriggerFlags.readLVL1configFromXML=False
#TriggerFlags.triggerMenuSetup="MC_pp_v7" 

#from RecExConfig.RecFlags import rec
#from AthenaCommon.AthenaCommonFlags import athenaCommonFlags as acf
#from InDetRecExample.InDetJobProperties import InDetFlags
#from TriggerJobOpts.TriggerFlags import TriggerFlags
#from JetRec.JetRecFlags import jetFlags,JetContentDetail
#
#InDetFlags.doSlimming.set_Value_and_Lock(False)
#TriggerFlags.AODEDMSet.set_Value_and_Lock("AODFULL");
#jetFlags.detailLevel=JetContentDetail.Full
#rec.doFTK.set_Value_and_Lock(True)

### run config ###
isGrid = False
#isGrid = True

AllDebug = True
#AllDebug = False

#SpecificDebug = True
SpecificDebug = False


### import ###
from TrigL2MuonSA.TrigL2MuonSAConf import *
from TrigL2MuonSA.TrigL2MuonSAConfig import TrigL2MuonSAConfig
from TrigmuComb.TrigmuCombConf import *
from TrigmuComb.TrigmuCombConfig import TrigmuCombConfig
from TrigMuonHypo.TrigMuonHypoConfig import *
from TrigMuonHypo.TrigL2MuonOverlapRemoverConfig import *
from AthenaCommon.AppMgr import ServiceMgr
from AthenaCommon.AppMgr import ToolSvc

if not isGrid:
    ServiceMgr.MessageSvc.OutputLevel = INFO
    ServiceMgr.MessageSvc.defaultLimit = 99999999

### Hypo ###
if AllDebug: MufastHypoConfig.OutputLevel = DEBUG
if AllDebug: MucombHypoConfig.OutputLevel = DEBUG
if AllDebug: TrigL2MuonOverlapRemoverConfig.OutputLevel = DEBUG
if SpecificDebug: MufastHypoConfig.OutputLevel = INFO
if SpecificDebug: MucombHypoConfig.OutputLevel = INFO
if SpecificDebug: TrigL2MuonOverlapRemoverConfig.OutputLevel = INFO

### TrigL2MuonSA ###
#MuFastSteering
if AllDebug: TrigL2MuonSAConfig.OutputLevel = DEBUG
if SpecificDebug: TrigL2MuonSAConfig.OutputLevel = INFO

#CscDataPreparator
preparator = TrigL2MuonSA__CscDataPreparator()
if AllDebug: preparator.OutputLevel = DEBUG
if SpecificDebug: preparator.OutputLevel = INFO
ToolSvc += preparator

#CscSegmentMaker
segment_maker = TrigL2MuonSA__CscSegmentMaker()
if AllDebug: segment_maker.OutputLevel= INFO
if SpecificDebug: segment_maker.OutputLevel= INFO


#AlphaBetaEstimate
alphabeta = TrigL2MuonSA__AlphaBetaEstimate()
if AllDebug: alphabeta.OutputLevel = DEBUG
if SpecificDebug: alphabeta.OutputLevel = INFO
ToolSvc += alphabeta

#PtFromAlphaBeta
ptfromangle = TrigL2MuonSA__PtFromAlphaBeta()
if AllDebug: ptfromangle.OutputLevel = DEBUG
if SpecificDebug: ptfromangle.OutputLevel = INFO
#ptfromangle.useCscPt = True
ToolSvc += ptfromangle


#PtEndcapLUT
pt_endcap_lut = TrigL2MuonSA__PtEndcapLUT()
if AllDebug: pt_endcap_lut.OutputLevel = DEBUG
if SpecificDebug: pt_endcap_lut.OutputLevel = INFO
ToolSvc += pt_endcap_lut

#AlignmentBarrelLUT
alignment_barrel = TrigL2MuonSA__AlignmentBarrelLUT()
if AllDebug: alignment_barrel.OutputLevel = DEBUG
if SpecificDebug: alignment_barrel.OutputLevel = INFO
ToolSvc += alignment_barrel

#AlignmentBarrelLUTSvc
#alignment_barrelSvc = TrigL2MuonSA__AlignmentBarrelLUTSvc()
#if AllDebug: alignment_barrelSvc.OutputLevel = DEBUG
#if SpecificDebug: alignment_barrelSvc.OutputLevel = INFO
#ToolSvc += alignment_barrelSvc


#MdtDataPreparator
mdt_prep = TrigL2MuonSA__MdtDataPreparator()
if AllDebug: mdt_prep.OutputLevel = DEBUG
if SpecificDebug: mdt_prep.OutputLevel = INFO
ToolSvc += mdt_prep

#MdtRegionDefiner
mdt_reg = TrigL2MuonSA__MdtRegionDefiner()
if AllDebug: mdt_reg.OutputLevel = DEBUG
if SpecificDebug: mdt_reg.OutputLevel = INFO
ToolSvc += mdt_reg


#MuCalStreamerTool
mucal_streamer = TrigL2MuonSA__MuCalStreamerTool()
if AllDebug: mucal_streamer.OutputLevel = DEBUG
if SpecificDebug: mucal_streamer.OutputLevel = INFO
ToolSvc += mucal_streamer


#MuFastDataPreparator
mufast_prep  = TrigL2MuonSA__MuFastDataPreparator()
if AllDebug: mufast_prep.OutputLevel = DEBUG
if SpecificDebug: mufast_prep.OutputLevel = INFO
ToolSvc += mufast_prep

#MuFastPatternFinder
mufast_patfinder  = TrigL2MuonSA__MuFastPatternFinder()
if AllDebug: mufast_patfinder.OutputLevel = DEBUG
if SpecificDebug: mufast_patfinder.OutputLevel = INFO
ToolSvc += mufast_patfinder


#MuFastStationFitter
mufast_stationfitter  = TrigL2MuonSA__MuFastStationFitter()
if AllDebug: mufast_stationfitter.OutputLevel = DEBUG
if SpecificDebug: mufast_stationfitter.OutputLevel = INFO
ToolSvc += mufast_stationfitter



#MuFastTrackExtrapolator
mufast_extrapolator  = TrigL2MuonSA__MuFastTrackExtrapolator()
if AllDebug: mufast_extrapolator.OutputLevel = DEBUG
if SpecificDebug: mufast_extrapolator.OutputLevel = INFO
ToolSvc += mufast_extrapolator



#MuFastTrackFitter
mufast_fitter = TrigL2MuonSA__MuFastTrackFitter()
if AllDebug: mufast_fitter.OutputLevel = DEBUG
if SpecificDebug: mufast_fitter.OutputLevel = INFO
ToolSvc += mufast_fitter



#PtBarrelLUT
ptbarrel  = TrigL2MuonSA__PtBarrelLUT()
if AllDebug: ptbarrel.OutputLevel = DEBUG
if SpecificDebug: ptbarrel.OutputLevel = INFO
ToolSvc += ptbarrel




#PtBarrelLUTSvc
#ptbarrelSvc  = TrigL2MuonSA__PtBarrelLUTSvc()
#if AllDebug: ptbarrelSvc.OutputLevel = DEBUG
#if SpecificDebug: ptbarrelSvc.OutputLevel = INFO
#ToolSvc += ptbarrelSvc



#PtFromRadius
pt_radius = TrigL2MuonSA__PtBarrelLUT()
if AllDebug: pt_radius.OutputLevel = DEBUG
if SpecificDebug: pt_radius.OutputLevel = INFO
ToolSvc += pt_radius



#RpcDataPreparator
rpc_prep =TrigL2MuonSA__RpcDataPreparator()
if AllDebug: rpc_prep.OutputLevel = DEBUG
if SpecificDebug: rpc_prep.OutputLevel = INFO
ToolSvc += rpc_prep



rpc_patfinder = TrigL2MuonSA__RpcPatFinder()
if AllDebug: rpc_patfinder.OutputLevel = DEBUG
if SpecificDebug: rpc_patfinder.OutputLevel = INFO
ToolSvc += rpc_patfinder



rpc_roaddefiner = TrigL2MuonSA__RpcRoadDefiner()
if AllDebug: rpc_roaddefiner.OutputLevel = DEBUG
if SpecificDebug: rpc_roaddefiner.OutputLevel = INFO
ToolSvc += rpc_roaddefiner


sagitta_radius = TrigL2MuonSA__SagittaRadiusEstimate()
if AllDebug: sagitta_radius.OutputLevel = DEBUG
if SpecificDebug: sagitta_radius.OutputLevel = INFO
ToolSvc += sagitta_radius


tgc_prep = TrigL2MuonSA__TgcDataPreparator()
if AllDebug: tgc_prep.OutputLevel = DEBUG
if SpecificDebug: tgc_prep.OutputLevel = INFO
ToolSvc += tgc_prep


tgc_fit = TrigL2MuonSA__TgcFit()
if AllDebug: tgc_fit.OutputLevel = DEBUG
if SpecificDebug: tgc_fit.OutputLevel = INFO
ToolSvc += tgc_fit


tgc_roaddefiner = TrigL2MuonSA__TgcRoadDefiner()
if AllDebug: tgc_roaddefiner.OutputLevel = DEBUG
if SpecificDebug: tgc_roaddefiner.OutputLevel = INFO
ToolSvc += tgc_roaddefiner


### TrigmuComb ###
#muComb
if AllDebug: TrigmuCombConfig.OutputLevel = DEBUG
if SpecificDebug: TrigmuCombConfig.OutputLevel = INFO

#muComb = TrigmuComb__muComb()
if AllDebug: muComb.OutputLevel = DEBUG
if SpecificDebug: muComb.OutputLevel = INFO
#ToolSvc += muComb

#from TrkExTools.AtlasExtrapolator import AtlasExtrapolator
#ToolSvc += AtlasExtrapolator()
#ToolSvc.AtlasExtrapolator.OutputLevel = VERBOSE

from TrigFastTrackFinder.TrigFastTrackFinderConf import TrigFastTrackFinder
TrigFastTrackFinder.OutputLevel = DEBUG

print "Setup my configuration..."

