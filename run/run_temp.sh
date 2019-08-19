#!/bin/sh

source $TestArea/../build/x86_64-centos7-gcc62-opt/setup.sh

pwd
cp ../*.sh ./
cp ../*.py ./
#INPUT_FILE="/gpfs/fs6001/yfukuhar/data/mc16_FTK/mc16_13TeV.424103.Pythia8B_A14NNPDF23LO_pp_Jpsi_mu15mu2p5.digit.RDO_FTK.e4397_s3126_r11234_d1505/RDO_FTK.17123182._000002.pool.root.1"
#INPUT_FILE="/gpfs/fs6001/yfukuhar/data/mc16_FTK/mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.digit.RDO_FTK.e5441_e5984_s3126_r11234_d1505/RDO_FTK17123282._000002.pool.root.1"
INPUT_FILE="/gpfs/fs6001/yfukuhar/data/mc16_FTK/mc16_13TeV.424100.Pythia8B_A14_CTEQ6L1_Jpsimu4mu4.digit.RDO_FTK.e3735_s3126_r11234_d1505/RDO_FTK.17123150._000001.pool.root.1"
#INPUT_FILE="/gpfs/fs6001/yfukuhar/data/mc16_FTK/mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.digit.RDO_FTK.e3601_s3126_d1480_d1471/RDO_FTK.15260455._000015.pool.root.1"

Reco_tf.py --outputAODFile AOD.pool.root  --inputRDO_FTKFile $INPUT_FILE --digiSteeringConf 'StandardSignalOnlyTruth'  --maxEvents 90 --conditionsTag 'default:OFLCOND-MC16-SDR-17' --valid 'True' --postInclude 'default:PyJobTransforms/UseFrontier.py' --autoConfiguration 'everything' --triggerConfig 'RDOtoRDOTrigger=MCRECO:MenuName:MC_pp_v7_tight_mc_prescale' --postExec 'all:CfgMgr.MessageSvc().setError+=["HepMcParticleLink"]' 'ESDtoAOD:fixedAttrib=[s if "CONTAINER_SPLITLEVEL = \"99\"" not in s else "" for s in svcMgr.AthenaPoolCnvSvc.PoolAttributes];svcMgr.AthenaPoolCnvSvc.PoolAttributes=fixedAttrib' 'RAWtoESD:topSequence.TrigSteer_HLT.enableRobRequestPreparation=False' --preExec 'all:rec.Commissioning.set_Value_and_Lock(True);from AthenaCommon.BeamFlags import jobproperties;jobproperties.Beam.numberOfCollisions.set_Value_and_Lock(20.0);from LArROD.LArRODFlags import larRODFlags;larRODFlags.NumberOfCollisions.set_Value_and_Lock(20);larRODFlags.nSamples.set_Value_and_Lock(4);larRODFlags.doOFCPileupOptimization.set_Value_and_Lock(True);larRODFlags.firstSample.set_Value_and_Lock(0);larRODFlags.useHighestGainAutoCorr.set_Value_and_Lock(True)' 'all:from InDetRecExample.InDetJobProperties import InDetFlags; InDetFlags.doSlimming.set_Value_and_Lock(False)' 'all:from TriggerJobOpts.TriggerFlags import TriggerFlags;TriggerFlags.AODEDMSet.set_Value_and_Lock("AODFULL");' 'all:from JetRec.JetRecFlags import jetFlags,JetContentDetail; jetFlags.detailLevel=JetContentDetail.Full' 'RAWtoESD:rec.doFTK.set_Value_and_Lock(True)' --geometryVersion 'default:ATLAS-R2-2016-01-00-01' --postInclude ../outputLevelConfig.py
