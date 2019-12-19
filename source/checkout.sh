BRANCH_NAME=$1
if [ "$BRANCH_NAME" = "" ] ; then
  BRANCH_NAME="Temp"
fi
echo "Branch name is $BRANCH_NAME"

#git atlas init-workdir ssh://git@gitlab.cern.ch:7999/yoyamagu/athena.git
git atlas init-workdir ssh://git@gitlab.cern.ch:7999/yfukuhar/athena.git
#git clone ssh://git@gitlab.cern.ch:7999/yyazawa/athena.git
cd athena
git atlas addpkg TrigL2MuonSA
git atlas addpkg TrigmuComb
git atlas addpkg TrigMuonEF
git atlas addpkg TriggerMenu
git atlas addpkg TrigEDMConfig
git atlas addpkg TrigMuonHypo
git atlas addpkg TriggerMenuMT
git atlas addpkg TrigMuonHypoMT

git fetch upstream
#git checkout TrigL2MuonSA-00-03-63-master-fixRpcPadMsg1
#git checkout -b $BRANCH_NAME origin/21.5.6_InOutL2Muon --no-track
git checkout -b tmp origin/FTFInsideOut1Nv2 --no-track

#svnpull.py TrigL2MuonSA-00-03-63
cd -

#21.0
#Tier0 and MC16c
#Athena, 21.0.X
#
#21.0-mc16a
#MC16a bugfixes
#Athena, 21.0.20.Y
#
#21.1
#HLT and P1 monitoring
#AthenaP1, 21.1.X
#
#21.1-dev
#HLT and P1 Run 2 Developments
#AthenaP1, No production releases
#
#21.2
#Derivations and Analysis
#AthDerivations + (Ath)AnalysisBase, 21.2.X
#
#21.3
#MC18 Simulation
#Athena + AthSimulation, 21.3.X
#
#master
#AthenaMT development
#Athena, 22.0.X
