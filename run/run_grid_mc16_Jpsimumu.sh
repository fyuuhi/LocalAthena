#!/bin/sh
DATE=$(date '+%Y%m%d%H%M')

source $TestArea/../build/$CMTCONFIG/setup.sh

INPUT_RDO="user.yfukuhar.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.recon.RDO.e3601_s3126_r9781.201806291853_EXT0"
INPUT_RDO_FOR_OUTPUT="user.yfukuhar.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.recon.RDO.e3601_s3126_r9781"
OUTPUT_AOD=$INPUT_RDO_FOR_OUTPUT".AOD."$DATE""
CMT_CONFIG="x86_64-slc6-gcc62-opt"
AMI_TAG="r9781"
NFILES_PER_JOB='5'
NEVENTS_PER_FILE='500'
NEVENTS_PER_JOB='100'
DESTSE='TOKYO-LCG2_SCRATCHDISK'

echo ""
echo "INPUT_RDO =    "$INPUT_RDO
echo "OUTPUT_AOD =   "$OUTPUT_AOD
echo ""

#LOG_NAME="log_"$LOG_INDEX"_"$DATE

SUBMIT_AREA="Grid_SubmitArea_"$DATE
mkdir $SUBMIT_AREA
cd $SUBMIT_AREA
cp ../*.sh ./

echo ${0} 2>&1 | tee -a jediTaskID.info
echo ""

COMMAND="pathena --trf \"Reco_tf.py --outputAODFile %OUT.AOD.pool.root --inputRDOFile %IN\" --nFilesPerJob=5 --inDS ${INPUT_RDO} --outDS ${OUTPUT_AOD} --cmtConfig ${CMTCONFIG}"

echo "COMMAND:"
echo $COMMAND

echo "" | tee -a jediTaskID.info
echo "LOG: " | tee -a jediTaskID.info

eval $COMMAND 2>&1 | tee -a jediTaskID.info

JEDITASKID=$(sed -n '/new jediTaskID=/s/INFO : succeeded. new jediTaskID=//p' jediTaskID.info)
#sed -n '/new jediTaskID/p' tmp_test.txt
echo ""
echo ""
echo "JEDITASKID: ${JEDITASKID}" 2>&1 | tee -a jediTaskID.info
echo "JEDI TASK URL: https://bigpanda.cern.ch/task/${JEDITASKID}/" 2>&1 | tee -a jediTaskID.info
echo ""
echo ""

echo "INFO: Create jediTaskID.info"
echo ""

cd -
