#!/bin/sh
DATE=$(date '+%Y%m%d%H%M')

#source $TestArea/../build/$CMTCONFIG/setup.sh

#INPUT_LIST="../ZmumuRDO_FTK_long.list"
#INPUT_LIST="../ZmumuRDO_FTK.list"
#INPUT_LIST="../JpsiRDO_FTK.list"
#INPUT_LIST="../shortZmumuRDO_FTK.list"
#INPUT_LIST="../Jpsi_424130.list"
#INPUT_LIST="../Jpsi_424103.list"
#INPUT_LIST="../Jpsi_424108.list"
#INPUT_LIST="../Jpsi_424100.list"
INPUT_LIST="../FTK_361107.list"
#INPUT_RDO_FOR_OUTPUT="user.yfukuhar.mc16_13TeV.Zmumu"
#OUTPUT_AOD=$INPUT_RDO_FOR_OUTPUT".AOD."$DATE".root"
TEMP_RUN_SCRIPT="bsub_run_temp.sh"
#MAXEVENTS='100'
LABEL=DIR_RDO_FTK_Zmumu
#LABEL=DIR_RDO_FTK_Jmumu

echo ""
echo "INPUT_LIST =   "$INPUT_LIST
echo "OUTPUT_AOD =   "$OUTPUT_AOD
echo ""


#echo "" | tee -a jediTaskID.info
#echo "LOG: " | tee -a jediTaskID.info
#
#eval $COMMAND 2>&1 | tee -a jediTaskID.info
#LOG_NAME="log_"$LOG_INDEX"_"$DATE

echo $PWD
COUNT=0

if [ -f "$INPUT_LIST" ] ; then
  echo $INPUT_LIST
  for INPUT_FILE in `cat $INPUT_LIST`
  do
    echo ""
    echo "INPUT_FILE =   "$INPUT_FILE
    TEMP_DIR=$LABEL"_"$COUNT
    DIR="$TEMP_DIR"
    mkdir -p $DIR
    echo "Submit Directory: "$DIR
    cp ../*.sh $DIR
    cp ../*.py $DIR
    cd $DIR
    sed -i -e "10i  \ \ --inputRDO_FTKFile $INPUT_FILE \ \\\\" $TEMP_RUN_SCRIPT
    echo "=========================="
    bsub -n 3 -q 1d -o log_out ./bsub_run_temp.sh
    echo "=========================="
    echo ""

    #COMMAND="bsub -q 2d -o log.out 'source $TestArea/../build/$CMTCONFIG/setup.sh && Reco_tf.py --maxEvents ${MAXEVENTS} --outputAODFile ${OUTPUT_AOD} --inputRDOFile ${INPUT_FILE}' --postInclude ./outputLevelConfig.py"
    #echo $COMMAND > tmp_athena_run.sh
    #sh tmp_athena_run.sh
    #echo $PWD
    cd -
    echo ""
    echo ""
    #sleep 1
    COUNT=`expr $COUNT + 1`
  done
fi
