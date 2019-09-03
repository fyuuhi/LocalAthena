#!/usr/bin/env python
import os, sys
from datetime import datetime

date = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
time =  datetime.now().strftime('%Y%m%d%H%M')
#testarea = os.environ.get('TestArea')
#testarea = '%s/..' %os.path(__file__)


# sample = 'tau3mu'
#sample = 'Jpsi_mu15mu2'
# sample = 'Jpsi_mu15_mu2_rel21'
# sample = 'Jpsi_mu4_mu4_rel21'
#sample = 'Zmumu_rel21'
# sample = 'Zmumu_21.0.20'
#sample = 'Zmumu_21.0.15'
#sample = 'Zmumu_21.0.32'
#sample = 'Jpsimumu_21.0.53'
sample = 'Jpsimumu_21.3.10'
#sample = 'Jpsimu4mu4_21.0.24'
maxEvents = 100
nFilesPerJob = 1
nEventsPerJob = 100
nEventsPerFile = 500
nFiles = 500
#site = "ANALY_RAL_SL6"
#ExcludedSite = 'ANALY_BU_ATLAS_Tier2_SL6,ANALY_CONNECT,ANALY_OU_OCHEP_SWT2,ANALY_AUSTRALIA,ANALY_SLAC,INFN,ANALY_RHUL_SL6,ANALY_TAIWAN'
#destSE = 'TOKYO-LCG2_SCRATCHDISK'
#destSE = 'SLACXRD_SCRATCHDISK'
#destSE = 'MWT2_DATADISK'
extOutFile = 'expert-monitoring.root'

from digiconfigftk import digiconfigftk
params = digiconfigftk[sample]
AMItag = params['AMItag']
athenaTag = params['athenaTag']
cmtConfig = params['cmtConfig']
InDS = params['InDS']
#highMinDS = params['highMinDS']
#nHighMin = params['nHighMin']
#lowMinDS = params['lowMinDS']
#nLowMin = params['nLowMin']
OutDS = params['OutDS']
nEventsPerFile = params['nEventsPerFile']



print ''
print date
print ''
print 'InputDS   = %s' %InDS
#print 'highMinDS   = %s' %highMinDS
#print 'lowMinDS   = %s' %lowMinDS
print 'OutputDS  = %s\033[32;1m' %OutDS
print 'nFiles    = %d\033[0m' %nFiles
print ''



#parameters = 'Reco_tf.py'
parameters = 'TrigFTKTM64SM1Un_tf.py'
parameters += ' --AMIConfig %s' %AMItag
#parameters += ' --inputHITSFile %IN'
parameters += ' --inputRDOFile %IN'
#parameters += ' --inputLowPtMinbiasHitsFile %LOMBIN'
#parameters += ' --inputHighPtMinbiasHitsFile %HIMBIN'
parameters += ' --outputRDO_FTKFile %OUT.RDO_FTK.pool.root'
#parameters += ' --jobNumber %RNDM:0'
#parameters += ' --digiSeedOffset1 %RNDM:0'
#parameters += ' --digiSeedOffset2 %RNDM:0'
parameters += ' --maxEvents %MAXEVENTS'
parameters += ' --skipEvents %SKIPEVENTS'


command = 'pathena '
#command = '/gpfs/fs2001/ynoguchi/PandaRun/pathena'
#command = './pathena'
command += ' --trf "%s"' %parameters
command += ' --inDS %s' %InDS
#command += ' --lowMinDS %s' %lowMinDS
#command += ' --nLowMin %d' %nLowMin 
#command += ' --highMinDS %s' %highMinDS
#command += ' --nHighMin %d' %nHighMin
command += ' --outDS %s' %OutDS
#command += ' --excludedSite %s' %ExcludedSite
#command += ' --site %s' %site
#command += ' --destSE %s' %destSE
command += ' --nEventsPerFile %d' %nEventsPerFile
command += ' --nEventsPerJob %d' %nEventsPerJob
command += ' --nFilesPerJob %d' %nFilesPerJob
command += ' --individualOutDS'
#command += ' --athenaTag %s' %athenaTag
command += ' --cmtConfig %s' %cmtConfig

## Conditional                                                                                                             
if nFiles > -1:
    command += ' --nFiles %d' %nFiles
#if supStream!='':
#    command += ' --supStream %s' %supStream

## Argv                                                                                                                    
for arg in sys.argv:
    if arg.find('submitftk.py') == -1:
        command += ' %s' %arg

com_display = command.replace('%s' %InDS, '[Inputs]')
#com_display = com_display.replace('%s' %highMinDS, '[InputHighMinBias]')
#com_display = com_display.replace('%s' %lowMinDS, '[InputLowMinBias]')
com_display = com_display.replace('%s' %OutDS, '[Outputs]')
#com_display = com_display.replace('%s' %ExcludedSite, '[Excludes]')


print '%s\n' %command
#print '%s\n' %com_display

#os.system('rm -r %s/../SubmitArea' %testarea)
#os.system('mkdir %s/../SubmitArea' %testarea)
#os.chdir('%s/../SubmitArea' %testarea)
#os.system('cp %s/../*.py ./' %testarea)
#os.system('cp %s/../*.sh ./' %testarea)

#os.system('rm -r SubmitArea_XFAhifhasdfIS')
os.system('mkdir Submit_%s' %(time))
os.chdir('Submit_%s' %(time))
os.system('cp ../../*.py ./')
os.system('cp ../../*.sh ./')
os.system('cp ../../pathena ./')

os.system(command)

#os.chdir('%s/scripts' %testarea)



print ''

