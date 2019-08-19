#!/usr/bin/env python
import os, sys
from datetime import datetime
time =  datetime.now().strftime('%Y%m%d%H%M')
user = os.environ.get('USER')

digiconfig = {}
dic={}
#'Jpsi_mu15mu2
    #http://bigpanda.cern.ch/task/?jeditaskid=8884673&display_limit=200


#dic = { 
#        'sample': 'Jpsi_mu15mu2', 
#        'AMItag': 'r7772',
#        'athenaTag': 'AtlasProd1,20.7.5.1.1',
#        'cmtConfig': 'x86_64-slc6-gcc49-opt',
#        'InDS': 'mc15_13TeV.424103.Pythia8B_A14NNPDF23LO_pp_Jpsi_mu15mu2p5.merge.HITS.e4311_s2608_s2183_tid06375809_00',
#        'highMinDS': 'mc15_13TeV.361035.Pythia8EvtGen_A2MSTW2008LO_minbias_inelastic_high.merge.HITS.e3581_s2578_s2195',
#        'nHighMin': 1,
#        'lowMinDS': 'mc15_13TeV.361034.Pythia8EvtGen_A2MSTW2008LO_minbias_inelastic_low.merge.HITS.e3581_s2578_s2195',
#        'nLowMin': 3,
#        'outDS': 'user.%s.mc15_13TeV.424103.Pythia8B_A14NNPDF23LO_pp_Jpsi_mu15mu2p5.recon.RDO.e4311_s2608_s2183_r7772.%s' %(user, time),
#        'EventsPerFile': 1000
#        }
#
#digiconfig[dic['sample']] = dic
#
#
#
#
#dic = {
#        'sample': 'tau3mu',
#        #http://bigpanda.cern.ch/task/8542355/
#        'AMItag': 'r8012',
#        'athenaTag': 'TrigMC,20.7.6.1.1',
#        'cmtConfig': 'x86_64-slc6-gcc49-opt',
#        'InDS': 'mc15_13TeV.180593.Pythia8_AUMSTW2008LO_Wtaunu_3mu_noPhotos.merge.HITS.e3802_s2608_s2183_tid05370755_00',
#        'highMinDS': 'mc15_valid.361035.Pythia8EvtGen_A2MSTW2008LO_minbias_inelastic_high.merge.HITS.e3581_s2578_s2169_tid05098387_00',
#        'nHighMin': 1,
#        'lowMinDS': 'mc15_valid.361034.Pythia8EvtGen_A2MSTW2008LO_minbias_inelastic_low.merge.HITS.e3581_s2578_s2169_tid05098374_00',
#        'nLowMin': 3,
#        'OutDS': 'user.%s.mc15_13TeV.180593.Pythia8_AUMSTW2008LO_Wtaunu_3mu_noPhotos.recon.RDO.e3802_s2608_s2183_r8012.%s' %(user, time),
#        'nEventsPerFile': 1000
#        }
#digiconfig[dic['sample']] = dic
#
#
#dic = {
#        'sample': 'Jpsi_mu15_mu2_rel21',
#        #http://bigpanda.cern.ch/task/?jeditaskid=10885974&display_limit=200
#        'AMItag': 'r9191',
#        'athenaTag': 'AtlasOffline,21.0.16',
#        'cmtConfig': 'x86_64-slc6-gcc49-opt',
#        'InDS': 'mc15_13TeV.424103.Pythia8B_A14NNPDF23LO_pp_Jpsi_mu15mu2p5.merge.HITS.e4311_s2608_s2183_tid06375809_00',
#        'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.merge.HITS.e4981_s3035_s3039',
#        'nHighMin': 5,
#        'lowMinDS': 'mc16_13TeV.361238.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_low.merge.HITS.e4981_s3035_s3039',
#        'nLowMin': 5,
#        'OutDS': 'user.%s.mc15_13TeV.424103.Pythia8B_A14NNPDF23LO_pp_Jpsi_mu15mu2p5.recon.RDO.e4311_s2608_s2183_r9191.%s' %(user, time),
#        'nEventsPerFile': 1000
#        }
#digiconfig[dic['sample']] = dic
#
#
#dic = {
#        'sample': 'Jpsi_mu4_mu4_rel21',
#        #http://bigpanda.cern.ch/task/?jeditaskid=10885974&display_limit=200
#        'AMItag': 'r9191',
#        'athenaTag': 'AtlasOffline,21.0.16',
#        'cmtConfig': 'x86_64-slc6-gcc49-opt',
#        'InDS': 'mc16_13TeV.300201.Pythia8BPhotospp_A14_CTEQ6L1_bb_Jpsimu4mu4.simul.HITS.e4397_s2997_tid10232564_00',
#        'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.merge.HITS.e4981_s3035_s3039',
#        'nHighMin': 5,
#        'lowMinDS': 'mc16_13TeV.361238.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_low.merge.HITS.e4981_s3035_s3039',
#        'nLowMin': 5,
#        'OutDS': 'user.%s.mc16_13TeV.300201.Pythia8BPhotospp_A14_CTEQ6L1_bb_Jpsimu4mu4.recon.RDO.e4397_s2997_r9191.%s' %(user, time),
#        'nEventsPerFile': 1000
#        }
#digiconfig[dic['sample']] = dic
#
#
#
#dic = { 
#        'sample': 'Zmumu_rel21',
#        #http://bigpanda.cern.ch/task/?jeditaskid=10872772&display_limit=200
#        'AMItag': 'r9191',
#        'athenaTag': 'AtlasOffline,21.0.16',
#        'cmtConfig': 'x86_64-slc6-gcc49-opt',
#        'InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s2997_tid09794988_00',
#        'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.merge.HITS.e4981_s3035_s3039',
#        'nHighMin': 5 ,
#        'lowMinDS': 'mc16_13TeV.361238.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_low.merge.HITS.e4981_s3035_s3039',
#        'nLowMin': 5 ,
#        'OutDS': 'user.%s.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.recon.RDO.e3601_s2997_r9191.%s' %(user, time),
#        'nEventsPerFile': 1000
#        }
#digiconfig[dic['sample']] = dic
#
#
#dic = {
#        'sample': 'Zmumu_21.0.20',
#        #http://bigpanda.cern.ch/task/?jeditaskid=11182595&display_limit=200
#        'AMItag': 'r9364',
#        'athenaTag': 'AtlasOffline,21.0.20',
#        'cmtConfig': 'x86_64-slc6-gcc62-opt',
#        '#InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00',
#        'InDS': 'user.%s.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00_201705141555',
#        'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.simul.HITS.e4981_s3087_s3111',
#        'nHighMin': 5,
#        'lowMinDS': 'mc16_13TeV.361238.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_low.simul.HITS.e4981_s3087_s3111',
#        'nLowMin': 5,
#        'OutDS': 'user.%s.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.recon.RDO.e3601_s2997_r9364.%s' %(user, time),
#        'nEventsPerFile': 1000,
#        }
#digiconfig[dic['sample']] = dic
#
#
#dic = {
#        'sample': 'Jpsimu4mu4_21.0.24',
#        #http://bigpanda.cern.ch/task/?jeditaskid=11182595&display_limit=200
#        'AMItag': 'r9573',
#        'athenaTag': 'Athena,21.0.24',
#        'cmtConfig': 'x86_64-slc6-gcc62-opt',
#        'InDS': 'valid1.424100.Pythia8B_A14_CTEQ6L1_Jpsimu4mu4.simul.HITS.e5112_s3091',
#        #'InDS': 'user.%s.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00_201705141555',
#        'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.simul.HITS.e4981_s3087_s3111',
#        'nHighMin': 5,
#        'lowMinDS': 'mc16_13TeV.361238.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_low.simul.HITS.e4981_s3087_s3111',
#        'nLowMin': 5,
#        'OutDS': 'user.%s.valid1.424100.Pythia8B_A14_CTEQ6L1_Jpsimu4mu4.recon.RDO.e5112_s3091_r9573._%s' %(user, time),
#        'nEventsPerFile': 100,
#        }
#digiconfig[dic['sample']] = dic


# Zmumu
#dic = {
#  'sample': 'Zmumu_21.0.32',
#  #http://bigpanda.cern.ch/task/?jeditaskid=11182595&display_limit=200
#  'AMItag': 'r9781',
#  'athenaTag': 'Athena,21.0.32',
#  'cmtConfig': 'x86_64-slc6-gcc62-opt',
#  # '#InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00',
#  # 'InDS': 'user.%s.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00_201705141555',
#  #'InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126',
#  #mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00',
#  #mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid11363382_00',
#  #mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid11012007_00',
#  #mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid13642339_00',
#  #'#InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00',
#  'InDS': 'user.yfukuhar.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_der1530232441',
#  #'InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid11363382_00',
#  #'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.simul.HITS.e4981_s3087_s3111_tid10701335_00',
#  #'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.simul.HITS.e4981_s3087_s3111',
#  'highMinDS': 'user.yfukuhar.mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.simul.HITS.e4981_s3087_s3111_der1530163766',
#  'nHighMin': 5,
#  'lowMinDS': 'user.yfukuhar.mc16_13TeV.361238.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_low.simul.HITS.e4981_s3087_s3111_der1530095342',
#  'nLowMin': 5,
#  'OutDS': 'user.%s.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.recon.RDO.e3601_s3126_r9781.%s' %(user, time),
#  'nEventsPerFile': 1000,
#}

# Zmumu
dic = {
  'sample': 'Jpsimumu_21.0.53',
  #http://bigpanda.cern.ch/task/?jeditaskid=11182595&display_limit=200
  'AMItag': 'r9781',
  'athenaTag': 'Athena,21.0.53',
  'cmtConfig': 'x86_64-slc6-gcc62-opt',
  # '#InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00',
  # 'InDS': 'user.%s.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00_201705141555',
  #'InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126',
  #mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00',
  #mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid11363382_00',
  #mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid11012007_00',
  #mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid13642339_00',
  #'#InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid10730601_00',
  #'InDS': 'user.yfukuhar.mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.simul.HITS.e5441_e5984_s3126_der1531370800',
  #'InDS': 'mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.simul.HITS.e5441_s3126_tid11330326_00',
  #'InDS': 'mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.simul.HITS.e5441_s3126_tid10730179_00',
  'InDS': 'mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.simul.HITS.e5441_e5984_s3126_tid14253740_00',
  #'InDS': 'mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.simul.HITS.e3601_s3126_tid11363382_00',
  #'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.simul.HITS.e4981_s3087_s3111_tid10701335_00',
  #'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.simul.HITS.e4981_s3087_s3111',
  'highMinDS': 'mc16_13TeV.361239.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_high.simul.HITS.e4981_s3087_s3111_tid10701335_00',
  'nHighMin': 5,
  'lowMinDS': 'mc16_13TeV.361238.Pythia8EvtGen_A3NNPDF23LO_minbias_inelastic_low.simul.HITS.e4981_s3087_s3111_tid10701323_00',
  'nLowMin': 5,
  'OutDS': 'user.%s.mc16_13TeV.424108.Pythia8B_A14_CTEQ6L1_Jpsimu6.simul.RDO.e5441_e5984_s3126.%s' %(user, time),
  #'OutDS': 'user.%s.mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.recon.RDO.e3601_s3126_r9781.%s' %(user, time),
  'nEventsPerFile': 1000,
}
digiconfig[dic['sample']] = dic
