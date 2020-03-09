# LocalAthena

## Usage
```cp
$ git clone "url of this repositry"
$ ls 
>> run/ source/
$ cd source/
$ source setup.sh
$ ./checkout.sh your_branch_name
$ ./compile.sh cmake
>> compiling athena with cmake
>> After compile, run athena in run directory
$ cd ../run/
$ bsub -q 30m -o log.out Reco_tf.py --AMI=q221 --imf=True --athenaopts=--threads=1 --maxEvents=20 --inputRDOFile=/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/TriggerTest/valid1.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.RDO.e4993_s3214_r11315/RDO.17533168._000001.pool.root.1  --outputAODFile=AOD.pool.root --steering=doRDO_TRIG

or

$ Reco_tf.py --AMI=q221 --imf=True --athenaopts=--threads=1 --maxEvents=10 --inputRDOFile=/gpfs/fs7001/yfukuhar/data/mc16_FTK/mc16_13TeV.300901.ParticleGunEvtGen_Jpsi_mu3p5mu3p5_prompt.digit.RDO_FTK.e7406_e5984_a875_r11558_d1534_tid18953291_00/RDO_FTK.18953291._000001.pool.root.1 --outputAODFile=AOD.pool.root --steering=doRDO_TRIG 2>&1 | tee log.out
```


## Each script

### setup.sh
When you login again, you must do `source setup.sh`.
If you want to use different asetup version, you must edit setup.sh.
Like...
```sh
asetup r2018-03-09T2222,Athena,21.0
asetup r2018-03-29T2156,Athena,21.3
asetup 21.0,Athena,r31
asetup AtlasOffline,21.0.20,slc6,gcc62,64
asetup Athena,master,2020-01-25T2132,here
asetup Athena,21.5.6,here
```

FYI:
https://gitlab.cern.ch/atlas/athena#branches
https://indico.cern.ch/event/646945/contributions/3156795/attachments/1727208/2790438/20181003_HLT.pdf


### checkout.sh
```sh
$ ./checkout.sh your_branch_name
$ git branch
master
* your_branch_name
```

### compile.sh
Once you edit code, you must compile.
Like...
```sh
$ ./compile.sh cmake
```
Then build directory is made.
And once you comile, you must do
```sh
$ source $TestArea/../build/$CMTCONFIG/setup.sh

```


## How to edit code
```sh
$ git fetch upstream
$ git checkout -b my_branch_name upstream/[project-branch] --no-track
```
[project-branch] is master, 21.0, 21.3 etc..

## How to merge code from master branch to your own branch under development.
```sh
$ git merge nightly/master/2020-01-25T2132
```



