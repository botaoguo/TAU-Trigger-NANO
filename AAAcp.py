import os
import sys
for i in range(3):
   _cmd = 'nohup xrdcp root://cmsxrootd.fnal.gov//store/user/boguo/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/221108_023750/0000/tree_{}.root . >> nohup.log 2>&1 &'.format(i+51)
   os.system(_cmd)
