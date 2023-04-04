import os
import sys

def run_TriggerTool(_dir, _isMC, _start, _end):
    for _jobid in range(_start, _end + 1):
        _cmd = "nohup python postproc.py --input /eos/cms/store/group/phys_tau/TauFW/nanoV10/Run2_2018/{0}/nano_{1}.root --isMC {2} --era 2018 > nohup_{1}.log 2>&1 &".format(_dir, _jobid, _isMC)
        print(_cmd)
        print("Running now!")
        os.system(_cmd)
        # exit(0)

if __name__ == '__main__':
    _dir = "SingleMuon_Run2018A"
    _isMC = int(sys.argv[1])
    _start = int(sys.argv[2])
    _end = int(sys.argv[3])
    run_TriggerTool(_dir, _isMC, _start, _end)
    print("End of the script...")