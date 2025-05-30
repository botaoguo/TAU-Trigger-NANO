import os,sys

era = "Run2024C"

# python3 submit_condor.py zvbftest/Muon0_Run2023D_v1/output_Muon*
# python3 submit_condor.py zvbftest/Muon0_Run2023D_v2/output_Muon*
# python3 submit_condor.py zvbftest/Muon1_Run2023D_v1/output_Muon*
# python3 submit_condor.py zvbftest/Muon1_Run2023D_v2/output_Muon*

# python3 submit_condor.py zvbftest/Muon0_Run2023C_v1/output_Muon*
# python3 submit_condor.py zvbftest/Muon0_Run2023C_v2/output_Muon*
# python3 submit_condor.py zvbftest/Muon0_Run2023C_v3/output_Muon*
# python3 submit_condor.py zvbftest/Muon0_Run2023C_v4/output_Muon*

# python3 submit_condor.py zvbftest/Muon1_Run2023C_v1/output_Muon*
# python3 submit_condor.py zvbftest/Muon1_Run2023C_v2/output_Muon*
# python3 submit_condor.py zvbftest/Muon1_Run2023C_v3/output_Muon*
# python3 submit_condor.py zvbftest/Muon1_Run2023C_v4/output_Muon*


# python3 submit_condor.py zvbftest/Run2024C/output_Muon*
# python3 submit_condor.py zvbftest/Run2024D/output_Muon*
# python3 submit_condor.py zvbftest/Run2024E/output_Muon*
# python3 submit_condor.py zvbftest/Run2024F/output_Muon*
# python3 submit_condor.py zvbftest/Run2024G/output_Muon*
# python3 submit_condor.py zvbftest/Run2024H/output_Muon*
# python3 submit_condor.py zvbftest/Run2024I/output_Muon*

# python3 submit_condor.py zvbftest/Run2025C/output_Muon*
count = 0
file_names = sys.argv[1:]
for file_name in file_names:
    fileprefix = os.path.basename(file_name)
    fileprefix = os.path.splitext(fileprefix)[0]
    fileprefix = fileprefix.strip("output_")
    idx = fileprefix.strip("Muon")
    idx = idx[2:]

    # Define the shell script content as a multi-line string
    _wrapper = '''#!/bin/bash
chmod +x sk2.py
# Define the variables
USER=$(whoami)
first_letter=${{USER:0:1}}
# Construct the EOS home directory path
eos_home_dir="/eos/home-${{first_letter}}/${{USER}}/"
eos_save_dir="${{eos_home_dir}}/botao/CMSSW_10_6_29/src/PhysicsTools/NanoAODTools/TAU-Trigger-NANO/zvbftest/jobs"
export eos_save_dir="${{eos_home_dir}}/botao/CMSSW_10_6_29/src/PhysicsTools/NanoAODTools/TAU-Trigger-NANO/zvbftest/jobs"
# Ensure the target directory exists
mkdir -p "${{eos_save_dir}}"
echo " output will be in ${{eos_save_dir}}/skimtuple_veto_sec_mu_{1}_{0}.root"
python3 sk2.py --inputlist ./output_{0}.txt --output ./ --version veto_sec_mu_{1}_{0}
# Now, transfer the output file
mv skimtuple_veto_sec_mu_{1}_{0}.root "${{eos_save_dir}}/skimtuple_veto_sec_mu_{1}_{0}.root"
    '''.format(fileprefix,era)
    _wrapper_name = "{1}_wrapper_{0}.sh".format(fileprefix,era)
    _wrapper_f = open(_wrapper_name, "w")
    _wrapper_f.write(_wrapper)
    _wrapper_f.close()
    
    # Define the menu for job
    _menu = '''Universe = vanilla
Executable = {1}_wrapper_{0}.sh
Proxy_filename = x509up
arguments = $(Proxy_filename) $(Cluster) $(Process)
should_transfer_files = YES
transfer_input_files = $(Proxy_filename), sk2.py, 2024C_Golden.txt, interface/tau_ntupler.h, Common/python/AnalysisTools.py, Common/python/RootObjects.py, zvbftest/{1}/output_{0}.txt
when_to_transfer_output = ON_EXIT
+AccountingGroup = group_u_CMS.u_zh.users
universe = vanilla
use_x509userproxy = true
x509userproxy = /afs/cern.ch/user/b/boguo/TAU-Trigger-NANO/x509up
transfer_output_files = ""
request_cpus = 1
Error = job/job_$(Cluster)_$(Process)_{0}_{1}.err
Output = job/job_$(Cluster)_$(Process)_{0}_{1}.out
Log = job/job_$(Cluster)_$(Process)_{0}_{1}.log
+MaxRuntime= 86400
queue 1
    '''.format(fileprefix,era)
    _menu_name = "submit_{1}_{0}.sub".format(fileprefix,era)
    _menu_f = open(_menu_name, "w")
    _menu_f.write(_menu)
    _menu_f.close()
    
    # submit the job
    print("condor_submit {0}".format(_menu_name))
    os.system('condor_submit {0}'.format(_menu_name))
    count += 1
    # exit(0)

# os.system("mkdir -p job/script")
# os.system("mv submit_2024*.sub 2024*_wrapper_*.sh job/script/")
print("total submit {} jobs".format(count))