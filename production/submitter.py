'''
Submitter script to produce samples for
https://gitlab.com/simonknapen/dark_showers_tool/-/blob/master/cards/test_scenarioA.cmnd
'''
#!/usr/bin/env python3

from metis.CondorTask import CondorTask
from metis.Sample import DirectorySample,DummySample
from metis.Path import Path
from metis.StatsParser import StatsParser
import time
import os


if __name__ == '__main__':



    total_summary = {}

    # NOTE make this script more flexible again
    '''
    requests = [
        {'name': 'DarkShower_ScenarioA', 'fragment/a_fragment.py'}
    ]
    '''

    #reqname = 'DarkShower_ScenarioB1_default'  # Scouting For Scouting

    #fragment = "psets/darkshower_cfg.py"  # This defines the physics
    #

    requests = {
        "ScenB1_30_9p9_4p8_ctau_26_filter": "psets/ScenB1_30_9p9_4p8_ctau_26_filter_cfg.py",
        "ScenA_20_5p0_1p2_ctau_23_filter": "psets/ScenA_20_5p0_1p2_ctau_23_filter_cfg.py",
    }
    events_per_point = 20
    events_per_job = 10
    njobs = int(events_per_point)//events_per_job

    tag = "v5p0"  # v0p8 first one with compiling on worker, 14 switch to el8 from rhel8
    # 17 - trying to switch of pythia multithreading
    campaign = "Run3Summer22GS"

    tasks = []

    for reqname in requests:
        fragment = requests[reqname]

        task = CondorTask(
            sample = DummySample(
                dataset=f"/{reqname}/{campaign}_{tag}/AODSIM",
                N=njobs,
                nevents=int(events_per_point),
            ),
            output_name = "output.root",
            executable = "executables/test_signal_production.sh",
            tarfile = "package.tar.gz",  # contains precompiled CMSSW_12_4_12, pythia and RAWSIM/AODSIM configs
            additional_input_files = [fragment],  # this is signal point dependent
            scram_arch = "el8_amd64_gcc10",  # not directly used, but better define it
            # should we switch to el8_amd64_gcc10??
            cmssw_version = "CMSSW_12_4_12",  # not directly used, but better define it
            open_dataset = False,
            files_per_output = 1,
            arguments = fragment,  # tell the job the name of the fragment
            condor_submit_params = {
                "sites":"T2_US_UCSD", #
                "classads": [
                    ["param_nevents",events_per_job],
                    ["metis_extraargs",""],
                    ["JobBatchName",reqname],
                    ["SingularityImage", "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel8-m-m20230223"],
                    #["SingularityImage", "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/cc7:x86_64-latest"],
                    #["SingularityImage", "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64-d20230317"],
                    #["SingularityImage", "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m20221104"],
                    ],
                "requirements_line": 'Requirements = (HAS_SINGULARITY=?=True)'  # && (HAS_CVMFS_cms_cern_ch =?= true) && {extra_requirements})'.format(extra_requirements=extra_requirements),
                },
            tag = tag,
            min_completion_fraction = 0.90,
            )
        tasks.append(task)

    for task in tasks:
        task.process()
        total_summary[task.get_sample().get_datasetname()] = task.get_task_summary()

    StatsParser(data=total_summary, webdir="~/public_html/dump/scouting/").do()
