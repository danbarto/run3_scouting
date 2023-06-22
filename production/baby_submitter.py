#!/usr/bin/env python3

from metis.CondorTask import CondorTask
from metis.Sample import FilelistSample
from metis.Path import Path
from metis.StatsParser import StatsParser
import time
import os
import glob

datasets = {
    'Run2022D': '/ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/'
}



if __name__ == '__main__':

    total_summary = {}

    reqname = 'ScoutingPFRun3'
    babymaker = 'babymaker.py'

    tag = "v3"

    small = False
    tasks = []

    for dataset in datasets:

        input_files = glob.glob(f"{datasets[dataset]}/*.root")
        if small: input_files = input_files[1:5]

        sample = FilelistSample(
            #dataset = f"{req_name}_{dataset}",
            dataset = dataset,
            filelist=input_files,
            #use_xroot=True,  # FIXME maybe needed?
        )

        task = CondorTask(
                sample = sample,
                output_name = "output.root",
                executable = "executables/baby_production.sh",
                tarfile = "package.tar.gz",  # contains propagation_utils etc
                additional_input_files = [babymaker],
                scram_arch = "el8_amd64_gcc10",  # not directly used, but better define it
                cmssw_version = "CMSSW_12_4_12",  # not directly used, but better define it
                open_dataset = False,
                files_per_output = 5,
                arguments = babymaker,
                condor_submit_params = {
                "sites":"T2_US_UCSD", #
                "classads": [
                        ["metis_extraargs",""],
                        ["JobBatchName",reqname],
                        ["SingularityImage", "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel8-m-m20230223"],
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
