#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: <data>
#@author: <name>
#@contact: <email>

import sys
import argparse
from settings import settings, tools_settings
from inspect import getsourcefile
from scripts import aragorn
from scripts import barrnap

def input_checker(file):
    #checks for existing and size
    pass

def output_checker(file):
    #check file for not existing before run
    pass
    
def tool_runner(tool_command):
    # if input exists and output don't, run command
    print(tool_command)
    os.system(tool_command)


def main(settings, tools_settings):
    # checks for files 
    
    # aragorn
    
    # here should be run of function which creates comand for tool (create it in scripts folder)
    aragorn_command = aragorn_command(settings["aragorn_binary"], tools_settings["aragorn_options"]["aragorn_input"], ) 
    #barrnap
    barrnap_command = "" # the same as for aragorn
    
    tools_commands = [arogorn_command, barrnap_command]
    
    for command in tools_commands:
        tool_runner(command)
        #check_for_output


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Program description.')
    parser.add_argument('-i','--input', help='Input file', required=True)
    parser.add_argument('-o','--output', help='Output folder', required=True)
    
    args = vars(parser.parse_args())
    
    input_file = args["input"]
    output_folder = args["output"]
    
    execution_folder = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
    
    settings = {
    "threads": 32,
    "input_file": os.path.abspath(input_file),
    "binary_folder": os.path.join(execution_folder ,"prokka/binaries/linux/"),
    "work_folder": os.path.abspath(exucution_folder),
    "output_folder": os.path.abspath(output_folde),
    "aragorn_genetic_code": "11",
    "barrnap_kingdom": "bac",
    "blastp_is" : os.path.abspath(execution_folder, "prokka/db/kingdom/Bacteria/IS"),
    "blastp_amr" : os.path.abspath(execution_folder, "prokka/db/kingdom/Bacteria/AMR"),
    "blastp_sprot" : os.path.abspath("./prokka/db/kingdom/Bacteria/sprot"),
}

    tools_settings = {
        "aragorn_options": {
            "aragorn_binary": os.path.join(settings["binary_folder"], "aragorn"),
            "aragorn_input": settings["input_file"],
            "aragorn_output": os.path.join(settings["output_folder"], "aragorn.out"),
            "aragorn_options": "-l -gc%(aragorn_genetic_code)s -w" % settings,
        },
        "barrnap_options": {
            "barrnap_binary": os.path.join(settings["binary_folder"], "barrnap"),
            "barrnap_options": "--kingdom %(barrnap_kingdom)s --threads %(threads)s --quiet" % settings,
            "barrnap_output": os.path.join(settings["output_folder"], "barrnap.out"),
            "barrnap_input": settings["input_file"],
        },
        "parallel_options": {
            "parallel_binary": "parallel",
            "parallel_input": settings["input_file"],
            "parallel_options": "--gnu --plain -j %(threads)s --block 22707 --recstart '>' --pipe" % settings,
        },
        "prodigal_options": {
            "prodigal_binary": os.path.join(settings["binary_folder"], "prodigal"),
            "prodigal_input": settings["input_file"],
            "prodigal_output": os.path.join(settings["output_folder"], "prodigal_out.faa"),
            "prodigal_options": "-c -m -g 11 -p single -f sco -q -a" % settings,
        },
        "blastp_options": {
            "blastp_binary": os.path.join(settings["binary_folder"], "blastp"),
            "blastp_is_input": settings["input_file"],
            "blastp_is_output": os.path.join(settings["output_folder"], "is.blastp"),
            "blastp_amr_output": os.path.join(settings["output_folder"], "amr.blastp"),
            "blastp_sprot_output": os.path.join(settings["output_folder"], "sprot.blastp"),
            "blastp_is_options": "-query - -db %(blastp_is)s -evalue 1e-30 -qcov_hsp_perc 90 \
            -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
            "blastp_amr_options": "-query - -db %(blastp_amr)s -evalue 9.99999999999999e-301 -qcov_hsp_perc 90 \
            -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
            "blastp_sprot_options": "-query - -db %(blastp_sprot)s -evalue 1e-09 -qcov_hsp_perc 80 \
            -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no" % settings,
        },
        "makeblastdb_options": {
            "makeblastdb_binary": os.path.join(settings["binary_folder"], "makeblastdb"),
            "makeblastdb_is_input": settings["blastp_is"],
            "makeblastdb_amr_input": settings["blastp_amr"],
            "makeblastdb_sprot_input": settings["blastp_sprot"],
            "makeblastdb_options": "-hash_index -dbtype prot -in",    
        }
    }

    main(settings)