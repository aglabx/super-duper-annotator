#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: <data>
#@author: <name>
#@contact: <email>

import sys
import argparse
from settings import settings, tools_settings
from scripts import aragorn
from scripts import barrnap

def input_checker(file):
    #checks for existing and size
    pass

def output_checker(file):
    #check file for not existing before run
    pass
    
def tool_runner(tool_command):
    #check for input and output files - it should be separate function
    
    # if input exists and output don't, run command
    print(tool_command)
    os.system(tool_command)


def main(settings, tools_settings):
    
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
    output_file = args["output"]

    main(settings)