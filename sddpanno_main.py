#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import os

from os.path import abspath
from os.path import join
from settings import settings
from settings import tools_settings
from inspect import getsourcefile
from scripts import aragorn
from scripts import barrnap
from collections import namedtuple


Tool = namedtuple('Tool', ['name', 'command'])


def input_checker(file):
    # checks for existing and size
    pass


def output_checker(file):
    # check file for not existing before run
    pass


def tool_runner(tool):
    # if input exists and output don't, run command
    print(f'Running {tool.name}')
    errcode = os.system(tool.command)
    msg = 'SUCCESS' if errcode == 0 else 'ERROR'

    print(f"[{msg:5}][{errcode:3}] {tool.name}")


def main(settings, tools_settings):
    # checks for files

    # aragorn

    # here should be run of function which creates comand for tool (create it in scripts folder)
    aragorn_opts = tools_settings['aragorn']
    aragorn_command = aragorn.command(
        aragorn_opts['binary'],
        aragorn_opts['input'],
        aragorn_opts['output'],
        aragorn_opts['options'],
    )

    # barrnap
    barrnap_opts = tools_settings['barrnap']
    barrnap_command = barrnap.command(
        barrnap_opts['binary'],
        barrnap_opts['options'],
        barrnap_opts['output'],
        barrnap_opts['input'],
    )  

    pipeline = [
        Tool('aragorn', aragorn_command)
    ]

    for command in pipeline:
        tool_runner(command)
        # check_for_output


def parse_args():
    parser = argparse.ArgumentParser(description='Program description.')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-o', '--output', help='Output folder', required=True)

    return parser.parse_args()


def build_settings(args):
    input_file = args.input
    output_folder = args.output
    execution_folder = os.path.dirname(abspath(getsourcefile(lambda: 0)))

    # settings override
    return {
        'input_file': abspath(input_file),
        'output_folder': abspath(output_folder),
        'work_folder': execution_folder
    }


if __name__ == '__main__':
    args = parse_args()
    user_settings = build_settings(args)
    settings.update(user_settings)

    main(settings, tools_settings)
