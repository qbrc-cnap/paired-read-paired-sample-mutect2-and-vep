#!/usr/bin/python3

import subprocess as sp
import json
import argparse
import pandas as pd
import sys
import os
from jinja2 import Environment, FileSystemLoader

# some variables for common reference.
# many refer to keys set in the json file of
# config variables
TEMPLATE = 'template'
OUTPUT = 'output'
CFG = 'config_vars'
R1 = 'r1_files'
R2 = 'r2_files'
GIT_REPO = 'git_repo'
GIT_COMMIT = 'git_commit'
GENOME = 'genome'


class InputDisplay(object):
    '''
    A simple object to carry info to the markdown report.
    '''
    def __init__(self, sample_name, r1, r2):
        self.sample_name = sample_name
        self.r1 = r1
        self.r2 = r2

class AnnotationDisplay(object):
    '''
    A simple object to carry info to the markdown report.
    '''
    def __init__(self, sample_name, condition):
        self.name = sample_name
        self.condition = condition



def get_jinja_template(template_path):
    '''
    Returns a jinja template to be filled-in
    '''
    template_dir = os.path.dirname(template_path)
    env = Environment(loader=FileSystemLoader(template_dir), lstrip_blocks=True, trim_blocks=True)
    return env.get_template(
        os.path.basename(template_path)
    )


def run_cmd(cmd, return_stderr=False):
    '''
    Runs a command through the shell
    '''
    p = sp.Popen(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    stdout, stderr = p.communicate()
    if return_stderr:
        return stderr.decode('utf-8')
    return stdout.decode('utf-8')


def get_versions():
    '''
    Runs the command to get the versions
    '''
    bwa_version_str = run_cmd('bwa', return_stderr=True)
    samtools_version_str = run_cmd('samtools --version')
    multiqc_version_str = run_cmd('multiqc --version')
    fastqc_version_str = run_cmd('fastqc --version')
    picard_md_str = run_cmd('java -jar /opt/software/picard/picard.jar MarkDuplicates --version', return_stderr=True)
    picard_am_str = run_cmd('java -jar /opt/software/picard/picard.jar CollectAlignmentSummaryMetrics --version', return_stderr=True)
    
    # post-process these to extract just the version number:
    bwa_version = bwa_version_str.split('\n')[2].split(' ')[1]
    samtools_version = samtools_version_str.split('\n')[0].split(' ')[1]
    multiqc_version = multiqc_version_str.strip().split(',')[-1].split(' ')[-1]
    fastqc_version = fastqc_version_str.strip().split(" ")[-1]
    picard_md_version = picard_md_str.strip().split(' ')[-1]
    picard_am_version = picard_am_str.strip().split(' ')[-1]

    d = {}
    d['bwa_version'] = bwa_version
    d['samtools_version'] = samtools_version
    d['multiqc_version'] = multiqc_version
    d['fastqc_version'] = fastqc_version
    d['picard_mark_duplicates_version'] = picard_md_version
    d['picard_alignment_metrics_version'] = picard_am_version
    return d


def parse_input():
    '''
    Parses the commandline input, returns a dict
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', required=True, dest=TEMPLATE)
    parser.add_argument('-o', required=True, dest=OUTPUT)
    parser.add_argument('-j', required=True, dest=CFG)
    parser.add_argument('-r1', required=True, dest=R1, nargs='+')
    parser.add_argument('-r2', required=True, dest=R2, nargs='+')

    args = parser.parse_args()
    return vars(args)


def fill_template(context, template_path, output):
    if os.path.isfile(template_path):
        template = get_jinja_template(template_path)
        with open(output, 'w') as fout:
            fout.write(template.render(context))
    else:
        print('The report template was not valid: %s' % template_path)
        sys.exit(1)


if __name__ == '__main__':

    # parse commandline args and separate the in/output from the
    # context variables
    arg_dict = parse_input()
    output_file = arg_dict.pop(OUTPUT)
    input_template_path = arg_dict.pop(TEMPLATE)

    # get all the software versions:
    versions_dict = get_versions()

    # parse the json file which has additional variables
    j = json.load(open(arg_dict[CFG]))

    # alter how the files are displayed:
    r1_files = arg_dict[R1]
    r1_files = [r.split('/')[-1] for r in r1_files]
    r2_files = arg_dict[R2]
    r2_files = [r.split('/')[-1] for r in r2_files]
    samples = [os.path.basename(x)[:-len('_R1.fastq.gz')] for x in r1_files]
    file_display = []
    for r1, r2, s in zip(r1_files, r2_files, samples):
        ipd = InputDisplay(s, r1, r2)
        file_display.append(ipd)


    # make the context dictionary
    context = {}
    context.update(versions_dict)
    context.update(arg_dict)
    context.update(j)
    
    context.update({'file_display': file_display})

    # fill and write the completed report:
    fill_template(context, input_template_path, output_file)
