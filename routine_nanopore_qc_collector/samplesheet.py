import glob
import json
import logging
import os
import re


def find_samplesheet_for_run(run_id, sequencer_output_dirs):
    """
    """
    samplesheets = []
    samplesheet_path = None

    sequencer_run_dir = ""

    for run_dir in os.listdir(sequencer_output_dir):
        if os.path.basename(run_dir) == run_id:
            sequencer_run_dir = os.path.join(sequencer_output_dir, run_dir)
            samplesheet_glob = os.path.join(sequencer_run_dir, 'sample_sheet_*.csv')
            samplesheets = glob.glob(samplesheet_glob)

        if len(samplesheets) == 1:
            samplesheet_path = samplesheets[0]

    return samplesheet_path


def parse_samplesheet(samplesheet_path):
    """
    """
    samplesheet = []
    return samplesheet
