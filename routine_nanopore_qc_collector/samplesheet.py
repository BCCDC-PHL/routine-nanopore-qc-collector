import csv
import logging
import os

log = logging.getLogger(__name__)

def find_parsed_samplesheet(run_id: str, analysis_dir: os.PathLike):
    """
    Find the parsed samplesheet for a run.
    The parsed samplesheet has fields: 'sample_id', 'barcode', 'project_id'
    """
    samplesheet_path = os.path.join(analysis_dir, f"{run_id}_sample_sheet.parsed.csv")

    if not os.path.exists(samplesheet_path):
        log.error({
            "event_type": "failed_to_find_file",
            "sequencing_run_id": run_id,
            "file_type": "parsed_samplesheet",
            "path": samplesheet_path,
        })
        samplesheet_path = None

    return samplesheet_path


def parse_samplesheet(samplesheet_path: os.PathLike) -> list[dict]:
    """
    Parse the samplesheet file into a list of dicts.
    """
    samplesheet = []
    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samplesheet.append(row)

    return samplesheet
