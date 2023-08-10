# Routine Nanopore QC Collector

Collects QC metrics from the [BCCDC-PHL/routine-nanopore-qc](https://github.com/BCCDC-PHL/routine-nanopore-qc) pipeline.
Outputs are supplied to the [BCCDC-PHL Routine Nanopore QC site](https://github.com/BCCDC-PHL/routine-nanopore-qc-site)

## Usage

```
usage: routine-nanopore-qc-collector [-h] [-c CONFIG] [--log-level LOG_LEVEL]

options:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
  --log-level LOG_LEVEL
```

```
routine-nanopore-qc-collector -c config.json
```

## Configuration

The tool takes a single config file, in json format. A `config_template.json` is provided in this repo:

```json
{
    "analysis_by_run_dir": "/path/to/analysis_by_run",
    "excluded_runs_list": "/path/to/excluded_runs.csv",
    "known_species_list": "/path/to/known_species.csv",
    "scan_interval_seconds": 3600,
    "output_dir": "/path/to/routine-nanopore-qc-collector/data"
}
```