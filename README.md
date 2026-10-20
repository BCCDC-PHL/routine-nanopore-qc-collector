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

## Output

The tool writes a file named `runs.json` to the directory listed as `output_dir` in the config. 
That file has the following format:

```json
[
    {
        "run_id": "20261031_1336_X1_FBA28306_e2f5b1c3",
        "instrument_type": "gridion",
        "run_qc_check": {
            "checked_metrics": [
                {
                    "metric": "PercentReadsPassed",
                    "threshold": 70,
                    "pass_above_or_below": "above",
                    "value": 87.351,
                    "pass_fail": "PASS"
                },
                {
                    "metric": "ReadN50",
                    "threshold": 4000,
                    "pass_above_or_below": "above",
                    "value": 10248,
                    "pass_fail": "PASS"
                },
                {
                    "metric": "NumReadsPassed",
                    "threshold": 100000,
                    "pass_above_or_below": "above",
                    "value": 2032927,
                    "pass_fail": "PASS"
                }
            ],
            "overall_qc_pass_fail": "PASS"
        }
    },
    {
        "run_id: "2026...",
        ...
    }
]
```

It also writes a file named like `<RUN_ID>_library_qc.json` under a `library-qc` directory in the `output_dir`.
That file has the format:

```json
[
  {
    "sample_id": "BC10A123A",
    "barcode": "barcode01",
    "project_id": "test_project",
    "library_id": "BC10A123A",
    "num_reads": 39291,
    "num_bases": 278453360,
    "read_n50": 14409,
    "longest_read": 119339,
    "shortest_read": 20,
    "median_read_length": 3948,
    "median_quality": 16.5,
    "inferred_species_name": "Klebsiella pneumoniae",
    "inferred_genus_name": "Klebsiella",
    "inferred_species_percent": 73.29,
    "inferred_genus_percent": 74.99,
    "inferred_species_genome_size_mb": 5.7,
    "inferred_species_estimated_depth": 35.803,
    "inferred_genus_estimated_depth": 36.634
  },
  {
    "sample_id": "BC10A124A",
	...
  }
  ...
]
```

...and a file named like `<RUN_ID>_species_abundance.json` under a `species-abundance` directory in the `output_dir`.
That file has the format:

```json
[
  {
    "library_id": "BC10A123A",
    "unclassified_fraction_total_reads": 0.0417,
    "abundance_1_genus_name": "Klebsiella",
    "abundance_1_genus_taxid": "570",
    "abundance_1_name": "Klebsiella pneumoniae",
    "abundance_1_fraction_total_reads": 0.7329,
    "abundance_2_genus_name": "Escherichia",
    "abundance_2_genus_taxid": "561",
    "abundance_2_name": "Escherichia coli",
    "abundance_2_fraction_total_reads": 0.08,
    "abundance_3_genus_name": "Klebsiella",
    "abundance_3_genus_taxid": "570",
    "abundance_3_name": "Klebsiella quasipneumoniae",
    "abundance_3_fraction_total_reads": 0.0108,
    "abundance_4_genus_name": "Klebsiella",
    "abundance_4_genus_taxid": "570",
    "abundance_4_name": "Klebsiella variicola",
    "abundance_4_fraction_total_reads": 0.0062,
    "abundance_5_genus_name": "Yersinia",
    "abundance_5_genus_taxid": "629",
    "abundance_5_name": "Yersinia enterocolitica",
    "abundance_5_fraction_total_reads": 0.0044,
    "abundance_6_genus_name": "Shigella",
    "abundance_6_genus_taxid": "620",
    "abundance_6_name": "Shigella flexneri",
    "abundance_6_fraction_total_reads": 0.0042
  },
  {
    "library_id": "BC10A123A",
	...
  }
  ...
]
```
