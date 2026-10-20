#!/usr/bin/env python

import argparse
import datetime
import json
import logging
import os
import time

from pathlib import Path

from routine_nanopore_qc_collector.config import load_config
from routine_nanopore_qc_collector.logging_config import configure_logging
import routine_nanopore_qc_collector.core as core

log = logging.getLogger(__name__)

DEFAULT_SCAN_INTERVAL_SECONDS = 3600.0

def reload_config(config_path: os.PathLike, current_config: dict) -> dict:
    """
    Reload the config file, so that changes take effect without restarting.
    If the config can't be loaded, continue with the last valid config.

    :param config_path: Path to the config file.
    :param current_config: The last valid config.
    :return: The reloaded config, or the current config if reloading failed.
    """
    try:
        config = load_config(config_path)
        log.info({"event_type": "config_loaded", "config_file": os.path.abspath(config_path)})
        return config
    except Exception as e:
        log.error({"event_type": "load_config_failed", "config_file": os.path.abspath(config_path), "error": str(e)})
        return current_config


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', type=Path)
    parser.add_argument('--log-level', default='info', choices=['debug', 'info', 'warning', 'error'])
    parser.add_argument('--log-file', help='Append logs to this file, rotated daily at midnight (default: log to stdout)')
    parser.add_argument('--log-retention-days', type=int, default=90, help='Number of rotated log files to keep. 0 keeps all of them (default: 90)')
    args = parser.parse_args()

    configure_logging(args.log_level, args.log_file, args.log_retention_days)
    log.debug({"event_type": "debug_logging_enabled"})

    config = load_config(args.config)
    log.info({
        "event_type": "config_loaded",
        "config_file": os.path.abspath(args.config)
    })

    core.create_output_dirs(config)

    try:
        while(True):
            scan_start_timestamp = datetime.datetime.now()

            runs = core.find_runs(config)
            runs_output_file = os.path.join(config['output_dir'], 'runs.json')
            with open(runs_output_file, 'w') as f:
                json.dump(runs, f, indent=2)
            log.info({"event_type": "write_runs_file_complete", "runs_file": runs_output_file})

            for run in core.scan(config):
                reload_config(args.config, config)
                try:
                    log.debug({
                        "event_type": "collection_starting",
                        "sequencing_run_id": run['sequencing_run_id'],
                    })
                    core.collect_outputs(config, run)
                    log.debug({
                        "event_type": "collection_complete",
                        "sequencing_run_id": run['sequencing_run_id'],
                    })
                # A problem with one run shouldn't stop other runs from being checked.
                # No 'qc_check_complete.json' is written, so the run will be retried on the next scan.
                except OSError as e:
                    # Expected, environmental problems (missing files, permissions, network storage).
                    # The message says what went wrong, so a traceback wouldn't add anything.
                    log.error({
                        "event_type": "collection_failed",
                        "sequencing_run_id": run['sequencing_run_id'],
                        "error_type": type(e).__name__, "error": str(e)
                    })
                except Exception as e:
                    # Unexpected problems are probably bugs, so include the traceback.
                    log.error({
                        "event_type": "collection_failed",
                        "sequencing_run_id": run['sequencing_run_id'],
                        "error_type": type(e).__name__,
                        "error": str(e)
                    }, exc_info=True)

            scan_complete_timestamp = datetime.datetime.now()
            scan_duration_delta = scan_complete_timestamp - scan_start_timestamp
            scan_duration_seconds = scan_duration_delta.total_seconds()
            next_scan_timestamp = scan_start_timestamp + datetime.timedelta(seconds=config['scan_interval_seconds'])
            log.info({
                "event_type": "scan_complete",
                "scan_duration_seconds": scan_duration_seconds,
                "scan_interval_seconds": config['scan_interval_seconds'],
                "timestamp_next_scan_start": next_scan_timestamp.isoformat(),
            })
            time.sleep(max(0, (next_scan_timestamp - datetime.datetime.now()).total_seconds()))
            config = reload_config(args.config, config)
    except KeyboardInterrupt:
        log.info({"event_type": "quit"})

if __name__ == '__main__':
    main()
