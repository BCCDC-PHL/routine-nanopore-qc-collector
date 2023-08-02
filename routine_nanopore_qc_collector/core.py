import collections
import csv
import glob
import json
import logging
import os
import re
import shutil
import subprocess

from typing import Iterator, Optional

import routine_nanopore_qc_collector.parsers as parsers
import routine_nanopore_qc_collector.samplesheet as samplesheet


def create_output_dirs(config):
    """
    """
    base_outdir = config['output_dir']
    output_dirs = [
        base_outdir,
        os.path.join(base_outdir, 'library-qc'),
        os.path.join(base_outdir, 'species-abundance'),
    ]
    for output_dir in output_dirs:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)    


def find_latest_routine_nanopore_qc_output(analysis_dir):
    """
    """
    routine_nanopore_qc_output_dir_glob = "routine-nanopore-qc-v*-output"
    routine_nanopore_qc_output_dirs = glob.glob(os.path.join(analysis_dir, routine_nanopore_qc_output_dir_glob))
    latest_routine_nanopore_qc_output_dir = None
    if len(routine_nanopore_qc_output_dirs) > 0:
        latest_routine_nanopore_qc_output_dir = os.path.abspath(routine_nanopore_qc_output_dirs[-1])

    return latest_routine_nanopore_qc_output_dir


def find_analysis_dirs(config, check_complete=True):
    """
    """
    gridion_run_id_regex = "\d{8}_\d{4}_X\d_[A-Z0-9]{8}_[a-z0-9]{8}$"
    promethion_run_id_regex = "\d{8}_\d{4}_P2S_\d{4}-\d_[A-Z0-9]{8}_[a-z0-9]{8}$"
    analysis_by_run_dir = config['analysis_by_run_dir']
    subdirs = os.scandir(analysis_by_run_dir)
    
    for subdir in subdirs:
        run_id = subdir.name
        matches_gridion_regex = re.match(gridion_run_id_regex, run_id)
        matches_promethion_regex = re.match(promethion_run_id_regex, run_id)
        sequencer_type = None
        if matches_gridion_regex:
            sequencer_type = 'miseq'
        elif matches_promethion_regex:
            sequencer_type = 'promethion'
        not_excluded = run_id not in config['excluded_runs']
        ready_to_collect = False
        if check_complete:
            latest_routine_nanopore_qc_output = find_latest_routine_nanopore_qc_output(subdir)
            if latest_routine_nanopore_qc_output is not None and os.path.exists(latest_routine_nanopore_qc_output):
                routine_nanopore_qc_analysis_complete = os.path.exists(os.path.join(latest_routine_nanopore_qc_output, 'analysis_complete.json'))
                ready_to_collect = routine_nanopore_qc_analysis_complete
        else:
            ready_to_collect = True

        conditions_checked = {
            "is_directory": subdir.is_dir(),
            "matches_nanopore_run_id_format": ((matches_gridion_regex is not None) or (matches_promethion_regex is not None)),
            "not_excluded": not_excluded,
            "ready_to_collect": ready_to_collect,
        }
        conditions_met = list(conditions_checked.values())

        analysis_directory_path = os.path.abspath(subdir.path)
        analysis_dir = {
            "path": analysis_directory_path,
            "sequencer_type": sequencer_type,
        }
        if all(conditions_met):
            logging.info(json.dumps({
                "event_type": "analysis_directory_found",
                "sequencing_run_id": run_id,
                "analysis_directory_path": analysis_directory_path
            }))

            yield analysis_dir
        else:
            logging.debug(json.dumps({
                "event_type": "directory_skipped",
                "analysis_directory_path": os.path.abspath(subdir.path),
                "conditions_checked": conditions_checked
            }))
            yield None

            
def find_runs(config):
    """
    Finda all runs that have routine sequence QC data.

    :param config: Application config.
    :type config: dict[str, object]
    :return: List of runs. Keys: ['run_id', 'sequencer_type']
    :rtype: list[dict[str, str]]
    """
    logging.info(json.dumps({"event_type": "find_runs_start"}))
    runs = []
    all_analysis_dirs = sorted(list(os.listdir(config['analysis_by_run_dir'])))
    all_run_ids = filter(lambda x: re.match('\d{8}_\d{4}_', x) != None, all_analysis_dirs)
    for run_id in all_run_ids:
        if run_id in config['excluded_runs']:
            continue

        sequencer_type = None
        if re.match('\d{8}_\d{4}_X', run_id):
            sequencer_type = 'gridion'
        elif re.match('\d{8}_\d{4}_P', run_id):
            sequencer_type = 'promethion'

        analysis_dir = os.path.join(config['analysis_by_run_dir'], run_id)
        latest_routine_nanopore_qc_output_dir = find_latest_routine_nanopore_qc_output(analysis_dir)
        
        if os.path.exists(os.path.join(latest_routine_nanopore_qc_output_dir, 'analysis_complete.json')):
            run = {
                'run_id': run_id,
                'sequencer_type': sequencer_type,
            }
            runs.append(run)

    logging.info(json.dumps({
        "event_type": "find_runs_complete"
    }))

    return runs


def scan(config: dict[str, object]) -> Iterator[Optional[dict[str, str]]]:
    """
    Scanning involves looking for all existing runs and...

    :param config: Application config.
    :type config: dict[str, object]
    :return: A run directory to analyze, or None
    :rtype: Iterator[Optional[dict[str, object]]]
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    for analysis_dir in find_analysis_dirs(config):    
        yield analysis_dir


def infer_species(config: dict[str, object], species_abundance):
    """
    :return: Name of inferred species.
    :rtype: str
    """
    species_abundance_keys = [
        'abundance_5_',
        'abundance_4_',
        'abundance_3_',
        'abundance_2_',
        'abundance_1_',
    ]
    inferred_species = None
    greatest_fraction_total_reads = 0.0
    for k in species_abundance_keys:
        if k + 'name' in species_abundance and species_abundance[k + 'name'] != 'Homo sapiens' and species_abundance[k + 'fraction_total_reads'] > greatest_fraction_total_reads:
            inferred_species = species_abundance[k + 'name']
            greatest_fraction_total_reads = species_abundance[k + 'fraction_total_reads']

    return inferred_species


def infer_genus(config: dict[str, object], species_abundance, inferred_species):
    """
    :return: Name of inferred genus.
    :rtype: str
    """
    species_abundance_keys = [
        'abundance_5_',
        'abundance_4_',
        'abundance_3_',
        'abundance_2_',
        'abundance_1_',
    ]
    inferred_genus = None
    for k in species_abundance_keys:
        if k + 'name' in species_abundance and species_abundance[k + 'name'] == inferred_species:
            inferred_genus = species_abundance[k + 'genus_name']

    return inferred_genus



def get_percent_reads_by_species_name(species_abundance, species_name):
    """
    """
    percent_reads = None
    species_abundance_keys = [
        'abundance_1_',
        'abundance_2_',
        'abundance_3_',
        'abundance_4_',
        'abundance_5_',
    ]

    for k in species_abundance_keys:
        try:
            if species_abundance[k + 'name'] == species_name:
                percent_reads = round(100 * species_abundance[k + 'fraction_total_reads'], 3)
        except KeyError as e:
            pass

    return percent_reads        


def get_percent_reads_by_genus_name(species_abundance, genus_name):
    """
    """
    percent_reads = 0.0
    species_abundance_keys = [
        'abundance_1_',
        'abundance_2_',
        'abundance_3_',
        'abundance_4_',
        'abundance_5_',
    ]

    for k in species_abundance_keys:
        try:
            if species_abundance[k + 'genus_name'] == genus_name:
                percent_reads += 100 * species_abundance[k + 'fraction_total_reads']
        except KeyError as e:
            pass

    percent_reads = round(percent_reads, 3)

    return percent_reads


def add_genus(kraken_species_record):
    """
    """
    taxid = kraken_species_record['ncbi_taxonomy_id']
    if taxid != '0':
        echo_cmd = [
            'echo',
            taxid,
        ]
        taxonkit_cmd = [
            'taxonkit',
            'reformat',
            '-F',
            '-I', '1',
            '-f', '"{g}|{s}"',
            '-t',
        ]
        try:
            echo_result = subprocess.run(echo_cmd, text=True, capture_output=True, check=True)
            taxonkit_result = subprocess.run(taxonkit_cmd, input=echo_result.stdout, text=True, capture_output=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(json.dumps({"event_type": "add_genus_failed", "ncbi_taxonomy_id": taxid}))
        taxonkit_output = list(map(lambda x: x.strip('"'), taxonkit_result.stdout.strip().split('\t')[1:]))
        genus_name = taxonkit_output[0].split('|')[0]
        genus_taxid = taxonkit_output[1].split('|')[0]
        kraken_species_record['genus_taxon_name'] = genus_name
        kraken_species_record['genus_ncbi_taxonomy_id'] = genus_taxid
    
    return kraken_species_record
    
    
def collect_outputs(config: dict[str, object], analysis_dir: Optional[dict[str, str]]):
    """
    Collect all routine sequence QC outputs for a specific analysis dir.

    :param config: Application config.
    :type config: dict[str, object]
    :param analysis_dir: Analysis dir. Keys: ['path', 'sequencer_type']
    :type analysis_dir: dict[str, str]
    :return: 
    :rtype: 
    """
    run_id = os.path.basename(analysis_dir['path'])
    logging.info(json.dumps({"event_type": "collect_outputs_start", "sequencing_run_id": run_id, "analysis_dir_path": analysis_dir['path']}))

    libraries_by_library_id = {}
    latest_routine_nanopore_qc_output_path = find_latest_routine_nanopore_qc_output(analysis_dir['path'])
    routine_nanopore_qc_output_dir_contents = os.listdir(latest_routine_nanopore_qc_output_path)
    for library_output_dir in routine_nanopore_qc_output_dir_contents:
        if os.path.isdir(os.path.join(latest_routine_nanopore_qc_output_path, library_output_dir)):
            library_id = os.path.basename(library_output_dir)
            libraries_by_library_id[library_id] = {}
    
    # species-abundance
    species_abundance_by_library_id = {library_id: {'library_id': library_id} for library_id in libraries_by_library_id.keys()}
    species_abundance_dst_file = os.path.join(config['output_dir'], "species-abundance", run_id + "_species_abundance.json")
    if not os.path.exists(species_abundance_dst_file):
        for library_id in species_abundance_by_library_id.keys():
            kraken_species_src_file = os.path.join(latest_routine_nanopore_qc_output_path, library_id, library_id + '_kraken2_species.csv')
            if os.path.exists(kraken_species_src_file):
                kraken_species = parsers.parse_kraken_species(kraken_species_src_file)
                library_species_abundance = {'library_id': library_id}
                abundance_num = 1
                for kraken_species_record in kraken_species[0:7]:
                    if kraken_species_record['rank_code'] == 'U':
                        library_species_abundance['unclassified_fraction_total_reads'] = round(kraken_species_record['percent_seqs_in_clade'] / 100, 6)
                    else:
                        kraken_species_record = add_genus(kraken_species_record)
                        if 'genus_taxon_name' in kraken_species_record:
                            logging.info(json.dumps({"event_type": "add_genus_complete", "sequencing_run_id": run_id, "library_id": library_id, "species": kraken_species_record['taxon_name'], "genus": kraken_species_record['genus_taxon_name']}))
                        library_species_abundance['abundance_' + str(abundance_num) + '_name'] = kraken_species_record['taxon_name']
                        library_species_abundance['abundance_' + str(abundance_num) + '_genus_name'] = kraken_species_record['genus_taxon_name']
                        library_species_abundance['abundance_' + str(abundance_num) + '_genus_taxid'] = kraken_species_record['genus_ncbi_taxonomy_id']
                        library_species_abundance['abundance_' + str(abundance_num) + '_fraction_total_reads'] = round(kraken_species_record['percent_seqs_in_clade'] / 100, 6)
                        abundance_num += 1
                species_abundance_by_library_id[library_id] = library_species_abundance

        with open(species_abundance_dst_file, 'w') as f:
            json.dump(list(species_abundance_by_library_id.values()), f, indent=2)

        logging.info(json.dumps({
            "event_type": "write_species_abundance_complete",
            "run_id": run_id,
            "dst_file": species_abundance_dst_file
        }))

    # library-qc
    library_qc_dst_file = os.path.join(config['output_dir'], "library-qc", run_id + "_library_qc.json")
    if not os.path.exists(library_qc_dst_file):
        for library_id in libraries_by_library_id:
            nanoq_path = os.path.join(latest_routine_nanopore_qc_output_path, library_id, library_id + '_nanoq.csv')
            if os.path.exists(nanoq_path):
                nanoq_report = parsers.parse_nanoq(nanoq_path)
                if len(nanoq_report) == 1:
                    nanoq = nanoq_report[0]
                    libraries_by_library_id[library_id]['library_id'] = library_id
                    libraries_by_library_id[library_id]['num_reads'] = nanoq['reads']
                    libraries_by_library_id[library_id]['num_bases'] = nanoq['bases']
                    libraries_by_library_id[library_id]['read_n50'] = nanoq['n50']
                    libraries_by_library_id[library_id]['longest_read'] = nanoq['longest']
                    libraries_by_library_id[library_id]['shortest_read'] = nanoq['shortest']
                    libraries_by_library_id[library_id]['median_read_length'] = nanoq['median_length']
                    libraries_by_library_id[library_id]['median_quality'] = nanoq['median_quality']

                inferred_species = None
                inferred_species = infer_species(config, species_abundance_by_library_id[library_id])
                inferred_genus = None
                inferred_genus = infer_genus(config, species_abundance_by_library_id[library_id], inferred_species)
                if inferred_species is not None:
                    logging.debug(json.dumps({'event_type': 'library_species_inferred', 'sequencing_run_id': run_id, 'library_id': library_id, 'inferred_species': inferred_species}))
                    libraries_by_library_id[library_id]['inferred_species_name'] = inferred_species
                    libraries_by_library_id[library_id]['inferred_genus_name'] = inferred_genus
                    percent_inferred_species = get_percent_reads_by_species_name(species_abundance_by_library_id[library_id], inferred_species)
                    if percent_inferred_species is None:
                        logging.error(json.dumps({"event_type": "collect_library_qc_metric_failed", "metric": "inferred_species_percent", 'library_id': library_id, 'inferred_species': inferred_species}))
                    libraries_by_library_id[library_id]['inferred_species_percent'] = percent_inferred_species
                    percent_inferred_genus = get_percent_reads_by_genus_name(species_abundance_by_library_id[library_id], inferred_genus)
                    if percent_inferred_genus is None:
                        logging.error(json.dumps({"event_type": "collect_library_qc_metric_failed", "metric": "inferred_genus_percent", 'library_id': library_id, 'inferred_genus': inferred_genus}))
                    libraries_by_library_id[library_id]['inferred_genus_percent'] = percent_inferred_genus
                    if 'known_species' in config and inferred_species in config['known_species']:
                        inferred_species_genome_size = config['known_species'][inferred_species]['genome_size_mb']
                        libraries_by_library_id[library_id]['inferred_species_genome_size_mb'] = inferred_species_genome_size
                    else:
                        logging.debug(json.dumps({'event_type': 'library_species_inference_failed', 'sequencing_run_id': run_id, 'library_id': library_id, 'inferred_species': inferred_species}))

                    if all(k in libraries_by_library_id[library_id] for k in ['num_bases', 'inferred_species_genome_size_mb', 'inferred_species_percent']):
                        num_bases = libraries_by_library_id[library_id]['num_bases']
                        genome_size = libraries_by_library_id[library_id]['inferred_species_genome_size_mb'] * 1000000
                        species_percent = libraries_by_library_id[library_id]['inferred_species_percent']
                        if all([num_bases, genome_size, species_percent]):
                            libraries_by_library_id[library_id]['inferred_species_estimated_depth'] = round((num_bases * (species_percent / 100)) / genome_size, 3)
                        if 'inferred_genus_percent' in libraries_by_library_id[library_id]:
                            genus_percent = libraries_by_library_id[library_id]['inferred_genus_percent']
                        if all([num_bases, genome_size, genus_percent]):
                            libraries_by_library_id[library_id]['inferred_genus_estimated_depth'] = round((num_bases * (genus_percent / 100)) / genome_size, 3)

        with open(library_qc_dst_file, 'w') as f:
            json.dump(list(libraries_by_library_id.values()), f, indent=2)

        logging.info(json.dumps({
            "event_type": "write_library_qc_complete",
            "run_id": run_id,
            "dst_file": library_qc_dst_file
        }))




    logging.info(json.dumps({"event_type": "collect_outputs_complete", "sequencing_run_id": run_id, "analysis_dir_path": analysis_dir['path']}))
