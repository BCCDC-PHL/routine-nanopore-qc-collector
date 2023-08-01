import collections
import re
import csv


def parse_nanoq(nanoq_path):
    """
    """
    nanoq = []
 
    int_fields = [
        'reads',
        'bases',
        'n50',
        'longest',
        'shortest',
        'mean_length',
        'median_length',
    ]

    float_fields = [
        'mean_quality',
        'median_quality',
    ]

    with open(nanoq_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for field in int_fields:
                try:
                    row[field] = int(row[field])
                except ValueError as e:
                    row[field] = None
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    row[field] = None

            nanoq.append(row)

    return nanoq


def parse_kraken_species(kraken_species_path):
    """
    """
    kraken_species = []
    int_fields = [
        'num_seqs_in_clade',
        'num_seqs_this_taxon',
    ]

    float_fields = [
        'percent_seqs_in_clade',
    ]

    with open(kraken_species_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for field in int_fields:
                try:
                    row[field] = int(row[field])
                except ValueError as e:
                    row[field] = None
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    row[field] = None

            kraken_species.append(row)

    return kraken_species
