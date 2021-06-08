#!/usr/bin/env python
import ast

import numpy
import json
import time
import pandas as p
from Bio import Entrez

species_to_keep = \
    ['Acauda', 'Cnidaria', 'Neobipteria', 'Agarella', 'Neohenneguya', 'Alatospora', 'Neomyxobolus', 'Auerbachia',
     'Neoparvicapsula', 'Bipteria', 'Neothelohanellus', 'Buddenbrockia', 'Noblea', 'Cardimyxobolus', 'Octospina',
     'Caudomyxum', 'Ortholinea', 'Ceratomyxa', 'Palliatus', 'Ceratonova', 'Paramyxoproteus', 'Chloromyxum',
     'Parvicapsula', 'Coccomyxa', 'Phlogospora', 'Cystodiscus', 'Pseudalatospora', 'Dicauda', 'Renispora', 'Ellipsomyxa'
        , 'Schulmania', 'Enteromyxum', 'Sigmomyxa', 'Fabespora', 'Sinuolinea', 'Gadimyxa', 'Soricimyxum', 'Globospora',
     'Sphaeromyxa', 'Hennegoides', 'Sphaerospora', 'Henneguya', 'Spirosuturia', 'Hoferellus', 'Tetracapsuloides',
     'Kentmoseria', 'Tetrauronema', 'Kudoa', 'Thelohanellus', 'Laterocaudata', 'Triangula', 'Latyspora', 'Trigonosporus'
        , 'Meglitschia', 'Trilospora', 'Myxidium', 'Trilosporoides', 'Myxobilatus', 'Unicapsula', 'Myxobolus',
     'Unicauda', 'Myxodavisia', 'Wardia', 'Myxoproteus', 'Zschokkella']

minimal_match_length = 100
minimal_identity_percent = 80
entrez_calls_counter = 0


def classify_sequence(self, lineage, distribution_dictionary, row):
    identity_percent = row['pident']
    match_length = row['length']
    contig_length = row['qlen']
    lengths_ratio = match_length / contig_length
    if lineage == 'unknown':
        distribution_dictionary['no_match'] += 1
    else:
        if lineage not in distribution_dictionary.keys():
            # Dictionary initialization
            lineage_dictionary = {'greater identity': {'greater ratio': 0, 'smaller ratio': 0},
                                   'smaller identity': {'greater ratio': 0, 'smaller ratio': 0}}
            distribution_dictionary.update({lineage: lineage_dictionary})
        # Populating the lineages counts dictionary according to identity percent and lengths ratio
        if identity_percent >= self.identity_threshold:
            if lengths_ratio >= self.ratio_threshold:
                distribution_dictionary[lineage]['greater identity']['greater ratio'] += 1
            else:
                distribution_dictionary[lineage]['greater identity']['smaller ratio'] += 1
        elif lengths_ratio >= self.ratio_threshold:
            distribution_dictionary[lineage]['smaller identity']['greater ratio'] += 1
        else:
            distribution_dictionary[lineage]['smaller identity']['smaller ratio'] += 1
    return distribution_dictionary


def connect_entrez(item, entrez_calls_counter):
    Entrez.email = 'e@mail.com'
    entrez_calls_counter += 1
    try:
        return Entrez.read(Entrez.efetch(db='taxonomy', id=item))[0]
    except Exception as e:
        print(e)
        print(item)
        if entrez_calls_counter == 4:
            return
        time.sleep(4.0)
        print(f'{entrez_calls_counter + 1} attempt to get taxonomy from Entrez')
        connect_entrez(item, entrez_calls_counter)


def call_taxonomic_db(item):
    taxonomy = connect_entrez(item, entrez_calls_counter)
    try:
        lineage_list = taxonomy['Lineage'].split('; ')
        if 'Myxozoa' in lineage_list:
            # print('Myxozoa')
            return 'Myxozoa'
        else:
            return taxonomy['ScientificName']
    except Exception as e:
        print(taxonomy)
        return


def calculate_time_passed(time_stamp):
    passed_time = round(time.time() - time_stamp)
    time_stamp = time.time()
    return time_stamp, passed_time


def get_taxonomic_data(dataframe):
    """
    Getting taxids from the blast output and taxonomy from taxonomy.txt and from call_taxonomic_db(taxid)
    :param dataframe
    :return dataframe
    """
    # ###### finding cells with multiple taxids and keeping only one:
    dataframe.staxids = dataframe.apply(lambda x: x.staxids.split(';')[0] if ';' in x.staxids else x.staxids, axis=1)
    dataframe.staxids = dataframe.staxids.astype('int32')
    dataframe.lineage.fillna('nan', inplace=True)
    # ###### Getting taxonomic data ########
    # Getting taxonomy dictionary from taxonomy.txt
    with open('taxonomy.txt', 'r') as taxonomy_file:
        file_dictionary = {int(key): value for key, value in ast.literal_eval(taxonomy_file.read()).items()}
    original_dictionary_size = len(file_dictionary)
    # Inserting lineages from file_dictionary (by taxid) into dataframe.lineage. If a taxid in the dataframe is not
    # found in file_dictionary, the lineage of that taxid is retrieved from Entrez database, put inserted in the
    # dataframe, and add to file_dictionary.
    for index, row in dataframe.iterrows():
        if row.staxids in file_dictionary.keys():
            dataframe.at[index, 'lineage'] = file_dictionary[row.staxids]
        else:
            dataframe.at[index, 'lineage'] = call_taxonomic_db(row.staxids)
            file_dictionary.update({row.staxids: dataframe.at[index, 'lineage']})

    # Writing taxids and their respective lineages to file_dataframe if they are not already in it.
    if len(file_dictionary) > original_dictionary_size:
        with open('taxonomy.txt', 'w') as writing_variable:
            writing_variable.write(json.dumps(file_dictionary))
    # Returning dataframe rows in which the lineage is not Myxozoa:
    return dataframe[dataframe['lineage'] != 'Myxozoa']


def get_dataframe(blast_results_path):
    dataframe = p.read_csv(blast_results_path, index_col=False, sep='\t', header=None, low_memory=False)
    dataframe.columns = ['qseqid', 'sseqid', 'pident', 'staxids', 'lineage', 'qstart', 'qend', 'qlen', 'length',
                         'sstart', 'send', 'slen', 'evalue', 'mismatch', 'gapopen', 'bitscore', 'stitle']
    dataframe = dataframe.drop(['sseqid', 'sstart', 'send', 'mismatch', 'gapopen', 'bitscore'], axis=1)
    dataframe[['pident', 'evalue']] = dataframe[['pident', 'evalue']].astype('float16')
    dataframe[['length', 'slen', 'qlen', 'qend', 'qstart']] = dataframe[['length', 'slen', 'qlen', 'qend', 'qstart']].astype('int32')
    dataframe.lineage = numpy.nan
    dataframe.staxids = dataframe.staxids.astype('string')
    # Removing rows from the dataframe in which the hits are smaller than 100 bp, or with identity that's less than
    # 80%, and the substring 'PREDICTED: ' from the species' description.
    # Removed rows WON'T be deleted from the fasta file, since we don't know what species they belong to, and could be
    # from Thelohanellus
    dataframe = dataframe[dataframe.length >= minimal_match_length]
    dataframe = dataframe[dataframe.pident >= minimal_identity_percent]
    dataframe.stitle = dataframe.apply(lambda x: x.stitle.replace('PREDICTED: ', '') if 'PREDICTED: ' in x.stitle else
                                       x.stitle, axis=1)
    return dataframe
