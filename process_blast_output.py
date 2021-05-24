#!/usr/bin/env python
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


def classify_sequence(self, division, distribution_dictionary, row):
    identity_percent = row['pident']
    match_length = row['length']
    contig_length = row['qlen']
    lengths_ratio = match_length / contig_length
    if division == 'unknown':
        distribution_dictionary['no_match'] += 1
    else:
        if division not in distribution_dictionary.keys():
            # Dictionary initialization
            division_dictionary = {'greater identity': {'greater ratio': 0, 'smaller ratio': 0},
                                   'smaller identity': {'greater ratio': 0, 'smaller ratio': 0}}
            distribution_dictionary.update({division: division_dictionary})
        # Populating the divisions counts dictionary according to identity percent and lengths ratio
        if identity_percent >= self.identity_threshold:
            if lengths_ratio >= self.ratio_threshold:
                distribution_dictionary[division]['greater identity']['greater ratio'] += 1
            else:
                distribution_dictionary[division]['greater identity']['smaller ratio'] += 1
        elif lengths_ratio >= self.ratio_threshold:
            distribution_dictionary[division]['smaller identity']['greater ratio'] += 1
        else:
            distribution_dictionary[division]['smaller identity']['smaller ratio'] += 1
    return distribution_dictionary


def call_taxonomic_db(item):
    Entrez.email = 'e@mail.com'
    try:
        return Entrez.read(Entrez.efetch(db='taxonomy', id=item))[0]['Division']
    except Exception as e:
        print(e)
        print(item)
        time.sleep(2.0)
        try:
            print('Second attempt to get taxonomy from Entrez')
            return Entrez.read(Entrez.efetch(db='taxonomy', id=item))[0]['Division']
        except Exception as e:
            print('Third attempt to get taxonomy from Entrez')
            print(e)
            print(item)
            time.sleep(2.0)
            return Entrez.read(Entrez.efetch(db='taxonomy', id=item))[0]['Division']


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
    dataframe.division.fillna('nan', inplace=True)
    # ###### Getting taxonomic data ########
    # Getting taxonomy dictionary from taxonomy.txt
    with open('taxonomy.txt', 'r') as taxonomy_file:
        file_list = taxonomy_file.read().split(',')[:-1]
    file_dictionary = {int(item.split(':')[0]): item.split(':')[1] for item in file_list}
    original_dictionary_size = len(file_dictionary)
    # Todo: after running the code on B1, improve speed by:
    #  Change dtype to int16 from -32768 to 32767 or int32 from -2147483648 to 2147483647
    # Inserting divisions from file_dictionary (by taxid) into dataframe.division. If a taxid in the dataframe is not
    # found in file_dictionary, the division of that taxid is retrieved from Entrez database, put inserted in the
    # dataframe, and add to file_dictionary.
    for index, row in dataframe.iterrows():
        if row.staxids in file_dictionary.keys():
            dataframe.at[index, 'division'] = file_dictionary[row.staxids]
        else:
            dataframe.at[index, 'division'] = call_taxonomic_db(row.staxids)
            file_dictionary.update({row.staxids: dataframe.at[index, 'division']})

    # Writing taxids and their respective divisions to file_dataframe, in case they are not found in it.
    if len(file_dictionary) > original_dictionary_size:
        with open('taxonomy.txt', 'w') as writing_variable:
            for key, item in file_dictionary.items():
                writing_variable.write(f'{key}: {item}, ')

    # Returning dataframe rows in which the division is not Myxozoa:
    return dataframe[~dataframe['division'].isin(species_to_keep)]


def get_dataframe(blast_results_path):
    dataframe = p.read_csv(blast_results_path, index_col=False, sep='\t', header=None, low_memory=False)
    dataframe.columns = ['qseqid', 'sseqid', 'pident', 'staxids', 'division', 'qstart', 'qend', 'qlen', 'length',
                         'sstart', 'send', 'slen', 'evalue', 'mismatch', 'gapopen', 'bitscore', 'stitle']
    dataframe = dataframe.drop(['sseqid', 'sstart', 'send', 'mismatch', 'gapopen', 'bitscore'], axis=1)
    dataframe[['pident', 'evalue']] = dataframe[['pident', 'evalue']].astype('float16')
    dataframe[['length', 'slen', 'qlen', 'qend', 'qstart']] = dataframe[['length', 'slen', 'qlen', 'qend', 'qstart']].astype('int32')
    dataframe.division = numpy.nan
    # Removing rows from the dataframe in which the hits are smaller than 100 bp, or with identity that's less than
    # 80%, and the substring 'PREDICTED: ' from the species' description.
    # Removed rows WON'T be deleted from the fasta file, since we don't know to what species they belong, and could be
    # from Thelohanellus
    dataframe = dataframe[dataframe.length >= minimal_match_length]
    dataframe = dataframe[dataframe.pident >= minimal_identity_percent]
    dataframe.stitle = dataframe.apply(lambda x: x.stitle.replace('PREDICTED: ', '') if 'PREDICTED: ' in x.stitle else
                                       x.stitle, axis=1)
    return dataframe
