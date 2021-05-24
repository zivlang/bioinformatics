#!/usr/bin/env python
import blastn
from process_blast_output import *
from remove_false_sequences import *
import os

# directory = 'dorotheeh/zivlang/output/'


def remove_big_ratios(dataframe, ranges):
    """
    Calculates the portion of each scaffold which is occupied by identified contamination matches.
    :param dataframe:
    :param ranges: all species' matches ranges.
    :return:
    """
    ratio_dictionary = {}
    # Removing duplications in scaffold names ('qseqid')
    scaffolds = dataframe.drop_duplicates(subset='qseqid')['qseqid'].tolist()
    deleted = 0
    for scaffold in scaffolds:
        current_ranges = ranges[scaffold]
        ratio_dictionary[scaffold] = 0
        for match in current_ranges:
            ratio_dictionary[scaffold] += match[1] - match[0]
        # Getting the query's length. Each query may appear multiple times, so a tuple is returned, from which the first
        # item is taken.
        query_length = dataframe.loc[dataframe['qseqid'] == scaffold]['qlen'].tolist()[0]
        # Deleting the item from ranges if the ratio is greater than 90%.
        if ratio_dictionary[scaffold] / query_length >= .9:
            del ranges[scaffold]
            deleted += 1
        ratio_dictionary[scaffold] = ratio_dictionary[scaffold] / query_length
    print('Deleted', deleted, 'scaffolds out of', len(ranges) + deleted, 'scaffolds.')
    return ranges


def print_time_passed(process, time_stamp):
    time_stamp, passed_time = calculate_time_passed(time_stamp)
    if passed_time > 0:
        print(process, 'took', round(passed_time / 60, 2), 'minutes.')
    return time_stamp


def main():
    time_stamp = time.time()
    # for file in os.listdir('fasta_files'):
    #     if file in ['btcaA1.fa', 'btcaB1.fa']:
            # blastn.blast(f'fasta_files/{file}')
    print(file)
    dataframe = get_dataframe(f'blast/blastn/{file[:-2]}txt')
    time_stamp = print_time_passed('get_dataframe()', time_stamp)
    dataframe = get_taxonomic_data(dataframe)
    time_stamp = print_time_passed('get_taxonomic_data()', time_stamp)
    dictionary = get_ranges_dictionary(dataframe)
    time_stamp = print_time_passed('fasta_to_dictionary()', time_stamp)
    merged = merge_overlaps(dictionary)  # Merging each species' hits
    time_stamp = print_time_passed('merge_overlaps()', time_stamp)
    # Merging all the hits of all the species into a list of lists:
    scaffold_ranges = {}
    for key, value in merged.items():
        scaffold_ranges[key] = []
        for inner_key, inner_value in value.items():
            scaffold_ranges[key].append(inner_value[0])
    time_stamp = print_time_passed('scaffold_ranges population()', time_stamp)
    # Now merges can be performed if needed:
    merged_scaffold_ranges = merge_overlaps({'outer_dictionary': scaffold_ranges})
    time_stamp = print_time_passed('merge_overlaps()', time_stamp)
    small_ranges = remove_big_ratios(dataframe, merged_scaffold_ranges['outer_dictionary'])
    time_stamp = print_time_passed('remove_big_ratios()', time_stamp)
    delete_false_hits(small_ranges, f'fasta_files/{file[:-2]}fa')
    time_stamp = print_time_passed('delete_false_hits() and merged_scaffold_ranges()', time_stamp)


if __name__ == '__main__':
    main()
