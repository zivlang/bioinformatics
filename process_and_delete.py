#!/usr/bin/env python
# import blastn
from process_blast_output import *
from remove_false_sequences import *

# directory = 'dorotheeh/zivlang/output/'


def mark_high_contamination(dataframe, ranges):
    """
    Calculates the portion of each scaffold which is occupied by identified contamination matches.
    :param dataframe:
    :param ranges: all species' matches ranges.
    :return:
    """
    blastn_dictionary = {}
    # Removing duplications in scaffold names ('qseqid') when creating scaffolds' names list
    scaffolds = dataframe.drop_duplicates(subset='qseqid')['qseqid'].tolist()
    contaminated = 0
    for scaffold in scaffolds:
        current_ranges = ranges[scaffold]
        blastn_dictionary[scaffold] = 0
        smallest_hit_length = 0
        for match in current_ranges:
            match_length = match[1] - match[0]
            blastn_dictionary[scaffold] += match_length
            if match_length < smallest_hit_length or smallest_hit_length == 0:
                smallest_hit_length = match_length
        # Getting the query's length. Each query may appear multiple times, so a tuple is returned, from which the first
        # item is taken.
        # Saving the query's length and occupancy of contamination:
        query_length = dataframe.loc[dataframe['qseqid'] == scaffold]['qlen'].tolist()[0]
        # Deleting the item from ranges if the ratio is greater than 90%.
        # print(scaffold, blastn_dictionary[scaffold] / query_length, 'smallest_hit_length:', smallest_hit_length)
        if blastn_dictionary[scaffold] / query_length >= .9:
            ranges[scaffold] = 'contaminated'
            contaminated += 1
    print(contaminated, 'scaffolds out of', len(ranges), 'scaffolds are marked for deletion.')
    return ranges


def print_time_passed(process, time_stamp):
    time_stamp, passed_time = calculate_time_passed(time_stamp)
    if passed_time > 0:
        minutes, sec = divmod(passed_time, 60)
        hours, minutes = divmod(minutes, 60)
        print(process, 'duration:', "%02d:%02d" % (minutes, sec))
    return time_stamp


def merge_hits(merged):
    scaffold_ranges = {}
    for key, value in merged.items():
        scaffold_ranges[key] = []
        for inner_key, inner_value in value.items():
            for item in inner_value:
                scaffold_ranges[key].append(item)
    return scaffold_ranges


def update(file):
    pass


def run_blastn_if_needed(fasta_file, sequences_for_blastn):
    with open(fasta_file, 'r') as read:
        fasta_lines = read.readlines()
    sequences = []
    for header, sequence in zip(fasta_lines[::2], fasta_lines[1::2]):
        if header.count('_') == 1:
            if header.split('_')[1][:-1] in sequences_for_blastn:
                sequences.append(header)
                sequences.append(sequence)
        elif re.findall('_(.*)_', header)[0] in sequences_for_blastn:
            sequences.append(header)
            sequences.append(sequence)
        else:
            sequences.append(header)
            sequences.append('\n')


def main(file):
    time_stamp = time.time()
    print(file)
    dataframe = get_dataframe(f'blast/blastn/{file}.txt')
    time_stamp = print_time_passed('get_dataframe()', time_stamp)
    dataframe = get_taxonomic_data(dataframe)
    time_stamp = print_time_passed('get_taxonomic_data()', time_stamp)
    dictionary = get_ranges_dictionary(dataframe)
    time_stamp = print_time_passed('get_ranges_dictionary()', time_stamp)
    merged = merge_overlaps(dictionary)  # Merging each species' hits
    time_stamp = print_time_passed('merge_overlaps()', time_stamp)
    # Merging all the hits of all the species into a list of lists:
    scaffold_ranges = merge_hits(merged)
    time_stamp = print_time_passed('merge_hits()', time_stamp)
    # Now merges can be performed if needed:
    merged_scaffold_ranges = merge_overlaps({'outer_dictionary': scaffold_ranges})
    time_stamp = print_time_passed('merge_overlaps()', time_stamp)
    smaller_contamination = mark_high_contamination(dataframe, merged_scaffold_ranges['outer_dictionary'])
    time_stamp = print_time_passed('mark_high_contamination()', time_stamp)
    delete_false_hits(smaller_contamination, f'fasta_files/{file}.fa')
    time_stamp = print_time_passed('delete_false_hits()', time_stamp)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file_name')
    args = parser.parse_args()
    main(args.file_name)
