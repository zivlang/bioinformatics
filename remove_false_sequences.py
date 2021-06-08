#!/usr/bin/env python
import os
import re

directory = 'fasta_files/'


def merge_overlaps(dictionary):
    assembly_occupancy = {}
    for scaffold, hits_ranges in dictionary.items():
        assembly_occupancy.update({scaffold: {}})
        for species, hits in hits_ranges.items():
            assembly_occupancy[scaffold].update({species: []})
            hits = sorted(hits)
            # Removing hits that are within ranges of other hits
            filtered_hits = []
            hit_to_keep = []
            for index in range(len(hits)):
                if not hit_to_keep:
                    if len(hits) == 1:
                        filtered_hits = hits
                        continue
                    else:
                        hit_to_keep = hits[index]
                    continue
                # Check if hit_to_keep engulfs hits[index]:
                if hits[index][0] >= hit_to_keep[0] and hits[index][1] <= hit_to_keep[1]:
                    continue
                elif hits[index][1] > hit_to_keep[1] > hits[index][0]:  # check for partial overlap
                    hit_to_keep = [hit_to_keep[0], hits[index][1]]  # increasing the hit range accordingly
                else:  # in case of a non overlapping range in hits[index]
                    filtered_hits.append(hit_to_keep)
                    hit_to_keep = hits[index]
            if not filtered_hits:
                filtered_hits.append(hit_to_keep)
            elif hit_to_keep:
                filtered_hits.append(hit_to_keep)
            assembly_occupancy[scaffold][species] = filtered_hits
    return assembly_occupancy


def get_ranges_dictionary(dataframe):

    # Storing hits' positions lists in a dictionary
    # Creating contigs_dictionary: {contig name_1: {species_1: [a list of hits ranges], species_2...}, contig_name_2}
    contigs_dictionary = {}
    for index, row in dataframe.iterrows():
        species = ' '.join(row['stitle'].split()[:2])
        hit_range = [row['qstart'], row['qend']]
        contig = row['qseqid']
        if contig not in contigs_dictionary.keys():
            contigs_dictionary.update({contig: {}})
            contigs_dictionary[contig][species] = [hit_range]
        elif species not in contigs_dictionary[contig].keys():
            contigs_dictionary[contig].update({species: [hit_range]})
        # Appending only hit ranges that are not already in the dictionary
        elif hit_range not in contigs_dictionary[contig][species]:
            contigs_dictionary[contig][species].append(hit_range)
    return contigs_dictionary


def delete_false_hits(hits, fasta_file):
    """
    Deletes hits of species that are neither Myxozoa nor Cnidaria.
    :param fasta_file
    :param hits: ranges on contigs to delete from the contigs
    """
    with open(fasta_file, 'r') as read:
        fasta_lines = read.readlines()
    if len(hits) > 1:
        local_dictionary = hits
    else:
        local_dictionary = {key: value for key, value in hits.items()}
    sequences = []
    cleaned = False
    for header, sequence in zip(fasta_lines[::2], fasta_lines[1::2]):
        # In case the current header from the fasta file in not in the local_dictionary.keys(), write the header
        # and the sequence to the file, since blastn might not be accurate, and thus coding Thelohanellus sequences
        # are sometimes not found by blastn while in blastx they are.
        if header[1:-1] not in local_dictionary.keys():
            sequences.append(header)
            sequences.append(sequence)
            continue
        for scaffold, ranges in local_dictionary.items():
            if scaffold == header[1:-1] and ranges != 'contaminated':
                if len(ranges) > 1:
                    # If 'scaffold' in 'local dictionary' matches 'header' from the fasta file, the false
                    # sequences are assigned from 'sequence' by slicing it with the proper indexes in 'hit'.
                    hits_for_removal = [sequence[hit[0] - 1: hit[1]] for hit in ranges if len(hit)]
                else:
                    hits_for_removal = [sequence[ranges[0][0] - 1: ranges[0][1]]]
                # The false sequences in 'hits_for_removal' are joint by '|' so that each of them can be found as a
                # pattern by python's re.
                string_for_re = '|'.join([f'{hit}' for hit in hits_for_removal])
                removed_sequences = re.split(string_for_re, sequence.strip())
                kept_sequences = [sequence for sequence in removed_sequences if len(sequence) >= 100 if sequence]
                # Writing 'kept_sequences' to the fasta file:
                if len(kept_sequences) > 1:
                    for index, subsequence in enumerate(kept_sequences):
                        sequences.append(f'{header.strip()}_{index}\n')
                        sequences.append(f'{subsequence}\n')
                    cleaned = True
                elif len(kept_sequences) == 1:
                    sequences.append(header)
                    sequences.append(f'{kept_sequences[0]}\n')
                    cleaned = True
                # No need to continue iterating over the rest of the items in local_dictionary, since the matching
                # 'header' and 'scaffold' were already found.
                break
        # Removing the original header and sequence from which false sequences were excised and correct ones were
        # kept (it makes the loop run much faster):
        if cleaned:
            fasta_lines.remove(fasta_lines[0])
            fasta_lines.remove(fasta_lines[0])
            cleaned = False
    with open(f'{fasta_file[:-3]}_filtered.fa', 'w') as write:
        write.writelines(sequences)


def main():
    # for file in os.listdir('blast'):
    #     contigs_dictionary = get_ranges_dictionary('blast/' + file)
    #     merged = merge_overlaps(contigs_dictionary)
    #     # Merging all the hits of all the species into a list of lists:
    #     scaffold_ranges = {}
    #     for key, value in merged.items():
    #         scaffold_ranges[key] = []
    #         for inner_key, inner_value in value.items():
    #             scaffold_ranges[key].append(inner_value[0])
    #     # Now merges can be performed if needed:
    #     merged_scaffold_ranges = merge_overlaps({'outer_dictionary': scaffold_ranges})
    delete_false_hits()


if __name__ == '__main__':
    main()
