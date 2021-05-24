# Comparing every parameter in the N50 output files for each pair: before and after the filter.
# This requires 8 lists: one for every parameter. Each list will include the parameter value of the pairs:
# 16 pairs - one for every genome. Each pair will be in a dictionary: The key will be the name of the file,
# and the value will be the parameter's value.

import os
import numpy as np
import matplotlib.pyplot as plt

# I don't think it will create images when running on the server. Run locally after modifying the following paths
after_path = r"C:\Users\97254\PycharmProjects\workingProject\after_filter"
before_path = r"C:\Users\97254\PycharmProjects\workingProject\before_filter"

files_names = []

N25_length, N50_length, N75_length, Longest_contig, Shortest_contig,\
    Number_of_bases, Number_of_contigs, Average_number_of_bases_per_contig = [[], []], [[], []], [[], []], [[], []], \
                                                                             [[], []], [[], []], [[], []], [[], []]

lists_list = [N25_length, N50_length, N75_length, Longest_contig, Shortest_contig, Number_of_bases, Number_of_contigs,
              Average_number_of_bases_per_contig]


class File:
    def __init__(self, file_name, before_content, after_content):
        self.file_name = file_name
        self.before_content = before_content
        self.after_content = after_content

    def get_files_names(self):
        files_names.append((self.file_name.split('.'))[0])

    @staticmethod
    def get_values(self, content, before_filter):
        value = ''
        for line, _list in zip(content, lists_list):
            if '=' in line:
                value = int(line.strip('\n').split(' = ')[1])
            elif ':' in line:
                value = line.strip('\n').split(': ')[1]
                if '.' not in value:
                    value = int(value)
                else:
                    value = float(value)
            if value != '':
                if before_filter:
                    _list[0].append(value)
                else:
                    _list[1].append(value)


def create_histograms():

    for secondary_list in lists_list:

        parameter_name = [name for name in globals() if globals()[name] is secondary_list][0]
        before_list = secondary_list[0]
        after_list = secondary_list[1]

        labels_locations = np.arange(len(files_names))  # the labels positions
        width = 0.20  # bars width
        fig, ax = plt.subplots()

        before_histogram = ax.bar(labels_locations - width / 2, before_list, width, label='Before filter')
        after_histogram = ax.bar(labels_locations + width / 2, after_list, width, label='After filter')

        ax.set_ylabel(str(parameter_name[0]))
        ax.set_title(str(parameter_name[0]) + '\n\n')
        ax.set_xticks(labels_locations)
        ax.set_xticklabels(files_names)
        plt.xticks(rotation=45, ha='right')
        ax.legend()
        # fig.tight_layout()

        def autolabel(values):
            """Attach a text label above each bar in *rects*, displaying its height."""
            for value in values:
                height = value.get_height()
                ax.annotate('{}'.format(height),
                            xy=(value.get_x() + value.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', rotation=45,
                            style='oblique',
                            fontsize=10,
                            va='bottom')

        autolabel(before_histogram)
        autolabel(after_histogram)
        plt.show()
        # plt.savefig(parameter_name)

        plt.clf()


if __name__ == '__main__':
    for path in os.listdir(before_path):
        file_name = os.path.basename(path)
        after_file = after_path + '\\' + file_name
        with open(before_path + '\\' + path, 'r') as before_file, open(after_file, 'r') as after_file:
            before_values = before_file.readlines()
            after_values = after_file.readlines()
            file = File(file_name, before_values, after_values)
            file.get_values(file, file.before_content[3:], before_filter=True)
            file.get_values(file, file.after_content[3:], before_filter=False)
            file.get_files_names()
    create_histograms()            
