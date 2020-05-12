#!/usr/bin/env python3

import os
import argparse
import subprocess
import numpy as np


def parse_command_line():
    """
    Parses the command line. argparse is used to declare the argument types the script
    requires to run.

    :return args (namespace object): The arguments passed by the user parsed into a single
                                     namespace object.
    """
    parser = argparse.ArgumentParser(description='Average Phred score per base calculator')
    # Add a positional argument for the FASTQ file path.
    parser.add_argument('input_fastq', help='input FASTQ file path', type=str, nargs=1)
    # Add an optional argument for assigning the number of threads to use.
    parser.add_argument('-n', '--processes', help='number of processes to assign', type=int, default=2)
    # Add an optional argument for outputting the results to a CSV file.
    parser.add_argument('-o', '--output', help='output CSV file path', type=str)

    args = parser.parse_args()

    return args


def get_phred_score_lines(file_path):
    """
    Retrieves all lines containing Phred scores from a FastQ file and stores them in a list.
    Has sed as a dependency to parse every fourth line.

    :param file_path: path pointing to a FastQ file.
    :return phred_scores: a list containing strings of Phred scores per read.
    """
    # Create a sed process that parses every 4th line (0~4p), and strips newline (s/\\n$//).
    get_4th_lines = subprocess.Popen("sed -n '0~4p' " + file_path,
                                     shell=True,
                                     stdout=subprocess.PIPE,
                                     universal_newlines=True)
    # Subprocess will then pipe the output.
    phred_score_lines = get_4th_lines.stdout.readlines()

    return phred_score_lines


def preallocate_numpy_matrix(phred_score_lines):
    """
    Preallocate a numpy matrix containing NaN values given a list of Phred score strings.
    Its dimensions are (number of reads/strings, longest read length).

    :param phred_score_lines: list of Phred score strings.
    :return phred_score_matrix: numpy matrix filled with NaN values for every possible base location.
    """
    rows = len(phred_score_lines)
    columns = len(max(phred_score_lines, key=len))

    phred_score_matrix = np.full(shape=(rows, columns), fill_value=np.nan)

    return phred_score_matrix


def phred_lines_to_chars(phred_score_lines):
    """
    Given a list of Phred score strings, unpack each string into a list of its individual characters.

    :param phred_score_lines: list of Phred score strings.
    :return phred_score_chars: list of lists containing Phred score characters.
    """
    phred_score_chars = [[*line.strip()] for line in phred_score_lines]

    return phred_score_chars


def main():
    args = parse_command_line()

    # Retrieve absolute path in case data is stored in working directory.
    file_path = os.path.abspath(args.input_fastq[0])

    phred_score_lines = get_phred_score_lines(file_path)

    phred_score_matrix = preallocate_numpy_matrix(phred_score_lines)
    phred_score_chars = phred_lines_to_chars(phred_score_lines)


if __name__ == "__main__":
    main()
