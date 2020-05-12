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
    # Create a sed process that parses every 4th line (0~4p).
    get_4th_lines = subprocess.Popen("sed -n '0~4p' " + file_path,
                                     shell=True,
                                     stdout=subprocess.PIPE,
                                     universal_newlines=True)
    # Subprocess will pipe the output, lines are stripped of whitespace / newline characters.
    phred_score_lines = [line.strip() for line in get_4th_lines.stdout.readlines()]

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
    # Using data-type 8-bit unsigned integer to save memory (Phred scores can't be negative nor > 104).
    phred_score_matrix = np.full(shape=(rows, columns), dtype=np.dtype('u1'), fill_value=np.nan)

    return phred_score_matrix


def phred_lines_to_matrix(phred_score_lines, phred_score_matrix):
    """
    Iterates over a list of Phred score strings, calculating the ASCII value of its characters
    and placing those values in the Phred score matrix.

    Indices are based on (read, base index), where reads shorter than the longest read will
    have trailing NaNs.

    :param phred_score_lines: list of Phred score strings.
    :param phred_score_matrix: numpy matrix filled with NaN values for every possible base location.
    :return phred_score_matrix: numpy matrix filled with integers (where applicable) or NaN values
                                for every possible base location.
    """
    for i, phred_line in enumerate(phred_score_lines):
        for j, char in enumerate(phred_line):
            phred_score_matrix[i, j] = ord(char)

    return phred_score_matrix


def main():
    args = parse_command_line()

    # Retrieve absolute path in case data is stored in working directory.
    file_path = os.path.abspath(args.input_fastq[0])

    # Preparing data.
    print("> Retrieving Phred score lines from FastQ file.")
    phred_score_lines = get_phred_score_lines(file_path)
    print("> Preallocating NumPy matrix.")
    phred_score_matrix = preallocate_numpy_matrix(phred_score_lines)

    print("> Filling NumPy matrix with Phred score integer values.")
    phred_score_matrix = phred_lines_to_matrix(phred_score_lines, phred_score_matrix)


if __name__ == "__main__":
    main()
