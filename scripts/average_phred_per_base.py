#!/usr/bin/env python3

import os
import argparse
import subprocess


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


def get_phred_scores(file_path):
    """
    Retrieves all lines containing Phred scores from a FastQ file and stores them in a list.
    Has sed as a dependency to parse every fourth line.

    :param file_path: path pointing to a FastQ file.
    :return phred_scores: a list containing strings of Phred scores per read.
    """
    # Retrieve absolute path in case data is stored in working directory.
    absolute_file_path = os.path.abspath(file_path)

    # Create a sed process that parses every 4th line, subprocess will pipe the output.
    get_4th_lines = subprocess.Popen("sed -n '0~4p' " + absolute_file_path, shell=True,
                                     stdout=subprocess.PIPE, universal_newlines=True)
    phred_scores = get_4th_lines.stdout.readlines()

    return phred_scores


def main():
    args = parse_command_line()
    phred_scores = get_phred_scores(args.input_fastq[0])


if __name__ == "__main__":
    main()
