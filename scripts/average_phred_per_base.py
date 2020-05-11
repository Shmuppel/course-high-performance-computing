#!/usr/bin/python3

import argparse


def parse_command_line():
    """
    Parses the command line. argparse is used to declare the argument types the script
    requires to run.

    --- Returns ---
    args (namespace object) : The given arguments by the user parsed into a single
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


def main():
    args = parse_command_line()


if __name__ == "__main__":
    main()
