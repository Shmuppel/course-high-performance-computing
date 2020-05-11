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
    parser.add_argument('fastq_file', help='path to the file containing FASTQ reads', type=str, nargs=1)
    # Add an argument for assigning the number of threads to use.
    parser.add_argument('-t', '--threads', help='number of threads to assign', type=int, default=2)

    args = parser.parse_args()

    return args


def main():
    args = parse_command_line()


if __name__ == "__main__":
    main()
