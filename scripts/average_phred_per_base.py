#!/usr/bin/env python3

import csv
import os
import sys
import math
import argparse
import subprocess
import numpy as np
import multiprocessing as mp


class FastQFileHandler:
    """
    """

    def __init__(self, file_path):
        print("> Retrieving Phred score lines from FastQ file...")
        self.file_path = file_path

    def get_phred_score_lines(self):
        """
        Retrieves all lines containing Phred scores from a FastQ file and stores them in a list.
        Has sed as a dependency to parse every fourth line.

        :return phred_scores: a list containing strings of Phred scores per read.
        """
        # Create a sed process that parses every 4th line (0~4p).
        get_4th_lines = subprocess.Popen("sed -n '0~4p' " + self.file_path,
                                         shell=True,
                                         stdout=subprocess.PIPE,
                                         universal_newlines=True)
        # Subprocess will pipe the output, lines are stripped of whitespace / newline characters.
        phred_score_lines = [line.strip() for line in get_4th_lines.stdout.readlines()]

        return phred_score_lines


class AveragePhredCalculator:
    """
    """

    def __init__(self, phred_score_lines, no_processes):
        # Start multiprocessing with passed amount of processes.
        phred_score_matrices = self.multiprocess(phred_score_lines, no_processes)
        # Combine results of processes.
        print("> Combining process results into matrix...")
        self.phred_score_matrix = np.vstack(phred_score_matrices)

    def multiprocess(self, phred_score_lines, no_processes):
        """
        :param no_processes:
        :param phred_score_lines:
        :return:
        """
        # Initialize multiprocessing variables.
        processes = []
        output_queue = mp.Queue()

        # Initialize data separation related variables.
        phred_score_lines_count = len(phred_score_lines)
        max_read_length = len(max(phred_score_lines, key=len))
        chunk_size = int(math.ceil(phred_score_lines_count / no_processes))

        print("> Starting {n} processes...".format(n=no_processes))
        for p in range(no_processes):
            process_no = p + 1
            phred_score_lines_subset = phred_score_lines[chunk_size * p: chunk_size * (p + 1)]

            process = mp.Process(target=self.__process_function,
                                 args=(process_no,
                                       phred_score_lines_subset,
                                       max_read_length,
                                       output_queue))

            processes.append(process)
            process.start()

        # Retrieve output from processes.
        output = [output_queue.get() for p in processes]

        # Wait for processes to complete.
        for p in processes:
            p.join()

        return output

    def __process_function(self, process_no, phred_score_lines_subset, max_read_length, output_queue):
        """
        :return:
        """
        print("{process} - Preallocating NumPy matrix...".format(process=process_no))
        phred_score_matrix = self.__preallocate_numpy_matrix(phred_score_lines_subset, max_read_length)
        print("{process} - Filling matrix with numerical Phred scores...".format(process=process_no))
        phred_score_matrix = self.__phred_lines_to_matrix(phred_score_lines_subset, phred_score_matrix)

        # Add process results to output queue.
        output_queue.put(phred_score_matrix)

        return 0

    def __preallocate_numpy_matrix(self, phred_score_lines, max_read_length):
        """
        Preallocate a numpy matrix containing NaN values given a list of Phred score strings.
        Its dimensions are (number of reads/strings, longest read length).

        :param phred_score_lines: list of Phred score strings.
        :param max_read_length:
        :return phred_score_matrix: numpy matrix filled with NaN values for every possible base location.
        """
        rows = len(phred_score_lines)
        columns = max_read_length

        # Using data-type 8-bit unsigned integer to save memory (Phred scores can't be negative nor > 104).
        phred_score_matrix = np.full(shape=(rows, columns), dtype=np.dtype('u1'), fill_value=np.nan)

        return phred_score_matrix

    def __phred_lines_to_matrix(self, phred_score_lines, phred_score_matrix):
        """
        Iterates over a list of Phred score strings, calculating the ASCII value of its characters
        and placing those values in the Phred score matrix.

        Indices of the matrix are (read, base index), where reads shorter than the longest read will
        have trailing NaNs.

        :param phred_score_lines: list of Phred score strings.
        :param phred_score_matrix: numpy matrix filled with NaN values for every possible base location.
        :return phred_score_matrix: numpy matrix filled with integers (where applicable) or NaN values
                                    for every possible base location.
        """
        for i, phred_line in enumerate(phred_score_lines):
            phred_score_matrix[i, 0:len(phred_line)] = [ord(char) for char in phred_line]

        return phred_score_matrix

    def calculate_average_phred_per_base(self):
        """
        Calculates the average value along each Phred score matrix column.

        :return average_phred_per_base: numpy array containing average Phred score per base location.
        """
        print("> Calculating average Phred score per base location...")
        # matrix.nanmean to support the trailing NaN values for bases shorter than the longest read.
        average_phred_per_base = np.round(np.nanmean(self.phred_score_matrix, axis=0), decimals=2)

        return average_phred_per_base


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
    parser.add_argument('-n', '--processes', help='number of processes to assign', type=int, default=4)
    # Add an optional argument for outputting the results to a CSV file.
    parser.add_argument('-o', '--output', help='output CSV file path', type=str)

    args = parser.parse_args()

    return args


def main():
    """
    :return:
    """
    args = parse_command_line()

    # Retrieve absolute path in case data is stored in working directory.
    input_file_path = os.path.abspath(args.input_fastq[0])
    output_file_path = os.path.abspath(args.output)
    no_processes = args.processes

    # FastQ file handling.
    fastq_file_handler = FastQFileHandler(input_file_path)
    phred_score_lines = fastq_file_handler.get_phred_score_lines()

    # Calculate Phred score average per base location using multiprocessing.
    average_phred_calculator = AveragePhredCalculator(phred_score_lines, no_processes)
    average_phred_per_base = average_phred_calculator.calculate_average_phred_per_base()

    # If output path is given, save results to file.
    if output_file_path is not None:
        print("> Saving results to: {output}".format(output=output_file_path))
        with open(output_file_path, 'w') as output_file:
            writer = csv.writer(output_file, delimiter=',')
            writer.writerow(['base index', 'average phred score'])
            writer.writerows([index + 1, average] for index, average in enumerate(average_phred_per_base))

    return 0


if __name__ == "__main__":
    sys.exit(main())
