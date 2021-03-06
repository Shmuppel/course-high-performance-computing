#!/usr/bin/env python3
"""
Description: Script that calculates the average Phred score per base location of a FastQ file.
Usage: average_phred_per_base.py <FastQ file path> [-n <number of processes to use>] [-o CSV output file path]
Author: Niels van der Vegt (#365285)
Dependencies:
    - numpy
    - sed (UNIX)
"""

import os
import sys
import csv
import math
import argparse
import subprocess
import numpy as np
import multiprocessing as mp


class FileHandler:
    def __init__(self, file_path):
        """
        Object containing methods for handling FastQ files.

        :param file_path: file path to the FastQ file to be used.
        """
        self.file_path = file_path

    def get_phred_score_lines(self):
        """
        Retrieves all lines containing Phred scores from a FastQ file and stores them in a list.
        Has sed as a dependency to parse every fourth line.

        :return phred_scores: a list containing strings of Phred scores per read.
        """
        print("> Retrieving Phred score lines from FastQ file...")
        # Create a sed process that parses every 4th line (0~4p).
        get_4th_lines = subprocess.Popen("sed -n '0~4p' " + self.file_path,
                                         shell=True,
                                         stdout=subprocess.PIPE,
                                         universal_newlines=True
                                         )

        # Subprocess will pipe the output, lines are stripped of whitespace / newline characters.
        phred_score_lines = [line.strip() for line in get_4th_lines.stdout.readlines()]

        return phred_score_lines

    def output_results_to_csv(self, average_phred_per_base):
        """
        Outputs the average Phred score per base to a given CSV file path.

        :param output_file_path: output file path to generate CSV file at.
        :param average_phred_per_base: numpy array containing the average phred per base location.
        """
        print("> Saving results to: {output}".format(output=self.file_path))
        try:
            with open(self.file_path, "w") as output_file:
                writer = csv.writer(output_file, delimiter=",")
                writer.writerow(["base index", "average phred score"])

                # Write average per base index, round average to 3 decimals.
                writer.writerows([index + 1, round(average, 3)] for index, average in enumerate(average_phred_per_base))

        except EnvironmentError:
            print("Could not create output file, please check if the given output file path is valid.")
            sys.exit(1)


class AveragePhredWorker:
    def __init__(self, output_queue, phred_score_lines_subset, max_read_length):
        """
        """
        phred_score_matrix = self.__preallocate_numpy_matrix(phred_score_lines_subset, max_read_length)
        phred_score_matrix = self.__phred_lines_to_matrix(phred_score_lines_subset, phred_score_matrix)
        # Add process results to output queue.
        output_queue.put(phred_score_matrix)

    def __preallocate_numpy_matrix(self, phred_score_lines, max_read_length):
        """
        Preallocate a numpy matrix containing NaN values given a list of Phred score strings.
        Its dimensions are (number of phred strings in subset, overall longest read length).

        :param phred_score_lines: list of Phred score strings.
        :param max_read_length: the longest read length found across all processes.
        :return phred_score_matrix: numpy matrix filled with NaN values for every possible base location.
        """
        rows = len(phred_score_lines)
        columns = max_read_length

        # Using data-type 8-bit unsigned integer (u1) to save memory (Phred scores can't be negative nor > 104).
        phred_score_matrix = np.full(shape=(rows, columns), dtype=np.dtype("u1"), fill_value=np.nan)

        return phred_score_matrix

    def __phred_lines_to_matrix(self, phred_score_lines, phred_score_matrix):
        """
        Iterates over a list of Phred score strings, calculating the ASCII value of its characters
        and placing those values in a given Phred score matrix. Indices of the matrix are (read, base index),
        where reads shorter than the overall longest read will have trailing NaNs.

        :param phred_score_lines: list of Phred score strings.
        :param phred_score_matrix: numpy matrix filled with NaN values.
        :return phred_score_matrix: numpy matrix filled with integers (where applicable) or NaN values.
        """
        for i, phred_line in enumerate(phred_score_lines):
            phred_score_matrix[i, 0: len(phred_line)] = [ord(char) for char in phred_line]

        return phred_score_matrix


class ProcessManager:
    def __init__(self):
        """
        Composite object that collects a list of processes and can start them on command.
        Will return output when started.
        """
        self.processes = []
        self.no_processes = 0
        self.output_queue = mp.Queue()

    def add(self, worker_function, *args):
        """
        Adds a process to the composite object.
        """
        process = mp.Process(target=worker_function, args=(self.output_queue, *args))
        self.processes.append(process)
        self.no_processes += 1

    def run(self):
        """
        Iterates over all processes stored in the composite object and starts them, waits for them
        to finish and collects their results.
        :return:
        """
        for process in self.processes:
            process.start()

        # Retrieve output from processes.
        output = [self.output_queue.get() for process in self.processes]

        # Wait for processes to complete.
        for process in self.processes:
            process.join()

        return output


def calculate_average_phred(phred_score_lines, no_processes):
    """
    :return:
    """
    phred_score_lines_count = len(phred_score_lines)
    max_read_length = len(max(phred_score_lines, key=len))
    chunk_size = int(math.ceil(phred_score_lines_count / no_processes))

    # Make and start processes for calculating average phred.
    client = ProcessManager()
    for p in range(no_processes):
        phred_score_lines_subset = phred_score_lines[chunk_size * p: chunk_size * (p + 1)]
        client.add(AveragePhredWorker, phred_score_lines_subset, max_read_length)
    client_output = client.run()

    # Combine and average results.
    print("> Combining process results into matrix...")
    phred_score_matrix = np.vstack(client_output)
    average_phred_per_base = np.nanmean(phred_score_matrix, axis=0)

    return average_phred_per_base


def parse_command_line():
    """
    Parses the command line. argparse is used to declare the argument types the script requires to run.

    :return args: The arguments passed by the user parsed into a single namespace object.
    """
    parser = argparse.ArgumentParser(description="Average Phred score per base calculator")
    parser.add_argument("input_fastq", help="input FASTQ file path", type=str, nargs=1)
    parser.add_argument("-n", "--processes", help="number of processes to assign", type=int, default=4)
    parser.add_argument("-o", "--output", help="output CSV file path", type=str)

    args = parser.parse_args()

    return args


def main():
    """
    Function called at script start up.

    :return: 0 when the program finishes correctly.
    """
    args = parse_command_line()
    # Retrieve absolute path in case data is stored in working directory.
    input_file_path = os.path.abspath(args.input_fastq[0])
    no_processes = args.processes if args.processes else 1

    # FastQ file handling.
    fastq_file_handler = FileHandler(input_file_path)
    phred_score_lines = fastq_file_handler.get_phred_score_lines()
    # Calculate average phred per base.
    average_phred_per_base = calculate_average_phred(phred_score_lines, no_processes)

    # If output path is given, save results to file.
    if args.output is not None:
        output_file_path = os.path.abspath(args.output)
        output_file_handler = FileHandler(output_file_path)
        output_file_handler.output_results_to_csv(average_phred_per_base)
    else:
        print(np.round(average_phred_per_base, decimals=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())
