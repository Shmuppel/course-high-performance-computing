#!/usr/bin/env python3
"""
Dependencies:
    - numpy
    - sed (UNIX)
"""

import os
import queue
import sys
import time
from multiprocessing.managers import BaseManager
import math
import argparse
import subprocess
import numpy as np
import multiprocessing as mp


class AveragePhredCalculator:

    def __init__(self, process_no, phred_score_lines_subset, output_queue):
        """
        Worker function to be called for every process. Given a subset of Phred score lines to process,
        allocate a NumPy matrix and fill it with the strings' character ASCII values. Places resulting
        NumPy matrix in the shared queue when finished.
        """
        print("{process} - Preallocating NumPy matrix...".format(process=process_no))
        phred_score_matrix = self.__preallocate_numpy_matrix(phred_score_lines_subset)
        print("{process} - Filling matrix with numerical Phred scores...".format(process=process_no))
        phred_score_matrix = self.__phred_lines_to_matrix(phred_score_lines_subset, phred_score_matrix)
        # TODO - Average results, send those results back

        # Add process results to output queue.
        output_queue.put(phred_score_matrix)

    def __preallocate_numpy_matrix(self, phred_score_lines):
        """
        Preallocate a numpy matrix containing NaN values given a list of Phred score strings.
        Its dimensions are (number of phred strings in subset, overall longest read length).

        :param phred_score_lines: list of Phred score strings.
        :return phred_score_matrix: numpy matrix filled with NaN values for every possible base location.
        """
        rows = len(phred_score_lines)
        columns = len(max(phred_score_lines, key=len))

        # Using data-type 8-bit unsigned integer (u1) to save memory (Phred scores can't be negative nor > 104).
        phred_score_matrix = np.full(shape=(rows, columns),
                                     dtype=np.dtype('u1'),
                                     fill_value=np.nan)

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
            phred_score_matrix[i, 0:len(phred_line)] = [ord(char) - 33 for char in phred_line]

        return phred_score_matrix

    def __calculate_average_phred_per_base(self, phred_score_matrix):
        """
        Calculates the average value along each Phred score matrix column.

        :return average_phred_per_base: numpy array containing average Phred score per base location.
        """
        print("> Calculating average Phred score per base location...")
        # matrix.nanmean to support the trailing NaN values for bases shorter than the longest read.
        average_phred_per_base = np.nanmean(phred_score_matrix, axis=0)

        return average_phred_per_base


class Client:

    def __init__(self, file_path, no_workers):
        """
        Object containing methods using multiprocessing to calculate the average phred per base location,
        given a list of Phred score strings and the no_processes to use. Initializes a Phred score matrix
        with integer values representing every Phred score character in every Phred score line.

        :param no_workers: number of process to assign calculations.
        """
        phred_score_lines = self.get_phred_score_lines(file_path)
        phred_score_matrices = self.run_workers(phred_score_lines, no_workers)
        print(phred_score_matrices)

    def get_phred_score_lines(self, file_path):
        """
        Retrieves all lines containing Phred scores from a FastQ file and stores them in a list.
        Has sed as a dependency to parse every fourth line.

        :return phred_scores: a list containing strings of Phred scores per read.
        """
        print("> Retrieving Phred score lines from FastQ file...")
        # Create a sed process that parses every 4th line (0~4p).
        get_4th_lines = subprocess.Popen("sed -n '0~4p' " + file_path,
                                         shell=True,
                                         stdout=subprocess.PIPE,
                                         universal_newlines=True)

        # Subprocess will pipe the output, lines are stripped of whitespace / newline characters.
        phred_score_lines = [line.strip() for line in get_4th_lines.stdout.readlines()]

        return phred_score_lines

    def run_workers(self, phred_score_lines, no_processes):
        """
        Starts a specified number of processes, with each process receiving a subset of Phred score
        string to perform conversion to ASCII values upon.

        :param phred_score_lines: list of Phred score strings.
        :param no_processes: number of process to assign calculations.
        :return output: list of NumPy matrices containing each process' results.
        """
        # Initialize multiprocessing variables.
        processes = []
        output_queue = mp.Queue()

        # Initialize data separation related variables.
        phred_score_lines_count = len(phred_score_lines)
        chunk_size = int(math.ceil(phred_score_lines_count / no_processes))

        print("> Starting {n} processes...".format(n=no_processes))
        for p in range(no_processes):
            process_no = p + 1
            phred_score_lines_subset = phred_score_lines[chunk_size * p: chunk_size * (p + 1)]

            process = mp.Process(target=AveragePhredCalculator,
                                 args=(process_no,
                                       phred_score_lines_subset,
                                       output_queue))

            processes.append(process)
            process.start()

        # Retrieve output from processes.
        output = [output_queue.get() for p in processes]

        # Wait for processes to complete.
        for p in processes:
            p.join()

        return output


class Server:

    def __init__(self, function, file_path, port):
        """

        """
        self.file_path = file_path
        self.function = function

        self.manager = self.initialize_server_manager(port, b'authentication')
        self.shared_job_queue = self.manager.get_job_queue()
        self.shared_result_queue = self.manager.get_result_queue()

    def initialize_server_manager(self, port, authkey):
        """
        :param port:
        :param authkey:
        :return:
        """
        job_queue = queue.Queue()
        result_queue = queue.Queue()

        class ServerManager(BaseManager):
            pass

        ServerManager.register('get_job_queue', callable=lambda: job_queue)
        ServerManager.register('get_result_queue', callable=lambda: result_queue)

        manager = ServerManager(address=('', port), authkey=authkey)
        manager.start()

        return manager

    def run_server(self):
        """
        :return:
        """
        for item in self.data:
            self.shared_job_queue.put({'function': self.function, 'data': item})

        time.sleep(2)

        results = []
        while True:
            try:
                result = self.shared_result_queue.get_nowait()
                results.append(result)

                if len(results) == len(self.data):
                    break

            except queue.Empty:
                time.sleep(1)
                continue

        self.shared_job_queue.put("KILL")
        time.sleep(5)
        self.manager.shutdown()
        print(results)


def parse_command_line():
    """
    Parses the command line. argparse is used to declare the argument types the script requires to run.

    :return args: The arguments passed by the user parsed into a single namespace object.
    """
    parser = argparse.ArgumentParser(description='Average Phred score per base calculator')
    parser.add_argument('-f', '--input_fastq', help='input FASTQ file path', type=str, required=True)
    parser.add_argument('-o', '--output', help='output CSV file path', type=str)
    parser.add_argument('-n', '--processes', help='number of processes to assign per client', type=int, default=4)
    parser.add_argument('-p', '--port', help='port to use for server-host communication', type=int, required=True)
    parser.add_argument('--hosts', help='hosts, first host should refer to the server', type=str, required=True)

    execution_type = parser.add_mutually_exclusive_group(required=True)
    execution_type.add_argument('-c', '--client', help='run script as client', action='store_false')
    execution_type.add_argument('-s', '--server',  help='run script as server', action='store_true')

    args = parser.parse_args()

    return args


def main():
    """
    Function called at script start up.

    :return: 0 when the program finishes correctly.
    """
    args = parse_command_line()
    input_file_path = os.path.abspath(args.input_fastq[0])
    no_processes = args.processes
    port = args.port
    hosts = args.hosts
    execution_type = args.server

    if execution_type:
        #execute server
        server = mp.Process(target=Server, args=(AveragePhredCalculator, input_file_path))
        #server.start()
        #server.join()
        pass
    else:
        #execute client
        server = mp.Process(target=Client, args=(no_processes))
        #client.start()
        #client.join()
        pass


if __name__ == "__main__":
    sys.exit(main())
