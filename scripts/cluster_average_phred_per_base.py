#!/usr/bin/env python3
"""
Dependencies:
    - numpy
    - sed (UNIX)
"""
import csv
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

    def __init__(self, process_no, phred_score_lines_subset, max_read_length, output_queue):
        """
        Worker function to be called for every process. Given a subset of Phred score lines to process,
        allocate a NumPy matrix and fill it with the strings' character ASCII values. Places resulting
        NumPy matrix in the shared queue when finished.
        """
        print("{process} - Preallocating NumPy matrix...".format(process=process_no))
        phred_score_matrix = self.__preallocate_numpy_matrix(phred_score_lines_subset, max_read_length)
        print("{process} - Filling matrix with numerical Phred scores...".format(process=process_no))
        phred_score_matrix = self.__phred_lines_to_matrix(phred_score_lines_subset, phred_score_matrix)
        # TODO: Cleanup average declaration
        # Add process average results to output queue.
        output_queue.put(self.__calculate_average_phred_per_base(phred_score_matrix))

    def __preallocate_numpy_matrix(self, phred_score_lines, max_read_length):
        """
        Preallocate a numpy matrix containing NaN values given a list of Phred score strings.
        Its dimensions are (number of phred strings in subset, overall longest read length).

        :param phred_score_lines: list of Phred score strings.
        :return phred_score_matrix: numpy matrix filled with NaN values for every possible base location.
        """
        rows = len(phred_score_lines)
        columns = max_read_length

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

    def __init__(self, file_path, no_processes, server_ip, port, authkey):
        """
        :param file_path:
        :param no_processes:
        :param server_ip:
        :param port:
        :param authkey:
        """
        manager = self.__initialize_client_manager(server_ip, port, authkey)
        job_queue = manager.get_job_queue()
        result_queue = manager.get_result_queue()

        self.__run_client(job_queue, result_queue, file_path, no_processes)

    def get_phred_score_lines(self, file_path, start, stop):
        """
        :param file_path:
        :param start:
        :param stop:
        :return:
        """
        print("> Retrieving Phred score lines from FastQ file...")
        # Create a sed process that parses every 4th line (0~4p).
        get_4th_lines = subprocess.Popen("sed -n '{start}~4p;{stop}q' {file_path}".format(start=start,
                                                                                          stop=stop,
                                                                                          file_path=file_path),
                                         shell=True,
                                         stdout=subprocess.PIPE,
                                         universal_newlines=True)

        # Subprocess will pipe the output, lines are stripped of whitespace / newline characters.
        phred_score_lines = [line.strip() for line in get_4th_lines.stdout.readlines()]

        return phred_score_lines

    def __initialize_client_manager(self, server_ip, port, authkey):
        """
        :param server_ip:
        :param port:
        :param authkey:
        :return:
        """

        class ClientManager(BaseManager):
            pass

        ClientManager.register('get_job_queue')
        ClientManager.register('get_result_queue')

        manager = ClientManager(address=(server_ip, port), authkey=authkey)
        manager.connect()

        print('Client connected to %s:%s' % (server_ip, port))

        return manager

    def __run_client(self, job_queue, result_queue, file_path, no_processes):
        """
        :param job_queue:
        :param result_queue:
        :param file_path:
        :param no_processes:
        :return:
        """
        while True:
            try:
                job = job_queue.get_nowait()

                if job == "KILL":
                    job_queue.put("KILL")
                    break
                else:
                    function = job['function']
                    chunk_indices = job['chunk_indices']

                    phred_score_lines = self.get_phred_score_lines(file_path, chunk_indices[0], chunk_indices[1])

                    results = self.__start_processes(phred_score_lines, no_processes)
                    result_queue.put(results)

            except queue.Empty:
                print("sleep")
                time.sleep(1)

    def __start_processes(self, phred_score_lines, no_processes):
        """
        :param job_queue:
        :param result_queue:
        :param no_processes:
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

            process = mp.Process(target=AveragePhredCalculator,
                                 args=(process_no,
                                       phred_score_lines_subset,
                                       max_read_length,
                                       output_queue))

            processes.append(process)
            process.start()

        # Retrieve output from processes and average it.
        output = np.vstack([output_queue.get() for p in processes])
        output = np.sum(output, axis=0) / no_processes

        # Wait for processes to complete.
        for p in processes:
            p.join()

        return output


class Server:

    def __init__(self, file_path, function, no_clients, port, authkey):
        """
        :param file_path:
        :param function:
        :param port:
        :param no_clients:
        """
        self.file_path = file_path
        self.function = function
        self.port = port
        self.no_clients = no_clients

        self.chunks = self.get_chunk_indices(no_clients)
        self.manager = self.__initialize_server_manager(port, authkey)

    def get_file_length(self):
        """
        :return:
        """
        output = subprocess.run(['wc', '-l', self.file_path], capture_output=True, text=True)
        # get first element from wc -l output (line count of file)
        output = int(output.stdout.split()[0])

        return output

    def get_chunk_indices(self, no_clients):
        """
        :param no_clients:
        :return:
        """
        file_length = self.get_file_length()
        chunks = []

        chunk_size = math.ceil(file_length / no_clients)
        chunk_size += chunk_size % 4

        for client in range(no_clients):
            chunks.append((chunk_size * client,
                           chunk_size * (client + 1) - 1))

        return chunks

    def __initialize_server_manager(self, port, authkey):
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
        shared_job_queue = self.manager.get_job_queue()
        shared_result_queue = self.manager.get_result_queue()

        for chunk_indices in self.chunks:
            shared_job_queue.put({'function': self.function, 'chunk_indices': chunk_indices})

        time.sleep(5)

        results = []
        while True:
            try:
                result = shared_result_queue.get_nowait()
                results.append(result)

                if len(results) == len(self.chunks):
                    break

            except queue.Empty:
                time.sleep(1)
                continue

        # TODO - Clean this up
        concatenated_results = np.full(shape=(len(results), max([len(result) for result in results])),
                                       fill_value=np.NaN)
        for i, result in enumerate(results):
            concatenated_results[i, 0:len(result)] = result

        output = np.sum(concatenated_results, axis=0) / self.no_clients

        shared_job_queue.put("KILL")
        time.sleep(5)

        self.manager.shutdown()

        return output


def output_results_to_csv(output_file_path, average_phred_per_base):
    """
    Outputs the average Phred score per base to a given CSV file path.

    :param output_file_path: output file path to generate CSV file at.
    :param average_phred_per_base: numpy array containing the average phred per base location.
    """
    print("> Saving results to: {output}".format(output=output_file_path))
    try:
        with open(output_file_path, 'w') as output_file:
            writer = csv.writer(output_file, delimiter=',')
            writer.writerow(['base index', 'average phred score'])

            # Write average per base index, round average to 3 decimals.
            writer.writerows([index + 1, round(average, 3)]
                             for index, average in enumerate(average_phred_per_base))

    except EnvironmentError:
        print("Could not create output file, please check if the given output file path is valid.")
        sys.exit(1)


def initiate_clients(client_list, args):
    """
    :param client_list:
    :param server_ip:
    :return:
    """
    script_path = os.path.realpath(__file__)

    for client in client_list:
        subprocess.Popen("ssh {client} " \
                         "{script_path} -f {input_file} " \
                         "-n {processes} -p {port} --hosts {server_ip} -c"
                         .format(client=client,
                                 script_path=script_path,
                                 input_file=args.input_fastq[0],
                                 processes=args.processes,
                                 port=args.port,
                                 server_ip=args.hosts[0]),
                         shell=True)


def parse_command_line():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Average Phred score per base calculator')
    parser.add_argument('-f', '--input_fastq', help='input FASTQ file path', type=str, required=True, nargs=1)
    parser.add_argument('-o', '--output', help='output CSV file path', type=str)
    parser.add_argument('-n', '--processes', help='number of processes to assign per client', type=int, default=4)
    parser.add_argument('-p', '--port', help='port to use for server-host communication', type=int, required=True)
    parser.add_argument('--hosts', help='server and client ips respectively', type=str, required=True, nargs='+')

    execution_type = parser.add_mutually_exclusive_group(required=True)
    execution_type.add_argument('-c', '--client', help='run script as client', action='store_false')
    execution_type.add_argument('-s', '--server', help='run script as server', action='store_true')

    args = parser.parse_args()

    return args


def main():
    """
    Function called at script start up.

    :return: 0 when the program finishes correctly.
    """
    args = parse_command_line()

    # get absolute filepath in case a relative path is given.
    input_file_path = os.path.abspath(args.input_fastq[0])
    no_processes = args.processes

    port = args.port
    hosts = args.hosts
    authkey = b'authentication'
    execution_type = args.server

    # server & client arguments set to true if server, false if client.
    if execution_type:
        server = Server(input_file_path, AveragePhredCalculator, len(hosts) - 1, port, authkey)

        # launch remote clients
        if len(hosts) > 1:
            client_list = hosts[1:]
            initiate_clients(client_list, args)

        results = server.run_server()

        # output results to file if -o argument is passed.
        if args.output is not None:
            output_file_path = os.path.abspath(args.output)
            output_results_to_csv(output_file_path, results)
        else:
            print(np.round(results, decimals=2))

    else:
        server_ip = hosts[0]
        client = Client(input_file_path, no_processes, server_ip, port, authkey)

    return 0


if __name__ == "__main__":
    sys.exit(main())
