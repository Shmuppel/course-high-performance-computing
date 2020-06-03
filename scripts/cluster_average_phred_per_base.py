#!/usr/bin/env python3
"""
TODO
"""
import csv
import os
import queue
import sys
import time
import math
import argparse
import subprocess
import numpy as np
import multiprocessing as mp
from threading import Timer
from functools import lru_cache
from multiprocessing.managers import BaseManager


class AveragePhredCalculator:
    def __init__(self, phred_score_lines_subset, max_read_length, output_queue):
        """
        Worker function to be called for every process. Given a subset of Phred score lines to process,
        allocate a NumPy matrix and fill it with the strings' character ASCII values. Places resulting
        NumPy matrix in the shared queue when finished.
        """
        phred_score_matrix = self.__preallocate_numpy_matrix(phred_score_lines_subset, max_read_length)
        phred_score_matrix = self.__phred_lines_to_matrix(phred_score_lines_subset, phred_score_matrix)
        output_queue.put(self.__calculate_average_phred_per_base(phred_score_matrix))

    @lru_cache(maxsize=None)
    def calculate_phred_score(self, char):
        """
        TODO
        """
        return ord(char) - 33

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
            phred_score_matrix[i, 0: len(phred_line)] = [self.calculate_phred_score(char) for char in phred_line]

        return phred_score_matrix

    def __calculate_average_phred_per_base(self, phred_score_matrix):
        """
        Calculates the average value along each Phred score matrix column.

        :return average_phred_per_base: numpy array containing average Phred score per base location.
        """
        # np.nanmean to support the trailing NaN values for bases shorter than the longest read.
        average_phred_per_base = np.nanmean(phred_score_matrix, axis=0)

        return average_phred_per_base


class Client:
    def __init__(self, file_path, no_processes, server_ip, port, authkey):
        """
        TODO
        """
        self.manager = self.__make_client_manager(server_ip, port, authkey)
        self.run_client(file_path, no_processes)

    def get_phred_score_lines(self, file_path, start, stop):
        """
        TODO
        """
        # Create a sed process that parses every 4th line (0~4p).
        get_4th_lines = subprocess.Popen("sed -n '{start}~4p;{stop}q' {file_path}"
                                         .format(start=start, stop=stop, file_path=file_path),
                                         shell=True,
                                         stdout=subprocess.PIPE,
                                         universal_newlines=True
                                         )
        # Subprocess will pipe the output, lines are stripped of whitespace / newline characters.
        phred_score_lines = [line.strip() for line in get_4th_lines.stdout.readlines()]

        return phred_score_lines

    def __make_client_manager(self, server_ip, port, authkey):
        """
        TODO
        """
        # Initialize ClientManager class and methods.
        class ClientManager(BaseManager):
            pass
        ClientManager.register("get_job_queue")
        ClientManager.register("get_result_queue")

        print("> Attempting client connection on {server_ip}:{port}".format(server_ip=server_ip, port=port))
        manager = ClientManager(address=(server_ip, port), authkey=authkey)
        manager.connect()
        print("> Client connected to {server}:{port}".format(server=server_ip, port=port))

        # Start sending heartbeat signals to the server.
        # Daemon such that the timer process ends when the client is finished.
        result_queue = manager.get_result_queue()
        timer = mp.Process(target=self.__send_heartbeat, args=(result_queue,), daemon=True)
        timer.start()

        return manager

    def __send_heartbeat(self, result_queue):
        """
        TODO
        """
        result_queue.put("HEARTBEAT")
        time.sleep(5)
        self.__send_heartbeat(result_queue)

    def run_client(self, file_path, no_processes):
        """
        TODO
        """
        job_queue = self.manager.get_job_queue()
        result_queue = self.manager.get_result_queue()

        while True:
            try:
                job = job_queue.get_nowait()

                if job == "KILL":
                    break
                else:
                    # Get assigned job and start processes.
                    chunk_indices = job["chunk"]
                    phred_score_lines = self.get_phred_score_lines(file_path, chunk_indices[0], chunk_indices[1])

                    results = self.__start_processes(phred_score_lines, no_processes)
                    result_queue.put(results)

            except queue.Empty:
                time.sleep(1)

    def __start_processes(self, phred_score_lines, no_processes):
        """
        TODO
        """
        processes = []
        output_queue = mp.Queue()

        # Initialize data division related variables.
        phred_score_lines_count = len(phred_score_lines)
        max_read_length = len(max(phred_score_lines, key=len))
        chunk_size = int(math.ceil(phred_score_lines_count / no_processes))

        for p in range(no_processes):
            # Get a chunk of data to process.
            phred_score_lines_subset = phred_score_lines[chunk_size * p: chunk_size * (p + 1)]
            # Declare and start process.
            process = mp.Process(target=AveragePhredCalculator,
                                 args=(phred_score_lines_subset, max_read_length, output_queue))
            processes.append(process)
            process.start()

        # Retrieve output from processes and average it.
        output = np.vstack([output_queue.get() for p in processes])
        average = np.sum(output, axis=0) / no_processes

        # Wait for processes to complete.
        for p in processes:
            p.join()

        return average


class Server:
    def __init__(self, server_ip, port, input_file_path, clients, no_processes, chunks, authkey):
        """
        TODO
        """
        # Server arguments.
        self.server_ip = server_ip
        self.port = port
        # Client arguments
        self.input_file_path = input_file_path
        self.no_processes = no_processes
        self.chunks = self.get_chunk_indices(chunks) if chunks else self.get_chunk_indices(len(clients))
        # Sockets
        self.sockets = self.__start_sockets(clients, port, authkey)

    def get_file_length(self):
        """
        TODO
        """
        output = subprocess.run(["wc", "-l", self.input_file_path], capture_output=True, text=True)
        # Get first element from wc -l output (line count of file).
        output = int(output.stdout.split()[0])

        return output

    def get_chunk_indices(self, no_chunks):
        """
        TODO
        """
        file_length = self.get_file_length()
        chunks = []
        chunk_size = math.ceil(file_length / no_chunks)
        chunk_size += chunk_size % 4

        for i in range(no_chunks):
            chunks.append((chunk_size * i, chunk_size * (i + 1) - 1))

        return chunks

    def __start_sockets(self, clients, port, authkey):
        """
        TODO
        """
        sockets = []
        for client in clients:
            port += 1
            # Define a socket manager and corresponding client.
            socket_manager = self.__make_socket_manager(port, authkey)
            self.__make_client(client, port, self.no_processes)
            sockets.append({"client": client,
                            "port": port,
                            "socket_manager": socket_manager,
                            "timer": Timer(10, self.__make_client, args=[client, port, self.no_processes]),
                            "needs_job": True})

        return sockets

    def __make_socket_manager(self, port, authkey):
        """
        TODO
        """
        # Initialize ServerManager class and methods.
        job_queue = queue.Queue()
        result_queue = queue.Queue()

        class ServerManager(BaseManager):
            pass
        ServerManager.register("get_job_queue", callable=lambda: job_queue)
        ServerManager.register("get_result_queue", callable=lambda: result_queue)

        manager = ServerManager(address=("", port), authkey=authkey)
        manager.start()
        print("> Running socket on {server_ip}:{port}".format(server_ip=self.server_ip, port=port))

        return manager

    def __make_client(self, client, port, no_processes):
        """
        TODO
        """
        # Get path where this script is stored.
        script_path = os.path.realpath(__file__)
        # Connect to the client through SSH and run this script as a client.
        subprocess.Popen(list(map(str, ["ssh", client, script_path,
                                        "-f", self.input_file_path,
                                        "-n", no_processes,
                                        "-p", port,
                                        "--hosts", self.server_ip,
                                        "-c"])))

    def concatenate_results(self, results):
        """
        TODO
        """
        # Pre-allocate combined results array.
        concat_results = np.full(shape=(len(results), max([len(result) for result in results])), fill_value=np.NaN)
        # Fill array with results.
        for i, result in enumerate(results):
            concat_results[i, 0:len(result)] = result

        return concat_results

    def run_server(self):
        """
        TODO
        """
        results = []
        job_counter = 0

        while len(results) != len(self.chunks):
            for socket in self.sockets:
                client = socket["client"]
                port = socket["port"]
                socket_manager = socket["socket_manager"]

                job_queue = socket_manager.get_job_queue()
                result_queue = socket_manager.get_result_queue()

                # Check if client is connected and in need of a job.
                if socket["needs_job"] and socket["timer"].is_alive():
                    # If there are still unassigned jobs, give one to this client.
                    if job_counter != len(self.chunks):
                        print("> Assigned job number {job} to {client}".format(job=job_counter + 1, client=client))
                        job_queue.put({"chunk": self.chunks[job_counter]})
                        socket["needs_job"] = False
                        job_counter += 1

                for r in range(result_queue.qsize()):
                    result = result_queue.get()
                    # Handle heartbeat to know that server is still connected.
                    if isinstance(result, str) and result == "HEARTBEAT":
                        socket["timer"].cancel()
                        # Declare and start a new thread timer.
                        socket["timer"] = Timer(10, self.__make_client, args=[client, port, self.no_processes])
                        socket["timer"].start()
                    else:
                        results.append(result)
                        socket["needs_job"] = True

        # Clean up processes.
        print("> Cleaning up processes")
        for socket in self.sockets:
            socket["timer"].cancel()
            # Terminate manager.
            job_queue = socket["socket_manager"].get_job_queue()
            job_queue.put("KILL")
            time.sleep(1)
            socket["socket_manager"].shutdown()

        # Process results.
        print("> Processing results")
        concatenated_results = self.concatenate_results(results)
        weighted_average = np.sum(concatenated_results, axis=0) / len(self.chunks)

        return weighted_average


def output_results_to_csv(output_file_path, average_phred_per_base):
    """
    Outputs the average Phred score per base to a given CSV file path.

    :param output_file_path: output file path to generate CSV file at.
    :param average_phred_per_base: numpy array containing the average phred per base location.
    """
    print("> Saving results to: {output}".format(output=output_file_path))
    try:
        with open(output_file_path, "w") as output_file:
            writer = csv.writer(output_file, delimiter=",")
            writer.writerow(["base index", "average phred score"])

            # Write average per base index, round average to 3 decimals.
            writer.writerows([index + 1, round(average, 3)] for index, average in enumerate(average_phred_per_base))

    except EnvironmentError:
        print("Could not create output file, please check if the given output file path is valid.")
        sys.exit(1)


def parse_command_line():
    """
    TODO
    """
    parser = argparse.ArgumentParser(description="Average Phred score per base calculator.")
    parser.add_argument("-f", "--input_fastq", help="input FASTQ file path", type=str, required=True)
    parser.add_argument("-o", "--output", help="output CSV file path", type=str)
    parser.add_argument("-n", "--processes", help="number of processes to assign per client", type=int, default=4)
    parser.add_argument("-p", "--port", help="port to use for server-host communication", type=int, required=True)
    parser.add_argument("--hosts", help="server and client IPs respectively", type=str, required=True, nargs="+")
    parser.add_argument("--chunks", help="amount of chunks to divide data into", type=int)

    # Arguments to run script as Client or Server.
    execution_type = parser.add_mutually_exclusive_group(required=True)
    execution_type.add_argument("-c", "--client", help="run script as a client", action="store_false")
    execution_type.add_argument("-s", "--server", help="run script as a server", action="store_true")

    args = parser.parse_args()

    return args


def main():
    """
    Function called at script start up.
    :return 0: when the program finishes correctly.
    """
    args = parse_command_line()

    # Get absolute filepath in case a relative path is given.
    input_file_path = os.path.abspath(args.input_fastq)
    port = args.port
    server_ip = args.hosts[0]
    clients = args.hosts[1:] if len(args.hosts) >= 2 else []
    no_processes = args.processes
    chunks = args.chunks
    execution_type = args.server
    authkey = b"authentication"

    # Execution type argument group set to true if server.
    if execution_type:
        server = Server(server_ip, port, input_file_path, clients, no_processes, chunks, authkey)
        results = server.run_server()

        # Output results to file if -o argument is passed.
        if args.output is not None:
            output_file_path = os.path.abspath(args.output)
            output_results_to_csv(output_file_path, results)
        else:
            print(results)

    # Execution type argument group set to false if server.
    else:
        Client(input_file_path, no_processes, server_ip, port, authkey)

    return 0


if __name__ == "__main__":
    sys.exit(main())
