#!/usr/bin/env python3
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
    def __init__(self, file_path, chunk_indices, output_queue):
        """
        Worker function to be called for every process. Given a subset of Phred score lines to process,
        allocate a NumPy matrix and fill it with the strings' character ASCII values. Places resulting
        NumPy matrix in the shared queue when finished.
        """
        phred_score_lines = self.__get_phred_score_lines(file_path, chunk_indices[0], chunk_indices[1])
        max_read_length = len(max(phred_score_lines, key=len))
        # Convert Phred scores to P values for every Phred line and put the results in a matrix.
        matrix = self.__preallocate_numpy_matrix(phred_score_lines, max_read_length)
        p_value_matrix = self.__phred_lines_to_matrix(phred_score_lines, matrix)
        # Calculate average P value per base.
        average_p_value_per_base = np.nanmean(p_value_matrix, axis=0)
        output_queue.put(average_p_value_per_base)

    @staticmethod
    @lru_cache(maxsize=None)
    def calculate_p_value_from_phred(char):
        """
        TODO
        """
        return 10 ** (-(ord(char) - 33) / 10)

    @staticmethod
    def calculate_phred_from_p_value(p_value_array):
        """
        :param p_value_array:
        :return:
        """
        # Convert P values back to Phred values (Q Score).
        return -10 * np.log10(p_value_array)

    def __get_phred_score_lines(self, file_path, start, stop):
        """
        TODO
        """
        # Create a sed process that parses every 4th line (0~4p).
        get_4th_lines = subprocess.Popen("sed -n '{start}~4p;{stop}q' {file_path} "
                                         .format(start=start, stop=stop, file_path=file_path),
                                         stdout=subprocess.PIPE,
                                         shell=True,
                                         universal_newlines=True)
        # Subprocess will pipe the output, lines are stripped of whitespace / newline characters.
        phred_score_lines = [line.strip() for line in get_4th_lines.stdout.readlines()]
        return phred_score_lines

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
        phred_score_matrix = np.full(shape=(rows, columns), dtype=np.dtype("float32"), fill_value=np.nan)

        return phred_score_matrix

    def __phred_lines_to_matrix(self, phred_score_lines, matrix):
        """
        Iterates over a list of Phred score strings, calculating the ASCII value of its characters
        and placing those values in a given Phred score matrix. Indices of the matrix are (read, base index),
        where reads shorter than the overall longest read will have trailing NaNs.

        :param phred_score_lines: list of Phred score strings.
        :param phred_score_matrix: numpy matrix filled with NaN values.
        :return phred_score_matrix: numpy matrix filled with integers (where applicable) or NaN values.
        """
        for i, phred_line in enumerate(phred_score_lines):
            matrix[i, 0: len(phred_line)] = [AveragePhredCalculator.calculate_p_value_from_phred(char)
                                             for char in phred_line]

        return matrix


class Client:
    def __init__(self, file_path, no_processes, server_ip, client, port, authkey):
        """
        TODO
        """
        self.client = client
        self.file_path = file_path
        self.no_processes = no_processes
        self.manager = self.__make_client_manager(server_ip, port, authkey)
        # Run client if the client is connected and ready to start operations.
        if self.manager is not None:
            self.__run_client()

    def __make_client_manager(self, server_ip, port, authkey):
        """
        TODO
        """
        # Initialize ClientManager class and methods.
        class ClientManager(BaseManager):
            pass
        ClientManager.register("get_job_queue")
        ClientManager.register("get_result_queue")

        print("{client}> Attempting client connection on {server_ip}:{port}..."
              .format(client=self.client, server_ip=server_ip, port=port))
        try:
            manager = ClientManager(address=(server_ip, port), authkey=authkey)
            manager.connect()
        except ConnectionRefusedError:
            # Failed to connect, a timer of the client's socket manager will attempt a reconnect.
            print("{client}> Failed to connect to {server}:{port}, retrying..."
                  .format(client=self.client, server=server_ip, port=port))
        else:
            print("{client}> Client connected to {server}:{port}"
                  .format(client=self.client, server=server_ip, port=port))
            return manager

    def __send_heartbeat(self, result_queue):
        """
        TODO
        """
        result_queue.put("HEARTBEAT")
        time.sleep(5)
        self.__send_heartbeat(result_queue)

    def __run_client(self):
        """
        TODO
        """
        job_queue = self.manager.get_job_queue()
        result_queue = self.manager.get_result_queue()

        # Start sending heartbeat signals to the server.
        # Daemon such that the timer process ends when the client is finished.
        timer = mp.Process(target=self.__send_heartbeat, args=(result_queue,), daemon=True)
        timer.start()

        while True:
            try:
                job = job_queue.get_nowait()

                if job == "KILL":
                    break
                else:
                    # Get assigned job and start processes.
                    chunk_indices = job["chunk"]
                    results = self.__start_processes(chunk_indices)
                    result_queue.put(results)

            except queue.Empty:
                time.sleep(1)

    def __start_processes(self, chunk_indices):
        """
        TODO
        """
        processes = []
        output_queue = mp.Queue()

        # Initialize data division related variables.
        phred_line_count = chunk_indices[1] - chunk_indices[0]
        chunks = calculate_chunk_indices(phred_line_count, self.no_processes)

        for p in range(self.no_processes):
            # Declare and start process.
            process = mp.Process(target=AveragePhredCalculator,
                                 args=(self.file_path, chunks[p], output_queue))
            processes.append(process)
            process.start()

        # Wait for processes to complete.
        for p in processes:
            p.join()

        # Retrieve output from processes and average it over the amount of processes.
        concat_p_value_arrays = concatenate_uneven_numpy_arrays([output_queue.get() for p in processes])
        p_average = np.sum(concat_p_value_arrays, axis=0) / self.no_processes

        return p_average


class Server:
    def __init__(self, server_ip, port, input_file_path, clients, no_processes, chunks, authkey):
        """
        Object that sets up a number of sockets, launches clients specified by the user, and facilitates
        communication with these clients. Handles the assignment of jobs to client processes and parses
        the results returned by them.
        """
        # Get chunk indices to be used for jobs; the amount of chunks default to the amount of clients
        # unless specified by the user.
        file_length = get_file_length(input_file_path)
        self.chunk_indices = calculate_chunk_indices(file_length, chunks) if chunks \
            else calculate_chunk_indices(file_length, len(clients))

        # Client arguments
        self.input_file_path = input_file_path
        self.no_processes = no_processes

        # Server arguments.
        self.server_ip = server_ip
        self.port = port
        self.authkey = authkey
        self.sockets = self.__make_sockets(clients, authkey)

    def __make_sockets(self, clients, authkey):
        """
        Declares a socket manager for every client specified by the user, then launches the client which
        will attempt to connect to the socket.
        :param clients: clients that will need a connection to the server.
        :return sockets: a list of dictionaries wherein all necessary socket information is stored.
        """
        sockets = []
        dynamic_port = self.port
        for client in clients:
            # Assign a unique port for every client, based on the port specified by the user.
            dynamic_port += 1
            # Define a socket manager and launch the corresponding client.
            server_manager = self.__make_server_manager(dynamic_port, authkey)
            self.__make_client(client, dynamic_port, self.no_processes)
            sockets.append({"client": client,
                            "port": dynamic_port,
                            "server_manager": server_manager,
                            "job": None,
                            "timer": Timer(10, self.__make_client, args=[client, dynamic_port, self.no_processes])})

        return sockets

    def __make_server_manager(self, port, authkey):
        """
        Declares a multiprocessing ServerManager running on the server's IP, on a specified port.
        :param port: which port the ServerManager should use for communications.
        :return manager: the started ServerManager object.
        """
        job_queue = queue.Queue()
        result_queue = queue.Queue()

        # Initialize ServerManager class and methods.
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
        Launches a SSH subprocess to a client, where it will launch this script as a client and attempt
        to connect to the server.
        :param client: client IP address to navigate to through SSH.
        :param port: the port the client should use for server communications.
        :param no_processes: the number of processes that should be launched on the client during jobs.
        """
        # Get path where this script is stored.
        script_path = os.path.realpath(__file__)
        # Connect to the client through SSH and run this script as a client.
        subprocess.Popen(list(map(str, ["ssh", client, script_path,
                                        "-f", self.input_file_path,
                                        "-n", no_processes,
                                        "-p", port,
                                        "--hosts", self.server_ip + ' ' + client,
                                        "-c"])))

    def run_server(self):
        """
        TODO
        """
        results = []

        # Define queue containing unassigned jobs shared among all clients.
        unassigned_job_queue = queue.Queue()
        for i, chunk in enumerate(self.chunk_indices):
            unassigned_job_queue.put({"job_id": i + 1, "chunk": chunk})

        while len(results) != len(self.chunk_indices):
            for socket in self.sockets:
                client = socket["client"]
                port = socket["port"]
                server_manager = socket["server_manager"]
                job_queue = server_manager.get_job_queue()
                result_queue = server_manager.get_result_queue()

                # Check if client is connected and in need of a job.
                if socket["job"] is None and socket["timer"].is_alive():
                    # If there are still unassigned jobs, give one to this client.
                    if not unassigned_job_queue.empty():
                        assigned_job = unassigned_job_queue.get()
                        # Set which job this socket is currently working on.
                        socket["job"] = assigned_job["job_id"]
                        job_queue.put(assigned_job)
                        print("> Assigned job number {job} to {client}".format(job=assigned_job["job_id"],
                                                                               client=client))

                # Should a client disconnect whilst working on a job, put that job back in unassigned job queue.
                elif socket["job"] is not None and not socket["timer"].is_alive():
                    print("> {client} has disconnected whilst operating on job {job}, attempting to reconnect"
                          .format(client=socket["client"],
                                  job=socket["job"]))

                    unassigned_job_queue.put({"job_id": socket["job"], "chunk": self.chunk_indices[socket["job"] - 1]})
                    socket["job"] = None

                # Handle heartbeats and calculated results placed in results queue.
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
                        socket["job"] = None

        # Clean up processes.
        print("> Cleaning up processes")
        for socket in self.sockets:
            socket["timer"].cancel()
            # Terminate manager.
            job_queue = socket["server_manager"].get_job_queue()
            job_queue.put("KILL")

        time.sleep(1)

        # Terminate sockets.
        for socket in self.sockets:
            socket["server_manager"].shutdown()

        # Process results.
        print("> Processing results")
        concatenated_results = concatenate_uneven_numpy_arrays(results)
        mean_results = np.nanmean(concatenated_results, axis=0)
        # Convert P values back to Phred values (Q Score).
        phred_average = AveragePhredCalculator.calculate_phred_from_p_value(mean_results)

        return phred_average


# --- HELPER FUNCTIONS ---


def calculate_chunk_indices(number, no_chunks):
    """
    Calculates the indices to use when 'binning' a certain number, such that all bins, or chunks, are of
    equal size. In this implementation chunks are always divisible by four, to account for the no. lines
    per read in a FastQ file.
    :param number: the number to be partitioned into equal chunks.
    :param no_chunks: the amount of chunks to partition the number into.
    :return chunk_indices: the indices of chunks.
    """
    chunk_indices = []
    chunk_size = math.ceil(number / no_chunks)
    # Make sure chunks are always divisible by four.
    chunk_size += chunk_size % 4

    for i in range(no_chunks):
        # -1 to prevent overlapping indices between chunks.
        chunk_indices.append([chunk_size * i, chunk_size * (i + 1) - 1])
    # On the end index of the last chunk: disregard the modulo (and its additive effect on chunk size).
    chunk_indices[-1][1] = number

    return chunk_indices


def concatenate_uneven_numpy_arrays(arrays):
    """
    Concatenated two Numpy arrays of uneven length by taking the length of the longest array and forming
    a 2D array with dimensions (no. arrays, length of longest array).
    :param arrays: list of Numpy arrays to be concatenated.

    """
    # Pre-allocate 2D array.
    concat_array = np.full(shape=(len(arrays), max([len(array) for array in arrays])), fill_value=np.NaN)
    # Fill 2D array with arrays to be concatenated.
    for i, array in enumerate(arrays):
        concat_array[i, 0:len(array)] = array

    return concat_array


def get_file_length(file_path):
    """
    Returns the amount of lines present in a given file using a unix wc subprocess.
    :param file_path: path to file of interest.
    :return file_length: an integer denoting the line count of the file.
    """
    output = subprocess.run(["wc", "-l", file_path], capture_output=True, text=True)
    # Get first element from wc -l output (line count of file).
    file_length = int(output.stdout.split()[0])

    return file_length


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
    Parses all argument passed by the user on the command line to a single namespace object.
    :return args: namespace object containing arguments passed by user.
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


# --- MAIN ---


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

    # Execution type argument group set to false if client.
    else:
        Client(input_file_path, no_processes, server_ip, clients[0], port, authkey)

    return 0


if __name__ == "__main__":
    sys.exit(main())
