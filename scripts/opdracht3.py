#!/usr/bin/env python3
"""
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
from multiprocessing.managers import BaseManager
from functools import lru_cache
from threading import Timer


class FastQFileHandler:
    def __init__(self, file_path):
        """
        Object containing methods for handling FastQ files.

        :param file_path: file path to the FastQ file to be used.
        """
        self.file_path = file_path

    def get_file_length(self):
        """
        Returns file length by running a wc subprocess and retrieving its output.
        :return:
        """
        output = subprocess.run(['wc', '-l', self.file_path], capture_output=True, text=True)
        # get first element from wc -l output (line count of file)
        output = int(output.stdout.split()[0])

        return output

    def get_phred_score_lines(self, file_path, start=None, stop=None):
        """
        Retrieves the Phred score quality lines from a specified section of a given FastQ file.
        A section is denoted by a start and stop line count. Uses unix sed to stream Phred score
        lines into memory.
        :param file_path: file path to a FastQ file.
        :param start: line count from which Phred lines should be retrieved.
        :param stop: line count to which Phred lines should be retrieved.
        :return phred_score_lines: a list of stripped Phred score lines.
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


class FileWriter:
    def __init__(self, file_path):
        """
        Object that deals with writing content to a given file path.
        :param file_path: file path to write to.
        """
        self.file_path = file_path

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


class MultiProcessingManager:
    def __init__(self, ip, port, authkey):
        """
        Base class that creates a multiprocessing manager with a job and result queue for communications.
        Is assigned an IP, port and authentication key for starting or connecting to.
        """
        self.ip = ip
        self.port = port
        self.authkey = authkey
        self.manager = self.initialize_manager()

    def initialize_manager(self):
        job_queue = queue.Queue()
        result_queue = queue.Queue()

        class Manager(BaseManager):
            pass

        Manager.register('get_job_queue', callable=lambda: job_queue)
        Manager.register('get_result_queue', callable=lambda: result_queue)

        manager = Manager(address=(self.ip, self.port), authkey=self.authkey)

        return manager


class ServerProxy(MultiProcessingManager):
    def __init__(self, port, authkey, ip=''):
        super().__init__(ip, port, authkey)

    def start(self):
        print("> Running socket on port: {port}".format(port=self.port))
        self.manager.start()


class ClientProxy(MultiProcessingManager):
    def __init__(self, port, authkey, ip, client):
        super().__init__(ip, port, authkey)
        self.client = client

    def start(self):
        """
        Declares a multiprocessing manager which will attempt connection on a given server_IP and
        port. Returns this manager if a connection was established successfully.
        :return manager: multiprocessing manager with an active connection to the server.
        """
        print("{client}> Attempting client connection on {server_ip}:{port}..."
              .format(client=self.client, server_ip=self.ip, port=self.port))
        try:
            self.manager.connect()

        except ConnectionRefusedError:
            # Failed to connect, the timer of the client's server manager will attempt a reconnect.
            print("{client}> Failed to connect to {server}:{port}, retrying..."
                  .format(client=self.client, server=self.ip, port=self.port))
        else:
            print("{client}> Client connected to {server}:{port}"
                  .format(client=self.client, server=self.ip, port=self.port))


class AveragePhredCalculator:
    def __init__(self, output_queue, file_path, chunk_indices):
        """
        Worker function to be called for every process. Given a subset of Phred score lines to process,
        allocate a NumPy 2D array and fill it with the strings' character respective P-values. Places
        resulting NumPy 2D array in the shared queue when finished.
        """
        fastq_file_handler = FastQFileHandler(file_path)
        phred_score_lines = fastq_file_handler.get_phred_score_lines(file_path, chunk_indices[0], chunk_indices[1])
        max_read_length = len(max(phred_score_lines, key=len))
        # Convert Phred scores to P values for every Phred line and put the results in a matrix.
        matrix = self.__preallocate_numpy_array(phred_score_lines, max_read_length)
        p_value_matrix = self.__phred_lines_to_array(phred_score_lines, matrix)
        # Calculate average P value per base.
        average_p_value_per_base = np.nanmean(p_value_matrix, axis=0)
        output_queue.put(average_p_value_per_base)

    @staticmethod
    @lru_cache(maxsize=None)
    def calculate_p_value_from_phred(char):
        """
        Calculates the P value from a Phred score character. The specified character is first converted
        to its ASCII value and then converted to its respective P-value. Makes use of the lru_cache
        decorator (with maxsize=None, disregarding the 'least-recently-used' operations) allowing this
        function to be memoized.
        :param char: Phred score character to be converted to P-value.
        :return: P-value of a given Phred score character.
        """
        return 10 ** (-(ord(char) - 33) / 10)

    @staticmethod
    def calculate_phred_from_p_value(p_value_array):
        """
        Converts a Numpy array of P-values to one of Phred scores.
        :param p_value_array: a Numpy array of P-values.
        :return: a Numpy array of Phred scores.
        """
        return -10 * np.log10(p_value_array)

    def __preallocate_numpy_array(self, phred_score_lines, max_read_length):
        """
        Preallocate a numpy array containing NaN values given a list of Phred score strings.
        Its dimensions are (no. of strings in list, longest string  in list).

        :param phred_score_lines: list of Phred score strings.
        :return phred_score_matrix: Numpy 2D array filled with NaN values.
        """
        rows = len(phred_score_lines)
        columns = max_read_length
        # Using data-type 32-bit float for storage of P-values.
        phred_score_matrix = np.full(shape=(rows, columns), dtype=np.dtype("float32"), fill_value=np.nan)

        return phred_score_matrix

    def __phred_lines_to_array(self, phred_score_lines, matrix):
        """
        Iterates over a list of Phred score strings, calculating the P-value value of its characters
        and placing those values in a given 2D Numpy array. Indices of the array are (read, base index),
        where reads shorter than the overall longest read will have trailing NaNs.

        :param phred_score_lines: list of Phred score strings.
        :param matrix: A 2D Numpy array filled with NaN values.
        :return matrix: Numpy 2D array filled with P-values (where applicable) or NaN values.
        """
        for i, phred_line in enumerate(phred_score_lines):
            matrix[i, 0: len(phred_line)] = [AveragePhredCalculator.calculate_p_value_from_phred(char)
                                             for char in phred_line]

        return matrix


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


class Client:
    def __init__(self, client_name, server_ip, port, authkey, no_processes):
        """
        Object used when running this script as a Client. Attempts communication with a specified
        socket and assigns local processes to perform average Phred score calculations, returning
        their results to the server.
        """
        self.client_name = client_name
        self.no_processes = no_processes

        self.proxy = ClientProxy(port, authkey, server_ip, client_name)
        self.manager = self.proxy.manager
        self.proxy.start()

        # Run client if the client is connected and ready to start operations.
        if self.manager is not None:
            self.run()

    def __send_heartbeat(self, result_queue):
        """
        Adds a heartbeat signal (string) to a given queue, then after a short delay performs recursion.
        Must be called with a background (daemon) process as it is designed to not halt.
        :param result_queue: queue to put heartbeat signal into.
        """
        result_queue.put("HEARTBEAT")
        time.sleep(5)
        self.__send_heartbeat(result_queue)

    def run(self):
        """
        Starts operations on the client: retrieving jobs from a server, starting processes to perform
        this job, and returning the result to the server. Whilst active a heartbeat signal will be send
        from a background process running on the client.
        """
        job_queue = self.manager.get_job_queue()
        result_queue = self.manager.get_result_queue()

        # Start sending heartbeat signals to the server.
        # Daemon such that the heartbeat signal ends when the client is finished or crashes.
        timer = mp.Process(target=self.__send_heartbeat, args=(result_queue,), daemon=True)
        timer.start()

        while True:
            try:
                job = job_queue.get_nowait()

                if job == "KILL":
                    # End operations on this client.
                    break
                else:
                    # Retrieve assigned job and start processes.
                    results = self.perform_job(job)
                    result_queue.put(results)

            except queue.Empty:
                time.sleep(1)


class AveragePhredCalculatorClient(Client):
    def __init__(self, client_name, server_ip, port, authkey, no_processes):
        """
        Client implementation for assigning average phred calculation jobs.
        """
        super().__init__(client_name, server_ip, port, authkey, no_processes)

    def perform_job(self, job):
        """
        Assigns local threads to calculating average phred per base of the given FastQ line indices.
        Averages the output of these threads and returns this array.
        """
        chunk_indices = job["chunk"]
        file_path = job["file_path"]

        phred_line_count = chunk_indices[1] - chunk_indices[0]
        chunks = calculate_chunk_indices(phred_line_count, self.no_processes)

        process_manager = ProcessManager()
        for p in range(self.no_processes):
            process_manager.add(AveragePhredCalculator, file_path, chunks[p])

        process_output = process_manager.run()

        # Retrieve output from processes and average it.
        concat_p_value_arrays = concatenate_uneven_numpy_arrays(process_output)
        p_value_average = np.nanmean(concat_p_value_arrays, axis=0)

        return p_value_average


class Server:
    """
    Implement Server class as a singleton instance.
    """
    instance = None

    class __Server:
        def __init__(self, server_ip, port, clients, no_processes, authkey, input_file_path, chunks):
            """
            Object that sets up a number of sockets, launches clients specified by the user, and facilitates
            communication with these clients. Handles the assignment of jobs to client processes and parses
            the results returned by them.
            """
            # Get chunk indices to be used for jobs; the amount of chunks default to the amount of clients
            # unless specified by the user.
            fastq_file_handler = FastQFileHandler(input_file_path)
            file_length = fastq_file_handler.get_file_length()

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
                server_proxy = ServerProxy(dynamic_port, authkey)
                server_manager = server_proxy.manager
                server_proxy.start()

                self.__make_client(client, dynamic_port, self.no_processes)
                sockets.append({"client": client,
                                "port": dynamic_port,
                                "server_manager": server_manager,
                                "job": None,
                                "timer": Timer(10, self.__make_client, args=[client, dynamic_port, self.no_processes])})

            return sockets

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
            Runs the server to perform a number of job, executing the following actions for every socket,
            as long as not all jobs are finished:
                - Check if client is connected and in need of a job.
                - If not, check if client has disconnected whilst working on a job.
                - Check the results queue for heartbeats or results and handle them., if it's not empty.
            :return phred_score_average: Numpy array containing average phred score for bases in a FastQ file.
            """
            results = []

            # Define queue to contain jobs that are still unassigned.
            unassigned_job_queue = queue.Queue()
            for i, chunk in enumerate(self.chunk_indices):
                unassigned_job_queue.put({"job_id": i + 1,
                                          "chunk": chunk,
                                          "file_path": self.input_file_path})

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
                    elif not socket["timer"].is_alive() and socket["job"] is not None:
                        print("> {client} has disconnected whilst operating on job {job}, attempting to reconnect"
                              .format(client=socket["client"],
                                      job=socket["job"]))

                        unassigned_job_queue.put({"job_id": socket["job"], "chunk": self.chunk_indices[socket["job"] - 1]})
                        socket["job"] = None

                    # Handle heartbeats and calculated results placed in results queue.
                    while not result_queue.empty():
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
            phred_score_average = AveragePhredCalculator.calculate_phred_from_p_value(mean_results)

            return phred_score_average

    def __init__(self, *args):
        if not Server.instance:
            Server.instance = Server.__Server(*args)
        else:
            print(">Error: Server already running.")

    def __getattr__(self, name):
        return getattr(self.instance, name)


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
        server = Server(server_ip, port, clients, no_processes, authkey, input_file_path, chunks)
        results = server.run_server()

        # If output path is given, save results to file.
        if args.output is not None:
            output_file_path = os.path.abspath(args.output)
            output_file_handler = FileWriter(output_file_path)
            output_file_handler.output_results_to_csv(results)
        else:
            print(np.round(results, decimals=2))

    # Execution type argument group set to false if client.
    else:
        AveragePhredCalculatorClient(clients[0], server_ip, port, authkey, no_processes)

    return 0


if __name__ == "__main__":
    sys.exit(main())
