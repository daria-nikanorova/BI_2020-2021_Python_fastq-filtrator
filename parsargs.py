from os.path import exists

# PART 1. ARGUMENTS PARSING FUNCTIONS


def parse_args(arguments, argument_dict):
    """
    Looks through all input arguments and calls functions to parse them by name
    :param arguments: a list with all input arguments (may be a direct output from sys.argv())
    :param argument_dict: a dictionary with all possible arguments and their values by default
    :return: a dictionary of parsed arguments and their values
    """
    if len(arguments) == 1 or '--help' in arguments:  # call help page
        find_help()
        exit()
    used_args = []
    file = find_fastq_file(arguments)  # finds file for filtering (a positional argument)
    argument_dict['--output_base_name'] = file.replace('.fastq', '').split('/')[-1]
    used_args += [file]
    ind = 1
    while ind + 1 != len(arguments):  # loop over all arguments except fastq file and the 1st arg (path to a script)
        arg = arguments[ind]
        next_value = arguments[ind + 1]
        if arg in used_args:  # check if argument was assigned before
            print(f'\nArgument {arg} is given more than once.\n'
                  f'Please try again')
            exit()
        if arg.startswith('--'):  # check if arg is an argument, not a value
            pass
        else:
            print(f'\nUnknown value used: {arg}\n'
                  f'Please try again')
            exit()
        if arg not in argument_dict.keys():  # check if it is a known argument for this script
            print(f'\nUnknown argument used: {arg} \n'
                  f'List of possible arguments can be checked with --help\n'
                  f'Please try again')
            exit()
        if arg == '--keep_filtered':  # this flag is too small for special function
            argument_dict[arg] = True
        elif next_value.startswith('--') or ind + 2 == len(arguments):  # check if there is a value after an argument
            print(f'\nNo value(s) for {arg} is/are defined. Please try again')
            exit()
        else:
            parsing_result = globals()[arg.replace('--', 'find_')](arguments, ind)  # call a special function for each optional arg
            argument_dict[arg] = parsing_result[0]
            ind = parsing_result[1]  # get an index of the last value for an argument
        used_args += [arg]
        ind += 1  # make 1 step to the next arg
    return argument_dict, file


def find_gc_bounds(args, ind):
    """
    Parses argument --gc_bounds, finds its values, is called by parse_args
    :param args: a list with all input arguments (may be a direct output from sys.argv())
    :param ind: an index of --gc_bounds in a list of arguments
    :return: a list with min and max bounds for GC filtering and
    an index of the last value of this argument (for error catching)
    """
    start, stop = ind + 1, ind + 2
    bounds = [0, 100]
    try:
        float(args[start])
    except ValueError:
        print('\nMinimum GC content threshold for filtering must be a positive float number. \n'
              'Please check --min_length parameters and try again')
        exit()
    if float(args[start]) < 0:
        print('\nMinimum GC content threshold must be higher or equal to 0%. \n'
              'Please check --gc_bounds parameters and try again')
        exit()
    bounds[0] = float(args[start])
    last_value_ind = start
    try:
        float(args[stop])
    except ValueError:
        print(f'\nOnly minimum threshold ({bounds[0]}%) will be used for GC filtering.')
        return bounds, last_value_ind
    if float(args[stop]) > 100:
        print('\nMaximum GC content threshold must be lower or equal to 100%. \n'
              'Please check --gc_bounds parameters and try again')
        exit()
    bounds[1] = float(args[stop])
    last_value_ind = stop
    if bounds[1] <= bounds[0]:
        print('\nMaximum GC content threshold must be higher than minimum GC content threshold. \n'
              'Please check --gc_bounds parameters and try again')
        exit()
    return bounds, last_value_ind


def find_min_length(args, ind):
    """
    Parses argument --min_length, finds its value, is called by parse_args
    :param args: a list with all input arguments (may be a direct output from sys.argv())
    :param ind: an index of --min_length in a list of arguments
    :return: a list with min length for filtering and
    an index of the last value of this argument (for error catching)
    """
    try:
        int(args[ind + 1])
    except ValueError:
        print('\nLength for filtering must be a positive integer. \n'
              'Please check --min_length parameters and try again')
        exit()
    last_value_ind = ind + 1
    return int(args[ind + 1]), last_value_ind


def find_output_base_name(args, ind):
    """
    Parses argument --output_base_name, finds its value, is called by parse_args
    :param args: a list with all input arguments (may be a direct output from sys.argv())
    :param ind: an index of --output_base_name in a list of arguments
    :return: a list with min and max bounds for GC filtering and
    an index of the last value of this argument (for error catching)
    """
    name = args[ind + 1]
    last_value_ind = ind + 1
    return name, last_value_ind


def find_fastq_file(args):
    """
    Parses the last argument fastq file with reads for filtration, is called by parse_args
    :param args: a list with all input arguments (may be a direct output from sys.argv())
    :return: True
    """
    if args[-1].endswith('.fastq'):
        if exists(args[-1]):
            file = args[-1]
            return file
        print(f'{args[-1]} does not exist.')
        exit()
    print('\nNo fastq file for filtering is provided.\n'
          'Please try again and put fastq file as THE LAST argument')
    exit()


def find_help():
    """
    :return: help page
    """
    print('\n       filter_fastq.py filtrates reads in .fastq format by GC content and length.\n\n'
          'SYNOPSIS:\n\n'
          '         python filter_fastq.py --some_arguments some_values my_reads.fastq\n\n'
          'DESCRIPTION:\n\n'
          '         The .fastq file is provided as positional argument and must be put in the end of command\n'
          '         All other arguments are optional:\n\n'
          '--min_length             Specifies the length for read to pass filtration.\n'
          '                         Value must be > 0\n\n'
          '--gc_bounds              Specifies minimum (value1) and maximum (value2)\n'
          '                         thresholds of GC content to pass filtration.\n'
          '                         Values must be between 0% and 100% (these are default)\n'
          '                         Only minimum threshold may be used\n\n'
          '--keep_filtered          If this flag is used, non-filtered reads are written\n'
          '                         to output_base_name__failed.fastq\n\n'
          '--output_base_name       Specifies basename for output files\n'
          '                         Reads passed the filtration will be stored\n'
          '                         in output_base_name__passed.fastq\n'
          '                         Non-selected reads will be stored in\n'
          '                         output_base_name__failed.fastq (in case a flag\n'
          '                         --keep_filtered is used)\n\n'
          '--help                   Help page\n\n')
    return


# PART 2. READS FILTERING FUNCTIONS


def calculate_gc(read_seq):
    """
    Calculates GC-content as a percentage of G or C bases in DNA sequence: Count(G + C)/Count(A + T + G + C) * 100%
    :param read_seq: a read sequence from the FASTQ file
    :return: GC-content of a read sequence, %
    """

    return (read_seq.count('G') + read_seq.count('C')) * 100 / len(read_seq)


def pass_read_check(read_seq, min_length, min_gc_bound, max_gc_bound):
    """
    Check whether a read sequence passes the filtration by its length and GC-content
    :param read_seq: a read sequence from the FASTQ file
    :param min_length: minimum length for a read to pass the filtration
    :param min_gc_bound: minimum GC-content value of a read to pass the filtration.
    :param max_gc_bound: maximum GC-content of a read to pass the filtration.
    :return: True if a read passed the filtration parameters, otherwise False.
    """
    if len(read_seq) >= min_length:
        gc_content = calculate_gc(read_seq)
        if (gc_content >= min_gc_bound) and (gc_content <= max_gc_bound):
            return True
        else:
            return False
    else:
        return False


def test_read_seq_line(read_seq):
    """
    Check whether a read sequence line from the FASTQ file contains bases A, T, C, G, or N.
    :param read_seq: a read sequence line from the FASTQ file.
    :return: True if a read sequence consists of A, T, C, G, or N bases, otherwise False.
    """
    for base in read_seq:
        if base not in ['A', 'T', 'C', 'G', 'N']:
            return False
    else:
        return True


def print_error_message(fastq_file, number_line, error_header=False,
                        error_read=False, error_sep=False, error_qual=False):
    """
    Prints out the error message if an input FASTQ file occurs to be corrupted.
    :param fastq_file: name of an input FASTQ file provided for filtering read sequences.
    :param number_line: specifies the number of the line that does not meet check criteria.
    :param error_header: True if a header line does not start with the '@' symbol, otherwise False.
    :param error_read: True if a read sequence does not contain bases A, T, G, C, or N; otherwise False.
    :param error_sep: True if a separator line is not a plus (+) sign.
    :param error_qual: True if a quality line is not of the same length as a read sequence line.
    :return: None if there no error was detected in the line.
    """
    if error_header:
        print(f"\nError! FASTQ file {fastq_file} seems to be corrupted. "
              f"\nLine number {number_line + 1}: a header line should start with the '@' symbol. "
              f"\nPlease fix the file, then try again.")
        exit()
    elif error_read:
        print(f"\nError! FASTQ file {fastq_file} seems to be corrupted. "
              f"\nLine number {number_line + 1}: a read sequence should contain only the following bases: A, T, G, C, or N."
              f"\nPlease fix the file, then try again.")
        exit()
    elif error_sep:
        print(f"\nError! FASTQ file {fastq_file} seems to be corrupted. "
              f"\nLine number {number_line + 1}: a separator line should be a plus (+) sign. "
              f"\nPlease fix the file, then try again.")
        exit()
    elif error_qual:
        print(f"\nError! FASTQ file {fastq_file} seems to be corrupted. "
              f"\nLine number {number_line + 1}: a quality line should be of the same length as a read sequence line. "
              f"\nPlease fix the file, then try again.")
        exit()
    else:
        return None
