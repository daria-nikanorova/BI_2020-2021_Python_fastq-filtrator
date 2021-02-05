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


# PART 2. READS FILTRATION FUNCTIONS


def calculate_gc(read_seq):
    """
    Calculates GC-content as a percentage of G or C bases in DNA sequence: Count(G + C)/Count(A + T + G + C) * 100%
    :param read_seq: a read sequence from the FASTQ file
    :return: GC-content of a read sequence, %
    """

    return (read_seq.count('G') + read_seq.count('C')) * 100 / len(read_seq)


def is_read_filtered(read_seq, length, min_gc, max_gc):
    """
    Check whether a read sequence passes the filtration by its length and GC-content
    :param read_seq: a read sequence from the FASTQ file
    :param length: minimum length for a read to pass the filtration
    :param min_gc: minimum GC-content value of a read to pass the filtration.
    :param max_gc: maximum GC-content of a read to pass the filtration.
    :return: True if a read passed the filtration parameters, otherwise False.
    """
    if len(read_seq) >= length:
        gc_content = calculate_gc(read_seq)
        if (gc_content >= min_gc) and (gc_content <= max_gc):
            return True
        else:
            return False
    else:
        return False


def is_correct_seq(read_seq, number_of_line, file):
    """
    Check whether a read sequence line from the FASTQ file contains bases A, T, C, G, or N
    :param read_seq: a read sequence line from the FASTQ file
    :param number_of_line: number of line with seq in FASTQ file
    :param file: FASTQ file
    :return: True if a read sequence consists of A, T, C, G, V or N bases, otherwise False
    """
    for base in read_seq:
        if base in ['A', 'T', 'C', 'G', 'N', 'V']:
            return True
        else:
            print(f"\nError! FASTQ file {file} seems to be corrupted. "
                  f"\nLine number {number_of_line}: a read sequence contains {base}, but should contain only A, T, G, C, V or N."
                  f"\nPlease fix the file, then try again.")
            exit()


def is_correct_separator(separator, number_of_line, file):
    """
    Check whether a separator is "+"
    :param separator: a separator line from the FASTQ file
    :param number_of_line: number of line with separator in FASTQ file
    :param file: FASTQ file
    :return: True if separator is "+", otherwise False
    """
    if separator.rstrip() == '+':
        return True
    else:
        print(f"\nError! FASTQ file {file} seems to be corrupted. "
              f"\nLine number {number_of_line}: a separator line should be a plus (+) sign. "
              f"\nPlease fix the file, then try again.")
        exit()


def is_correct_header(header, number_of_line, file):
    """
    Check whether a header starts with "@"
    :param header: a separator line from the FASTQ file
    :param number_of_line: number of line with header in FASTQ file
    :param file: FASTQ file
    :return: True if a header starts with "@", otherwise False
    """
    if header.startswith('@'):
        return True
    else:
        print(f"\nError! FASTQ file {file} seems to be corrupted. "
              f"\nLine number {number_of_line}: a header line should start with the '@' symbol. "
              f"\nPlease fix the file, then try again.")
        exit()


def is_correct_quality_line(quality_line, read_len, number_of_line, file):
    """
    Check whether a quality_line has the same length as sequence
    :param quality_line: a separator line from the FASTQ file
    :param read_len: a length of the sequence
    :param number_of_line: number of line with quality line in FASTQ file
    :param file: FASTQ file
    :return: True if a quality_line has the same length as sequence, otherwise False
    """
    if read_len == len(quality_line.rstrip()):
        return True
    else:
        print(f"\nError! FASTQ file {file} seems to be corrupted. "
              f"\nLine number {number_of_line}: a quality line should be of the same length as a read sequence line. "
              f"\nPlease fix the file, then try again.")
        exit()


def filter_reads(file_fastq, output_basename, min_len, gc_bound_min, gc_bound_max, keep_filt):
    """
    Filter reads by its length and GC-content
    :param file_fastq: the FASTQ file
    :param output_basename: basename for output FASTQ file(s)
    :param min_len: minimum length for a read to pass the filtration
    :param gc_bound_min: minimum GC-content value of a read to pass the filtration
    :param gc_bound_max: maximum GC-content of a read to pass the filtration
    :param keep_filt: if used, reads that fail filtering will be written to file
    :return: True if filtration process was successfully finished
    """
    with open(file_fastq) as fastq_input:
        print(f"\nRead sequences that pass the filtration will be written to the file:\n "
              f"{output_basename}__passed.fastq\n")
        if keep_filt:
            fastq_failed = open(output_basename + '__failed.fastq', 'w')
            print(f"Read sequences that fail the filtration will be written to the file:\n "
                  f"{output_basename}__failed.fastq")
        with open(output_basename + '__passed.fastq', 'w') as fastq_passed:
            index = 0  # current index of the line (0-3)
            line_number = 0  # exact line number, as in input file (0-length(input file))
            number_passed_reads = 0
            number_all_reads = 0
            tmp_lines = []  # to save temporarily a block of 4 lines for a read passed the filtration
            for line in fastq_input:
                line_number += 1
                if (index == 0) and (is_correct_header(line.rstrip(), line_number, file_fastq)):  # a header line
                    tmp_lines.append(line.rstrip())
                    index += 1
                    continue
                if (index == 1) and (is_correct_seq(line.rstrip(), line_number, file_fastq)):  # a sequence line
                    number_all_reads += 1
                    read_passed_filtration = is_read_filtered(line.rstrip(), min_len, gc_bound_min, gc_bound_max)
                    tmp_lines.append(line.rstrip())
                    read_length = len(line.rstrip())
                    index += 1
                    continue
                if (index == 2) and (is_correct_separator(line.rstrip(), line_number, file_fastq)):  # a separator line
                    tmp_lines.append(line.rstrip())
                    index += 1
                    continue
                if (index == 3) and (is_correct_quality_line(line.rstrip(), read_length, line_number, file_fastq)):  # a quality line
                    tmp_lines.append(line.rstrip())
                if read_passed_filtration:
                    number_passed_reads += 1
                    fastq_passed.write('\n'.join(tmp_lines) + '\n')
                elif not read_passed_filtration and keep_filt:
                    fastq_failed.write('\n'.join(tmp_lines) + '\n')
                tmp_lines = []
                index = 0
            if keep_filt:
                fastq_failed.close()
    print(f"\n{number_passed_reads} out of {number_all_reads} read sequences passed the filtration "
          f"(about {round(number_passed_reads / number_all_reads * 100, 2)}%) ")
    pass



