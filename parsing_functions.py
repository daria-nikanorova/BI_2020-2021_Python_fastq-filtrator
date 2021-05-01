from os.path import exists


def get_gc_bounds(args, ind):
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
        raise ValueError('\nMinimum GC content threshold for filtering must be a positive float number. \n'
                         'Please check --min_length parameter and try again')
    if float(args[start]) < 0:
        raise ValueError('\nMinimum GC content threshold must be higher or equal to 0%. \n'
                         'Please check --gc_bounds parameters and try again')
    bounds[0] = float(args[start])
    last_value_ind = start
    try:
        float(args[stop])
    except ValueError:
        return bounds, last_value_ind
    if float(args[stop]) > 100:
        raise ValueError('\nMaximum GC content threshold must be lower or equal to 100%. \n'
                         'Please check --gc_bounds parameters and try again')
    bounds[1] = float(args[stop])
    last_value_ind = stop
    if bounds[1] <= bounds[0]:
        raise ValueError('\nMaximum GC content threshold must be higher than minimum GC content threshold. \n'
                         'Please check --gc_bounds parameters and try again')
    return bounds, last_value_ind


def get_min_length(args, ind):
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
        raise ValueError('\nLength for filtering must be of integer type. \n'
                         'Please check --min_length parameters and try again')
    if int(args[ind + 1]) < 0:
        raise ValueError('\nLength for filtering must be a positive integer. \n'
                         'Please check --min_length parameters and try again')

    last_value_ind = ind + 1
    return int(args[ind + 1]), last_value_ind


def get_output_base_name(args, ind):
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


def get_fastq_file(args):
    """
    Parses the last argument fastq file with reads for filtration, is called by parse_args
    :param args: a list with all input arguments (may be a direct output from sys.argv())
    :return: parsed argument fastq file
    """
    if args[-1].endswith('.fastq'):
        if exists(args[-1]):
            file = args[-1]
            return file
        else:
            raise ValueError(f'{args[-1]} does not exist.')
    else:
        raise ValueError('\nNo fastq file for filtering is provided.\n'
                         'Please try again and put fastq file as THE LAST argument')


def print_help():
    """
    :return: None
    """
    return ('\nfilter_fastq.py filtrates reads in .fastq format by GC content and length.\n\n'
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


functions = {'--gc_bounds': get_gc_bounds,
             '--min_length': get_min_length,
             '--output_base_name': get_output_base_name}


def parse_args(arguments):
    """
    Looks through all input arguments and calls functions to parse them by name
    :param arguments: a list with all input arguments (may be a direct output from sys.argv())
    :return: a dictionary of parsed arguments and their values
    """
    # a dict that contains all optional arguments and their values by default
    default_args = {
        '--min_length': 0,
        '--gc_bounds': [0.0, 100.0],
        '--keep_filtered': False,
        '--output_base_name': ''
    }
    if len(arguments) == 1 or '--help' in arguments:  # call help page
        return print_help()
    used_args = []
    file = get_fastq_file(arguments)  # finds file for filtering (a positional argument)
    default_args['--output_base_name'] = file.replace('.fastq', '').split('/')[-1]
    used_args += [file]
    ind = 1
    while ind + 1 != len(arguments):  # loop over all arguments except fastq file and the 1st arg (path to a script)
        arg = arguments[ind]
        next_value = arguments[ind + 1]
        if arg in used_args:  # check if argument was assigned before
            raise ValueError(f'\nArgument {arg} is given more than once.\n'
                             f'Please try again')
        if not arg.startswith('--'):  # check if arg is an argument, not a value
            raise ValueError(f'\nUnknown value used: {arg}\n'
                             f'Please try again')
        if arg not in default_args:
            raise ValueError(f'\nUnknown argument used: {arg} \n'
                             f'List of possible arguments can be checked with --help\n'
                             f'Please try again')
        if arg == '--keep_filtered':  # this flag is too small for special function
            default_args[arg] = True
        elif next_value.startswith('--') or ind + 2 == len(arguments):  # check if there is a value after an argument
            raise ValueError(f'\nNo value(s) for {arg} is/are defined. Please try again')
        else:
            parsing_result = functions[arg](arguments, ind)
            default_args[arg] = parsing_result[0]
            ind = parsing_result[1]  # get an index of the last value for an argument
        used_args += [arg]
        ind += 1  # make 1 step to the next arg
    return default_args, file
