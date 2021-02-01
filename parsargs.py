from os.path import exists


def parse_args(arguments, argument_dict):
    if len(arguments) == 1 or arguments[1] == '--help':
        find_help()
        exit()
    used_args = []
    file = find_fastq_file(arguments)
    argument_dict['--output_base_name'] = file.replace('.fastq', '').split('/')[-1]
    used_args += [file]
    for arg in arguments[1:-1]:
        if arg in used_args:
            print(f'\nArgument {arg} is given more than once.\n'
                  f'Please try again')
            exit()
        elif arg.startswith('--'):
            if arg not in argument_dict.keys():
                print(f'\nUnknown argument used: {arg} \n'
                      f'List of possible arguments can be checked with --help\n'
                      f'Please try again')
                exit()
            else:
                ind = arguments.index(arg)
                value = arguments[ind + 1]
                if arg != '--keep_filtered' and value.startswith('--'):
                    print(f'\nNo value(s) for {arg} is/are defined. Please try again')
                    exit()
                parsing_result = globals()[arg.replace('--', 'find_')](arguments, ind)
                argument_dict[arg] = parsing_result[0]
                last_value_ind = parsing_result[1]
                used_args += [arg]
                if arguments[last_value_ind + 1] not in argument_dict.keys() and arguments[last_value_ind + 1] not in used_args:
                    if arguments[last_value_ind + 1].startswith('--'):
                        print(f'\nUnknown argument used: {arguments[last_value_ind + 1]}\n'
                              f'List of possible arguments can be checked with --help\n'
                              f'Please try again')
                        exit()
                    else:
                        print(f'\nUnknown value {arguments[last_value_ind + 1]} is used for {arg}\n'
                              f'Please try again')
                        exit()
        elif arguments.index(arg) == 1:
            print(f'\nUnknown value {arg} is used without any argument\n'
                  f'Please try again')
            exit()
    return argument_dict, file


def find_gc_bounds(args, ind):
    start, stop = ind + 1, ind + 2
    bounds = []
    try:
        float(args[start])
        if float(args[start]) < 0:
            print('\nMinimum GC content threshold must be higher or equal to 0%. \n'
                  'Please check --gc_bounds parameters and try again')
            exit()
        bounds += [float(args[start])]
        bounds += [100.0]
        try:
            float(args[stop])
            if float(args[stop]) > 100:
                print('\nMaximum GC content threshold must be lower or equal to 100%. \n'
                      'Please check --gc_bounds parameters and try again')
                exit()
            bounds[1] = float(args[stop])
            last_value_ind = stop
        except ValueError:
            print(f'\nOnly minimum threshold ({bounds[0]}%) will be used for GC filtering.')
            last_value_ind = start
        if bounds[1] <= bounds[0]:
            print('\nMaximum GC content threshold must be higher than minimum GC content threshold. \n'
                  'Please check --gc_bounds parameters and try again')
            exit()
    except ValueError:
        last_value_ind = ind
    return bounds, last_value_ind


def find_min_length(args, ind):
    try:
        int(args[ind + 1])
    except ValueError:
        print('\nLength for filtering must be a positive integer. \n'
              'Please check --min_length parameters and try again')
        exit()
    last_value_ind = ind + 1
    return int(args[ind + 1]), last_value_ind


def find_output_base_name(args, ind):
    name = args[ind + 1]
    last_value_ind = ind + 1
    return name, last_value_ind


def find_keep_filtered(args, ind):
    last_value_ind = ind
    return True, last_value_ind


def find_fastq_file(args):
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
