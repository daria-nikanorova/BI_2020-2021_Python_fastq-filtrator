import sys


def is_a_number(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def parse_args(arguments, argument_dict):
    for arg in arguments:
        if arg == "--keep_filtered":
            argument_dict[arg] = True
        elif arg in argument_dict.keys():
            ind = arguments.index(arg)
            value = arguments[ind + 1]
            if value in argument_dict.keys():
                print(f'No value(s) for {arg} is/are defined. Please try again')
                exit()
            argument_dict[arg] = globals()[arg.replace('--', 'find_')](arguments, ind)
    return argument_dict


def find_gc_bounds(args, ind):
    start, stop = ind + 1, ind + 2
    bounds = []
    if is_a_number(args[start]):
        bounds += [float(args[start])]
        bounds += [100.0]
        if is_a_number(args[stop]):
            bounds[1] = float(args[stop])
        else:
            print('NB: only minimum threshold for filtering by GC count is defined')
    else:
        print('The first value for --gc_bounds is not numeric. Please try again')
        exit()
    return bounds


def find_min_length(args, ind):
    if is_a_number(args[ind + 1]):
        return int(args[ind + 1])
    else:
        print('The value for --min_length is not numeric. Please try again')
        exit()


def find_output_base_name(args, ind):
    return args[ind + 1]


# ARGUMENTS PARSING

# fastq file is always the positional argument (the last one):
if sys.argv[-1].endswith('.fastq'):
    fastq_file = sys.argv[-1]
    output_base_name = fastq_file.replace('.fastq', '')
else:
    print('No fastq file for filtering is provided. Please try again and put fastq file as THE LAST argument')
    exit()

# the others are optional (the 1st argument, which is a path, is excluded)
unparsed_args = sys.argv[1:-1]

# a dict that contains all optional arguments and their values by default
default_args = {
    '--min_length': 0,
    '--gc_bounds': [0.0, 100.0],
    '--keep_filtered': False,
    '--output_base_name': output_base_name
}

# parse optional arguments, values will be stored in parsed_args
parsed_args = parse_args(unparsed_args, default_args.copy())

# optional args can be found in parsed_args:
min_length = parsed_args['--min_length']
gc_bounds = parsed_args['--gc_bounds']
min_gc_bound = parsed_args['--gc_bounds'][0]
max_gc_bound = parsed_args['--gc_bounds'][1]
keep_filtered = parsed_args['--keep_filtered']
output_base_name = parsed_args['--output_base_name']

print(f'{fastq_file} will be filtered with further parameters: \n'
      f'--min_length = {min_length}, \n'
      f'--gc_bounds = {gc_bounds}, \n'
      f'--keep_filtered = {keep_filtered},\n'
      f'--output_base_name = {output_base_name}')
