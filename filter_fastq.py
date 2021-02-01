import parsargs as par
from sys import argv

# ARGUMENTS PARSING

# the arguments are stored in sys.argv
unparsed_args = argv

# a dict that contains all optional arguments and their values by default
default_args = {
    '--min_length': 0,
    '--gc_bounds': [0.0, 100.0],
    '--keep_filtered': False,
    '--output_base_name': ''
}

# parse optional arguments, values will be stored in parsed_args
parsed_args = par.parse_args(unparsed_args, default_args.copy())

# args after parsing can be found in parsed_args:
fastq_file = parsed_args[1]
optional_parsed_args = parsed_args[0]

min_length = optional_parsed_args['--min_length']
gc_bounds = optional_parsed_args['--gc_bounds']
min_gc_bound = optional_parsed_args['--gc_bounds'][0]
max_gc_bound = optional_parsed_args['--gc_bounds'][1]
keep_filtered = optional_parsed_args['--keep_filtered']
output_base_name = optional_parsed_args['--output_base_name']

print(f'\n{fastq_file} will be filtered with further parameters: \n'
      f'--min_length = {min_length}bp, \n'
      f'--gc_bounds = {min_gc_bound}% - {max_gc_bound}%, \n'
      f'--keep_filtered = {keep_filtered},\n'
      f'--output_base_name = {output_base_name}')
