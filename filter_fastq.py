import parsing_functions
import filtering_functions
from sys import argv

# ARGUMENTS PARSING

# the arguments are stored in sys.argv
unparsed_args = argv

# parse arguments, values will be stored in parsed_args
parsed_args = parsing_functions.parse_args(unparsed_args)

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

# FILTRATION

print("\nThe filtration process has started. Please be patient, it will take some time.")

number_passed_reads, number_all_reads = filtering_functions.filter_reads(fastq_file, output_base_name, min_length, min_gc_bound, max_gc_bound, keep_filtered)

print(f"\n{number_passed_reads} out of {number_all_reads} read sequences passed the filtration "
      f"(about {round(number_passed_reads / number_all_reads * 100, 2)}%) ")

print('\nThank you for using filter_fastq.py!')

exit()
