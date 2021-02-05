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

# parse arguments, values will be stored in parsed_args
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

# FILTRATION

print("\nThe filtration process has started. Please be patient, it will take some time.")

with open(fastq_file) as fastq_input:
    print(f"\nRead sequences that pass the filtration would be written to the file:\n "
          f"{output_base_name}__passed.fastq")
    if keep_filtered:  # create a file to save FASTQ 4-lines blocks where reads failed the filtration (if required)
        print(f"\nRead sequences that fail the filtration would be written to the file:\n "
              f"{output_base_name}__failed.fastq")
        fastq_failed = open(output_base_name + '__failed.fastq', 'w')
    with open(output_base_name + '__passed.fastq', 'w') as fastq_passed:
        count = 0  # count value specifies certain line in FASTQ file
        number_passed_reads = 0  # number of reads passed the filtration
        number_line = 0  # number of the line
        tmp_lines_passed = []  # to save temporarily a block of 4 lines for a read passed the filtration
        tmp_lines_failed = []  # to save temporarily a block of 4 lines for a read failed the filtration (if required)
        for line in fastq_input:
            if count == 0:  # a header line
                if line.rstrip().startswith('@'):  # test if a header line starts with the '@' symbol
                    tmp_lines_passed.append(line.rstrip())
                    count += 1
                else:
                    par.print_error_message(fastq_file, number_line, error_header=True)
            elif count == 1:  # a read sequence line; check whether its length and GC-content pass the filtration
                if par.test_read_seq_line(line.rstrip()):  # test if a read sequence line contains A, T, C, G, or N bases
                    read_check_flag = par.pass_read_check(line.rstrip(), min_length, min_gc_bound, max_gc_bound)
                    if read_check_flag:
                        tmp_lines_passed.append(line.rstrip())
                        count += 1
                        number_passed_reads += 1
                    elif (not read_check_flag) and keep_filtered:
                        tmp_lines_failed.append(tmp_lines_passed[0])
                        tmp_lines_failed.append(line.rstrip())
                        count += 3
                    else:
                        tmp_lines_failed.append(line.rstrip())  # save only a read to further check the quality line
                        count += 5
                else:
                    par.print_error_message(fastq_file, number_line, error_read=True)
            elif count == 2:  # "+" line, from a 4-lines block where a read passed the filtration
                if line.rstrip() == '+':  # test if a separator line is a plus (+) sign
                    count += 1
                    tmp_lines_passed.append(line.rstrip())
                else:
                    par.print_error_message(fastq_file, number_line, error_sep=True)
            elif count == 3:  # quality line, from a 4-lines block where a read passed the filtration
                if len(tmp_lines_passed[1]) == len(line.rstrip()):  # compare the length of quality and read lines
                    tmp_lines_passed.append(line.rstrip())
                    fastq_passed.write('\n'.join(tmp_lines_passed) + '\n')
                    tmp_lines_passed = []
                    count = 0
                else:
                    par.print_error_message(fastq_file, number_line, error_sep=True)
            elif count == 4:  # "+" line, from a 4-lines block where a read failed the filtration; keep failed reads
                if line.rstrip() == '+':  # test if a separator line is a plus (+) sign
                    tmp_lines_failed.append(line.rstrip())
                    count += 1
                else:
                    par.print_error_message(fastq_file, number_line, error_sep=True)
            elif count == 5:  # quality line, from a 4-lines block where a read failed the filtration; keep failed reads
                if len(tmp_lines_failed[1]) == len(line.rstrip()):  # compare the length of quality and read lines
                    tmp_lines_failed.append(line.rstrip())
                    fastq_failed.write('\n'.join(tmp_lines_failed) + '\n')
                    tmp_lines_failed = []
                    tmp_lines_passed = []
                    count = 0
                else:
                    par.print_error_message(fastq_file, number_line, error_qual=True)
            elif count == 6:  # "+" line, from a 4-lines block where a read failed the filtration; not save to a file
                if line.rstrip() == '+':  # test if a separator line is a plus (+) sign
                    count += 1
                else:
                    par.print_error_message(fastq_file, number_line, error_sep=True)
            elif count == 7:  # quality line, from a 4-line block where a read failed the filtration; not save to a file
                if len(tmp_lines_failed[0]) == len(line.rstrip()):  # compare the length of quality and read lines
                    tmp_lines_failed = []
                    tmp_lines_passed = []
                    count = 0
                else:
                    par.print_error_message(fastq_file, number_line, error_qual=True)
            number_line += 1
        if keep_filtered:
            fastq_failed.close()

if (number_line % 4) == 0:
    print(f"\n{number_passed_reads} out of {int(number_line/4)} read sequences passed the filtration "
          f"(about {round(number_passed_reads/(number_line/4) * 100, 2)}%) ")
