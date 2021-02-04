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

# FILTRATION


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
    wrong_base = 0
    for base in read_seq:
        if base not in ['A', 'T', 'C', 'G', 'N']:
            wrong_base+=1
    if wrong_base > 0:
        return False
    else:
        return True


def print_error_message(fastq_file, number_line, error_header = False,
                        error_read = False, error_sep = False, error_qual = False):
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
              f"\nLine number {number_line}: a header line should start with the '@' symbol. "
              f"\nPlease fix the file, then try again.")
    elif error_read:
        print(f"\nError! FASTQ file {fastq_file} seems to be corrupted. "
              f"\nLine number {number_line}: a read sequence should contain only the following bases: A, T, G, C, or N."
              f"\nPlease fix the file, then try again.")
    elif error_sep:
        print(f"\nError! FASTQ file {fastq_file} seems to be corrupted. "
              f"\nLine number {number_line}: a separator line should be a plus (+) sign. "
              f"\nPlease fix the file, then try again.")
    elif error_qual:
        print(f"\nError! FASTQ file {fastq_file} seems to be corrupted. "
              f"\nLine number {number_line}: a quality line should be of the same length as a read sequence line. "
              f"\nPlease fix the file, then try again.")
    else:
        return None


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
                    print_error_message(fastq_file, number_line, error_header = True)
                    break
            elif count == 1:  # a read sequence line; check whether its length and GC-content pass the filtration
                if test_read_seq_line(line.rstrip()):  # test if a read sequence line contains A, T, C, G, or N bases
                    read_check_flag = pass_read_check(line.rstrip(), min_length, min_gc_bound, max_gc_bound)
                    if read_check_flag:
                        tmp_lines_passed.append(line.rstrip())
                        count += 1
                        number_passed_reads += 1
                    elif (not read_check_flag) and keep_filtered:
                        tmp_lines_failed.append(tmp_lines_passed[0])
                        tmp_lines_failed.append(line.rstrip())
                        count += 3
                    else:
                        tmp_lines_failed.append(line.rstrip()) # save only a read to further check the quality line
                        count += 5
                else:
                    print_error_message(fastq_file, number_line, error_read = True)
                    break
            elif count == 2:  # "+" line, from a 4-lines block where a read passed the filtration
                if line.rstrip() == '+':  # test if a separator line is a plus (+) sign
                    count += 1
                    tmp_lines_passed.append(line.rstrip())
                else:
                    print_error_message(fastq_file, number_line, error_sep = True)
                    break
            elif count == 3:  # quality line, from a 4-lines block where a read passed the filtration
                if len(tmp_lines_passed[1]) == len(line.rstrip()):  # compare the length of quality and read lines
                    tmp_lines_passed.append(line.rstrip())
                    fastq_passed.write('\n'.join(tmp_lines_passed) + '\n')
                    tmp_lines_passed = []
                    count = 0
                else:
                    print_error_message(fastq_file, number_line, error_sep = True)
                    break
            elif count == 4:  # "+" line, from a 4-lines block where a read failed the filtration; keep failed reads
                if line.rstrip() == '+':  # test if a separator line is a plus (+) sign
                    tmp_lines_failed.append(line.rstrip())
                    count += 1
                else:
                    print_error_message(fastq_file, number_line, error_sep = True)
                    break
            elif count == 5:  # quality line, from a 4-lines block where a read failed the filtration; keep failed reads
                if len(tmp_lines_failed[1]) == len(line.rstrip()):  # compare the length of quality and read lines
                    tmp_lines_failed.append(line.rstrip())
                    fastq_failed.write('\n'.join(tmp_lines_failed) + '\n')
                    tmp_lines_failed = []
                    tmp_lines_passed = []
                    count = 0
                else:
                    print_error_message(fastq_file, number_line, error_qual = True)
                    break
            elif count == 6:  # "+" line, from a 4-lines block where a read failed the filtration; not save to a file
                if line.rstrip() == '+':  # test if a separator line is a plus (+) sign
                    count += 1
                else:
                    print_error_message(fastq_file, number_line, error_sep = True)
                    break
            elif count == 7:  # quality line, from a 4-line block where a read failed the filtration; not save to a file
                if len(tmp_lines_failed[0]) == len(line.rstrip()):  # compare the length of quality and read lines
                    tmp_lines_failed = []
                    tmp_lines_passed = []
                    count = 0
                else:
                    print_error_message(fastq_file, number_line, error_qual = True)
                    break
            number_line += 1
        if keep_filtered:
            fastq_failed.close()

if (number_line % 4) == 0:
    print(f"\n{number_passed_reads} out of {int(number_line/4)} read sequences passed the filtration "
          f"(about {round(number_passed_reads/(number_line/4) * 100, 2)} %) ")
