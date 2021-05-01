# PART 2. READS FILTRATION FUNCTIONS

def is_filtered(sequence, min_length, min_gc_bound, max_gc_bound):
    """
    Check whether a read sequence passes the filtration by its length and GC-content
    :param sequence: a read sequence from the FASTQ file
    :param min_length: minimum length for a read to pass the filtration
    :param min_gc_bound: minimum GC-content value of a read to pass the filtration.
    :param max_gc_bound: maximum GC-content of a read to pass the filtration.
    :return: True if a read passed the filtration parameters, otherwise False.
    """
    gc_content = calculate_gc(sequence)
    return (gc_content >= min_gc_bound) and (gc_content <= max_gc_bound) and (len(sequence) > min_length)


def check_read(header, header_line_number, sequence, sequence_line_number, separator, separator_line_number,
               quality, quality_line_number, file_fastq):
    """
    Check whether all 4 lines of a read do not have mistakes
    :param header: a header line from the FASTQ file
    :param header_line_number: a line number of a header
    :param sequence: a sequence line from the FASTQ file
    :param sequence_line_number: a line number of a sequence
    :param separator: a separator line from the FASTQ file
    :param separator_line_number: a line number of a separator
    :param quality: a quality line from the FASTQ file
    :param quality_line_number: a line number of a quality line
    :param file_fastq: the FASTQ file
    :return: True if all 4 read lines are not corrupted.
    """

    if not header.startswith('@'):
        raise ValueError(f"\nError! FASTQ file {file_fastq} seems to be corrupted."
                         f"\nLine number {header_line_number}: a header line must start with the '@' symbol."
                         f"\nPlease fix the file, then try again.")
    for base in sequence:
        if base not in ('A', 'T', 'C', 'G', 'N'):
            raise ValueError(f"\nError! FASTQ file {file_fastq} seems to be corrupted."
                             f"\nLine number {sequence_line_number}: a read sequence contains {base}, but must contain only A, T, G, C or N."
                             f"\nPlease fix the file, then try again.")
    if separator != '+':
        raise ValueError(f"\nError! FASTQ file {file_fastq} seems to be corrupted."
                         f"\nLine number {separator_line_number}: a separator line should be a plus (+) sign. "
                         f"\nPlease fix the file, then try again.")
    if len(quality) != len(sequence):
        raise ValueError(f"\nError! FASTQ file {file_fastq} seems to be corrupted."
                         f"\nLine number {quality_line_number}: a quality line should be of the same length as a read sequence line. "
                         f"\nPlease fix the file, then try again.")
    return True


def calculate_gc(read_seq):
    """
    Calculates GC-content as a percentage of G or C bases in DNA sequence: Count(G + C)/Count(A + T + G + C) * 100%
    :param read_seq: a read sequence from the FASTQ file
    :return: GC-content of a read sequence, %
    """
    if len(read_seq) > 0:
        return (read_seq.count('G') + read_seq.count('C')) * 100 / len(read_seq)
    else:
        raise ValueError('A sequence of zero length was found.\n'
                         'Please, check your fastq file and try again')


def filter_reads(file_fastq, output_basename, min_length, min_gc_bound, max_gc_bound, keep_filt):
    """
    Filter reads by its length and GC-content
    :param file_fastq: the FASTQ file
    :param output_basename: basename for output FASTQ file(s)
    :param min_length: minimum length for a read to pass the filtration
    :param min_gc_bound: minimum GC-content value of a read to pass the filtration
    :param max_gc_bound: maximum GC-content of a read to pass the filtration
    :param keep_filt: if used, reads that fail filtering will be written to file
    :return: number_passed_reads, number_all_reads if filtration process was successfully finished
    """
    number_passed_reads = 0
    number_all_reads = 0
    with open(file_fastq) as fastq_input, open(output_basename + '__passed.fastq', 'w') as fastq_passed:
        print(f"\nRead sequences that pass the filtration will be written to the file:\n "
              f"{output_basename}__passed.fastq")

        # Get only 4 lines as a read
        index = 0
        read = list()
        for line_number, line in enumerate(fastq_input, 1):
            read.append(line_number)
            read.append(line.rstrip())
            index += 1
            if index == 4:
                number_all_reads += 1

                # Check read quality and if it passes filtration
                if check_read(header_line_number=read[0],
                              header=read[1],
                              sequence_line_number=read[2],
                              sequence=read[3],
                              separator_line_number=read[4],
                              separator=read[5],
                              quality_line_number=read[6],
                              quality=read[7],
                              file_fastq=file_fastq) \
                    and is_filtered(sequence=read[3],
                                    min_length=min_length,
                                    min_gc_bound=min_gc_bound,
                                    max_gc_bound=max_gc_bound):

                    # Write down reads that passed filtration
                    fastq_passed.write('\n'.join(read[1:9:2]) + '\n')
                    number_passed_reads += 1

                # Write down reads that did not pass filtration (if needed)
                else:
                    if keep_filt:
                        fastq_failed = open(output_basename + '__failed.fastq', 'a')
                        fastq_failed.write('\n'.join(read[1:9:2]) + '\n')
                        fastq_failed.close()

                index = 0
                read = list()
    return number_passed_reads, number_all_reads

