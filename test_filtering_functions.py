import unittest
import os

from filtering_functions import check_read, is_filtered, \
    calculate_gc, filter_reads

class TestFilteringFunctions(unittest.TestCase):
    def setUp(self):
        self.header = '@SRRABCABCABC?'
        self.header_wrong = '$@SRR'
        self.header_line_number = 1
        self.sequence = "ATGCAACATCAGCT"
        self.sequence_wrong1 = "ATGCAACATCAGCTTTTTT"
        self.sequence_wrong2 = "ATGCAACATCAX12"
        self.sequence_line_number = 2
        self.separator = '+'
        self.separator_wrong = '^'
        self.separator_line_number = 3
        self.quality = "@#$AAAAZATI?!%"
        self.quality_line_number = 4
        self.file_fastq = 'test_file.fastq'
        self.output_basename = 'test_file'


    def test_check_read_true(self):
        self.assertTrue(check_read(self.header,
                                   self.header_line_number,
                                   self.sequence,
                                   self.sequence_line_number,
                                   self.separator,
                                   self.separator_line_number,
                                   self.quality,
                                   self.quality_line_number,
                                   self.file_fastq
                                   ))

    def test_check_read_header_error(self):
        self.assertRaises(ValueError, check_read, self.header_wrong,
                          self.header_line_number, self.sequence,
                          self.sequence_line_number, self.separator,
                          self.separator_line_number, self.quality,
                          self.quality_line_number, self.file_fastq
                          )

    def test_check_read_sequence_longer_quality_error(self):
        self.assertRaises(ValueError, check_read, self.header,
                          self.header_line_number, self.sequence_wrong1,
                          self.sequence_line_number, self.separator,
                          self.separator_line_number, self.quality,
                          self.quality_line_number, self.file_fastq
                          )

    def test_check_read_sequence_bases_error(self):
        self.assertRaises(ValueError, check_read, self.header,
                          self.header_line_number, self.sequence_wrong2,
                          self.sequence_line_number, self.separator,
                          self.separator_line_number, self.quality,
                          self.quality_line_number, self.file_fastq
                          )


    def test_check_read_separator_error(self):
        self.assertRaises(ValueError, check_read, self.header,
                          self.header_line_number, self.sequence,
                          self.sequence_line_number, self.separator_wrong,
                          self.separator_line_number, self.quality,
                          self.quality_line_number, self.file_fastq
                          )


    def test_calculate_gc(self):
        self.assertEqual(calculate_gc(self.sequence), 42.857142857142854)

    def test_calculate_gc_zero_length(self):
        self.assertRaises(ValueError, calculate_gc, '')

    def test_is_filtered_min_gc_included(self):
        self.assertTrue(is_filtered(self.sequence,
                                         13,
                                         42.857142857142854,
                                         50))

    def test_is_filtered_max_gc_included(self):
        self.assertTrue(is_filtered(self.sequence,
                                    13,
                                    40,
                                    42.857142857142854))

    def test_is_filtered_length_not_included(self):
        self.assertFalse(is_filtered(self.sequence,
                                    14,
                                    40,
                                    42.857142857142854))



    def test_filter_read_reads_number(self):
        self.assertTrue(filter_reads(self.file_fastq, self.output_basename,
                                     min_length=5, min_gc_bound=40,
                                     max_gc_bound=50, keep_filt=True), (2, 5))

    def test_filter_read_output_crated(self):
        filter_reads(self.file_fastq, self.output_basename,
                     min_length=5, min_gc_bound=40,
                     max_gc_bound=50, keep_filt=True)
        self.assertTrue(os.path.exists(self.output_basename+'__passed.fastq') and
                        os.path.exists(self.output_basename+'__failed.fastq'))

if __name__ == '__main__':
    unittest.main()
