import unittest
import mock
import parsing_functions

# To pass all tests you should provide a path to test_file.fastq:
path = "/home/daria/Documents/git_repos/BI_2020-2021_Python_fastq-filtrator/"


class TestGcBounds(unittest.TestCase):
    def setUp(self):
        self.ind = 1

    def test_gc_bound_with_min_and_max_thresholds(self):
        args = ["filter_fastq.py", "--gc_bounds", "30", "70", "test_file.fastq"]

        self.assertEqual(parsing_functions.get_gc_bounds(args, self.ind), ([30, 70], 3))

    def test_get_gc_bounds_with_only_min_threshold(self):
        args = ["filter_fastq.py", "--gc_bounds", "50", "some_word", "test_file.fastq"]

        self.assertEqual(parsing_functions.get_gc_bounds(args, self.ind), ([50, 100], 2))

    def test_get_gc_bounds_with_negative_min_threshold(self):
        args = ["filter_fastq.py", "--gc_bounds", "-5", "50", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.get_gc_bounds(args, self.ind)
        obtained_msg = ctx.exception.args[0]
        expected_msg = '\nMinimum GC content threshold must be higher or equal to 0%. \n' \
                       'Please check --gc_bounds parameters and try again'
        self.assertEqual(obtained_msg, expected_msg)

    def test_get_gc_bounds_with_not_int_value(self):
        args = ["filter_fastq.py", "--gc_bounds", "as", "50", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.get_gc_bounds(args, self.ind)
        obtained_msg = ctx.exception.args[0]
        expected_msg = '\nMinimum GC content threshold for filtering must be a positive float number. \n' \
                       'Please check --min_length parameter and try again'
        self.assertEqual(obtained_msg, expected_msg)

    def test_get_gc_bounds_with_max_threshold_higher_than_100(self):
        args = ["filter_fastq.py", "--gc_bounds", "10", "150", "test_file.fastq"]

        self.assertRaises(ValueError, parsing_functions.get_gc_bounds, args, self.ind)

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.get_gc_bounds(args, self.ind)
        obtained_msg = ctx.exception.args[0]
        expected_msg = '\nMaximum GC content threshold must be lower or equal to 100%. \n' \
                       'Please check --gc_bounds parameters and try again'
        self.assertEqual(obtained_msg, expected_msg)

    def test_get_gc_bounds_with_min_threshold_higher_than_max_threshold(self):
        args = ["filter_fastq.py", "--gc_bounds", "105", "15", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.get_gc_bounds(args, self.ind)
        obtained_msg = ctx.exception.args[0]
        expected_msg = '\nMaximum GC content threshold must be higher than minimum GC content threshold. \n' \
                       'Please check --gc_bounds parameters and try again'
        self.assertEqual(obtained_msg, expected_msg)


class TestGetMinLength(unittest.TestCase):
    def setUp(self):
        self.ind = 1

    def test_get_min_length_with_not_int_value(self):
        args = ["filter_fastq.py", "--min_length", "some_word", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.get_min_length(args, self.ind)
        obtained_msg = ctx.exception.args[0]
        expected_msg = '\nLength for filtering must be of integer type. \n' \
                       'Please check --min_length parameters and try again'
        self.assertEqual(obtained_msg, expected_msg)

    def test_get_min_length_with_negative_min_threshold(self):
        args = ["filter_fastq.py", "--min_length", "-50", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.get_min_length(args, self.ind)
        obtained_msg = ctx.exception.args[0]
        expected_msg = '\nLength for filtering must be a positive integer. \n' \
                       'Please check --min_length parameters and try again'
        self.assertEqual(obtained_msg, expected_msg)


class TestGetFastqFile(unittest.TestCase):
    def test_get_fastq_file_not_in_fastq_format(self):
        args = ["filter_fastq.py", "--min_length", "35", "some_file.not_fastq_format"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.get_fastq_file(args)
        obtained_msg = ctx.exception.args[0]
        expected_msg = '\nNo fastq file for filtering is provided.\n' \
                       'Please try again and put fastq file as THE LAST argument'
        self.assertEqual(obtained_msg, expected_msg)

    def test_get_fastq_file_which_does_not_exist(self):
        args = ["filter_fastq.py", "--min_length", "35", "non-existing_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.get_fastq_file(args)
        obtained_msg = ctx.exception.args[0]
        expected_msg = f'{args[-1]} does not exist.'
        self.assertEqual(obtained_msg, expected_msg)


class TestParseArgs(unittest.TestCase):
    def test_parse_args_can_extract_output_basename(self):

        args = ["filter_fastq.py", "--min_length", "35", f"{path}test_file.fastq"]
        parsed_args = parsing_functions.parse_args(args)
        self.assertEqual(parsed_args[0]['--output_base_name'], "test_file")

    def test_parse_args_when_argument_assigned_more_than_once(self):
        args = ["filter_fastq.py", "--min_length", "45", "--min_length", "35", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.parse_args(args)
        obtained_msg = ctx.exception.args[0]
        expected_msg = f'\nArgument --min_length is given more than once.\n' \
                       f'Please try again'
        self.assertEqual(obtained_msg, expected_msg)

    def test_parse_args_when_argument_does_not_have_double_slash(self):
        args = ["filter_fastq.py", "-min_length", "35", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.parse_args(args)
        obtained_msg = ctx.exception.args[0]
        expected_msg = f'\nUnknown value used: -min_length\n' \
                       f'Please try again'
        self.assertEqual(obtained_msg, expected_msg)

    def test_parse_args_with_unknown_argument(self):
        args = ["filter_fastq.py", "--hello_im_new_argument", "35", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.parse_args(args)
        obtained_msg = ctx.exception.args[0]
        expected_msg = f'\nUnknown argument used: --hello_im_new_argument \n' \
                       f'List of possible arguments can be checked with --help\n' \
                       f'Please try again'
        self.assertEqual(obtained_msg, expected_msg)

    def test_parse_args_can_remember_keep_filtered_flag(self):
        args = ["filter_fastq.py", "--min_length", "35", "--keep_filtered", "test_file.fastq"]
        parsed_args = parsing_functions.parse_args(args)
        self.assertTrue(parsed_args[0]['--keep_filtered'])

    def test_parse_args_can_do_without_keep_filtered_flag(self):
        args = ["filter_fastq.py", "--min_length", "35", "test_file.fastq"]
        parsed_args = parsing_functions.parse_args(args)
        self.assertFalse(parsed_args[0]['--keep_filtered'])

    def test_parse_args_when_argument_is_given_without_value(self):
        args_min_length = ["filter_fastq.py", "--min_length", "test_file.fastq"]
        args_min_gc_bound = ["filter_fastq.py", "--gc_bounds", "test_file.fastq"]
        args_output_base_name = ["filter_fastq.py", "--output_base_name", "test_file.fastq"]

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.parse_args(args_min_length)
        obtained_msg = ctx.exception.args[0]
        expected_msg = f'\nNo value(s) for --min_length is/are defined. Please try again'
        self.assertEqual(obtained_msg, expected_msg)

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.parse_args(args_min_gc_bound)
        obtained_msg = ctx.exception.args[0]
        expected_msg = f'\nNo value(s) for --gc_bounds is/are defined. Please try again'
        self.assertEqual(obtained_msg, expected_msg)

        with self.assertRaises(ValueError) as ctx:
            parsing_functions.parse_args(args_output_base_name)
        obtained_msg = ctx.exception.args[0]
        expected_msg = f'\nNo value(s) for --output_base_name is/are defined. Please try again'
        self.assertEqual(obtained_msg, expected_msg)

    def test_parse_args_offers_help(self):
        expected_output = ('\nfilter_fastq.py filtrates reads in .fastq format by GC content and length.\n\n'
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

        only_one_arg = ["filter_fastq.py"]
        args_with_help = ["filter_fastq.py", "--I_really_dont_know", "HoWto_use_it", "I_need", "--help"]

        with mock.patch('builtins.input'):
            assert parsing_functions.parse_args(only_one_arg) == expected_output

        with mock.patch('builtins.input'):
            assert parsing_functions.parse_args(args_with_help) == expected_output


if __name__ == '__main__':
    unittest.main()
