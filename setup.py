import sys
from setuptools import setup, find_packages, Command
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    """
    pytest class for runnin py.test test.
    check out: http://pytest.org/latest/goodpractises.html
    fro more info
    """
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


setup(
    name='qiime-tools',
    version='0.5.1',
    install_requires=[
        'Click >= 0.6.0',
        'Biopython >=1.6.5',
        'cytoolz',
        'pandas'
    ],
    packages=find_packages(),
    entry_points='''
        [console_scripts]
        debarcodepairedfastq = qiime_tools.debarcodepairedfastq:debarcodepairedfastq
        demultiplexfasta = qiime_tools.demultiplexfasta:demultiplex
        demultiplexfasta_parallel = qiime_tools.demultiplexfasta:demultiplex_parallel
        fastqconcat = qiime_tools.fastq_concat:fastqconcat
        fastq_demultiplex = qiime_tools.demultiplex:demultiplex
        filterfasta_by_length = qiime_tools.filter_by_length:filter_by_length
        merge_OTU_UC = qiime_tools.mergeOTU_UC:merge_OTU_UCfile
        parallel_split_library_fastq = qiime_tools.parallel_split_libraries_fastq:parallel_splitlibraries_fastq
        parallel_split_libraries = qiime_tools.parallel_split_library:parallel_split_library
        split_fastq_by_name = qiime_tools.split_fastq_by_name:split_fastq_by_name
        split_fasta_by_name = qiime_tools.split_fasta_by_name:split_fasta_by_name
        taxtable_to_otutable = qiime_tools.UCfiles_to_taxtable:taxtable_to_otutable
        UCfiles_to_taxtable = qiime_tools.UCfiles_to_taxtable:UCfiles_to_taxtable
        UCfile_to_taxtable = qiime_tools.mergeOTU_UC:UC_to_taxtable

    ''',
    test_requirements = ['pytest>=2.1'],
    tests_require=['pytest'],
    cmdclass = {'test': PyTest},

)
