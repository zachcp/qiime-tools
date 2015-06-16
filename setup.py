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
    version='0.2.10',
    install_requires=[
        'Click >= 0.4.0',
        'Biopython >=1.6.5',
        'qiime >=1.9.0'
    ],
    packages=find_packages(),
    entry_points='''
        [console_scripts]
        fastqconcat = qiime_tools.fastq_concat:fastqconcat
        parallel_split_library_fastq = qiime_tools.parallel_split_libraries_fastq:parallel_splitlibraries_fastq
        parallel_split_libraries = qiime_tools.parallel_split_library:parallel_split_library
        UCfiles_to_taxtable = qiime_tools.UCfiles_to_taxtable:UCfiles_to_taxtable
        taxtable_to_otutable = qiime_tools.UCfiles_to_taxtable:taxtable_to_otutable
    ''',
    test_requirements = ['pytest>=2.1'],
    tests_require=['pytest'],
    cmdclass = {'test': PyTest},

)
