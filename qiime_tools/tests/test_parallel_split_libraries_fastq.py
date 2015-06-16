# content of test_sample.py
"""
This file test the parallel_split_libraries_fastq script.
"""


from subprocess import call, check_call

def test_coreutils():
    """
    Checking that coreutils is on your system. Note: You must have coreutils >= 8.4

    Since we use the CLI of Split and require it to work as expected we require that 8.4 or greater is used.
    this will cause issues on MacOSX. This will fail if

    """
    check_call(['split', '--h'])
    check_call(['cat', '--h'])
    check_call(['zcat', '--h'])

def test_seqtk():
    """
    Checking that coreutils is on your system. Note: You must have coreutils >= 8.4

    Since we use the CLI of Split and require it to work as expected we require that 8.4 or greater is used.
    this will cause issues on MacOSX. This will fail if

    """
    assert(call('seqtk', shell=True) == 1)
