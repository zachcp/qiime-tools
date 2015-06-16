from setuptools import setup, find_packages

setup(
    name='qiime-tools',
    version='0.2.9',
    install_requires=[
        'Click >= 0.4.0',
        'Biopython >=1.6.5',
        'qiime >=1.9.0'
    ],
    packages=find_packages(),
    entry_points='''
        [console_scripts]
        fastqconcat = qiime_tools.fastq_concat:fastqconcat
        parallelconcat = qiime_tools.fastq_concat:parallel_concat
        parallel_split_library_fastq = qiime_tools.parallel_split_libraries_fastq:parallel_splitlibraries_fastq
        parallel_split_libraries = qiime_tools.parallel_split_library:parallel_split_library
        UCfiles_to_taxtable = qiime_tools.UCfiles_to_taxtable:UCfiles_to_taxtable
        taxtable_to_otutable = qiime_tools.UCfiles_to_taxtable:taxtable_to_otutable
    ''',
    test_requirements = ['pytest>=2.1']

)
