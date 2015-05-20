from setuptools import setup

setup(
    name='qiime-tools',
    version='0.1.13',
    py_modules=['fastq_concat',
                'parallel_split_libraries_fastq'],
    install_requires=[
        'Click',
        'Biopython'
    ],
    entry_points='''
        [console_scripts]
        fastqconcat=fastq_concat:fastqconcat
        parallelconcat=fastq_concat:parallel_concat
        parallel_split_library_fastq=parallel_split_libraries_fastq:parallel_splitlibraries_fastq
        parallel_split_libraries=parallel_split_library:parallel_split_library
        UCfiles_to_taxtable=UCfiles_to_taxtable:UCfiles_to_taxtable
        taxtable_to_otutable=UCfiles_to_taxtable:taxtable_to_otutable
    ''',
)
