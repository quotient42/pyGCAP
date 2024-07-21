from setuptools import setup, find_packages

setup(
    name='pygcap',
    version='0.0.7',
    description='Python package for probe-based gene cluster finding in large microbial genome database',
    author='jsrim',
    author_email='comfortindex@naver.com',
    url='https://github.com/jrim42/pyGCAP',
    install_requires=['pandas', 'matplotlib', 'requests'],
    packages=find_packages(),
    package_data={
        'pygcap': ['data/*.tsv'],
    },
    keywords=['gene', 'cluster', 'genomics', 'bioinformatics'],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'pygcap=pygcap.cli:main',
        ],
    },
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)
