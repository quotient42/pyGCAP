from setuptools import setup, find_packages

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()
setup(
    name='pygcap',
    version='1.2.5',
    description='Python package for probe-based gene cluster finding in large microbial genome database',
	long_description=long_description,
    long_description_content_type='text/markdown',
    author='jsrim',
    author_email='comfortindex@naver.com',
    url='https://github.com/jrim42/pyGCAP',
    install_requires=[
        'pandas >= 2.1', 
        'matplotlib >= 3.7, < 3.9',  
        'requests', 
        'setuptools'
        ],
    packages=find_packages(),
    package_data={
        'pygcap': ['data/*.tsv'],
    },
    keywords=['gene', 
              'cluster', 
              'genomics', 
              'bioinformatics'],
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
