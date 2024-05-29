from setuptools import setup, find_packages

setup(
    name='pgc-finder',
    version='0.0.1',
    description='python package for probe-based gene cluster finding in large microbial genome database',
    author='jsrim',
    author_email='comfortindex@naver.com',
    url='https://github.com/jrim42/gene-cluster',
    install_requires=['pandas', 'matplotlib'],
    packages=find_packages(exclude=[]),
    keywords=[],
    python_requires='>=3.11',
    package_data={},
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
