from setuptools import setup, find_packages

setup(
    name='dcw-finder',
    version='0.0.26',
    description='Package for detecting dcw gene cluster',
    author='jsrim',
    author_email='comfortindex@naver.com',
    url='https://github.com/jrim42/gene-cluster',
    install_requires=['pandas', 'matplotlib'],
    packages=find_packages(exclude=[]),
    keywords=[],
    python_requires='>=3.6',
    package_data={},
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
