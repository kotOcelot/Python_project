from setuptools import setup, find_packages

setup(
    name = 'rings',
    version = '0.0.1',
    author = 'Ekaterina Podzolkova',
    author_email = 'podzolkova.en14@physics.msu.ru',
    description = 'Modeling of Galactic resonance rings',
    license = 'MIT',
    packages = find_packages(),
    test_suite = 'test',
    entry_points = {'console_scripts': ['rings = rings.__main__:main']},
    classifiers = [
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Topic :: Education',
        'Programming Language :: Python :: 3.7',],
    keywords = 'sample education',
)

