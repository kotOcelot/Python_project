from setuptools import setup, find_packages

setup(
    name = 'rings',
    version = '0.0.1',
    author = 'Ekaterina Podzolkova',
    author_email = 'podzolkova.en14@physics.msu.ru',
    description = 'Modeling of Galactic resonance rings',
    license = 'MIT',
    packages = find_packages(),
    entry_points = {'console_scripts': ['rings = rings.__main__:main']},
    test_suite = 'test',
    install_requires = ['numpy >= 1.16', 'matplotlib >= 3.1', 'scipy >= 1.3', 'astropy >= 3.2'],
    classifiers = [
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Topic :: Education',
        'Programming Language :: Python :: 3.7',],
    keywords = 'sample education',
)

