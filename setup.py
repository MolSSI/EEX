import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name='EEX',
        version="alpha",
        description='The Energy Expression Exchange.',
        author='MolSSI',
        author_email='dgasmith@vt.edu',
        url="https://github.com/molssi/eex",
        license='BSD-3C',
        packages=setuptools.find_packages(),
        install_requires=[
            'numpy>=1.7',
            'pandas>=0.18',
            'tables>=3.2',
            'pint>=0.8',
            'numexpr>=2.4',
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',  # autodoc was broken in 1.3.1
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },

        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],

        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True,
    )
