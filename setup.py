from distutils.core import setup, Extension

setup(
  name='PyRONN',
  packages=['pyronn'],
  version='4.0rc0',
  description='A refactoring/reimplementation of RONN',
  author='Shyam Saladi',
  author_email='saladi@caltech.edu',
  url='https://github.com/smsaladi/pyronn',
  download_url='https://github.com/smsaladi/pyronn/tarball/4.0',
  keywords=['protein', 'disorder', 'sequence', 'bioinformatics'],
  license='Non-commercial Academic Use License',
  ext_modules=[Extension('pyronn.libronn',
                         sources=['pyronn/libRONN.cpp',
                                  'pyronn/callBBF.cpp'])],
  package_data={'pyronn': ['data/c*/model.rec', 'data/c*/pdfs.rec']},
  zip_safe=True,
  include_package_data=True,
  setup_requires=['pytest-runner'],
  tests_require=['pytest', 'pandas']
)
