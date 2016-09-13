from distutils.core import setup, Extension, Command

class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        #import here, outside the eggs aren't loaded
        import pytest
        pytest.main()

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
  cmdclass = {'test': PyTest}
)
