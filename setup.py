from setuptools import setup, find_packages

setup(name='TIES_MD',
      version='0.1',
      description='TIES protocol',
      url='https://github.com/UCL-CCS/TIES_MD',
      platforms=['linux64',
                 'ppc64le'],
      author='adw62',
      author_email='None',
      license='None',
      packages=find_packages(),
      include_package_data = True,
      entry_points = {'console_scripts':['TIES_MD = TIES_MD.cli:main']})
