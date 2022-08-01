from setuptools import setup, find_packages

setup(name='TIES_MD',
      version='0.1',
      description='TIES protocol',
      url='https://github.com/UCL-CCS/TIES_MD',
      platforms=['linux64',
                 'ppc64le'],
      author='Alexander David Wade',
      author_email='None',
      license='LGPL',
      packages=find_packages(),
      include_package_data = True,
      entry_points = {'console_scripts':['ties_md = TIES_MD.cli:main', 'ties_ana = TIES_MD.ties_analysis.ties_analysis:main']})
