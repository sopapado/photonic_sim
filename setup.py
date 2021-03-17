from setuptools import setup

setup(name='photonic_sim',
      version='1.0',
      description='A collection of methods simulating propagation of optical waves in dielectric or lossy materials with 1D geometry',
      url='http://github.com/sopapado/photonic_sim',
      author='Sotirios Papadopoulos',
      author_email='send@git.com',
      license='MIT',
      packages=['photonic_sim'],
      install_requires = [ 'tabulate','tqdm','addcopyfighandler','mpmath','pickle-mixin'],
      zip_safe=False)
      


