from setuptools import setup

setup(name='GeCCo',
      version='0.1.0',
      author=['Sergio Garcia', 'Caleb Walker'],
      packages=['geccoco'],
      entry_points={
          'console_scripts': [
              'geccoco = geccoco.__main__:cli'
          ]
      },
      )
