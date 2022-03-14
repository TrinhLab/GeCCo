from setuptools import setup

setup(name='GeCCo',
      version='0.1.0',
      author=['Sergio Garcia', 'Caleb Walker'],
      packages=['gecco'],
      entry_points={
          'console_scripts': [
              'gecco = gecco.__main__:cli'
          ]
      },
      )
