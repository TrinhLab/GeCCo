from setuptools import setup

setup(name='GEC',
      version='0.1.0',
      author=['Sergio Garcia', 'Caleb Walker'],
      packages=['gec'],
      entry_points={
          'console_scripts': [
              'gec = gec.__main__:cli'
          ]
      },
      )