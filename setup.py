from setuptools import setup

setup(name='GEC',
      version='0.1.0',
      packages=['src'],
      entry_points={
          'console_scripts': [
              'gec = src.__main__:main'
          ]
      },
      )