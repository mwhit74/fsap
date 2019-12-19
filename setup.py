from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='fsap',
    version='0.1',
    description='Free structural analysis package',
    long_description=readme(),
    license='MIT',
    author='mlw',
    url='https://github.com/mwhit74/fsap',
    packages=['fsap.anen', 'fsap.load', 'fsap.input', 'fsap.outp',
              'fsap.refs', 'fsap.sepr', 'fsap.stif', 'fsap.test',
              'fsap.utils']
)
