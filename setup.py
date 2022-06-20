from setuptools import find_packages, setup
import os

README_PATH = os.path.join(os.path.dirname(__file__), "README.md")
with open(README_PATH, 'r') as f:
    README = f.read()

setup(
    name='site_motif',
    packages=find_packages(exclude=["tests"]),
    version='0.1.0',
    license='Apache License 2.0',
    description='A graph based method for aligning protein binding sites in sequence-order independent fashion',
    long_description=README,
    long_description_content_type="text/markdown",
    author='Nanome',
    author_email='hello@nanome.ai',
    url='https://github.com/nanome-ai/site-motif',
    platforms="any",
    keywords=['alignment', 'structural-biology', 'chemistry', 'python'],
)
