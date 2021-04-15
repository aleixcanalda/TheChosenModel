#!/usr/bin/env python
from distutils.core import setup
setup(name='TheChosenModel',
        version='1.0',
        description='Description of my project',
        author='Aleix Canalda Baltrons & Maria Díaz Ros',
        author_email=['aleix.canalda01@estudiant.upf.edu', 'maria.diaz07@estudiant.upf.edu'],
        url='https://github.com/aleixcanalda/TheChosenModel',
        install_requires=['biopython>=1.73'],
        license='LICENSE.txt',
        classifiers=[
                "Programming Language :: Python :: 3",
	              "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent"],
        packages=['scripts'],
        scripts=["scripts/TheChosenModel.py", "scripts/IOinterface.py", "scripts/classes.py",
                "scripts/build_complex4.py", "scripts/energies.py"]
        )
