#!/usr/bin/env python

from setuptools import setup


from setuptools import find_packages, setup


def get_requirements():
    """Collects requirements"""
    with open("requirements.txt", "r") as fp:
        requirements = [line.strip() for line in fp if line.strip()]
        print(requirements)
        return requirements

def main():
    setup(name='AnimtoD3plot',
        version='1.0',
        description='AnimToD3plot conversion package',
        license="Mozilla Public License Version 2.0",
        packages=['animtod3plot'],
        install_requires=get_requirements(),)


if __name__ == "__main__":
    main()
