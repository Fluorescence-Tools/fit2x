#! /usr/bin/env python
import os
import platform
import subprocess
import pathlib

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


def read_version(
        header_file='./include/fits2x.h'
):
    version = "0.0.0"
    with open(header_file, "r") as fp:
        for line in fp.readlines():
            if "#define" in line and "FIT2X_VERSION" in line:
                version = line.split()[-1]
    return version.replace('"', '')


__name__ = "fit2x"
__version__ = read_version()
__license__ = 'Mozilla Public License 2.0 (MPL 2.0)'


class CMakeExtension(Extension):

    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):

    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(
                self.get_ext_fullpath(ext.name)
            )
        ).replace('\\', '/')

        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DCMAKE_SWIG_OUTDIR=' + extdir
        ]
        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg, '-j 8']
        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        if platform.system() == "Windows":
            cmake_args += [
                '-DBUILD_PYTHON_INTERFACE=ON',
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                '-GVisual Studio 14 2015 Win64'
            ]
        else:
            build_args += ['--', '-j8']
            # When using conda try to convince cmake to use
            # the conda boost
            CONDA_PREFIX = os.getenv('CONDA_PREFIX')
            if CONDA_PREFIX is not None:
                print("Conda prefix is: ", CONDA_PREFIX)
                print("Convincing cmake to use the conda boost")
                cmake_args += [
                    '-DCMAKE_PREFIX_PATH=' + CONDA_PREFIX,
                    '-DBOOST_ROOT=' + CONDA_PREFIX,
                    '-DBoost_NO_SYSTEM_PATHS=ON',
                    '-DBoost_DEBUG=OFF',
                    '-DBoost_DETAILED_FAILURE_MESSAGE=ON'
                ]
        env = os.environ.copy()
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        try:
            # build the documentation.i file using doxygen and doxy2swig
            working_directory = pathlib.Path(__file__).parent.absolute()
            subprocess.check_call(
                ["doxygen"],
                cwd=str(working_directory / "docs"),
                env=env
            )
            subprocess.check_call(
                ["python", "doxy2swig.py", "../docs/_build/xml/index.xml", "../ext/python/documentation.i"],
                cwd=str(working_directory / "utility"),
                env=env
            )
        except:
            print("Problem calling doxygen")
        print("cmake building: " + " ".join(cmake_args))
        subprocess.check_call(
            ['cmake', ext.sourcedir] + cmake_args,
            cwd=self.build_temp,
            env=env
        )
        subprocess.check_call(
            ['cmake', '--build', '.'] + build_args,
            cwd=self.build_temp
        )


setup(
    name=__name__,
    version=__version__,
    license=__license__,
    author='Thomas-Otavio Peulen',
    author_email='thomas.otavio.peulen@gmail.com',
    ext_modules=[
        CMakeExtension('fit2x')
    ],
    cmdclass={
        'build_ext': CMakeBuild
    },
    install_requires=[
        'numpy'
    ],
    setup_requires=[
        'setuptools'
    ],
    zip_safe=False,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ]
)
