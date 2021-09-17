

import os, sys


from distutils.core import setup, Extension
from distutils import sysconfig
import pybind11
cpp_args = ['-std=c++11']

sfc_module = Extension(
    'VSM', sources = ['module.cpp'],
    include_dirs=['pybind11/include'],
    language='c++',
    extra_compile_args = cpp_args,
    )

setup(
    name = 'VSM-LATI',
    version = '1.0',
    description = 'Códigos VSM - Desenvolvido por Marcus',
    ext_modules = [sfc_module],
)

