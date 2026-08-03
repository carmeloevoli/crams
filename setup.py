# run SWIG first with
# > swig -c++ -python -Iinclude extension/crams.i
# TODO: make platform-independent and integrate with CMake

from setuptools import Extension, setup


swig_extension = Extension(
    "crams._crams",
    sources=["extension/crams_wrap.cxx"],
    include_dirs=["include", "external/plog/include"],
    library_dirs=["build", "/opt/homebrew/lib"],
    libraries=["crams_core", "gsl"],
)


setup(
    name="crams",
    version="0.0.1",
    author="Igor Vaiman",
    author_email="igor.vaiman@gssi.it",
    description="Python interface for cosmic-ray propagation solver CRAMS",
    packages=["crams"],
    package_dir={"crams": "extension"},
    ext_modules=[swig_extension],
)
