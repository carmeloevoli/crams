# TODO: make platform-independent and integrate with CMake

from pathlib import Path

from setuptools import Extension, setup

PROJECT_ROOT = Path(__file__).parent.resolve()

swig_extension = Extension(
    "crams._crams",
    sources=["extension/crams_wrap.cxx"],
    include_dirs=[str(PROJECT_ROOT / "include"), str(PROJECT_ROOT / "external/plog/include")],
    library_dirs=[str(PROJECT_ROOT / "build"), "/opt/homebrew/lib"],
    runtime_library_dirs=[str(PROJECT_ROOT / "build")],
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
    install_requires=["numpy~=2.0"],
)
