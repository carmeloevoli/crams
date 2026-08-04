# TODO: proper test suite for Python wrapper

import shutil
import subprocess
from pathlib import Path

from crams import main

TESTS_DIR = Path(__file__).parent
ROOT_DIR = TESTS_DIR.parent.resolve()
crams_ini = ROOT_DIR / "examples/crams.ini"
BUILD_DIR = ROOT_DIR / "build"
TEST_OUTPUT_DIR = Path.cwd() / "output"
TEST_OUTPUT_DIR.mkdir(exist_ok=True)

# 1. Running CRAMS for benchmark output
subprocess.run(
    [BUILD_DIR / "crams", crams_ini, "-q"],
    check=True,
    cwd=BUILD_DIR,
)
reference_output = TEST_OUTPUT_DIR / "reference.txt"
shutil.move(BUILD_DIR / "output/crams_spectra_R_7.txt", reference_output)

# 2. Running Python CRAMS wrapper with both native and Python side parsing
main(
    ini_file=crams_ini,
    native_parsing=True,
    quiet=False,
)
native_output = TEST_OUTPUT_DIR / "native_parsing.txt"
shutil.move(TEST_OUTPUT_DIR / "test_spectra_R_7.txt", native_output)

main(
    ini_file=crams_ini,
    native_parsing=False,
    quiet=False,
)
py_side_output = TEST_OUTPUT_DIR / "python_side_parsing.txt"
shutil.move(TEST_OUTPUT_DIR / "test_spectra_R_0.txt", py_side_output)

assert native_output.read_text() == reference_output.read_text(), "Native parsing output differs from the reference"
assert py_side_output.read_text() == reference_output.read_text(), (
    "Python-side parsing output differs from the reference"
)

shutil.rmtree(TEST_OUTPUT_DIR)
