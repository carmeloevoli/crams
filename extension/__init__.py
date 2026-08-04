import argparse
import pprint
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from crams.crams import (
    FluxSolver_CrankNicolson,
    Input,
    ParticleList,
    Runner,
    get_version,
    git_sha1,
    parseFluxSolver,
    parseFragmentationModel,
    parseInelasticModel,
)

ABUNDANCE_KEYS = (
    "qh",
    "qhe",
    "qli",
    "qbe",
    "qb",
    "qc",
    "qn",
    "qo",
    "qf",
    "qne",
    "qna",
    "qmg",
    "qal",
    "qsi",
    "qp",
    "qs",
    "qcl",
    "qar",
    "qk",
    "qca",
    "qsc",
    "qti",
    "qv",
    "qcr",
    "qmn",
    "qfe",
    "qco",
    "qni",
)


@dataclass
class PropagationParams:
    H_kpc: float
    v_A_km_sec: float
    R_b_GV: float
    delta: float
    ddelta: float
    D_0_cm2_sec: float
    X_src: float
    phi: float

    @staticmethod
    def from_ini_params(params: dict[str, float]):
        return PropagationParams(
            H_kpc=params["h"],
            v_A_km_sec=params["va"],
            R_b_GV=params["rb"],
            delta=params["delta"],
            ddelta=params["ddelta"],
            D_0_cm2_sec=params["d0"] * 1e28,
            X_src=params["xs"],
            phi=params["phi"],
        )

    def to_input(self) -> Input:
        return Input(
            H_kpc=self.H_kpc,
            v_A_km_sec=self.v_A_km_sec,
            R_b_GV=self.R_b_GV,
            delta=self.delta,
            ddelta=self.ddelta,
            D_0_cm2_sec=self.D_0_cm2_sec,
            X_s=self.X_src,
            modulationPotential=self.phi,
        )


@dataclass
class InjectionParams:
    abundances: Sequence[float]
    slopes: Sequence[float]

    @staticmethod
    def from_ini_params(params: dict[str, float]):
        return InjectionParams(
            abundances=[params[key] for key in ABUNDANCE_KEYS],
            slopes=[params["hslope"], params["heslope"], params["slope"]],
        )


class CramsRunner:
    def __init__(
        self,
        inelastic_model: str,
        fragmentation_model: str,
        verbose: bool = False,
        file_output: bool = False,
        _preloaded_injection: ParticleList | None = None,  # used mainly for testing
    ):
        self._runner = Runner(
            inelasticModel=parseInelasticModel(inelastic_model),
            fragmentationModel=parseFragmentationModel(fragmentation_model),
            injection=_preloaded_injection if _preloaded_injection is not None else ParticleList(),
        )
        self._verbose = verbose
        self._file_output = file_output

    def compute(
        self,
        propagation: PropagationParams | Input,
        injection: InjectionParams | None,  # None = use stored injection; mainly for tests
    ) -> np.ndarray:
        if self._file_output:
            Path("output").mkdir(exist_ok=True)  # hard-coded CRAMS output path

        if injection is not None:
            self._runner.setInjectionParams(abundances=injection.abundances, slopes=injection.slopes)
        R_spectra = self._runner.compute(
            propagation.to_input() if isinstance(propagation, PropagationParams) else propagation,
            dumpToFile=self._file_output,
            verbose=self._verbose,
            ignoreInputInitParams=True,
        )
        return np.array(R_spectra)


def cli():
    """Analog of CRAMS CLI, for testing and cross-validation purposes"""

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="CRAMS input file; see examples/crams.ini for details")
    parser.add_argument(
        "--native-parsing",
        action="store_true",
        help="If set, the CLI will use CRAMS native file parsing; useful for testing and cross-validation",
    )
    parser.add_argument("--quiet", "-q", action="store_true", help="Avoid verbose output")
    args = parser.parse_args()
    main(args.input_file, args.native_parsing, args.quiet)


def main(ini_file: str | Path, native_parsing: bool, quiet: bool):
    ini_file = Path(ini_file)
    params: dict[str, float] = {}
    params_str: dict[str, str] = {}
    for line in ini_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, value = line.split()
        key = key.lower().replace("_", "")
        try:
            params[key] = float(value)
        except ValueError:
            params_str[key] = value

    print(f"Running CRAMS v{get_version()} via Python steering")
    print(f"SHA1: {git_sha1()}")

    print(f"Params from {ini_file}:")
    pprint.pprint(params, indent=4, sort_dicts=False)

    if parseFluxSolver(params_str["solver"]) != FluxSolver_CrankNicolson:
        raise ValueError("Non-default solvers are not supported with Python steering")

    if native_parsing:
        _preloaded_injection = ParticleList()
        _preloaded_injection.readParamsFromFile(str(ini_file))
    else:
        _preloaded_injection = None

    runner = CramsRunner(
        inelastic_model=params_str["inelasticmodel"],
        fragmentation_model=params_str["fragmentationmodel"],
        verbose=not quiet,
        file_output=True,
        _preloaded_injection=_preloaded_injection,
    )

    if native_parsing:
        propagation = Input()
        propagation.readParamsFromFile(str(ini_file))
        print(propagation.describe())
        injection: InjectionParams | None = None
    else:
        propagation = PropagationParams.from_ini_params(params)
        print(propagation.to_input().describe())
        injection = InjectionParams.from_ini_params(params)

    result = runner.compute(propagation, injection)
    print(len(result))
