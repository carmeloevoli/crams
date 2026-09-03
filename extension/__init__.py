"""
Python interface for CRAMS code.

This module is a pythonic wrapper over a SWIG-generated interface.
It provides high-level access and type hints, as well as cleanly separated
parametrization of injection (per-element abundance and slope, common R-scaled
features) and propagation parameters.
"""

import abc
import argparse
import pprint
from collections.abc import MutableSequence, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .crams import (
    FluxSolver_CrankNicolson,
    GeV,
    Input,
    ParticleList,
    Result,
    Runner,
    SourceSpectraLognormalDist,
    SourceSpectrumBreak,
    SourceSpectrumExpCutoff,
    get_version,
    git_sha1,
    parseFluxSolver,
    parseFragmentationModel,
    parseInelasticModel,
)

__version__ = get_version()


@dataclass
class PropagationParams:
    H_kpc: float = 7.0
    v_A_km_sec: float = 4.40940

    # D coefficient params
    D_0_cm2_sec: float = 2.48255e28
    R_b_GV: float = 290.0
    delta: float = 0.56132
    ddelta: float = 0.22
    smoothness: float = 0.1

    X_src: float = -1.0
    phi: float = 4.87754e-01

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
            smoothness=self.smoothness,
            modulationPotential=self.phi,
        )


ELEMENT_NAMES = (
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
)

ABUNDANCE_INI_KEYS = ("q" + element.lower() for element in ELEMENT_NAMES)


class InjectionFeature(abc.ABC):
    @abc.abstractmethod
    def add_to_input(self, input: Input) -> None: ...


@dataclass
class InjectionBreak(InjectionFeature):
    R_GV: float
    delta_slope: float
    omega: float

    def add_to_input(self, input: Input) -> None:
        input.addSourceSpectrumBreak(
            SourceSpectrumBreak(
                R_GV=self.R_GV,
                deltaSlope=self.delta_slope,
                omega=self.omega,
            )
        )


@dataclass
class InjectionLognormalDist(InjectionFeature):
    R_mean_GV: float
    sigma: float

    # PL index in the dependence of CR accelerator luminocity on the maximum energy.
    # when convolving the individual cut-offs with population weight, we have
    # W(Emax) \propto Emax^beta, beta~1 for standard models of SNR acceleration
    beta: float

    is_lower: bool = False

    def add_to_input(self, input: Input) -> None:
        input.addSourceSpectrumLognormal(
            SourceSpectraLognormalDist(
                R_GV=self.R_mean_GV,
                sigma=self.sigma,
                beta=self.beta,
                isLower=self.is_lower,
            )
        )


@dataclass
class InjectionExpCutoff(InjectionFeature):
    R_cut_GV: float
    Delta: float
    is_lower: bool

    def add_to_input(self, input: Input) -> None:
        input.addSourceSpectrumExpCutoff(
            SourceSpectrumExpCutoff(
                R_GV=self.R_cut_GV,
                Delta=self.Delta,
                isLower=self.is_lower,
            )
        )


@dataclass
class InjectionParams:
    abundances: MutableSequence[float]
    slopes: Sequence[float]
    features: list[InjectionFeature] = field(default_factory=list)

    def __post_init__(self) -> None:
        if len(self.abundances) > len(ELEMENT_NAMES):
            raise ValueError(f"Too many abundances specified, expected exactly {len(ELEMENT_NAMES)}")
        if len(self.abundances) < len(ELEMENT_NAMES):
            raise ValueError(f"Too few abundances specified, expected exactly {len(ELEMENT_NAMES)}")
        if len(self.slopes) > len(ELEMENT_NAMES):
            raise ValueError(f"Too many slopes specified, expected at most {len(ELEMENT_NAMES)}")

    def set_abundance(self, Z: int, q: float) -> None:
        self.abundances[Z - 1] = q

    @staticmethod
    def default() -> "InjectionParams":
        return InjectionParams(
            abundances=[
                5.06605e-02,  # H
                2.54369e-02,  # He
                0.0,  # Li
                0.0,  # Be
                0.0,  # B
                3.98879e-03,  # C
                3.36117e-04,  # N
                7.15129e-03,  # O
                0.0,  # F
                1.34031e-03,  # Ne
                0.5e-4,  # Na
                2.38948e-03,  # Mg
                2.7e-4,  # Al
                2.77911e-03,  # Si
                1e-4,  # P
                4.87000e-04,  # S
                0.0,  # Cl
                3e-4,  # Ar
                0.0,  # K
                4e-4,  # Ca
                0.0,  # Sc
                0.0,  # Ti
                0.0,  # V
                2.5e-4,  # Cr
                0.0,  # Mn
                6.80000e-03,  # Fe
                0.0,  # Co
                4e-4,  # Ni
            ],
            slopes=[4.37, 4.30, 4.36],
        )

    @staticmethod
    def from_ini_params(params: dict[str, float]):
        return InjectionParams(
            abundances=[params[key] for key in ABUNDANCE_INI_KEYS],
            slopes=[params["hslope"], params["heslope"], params["slope"]],
        )


@dataclass
class LogGrid:
    min: float  # GeV / GV
    max: float  # GeV / GV
    size: int

    def __post_init__(self) -> None:
        assert self.size >= 2, "Size must be at least 2"
        assert self.max > self.min, "Maximum value of the grid must be greater than the minimum"

    def __str__(self) -> str:
        return f"[{self.min:.1e}, {self.max:.1e}] GeV, {self.size} points"

    def to_numpy(self) -> np.ndarray:
        return np.geomspace(self.min, self.max, self.size)


class CramsError(Exception):
    pass


CRAMS_DEFAULT_INPUT = Input()
CRAMS_DEFAULT_T_SIM_GRID = LogGrid(
    min=CRAMS_DEFAULT_INPUT.TSimMin() / GeV,
    max=CRAMS_DEFAULT_INPUT.TSimMax() / GeV,
    size=CRAMS_DEFAULT_INPUT.TSimSize(),
)
CRAMS_DEFAULT_R_OUT_GRID = LogGrid(
    min=CRAMS_DEFAULT_INPUT.ROutputMin() / GeV,
    max=CRAMS_DEFAULT_INPUT.ROutputMax() / GeV,
    size=CRAMS_DEFAULT_INPUT.ROutputSize(),
)


class CramsRunner:
    def __init__(
        self,
        inelastic_model: str = CRAMS_DEFAULT_INPUT.inelasticModelName(),
        fragmentation_model: str = CRAMS_DEFAULT_INPUT.fragmentationModelName(),
        verbose: bool = False,
        file_output: bool = False,
        T_sim_grid: LogGrid = CRAMS_DEFAULT_T_SIM_GRID,
        R_out_grid: LogGrid = CRAMS_DEFAULT_R_OUT_GRID,
        _preloaded_injection: ParticleList | None = None,  # used mainly for testing
    ):
        self._runner = Runner(
            inelasticModel=parseInelasticModel(inelastic_model),
            fragmentationModel=parseFragmentationModel(fragmentation_model),
            injection=_preloaded_injection if _preloaded_injection is not None else ParticleList(),
        )
        self._inelastic_model = inelastic_model
        self._fragmentation_model = fragmentation_model
        self._verbose = verbose
        self._file_output = file_output
        self._T_sim_grid = T_sim_grid
        self._R_out_grid = R_out_grid

    def compute(
        self,
        propagation: PropagationParams | Input,
        injection: InjectionParams | None,  # None = use stored injection; mainly for tests
    ) -> np.ndarray:
        """
        Main computation method. Returns table as a table of (n_points, n_elements + 1),
        the first column gives rigidities, the last n_elements --- elemental fluxes, summed
        over izotopes. The rigidity is in GV, the spectra are in 1 / GeV m^2 sec
        """
        if self._file_output:
            Path("output").mkdir(exist_ok=True)  # hard-coded CRAMS output path

        input = propagation.to_input() if isinstance(propagation, PropagationParams) else propagation
        input.setTSim(self._T_sim_grid.min, self._T_sim_grid.max, self._T_sim_grid.size)
        input.setROutput(self._R_out_grid.min, self._R_out_grid.max, self._R_out_grid.size)
        if injection is not None:
            # per-element injection params do not live inside the "input" object, but in a ParticleList
            # container inside the runner; here we modify them through the dedicated method
            self._runner.setInjectionParams(abundances=injection.abundances, slopes=injection.slopes)
            for feature in injection.features:
                feature.add_to_input(input)
        result: Result = self._runner.computeSafe(
            input=input,
            dumpToFile=self._file_output,
            verbose=self._verbose,
            # input object contains default values, but we need to use those already configured in the runner object
            ignoreInputInitParams=True,
        )
        if result.is_error:
            raise CramsError(result.error)
        return np.array(result.spectra)


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
    features: list[InjectionFeature] = []
    for line in ini_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        key = tokens[0]
        key = key.lower().replace("_", "")
        if key == "sourcebreak" or key == "sourceerfccutoff":
            if len(tokens) != 4:
                raise ValueError(f"Expected three parameters for {tokens[0]}")
            if key == "sourcebreak":
                t1, t2, t3 = tokens[1:]
                features.append(InjectionBreak(float(t1), float(t2), float(t3)))
            else:
                t1, t2, t3, t4 = tokens[1:]
                features.append(InjectionLognormalDist(float(t1), float(t2), float(t3), bool(t4)))
            continue
        if len(tokens) != 2:
            raise ValueError(f"Expected one value for {tokens[0]}")
        value = tokens[1]
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
        injection.features.extend(features)

    result = runner.compute(propagation, injection)
    print(len(result))
