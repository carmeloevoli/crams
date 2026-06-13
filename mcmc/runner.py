"""Run the crams binary and read back rigidity spectra."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import numpy as np

# Maps column index (1-based) to element symbol — matches output.cpp order
_ELEMENTS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
]
CHARGE: dict[str, int] = {sym: i for i, sym in enumerate(_ELEMENTS) if sym}

# Fragmentation cross-section models accepted by crams (see src/core/input.cpp).
FRAGMENTATION_MODELS = [
    "fluka4dragon",
    "usine_galprop17_opt12",
    "usine_galprop17_opt22",
    "usine_webber03_coste12",
]


def to_ini_params(params: dict[str, float]) -> dict[str, float]:
    """Convert fit-space keys to the keys crams expects (returns a new dict).

    The fit samples reparametrised quantities; crams needs the physical ones:
      - ``d0_h`` (= d0/h)   -> ``d0 = d0_h * h``
      - ``rb_log`` (log10R) -> ``rb = 10**rb_log``
    """
    p = dict(params)
    if "d0_h" in p:
        try:
            h = p["h"]
        except KeyError:
            raise KeyError("'d0_h' parametrization requires 'h' in params")
        p["d0"] = p.pop("d0_h") * h
    if "rb_log" in p:
        p["rb"] = 10.0 ** p.pop("rb_log")
    return p


class CramsRunner:
    """Thin wrapper around the crams binary.

    Parameters
    ----------
    build_dir:
        Directory containing the crams binary and its data/ sub-directory.
        Defaults to ``<repo_root>/build``.
    binary:
        Name of the crams executable within build_dir.
    timeout:
        Maximum wall-clock seconds to wait for a single run.
    """

    def __init__(
        self,
        build_dir: Path | str | None = None,
        binary: str = "crams",
        timeout: int = 120,
        quiet: bool = True,
        inelastic_model: str = "tripathi99",
        fragmentation_model: str | None = None,
        read_isotopes: bool = False,
    ) -> None:
        if build_dir is None:
            build_dir = Path(__file__).parent.parent / "build"
        self.build_dir = Path(build_dir).resolve()
        self.binary_path = self.build_dir / binary
        self.timeout = timeout
        self.quiet = quiet
        self.inelastic_model = inelastic_model
        # None -> let crams use its built-in default fragmentation model.
        # Accepted: fluka4dragon, usine_galprop17_opt12, usine_galprop17_opt22,
        #           usine_webber03_coste12
        self.fragmentation_model = fragmentation_model
        self.read_isotopes = read_isotopes
        self._counter = 0

        if not self.binary_path.exists():
            raise FileNotFoundError(f"crams binary not found: {self.binary_path}")
        (self.build_dir / "output").mkdir(exist_ok=True)

    # ------------------------------------------------------------------
    def run(self, params: dict[str, float]) -> dict[str, np.ndarray] | None:
        """Run crams with *params* and return interpolated spectra.

        Parameters
        ----------
        params:
            Mapping of .ini key → value (in the units expected by crams, e.g.
            ``d0`` in units of 1e28 cm²/s, ``delta`` dimensionless, etc.).
            As special cases, the key ``d0_h`` (= d0/h) is converted to
            ``d0 = d0_h * h`` and ``rb_log`` to ``rb = 10**rb_log`` before the
            .ini is written.

        Returns
        -------
        dict with keys ``'R'`` (rigidity in GV) and element symbols ``'H'``,
        ``'He'``, …, ``'Ni'`` (flux in 1/(GeV m² s sr)), or ``None`` if the
        run failed.
        """
        self._counter += 1
        tag = f"_mcmc_{os.getpid()}_{self._counter}"
        ini_name = f"{tag}.ini"
        ini_path = self.build_dir / ini_name
        output_file = self.build_dir / "output" / f"{tag}_spectra_R_0.txt"
        isotope_file = self.build_dir / "output" / f"{tag}_isotopes_R_0.txt"

        # Reconstruct the physical crams keys (d0 from d0/h and the fixed halo
        # height; rb from log10 rb) just before writing the .ini.
        params = to_ini_params(params)

        try:
            with open(ini_path, "w") as f:
                for key, value in params.items():
                    f.write(f"{key} {value:.6e}\n")
                f.write(f"inelastic_model {self.inelastic_model}\n")
                if self.fragmentation_model is not None:
                    f.write(f"fragmentation_model {self.fragmentation_model}\n")
                f.write("id 0\n")

            cmd = [str(self.binary_path), ini_name]
            if self.quiet:
                cmd.append("-q")
            result = subprocess.run(
                cmd,
                cwd=self.build_dir,
                capture_output=True,
                timeout=self.timeout,
            )
            if result.returncode != 0 or not output_file.exists():
                return None

            data = np.loadtxt(output_file, comments="#")
            spectra: dict[str, np.ndarray] = {"R": data[:, 0]}
            for Z, sym in enumerate(_ELEMENTS[1:], start=1):
                spectra[sym] = data[:, Z]

            # Isotope-resolved Be (columns: R, Be9, Be10) for the Be10/Be9 ratio.
            if self.read_isotopes and isotope_file.exists():
                iso = np.loadtxt(isotope_file, comments="#")
                spectra["Be9"] = iso[:, 1]
                spectra["Be10"] = iso[:, 2]
            return spectra

        except Exception:
            return None

        finally:
            for path in (ini_path, output_file, isotope_file):
                if path.exists():
                    path.unlink(missing_ok=True)
