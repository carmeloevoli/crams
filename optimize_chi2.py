#!/usr/bin/env python3
"""
Chi-squared minimization script using iminuit.
Reads initial parameters from params.ini, optimizes them to minimize
the sum of chi2 values from chi2_results.txt.
"""

import os
import subprocess
from iminuit import Minuit
import time

# Configuration
CRAMS_DIR = "/Users/bschroer/Desktop/UChicago/crams_fd/build"
CRAMS_EXECUTABLE = os.path.join(CRAMS_DIR, "crams")
INITIAL_PARAMS_FILE = os.path.join(CRAMS_DIR, "params.ini")
CHI2_RESULTS_FILE = os.path.join(CRAMS_DIR, "chi2_results.txt")
WORKING_PARAMS_FILE = os.path.join(CRAMS_DIR, "params_working.ini")

# Parameters to optimize (modify as needed)
PARAMS_TO_OPTIMIZE = {
    "v_A": (4.4, 0.5, 10.0),           # (initial, min, max)
    "D_0": (2.5, 0.5, 5.0),
    "delta": (0.565, 0.3, 1.0),
}

# Chi2 weights (modify to emphasize certain measurements)
# Leave unspecified or set to 1.0 for equal weighting
CHI2_WEIGHTS = {
    "B/C": 4.0,           # More sensitive to diffusion parameter D_0
    "B": 2.0,             # Reference weight
    "He": 1.0,
    "BeB": 2.0,
    "BeC": 2.0,
    "BeO": 2.0,
    "BO": 4.0,
    "CO": 1.0,
    "HeO": 1.0,
    "NeMg": 1.0,
    "SiMg": 1.0,
    # Add other chi2 names as needed
    # Default weight is 1.0 if not specified
}


def read_params(filepath):
    """Read parameters from ini file into a dictionary."""
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                key = parts[0]
                try:
                    value = float(parts[1])
                    params[key] = value
                except ValueError:
                    pass
    return params


def write_params(filepath, params):
    """Write parameters to ini file, preserving format of original."""
    with open(INITIAL_PARAMS_FILE, 'r') as f:
        lines = f.readlines()
    
    with open(filepath, 'w') as f:
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                f.write(line)
                continue
            
            parts = stripped.split()
            if len(parts) >= 1:
                key = parts[0]
                if key in params:
                    # Preserve the original formatting (scientific notation, etc.)
                    comment = ""
                    if '#' in line:
                        comment = " " + line[line.index('#'):]
                    f.write(f"{key} {params[key]:.5e}{comment}\n")
                else:
                    f.write(line)
            else:
                f.write(line)


def run_crams(params_file):
    """Run crams executable with given params file. Optimized for speed."""
    try:
        elapsed = time.time()
        # Stream output to console in real-time, don't capture
        #result = subprocess.run(
        #    [CRAMS_EXECUTABLE, params_file],
        #    cwd=CRAMS_DIR,
        #    timeout=300,  # 5 minute timeout
        #)
        subprocess.run(["/Users/bschroer/Desktop/UChicago/crams_fd/build/crams", "params.ini"])

        elapsed = time.time() - elapsed
        print(elapsed)
        exit()
        if result.returncode != 0:
            print(f"CRAMS execution failed with code {result.returncode}")
            return None
        return True
    except subprocess.TimeoutExpired:
        print("CRAMS execution timed out")
        return None
    except Exception as e:
        print(f"Error running CRAMS: {e}")
        return None


def read_chi2_results(filepath):
    """Read chi2 results and return dictionary of chi2 values by name. Optimized for speed."""
    try:
        chi2_dict = {}
        with open(filepath, 'r') as f:
            # Read all at once instead of line-by-line iteration
            content = f.read()
        
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Split only once on comma
            parts = line.split(',', 1)
            if len(parts) == 2:
                try:
                    name = parts[0].strip()
                    chi2_value = float(parts[1].strip())
                    chi2_dict[name] = chi2_value
                except ValueError:
                    pass
        
        return chi2_dict if chi2_dict else {}
    except FileNotFoundError:
        print(f"Chi2 results file not found: {filepath}")
        return {}
    except Exception as e:
        print(f"Error reading chi2 results: {e}")
        return {}


# Global cache for initial parameters to avoid repeated file reads
_initial_params_cache = None
_param_template_lines = None
_param_line_keys = None  # Cache which lines have which keys
_call_count = 0
_last_print_count = 0
_start_time = None
_call_start_time = None

def _get_initial_params():
    """Get cached initial parameters."""
    global _initial_params_cache
    if _initial_params_cache is None:
        _initial_params_cache = read_params(INITIAL_PARAMS_FILE)
    return _initial_params_cache

def _get_param_template():
    """Cache the template lines of params file for faster writing."""
    global _param_template_lines, _param_line_keys
    if _param_template_lines is None:
        with open(INITIAL_PARAMS_FILE, 'r') as f:
            _param_template_lines = f.readlines()
        # Preprocess: identify which line has which key
        _param_line_keys = {}
        for i, line in enumerate(_param_template_lines):
            stripped = line.strip()
            if stripped and not stripped.startswith('#'):
                parts = stripped.split()
                if len(parts) >= 1:
                    _param_line_keys[parts[0]] = i
    return _param_template_lines, _param_line_keys

def write_params_fast(filepath, params):
    """Write parameters to ini file with minimal overhead."""
    template_lines, line_keys = _get_param_template()
    
    # Only modify the lines that changed - much faster than rewriting everything
    output_lines = template_lines.copy()
    
    for key, value in params.items():
        if key in line_keys:
            i = line_keys[key]
            line = template_lines[i]
            
            # Extract comment if present
            comment = ""
            if '#' in line:
                comment = " " + line[line.index('#'):]
            
            # Replace just this line
            output_lines[i] = f"{key} {value:.5e}{comment}\n"
    
    # Write all at once (faster than line-by-line)
    with open(filepath, 'w') as f:
        f.writelines(output_lines)

def objective_function(**kwargs):
    """
    Objective function to minimize.
    Takes optimized parameters, writes them, runs crams, reads chi2.
    Applies weights to individual chi2 values before summing.
    Aggressively optimized for speed.
    """
    global _call_count, _start_time
    
    _call_count += 1
    
    # Get cached initial parameters and update only changed values
    params = _get_initial_params().copy()
    params.update(kwargs)
    
    # Write to working file (minimal overhead)
    write_params_fast(WORKING_PARAMS_FILE, params)
    
    # Run crams
    if not run_crams(WORKING_PARAMS_FILE):
        return 1e10
    
    # Read chi2 results (now returns dict)
    chi2_dict = read_chi2_results(CHI2_RESULTS_FILE)
    
    if not chi2_dict:
        return 1e10
    
    # Calculate weighted sum - inline for speed
    weighted_chi2_sum = sum(
        chi2_value * CHI2_WEIGHTS.get(name, 1.0)
        for name, chi2_value in chi2_dict.items()
    )
    
    # Print progress every iteration
    param_str = "  ".join([f"{k}={v:.4f}" for k, v in sorted(kwargs.items())])
    print(f"Call {_call_count:4d} | χ²={weighted_chi2_sum:10.4f} | {param_str}", flush=True)
    
    return weighted_chi2_sum


def create_cost_function(param_names):
    """
    Factory function to create a cost function with explicit parameter names.
    This allows iminuit to recognize the parameters.
    """
    # Create function code dynamically
    param_str = ", ".join(param_names)
    func_code = f"""
def cost_func({param_str}):
    kwargs = {{{", ".join([f"'{p}': {p}" for p in param_names])}}}
    return objective_function(**kwargs)
"""
    # Execute the code to create the function
    local_namespace = {"objective_function": objective_function}
    exec(func_code, local_namespace)
    return local_namespace["cost_func"]


def minimize_chi2():
    """
    Main optimization function.
    """
    print("=" * 70)
    print("Chi-squared Minimization using iminuit")
    print("=" * 70)
    
    # Verify files exist
    if not os.path.exists(CRAMS_EXECUTABLE):
        print(f"Error: CRAMS executable not found at {CRAMS_EXECUTABLE}")
        return
    
    if not os.path.exists(INITIAL_PARAMS_FILE):
        print(f"Error: Initial params file not found at {INITIAL_PARAMS_FILE}")
        return
    
    print(f"\nCRAMS executable: {CRAMS_EXECUTABLE}")
    print(f"Initial params file: {INITIAL_PARAMS_FILE}")
    print(f"Working params file: {WORKING_PARAMS_FILE}")
    print(f"Chi2 results file: {CHI2_RESULTS_FILE}")
    
    # Read and cache initial parameters
    initial_params = _get_initial_params()
    print(f"\nInitial parameters read: {len(initial_params)} parameters")
    
    # Setup Minuit
    print(f"\nParameters to optimize: {len(PARAMS_TO_OPTIMIZE)}")
    param_names = list(PARAMS_TO_OPTIMIZE.keys())
    
    for param_name, (init_val, min_val, max_val) in PARAMS_TO_OPTIMIZE.items():
        print(f"  {param_name}: initial={init_val}, range=[{min_val}, {max_val}]")
    
    # Build initial values dictionary for Minuit
    init_values = {name: value for name, (value, _, _) in PARAMS_TO_OPTIMIZE.items()}
    
    # Create a cost function with explicit parameter names
    cost_func = create_cost_function(param_names)
    
    # Create Minuit instance with named parameters
    m = Minuit(cost_func, **init_values)
    
    # Set bounds and errors using parameter names
    for param_name, (init_val, min_val, max_val) in PARAMS_TO_OPTIMIZE.items():
        m.limits[param_name] = (min_val, max_val)
        m.errors[param_name] = (max_val - min_val) / 10  # Initial step size
    
    print("\n" + "=" * 70)
    print("Starting minimization...")
    print("=" * 70)
    
    # Initialize timing
    global _start_time
    _start_time = time.time()
    
    # Run minimization with progress output
    m.migrad()
    
    elapsed = time.time() - _start_time
    
    print(m)
    
    # Print best-fit parameters
    print("\nBest-fit parameters:")
    for param_name in PARAMS_TO_OPTIMIZE.keys():
        value = m.values[param_name]
        error = m.errors[param_name]
        print(f"  {param_name} = {value:.6e} ± {error:.6e}")
    
    # Save best-fit parameters
    best_fit_file = os.path.join(CRAMS_DIR, "params_best_fit.ini")
    best_fit_params = _get_initial_params().copy()
    for param_name in PARAMS_TO_OPTIMIZE.keys():
        best_fit_params[param_name] = m.values[param_name]
    write_params_fast(best_fit_file, best_fit_params)
    print(f"\nBest-fit parameters saved to: {best_fit_file}")
    
    # Print final chi2
    print(f"\nFinal χ² sum = {m.fval:.6f}")
    
    # Create summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Total function calls: {m.nfcn}")
    print(f"Total elapsed time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print(f"Avg time per call: {elapsed/m.nfcn:.2f} seconds")
    print(f"Valid minimum: {m.valid}")
    print(f"Converged: {m.fval < 1e9}")


if __name__ == "__main__":
    minimize_chi2()
