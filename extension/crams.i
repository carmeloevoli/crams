%module crams
%include <std_common.i>
%include <std_vector.i>
%include <std_string.i>

%{
#include "crams/runner.h"
%}

%template(RigiditySpectraRow) std::vector<double>;
%template(RigiditySpectra) std::vector<std::vector<double>>;

%include "crams/utils/git_revision.h"

%include "crams/core/pid.h"
%include "crams/particle.h"
%include "crams/core/input.h"
%include "crams/core/output.h"
%include "crams/particlelist.h"
%include "crams/runner.h"

