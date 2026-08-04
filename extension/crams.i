%module crams
%include "std_string.i"

%{
#include "crams/runner.h"
%}

%include "include/crams/core/pid.h"
%include "crams/core/input.h"
%include "crams/particle.h"
%include "crams/particlelist.h"
%include "crams/runner.h"

%include "crams/utils/git_revision.h"
