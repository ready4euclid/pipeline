// SWIG Interface to Stat

%module cblStat

%{
#include "Chain.h"
#include "Prior.h"
#include "Parameter.h"
#include "Model.h"
#include "Data.h"
#include "Data1D.h"
#include "Data2D.h"
#include "Chi2.h"
%}

%include "Chain.h"
%include "Prior.h"
%include "Parameter.h"
%include "Model.h"
%include "Data.h"
%include "Data1D.h"
%include "Data2D.h"
%include "Chi2.h"

