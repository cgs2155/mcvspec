//
// Code auto-generated by initpackage (XSPEC12 local 
// model package code generator).  Do not edit.
// Package: mcvspec
// Function body: mcvspecFunctionMap.cxx

#include    "mcvspecFunctionMap.h"

#include    <XSFunctions/Utilities/XSModelFunction.h>

void 
createmcvspecFunctionMap()
{


	XSFunctionMap["mcvspec"]  = new XSCall<xsf77Call>(mcvspec_);

}