# Simple simulation of Langmuir turbulence with a buoyant tracer

The script `simple_buoyant_tracer.jl` implements a simple model for Langmuir turbulence
with a `BuoyancyTracer` as well as a second passive tracer that rises. We make sure to 
implement an `ImpenetrableBoundaryCondition` on the top for the rising velocity so that
total tracer is conserved.

