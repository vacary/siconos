// -*- c++ -*-
// SWIG interface for Siconos Mechanics/joints
%module(package="siconos.mechanics", directors="1", allprotected="1") czm

%include MechanicsBase.i


PY_FULL_REGISTER(BinaryCohesiveNSL, Mechanics); // Abstract

%inline
%{
  // For converting interaction.relations() to known Relations
  SP::BinaryCohesiveNSL cast_BinaryCohesiveNSL(SP::NonSmoothLaw nslaw)
    { return std::dynamic_pointer_cast<BinaryCohesiveNSL>(nslaw); }
%}
