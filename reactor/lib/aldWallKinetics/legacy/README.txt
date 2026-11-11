Original boundary-condition sources before the wall kinetics moved into lib/aldkin (work package 8).

patternedSurfaceFvPatchScalarField.{H,C}.orig  the multiscale wall: flat Langmuir surface plus a quasi-steady
                                               Knudsen/Langmuir hole per face, open-area blending, ARProfile and
                                               holeDensityProfile tables
langmuirStickingFvPatchScalarField.{H,C}.orig  the flat Langmuir wall with optional recombination (r) and
                                               desorption (kdes)

The files carry the suffix .orig so that wmake's lnInclude step does not link them next to the live sources of the
same name; they are not compiled. The live sources in the parent directory are thin adapters around libaldkin
(lib/aldkin/src/aldkin.cpp holds the transcribed arithmetic, expression by expression) and reproduce these files'
results bit for bit (tests/test_regression_reference.py; paper/wp/wp8_library.md).
