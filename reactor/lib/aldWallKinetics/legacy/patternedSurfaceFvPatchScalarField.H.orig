/*---------------------------------------------------------------------------*\
    patternedSurface: concurrent multiscale ALD wall for a passive scalar T = c/c0.

    Every patch face carries (a) a flat Langmuir surface with coverage thetaFlat and
    (b) a 1D Knudsen-diffusion / Langmuir model of a cylindrical hole (diameter d,
    aspect ratio AR, nCells cells) whose mouth sees the face gas density c0*T_face.
    holeDensity holes per m2 of wafer; open-area fraction phi = holeDensity*pi*d^2/4.

    Net uptake per unit wafer area:  j = K_eff c_face,
        K_eff = (1-phi) (vbar/4)[s(1-thetaFlat) + r]  +  phi * j_hole/c_face,
    where j_hole is the Knudsen flux into the hole mouth from the quasi-steady hole
    gas solution. Implemented as a mixed BC: refValue 0, refGrad 0,
    valueFraction = K_eff/(K_eff + D/delta).

    Spatial variation: AR and holeDensity are the wafer-wide defaults. Optional
    entries ARProfile and holeDensityProfile, each a table of (radius_m value)
    pairs, e.g.  ARProfile table ((0 100) (0.075 100) (0.1 60) (0.15 60));
    override them face by face, linearly interpolated at the face-centre radius
    sqrt(y^2 + z^2) (the wedge axis is x) and clamped outside the table range.
    A face with holeDensity 0 (or AR 0) is a flat Langmuir surface: no hole is
    solved there, its hole coverages stay 0 and phi = 0.

    Coverages advance once per time step. Hole gas is quasi-steady (Ylilammi 2018).
    Writes thetaFlat, thetaBottom (hole bottom), thetaMid (half depth), Keff,
    ARField and holeDensityField (per face), and thetaHole (all cells, for restart).
\*---------------------------------------------------------------------------*/
#ifndef patternedSurfaceFvPatchScalarField_H
#define patternedSurfaceFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Tuple2.H"

namespace Foam
{

class patternedSurfaceFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    typedef List<Tuple2<scalar, scalar>> profileTable;

    scalar s_, r_, vbar_, N0_, c0_, D_;
    scalar d_, AR_, holeDensity_;
    label nCells_;
    profileTable ARProfile_;            // empty when AR is uniform
    profileTable holeDensityProfile_;   // empty when holeDensity is uniform
    scalarField ARField_;               // per-face aspect ratio
    scalarField holeDensityField_;      // per-face holes per m2
    scalarField thetaFlat_;
    scalarField thetaHole_;      // nFaces * nCells, cell index fastest
    scalarField Keff_;
    label curTimeIndex_;

    // read an optional "name table ((r v) ...)" entry; returns an empty list when absent
    static profileTable readProfile(const dictionary& dict, const word& name);

    // piecewise-linear interpolation of a profile at radius r, clamped at the ends
    static scalar interpolateProfile(const profileTable& tab, const scalar r);

    // fill ARField_ and holeDensityField_ from the scalars and the profiles at the face-centre radii
    void setPatternFields();

    // true when the face carries holes that are solved
    bool hasHoles(const label facei) const
    {
        return holeDensityField_[facei] > 0 && ARField_[facei] > 0;
    }

    // quasi-steady hole gas for one face: returns mouth flux / (c0 T_face) and fills n (relative to mouth density = 1)
    scalar solveHole(const label facei, scalarField& nrel) const;

public:

    TypeName("patternedSurface");

    patternedSurfaceFvPatchScalarField(const fvPatch&, const DimensionedField<scalar, volMesh>&);
    patternedSurfaceFvPatchScalarField(const fvPatch&, const DimensionedField<scalar, volMesh>&, const dictionary&);
    patternedSurfaceFvPatchScalarField(const patternedSurfaceFvPatchScalarField&, const fvPatch&,
                                       const DimensionedField<scalar, volMesh>&, const fvPatchFieldMapper&);
    patternedSurfaceFvPatchScalarField(const patternedSurfaceFvPatchScalarField&, const DimensionedField<scalar, volMesh>&);

    virtual tmp<fvPatchScalarField> clone() const
    { return tmp<fvPatchScalarField>(new patternedSurfaceFvPatchScalarField(*this)); }
    virtual tmp<fvPatchScalarField> clone(const DimensionedField<scalar, volMesh>& iF) const
    { return tmp<fvPatchScalarField>(new patternedSurfaceFvPatchScalarField(*this, iF)); }

    virtual void autoMap(const fvPatchFieldMapper&);
    virtual void rmap(const fvPatchScalarField&, const labelList&);
    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
