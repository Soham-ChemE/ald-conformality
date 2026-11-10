/*---------------------------------------------------------------------------*\
    langmuirSticking: self-limiting (ALD) wall for a passive scalar T = c/c0.

    Kinetic uptake flux  j = (vbar/4) c_face [ s (1 - theta) + r ]   [m/s * c]
    balanced against diffusion  D (c_cell - c_face)/delta,  i.e. a mixed BC
    with refValue 0, refGrad 0 and valueFraction  f = K/(K + D/delta),
    K = (vbar/4)[s(1-theta) + r].

    Coverage per face advances once per time step:
        dtheta/dt = (vbar/4) s (1 - theta) c0 T_face / N0
    s   sticking coefficient on an open site      [-]
    r   loss probability on every collision (radical recombination; 0 for thermal ALD) [-]
    vbar mean thermal speed of the precursor      [m/s]
    N0  adsorption sites filled by one saturated half-cycle [1/m2]
    c0  precursor number density at T = 1          [1/m3]
    D   read from constant/transportProperties (DT) unless given as "D".
\*---------------------------------------------------------------------------*/
#ifndef langmuirStickingFvPatchScalarField_H
#define langmuirStickingFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{

class langmuirStickingFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    scalar s_;
    scalar r_;
    scalar vbar_;
    scalar N0_;
    scalar c0_;
    scalar D_;
    scalarField theta_;
    label curTimeIndex_;

public:

    TypeName("langmuirSticking");

    langmuirStickingFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    langmuirStickingFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    langmuirStickingFvPatchScalarField
    (
        const langmuirStickingFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    langmuirStickingFvPatchScalarField
    (
        const langmuirStickingFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new langmuirStickingFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new langmuirStickingFvPatchScalarField(*this, iF)
        );
    }

    const scalarField& theta() const { return theta_; }

    virtual void autoMap(const fvPatchFieldMapper&);
    virtual void rmap(const fvPatchScalarField&, const labelList&);
    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
