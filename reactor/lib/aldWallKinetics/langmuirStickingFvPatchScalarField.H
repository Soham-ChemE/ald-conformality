/*---------------------------------------------------------------------------*\
    langmuirSticking: self-limiting (ALD) wall for a passive scalar T = c/c0.

    Kinetic uptake flux  j = (vbar/4) c_face [ s (1 - theta) + r ]   [m/s * c]
    balanced against diffusion  D (c_cell - c_face)/delta,  i.e. a mixed BC
    with refValue 0, refGrad 0 and valueFraction  f = K/(K + D/delta),
    K = (vbar/4)[s(1-theta) + r].

    Coverage per face advances once per time step:
        dtheta/dt = (vbar/4) s (1 - theta) c0 T_face / N0  - kdes theta
    s   sticking coefficient on an open site      [-]
    kdes first-order desorption rate of the adsorbed layer [1/s]; 0 = irreversible Langmuir
    r   loss probability on every collision (radical recombination; 0 for thermal ALD) [-]
    vbar mean thermal speed of the precursor      [m/s]
    N0  adsorption sites filled by one saturated half-cycle [1/m2]
    c0  precursor number density at T = 1          [1/m3]
    D   read from constant/transportProperties (DT) unless given as "D".

    With desorption the wall balance K c_face - J_des = D (c_cell - c_face)/delta,
    J_des = kdes theta N0 / c0 (in T units), gives refValue = J_des/K.

    The physics lives in lib/aldkin (libaldkin, flatOnly mode); this class is the
    OpenFOAM adapter and is bit-identical to the pre-library implementation in legacy/.
\*---------------------------------------------------------------------------*/
#ifndef langmuirStickingFvPatchScalarField_H
#define langmuirStickingFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "aldkinHandle.H"

namespace Foam
{

class langmuirStickingFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    scalar s_;
    scalar r_;
    scalar kdes_;   // first-order desorption rate of adsorbed species, 1/s (0 = irreversible)
    scalar vbar_;
    scalar N0_;
    scalar c0_;
    scalar D_;
    aldkinHandle kin_;      // the coverage state of every face (libaldkin, flatOnly)
    scalarField theta_;     // mirror of the library's coverage for theta() and write()
    label curTimeIndex_;

    aldkin_params params() const;
    void resetKinetics();   // (re)create the handle for the current patch size with theta_ as its coverage

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
