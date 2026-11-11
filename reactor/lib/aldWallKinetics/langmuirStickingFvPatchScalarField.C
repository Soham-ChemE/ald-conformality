#include "langmuirStickingFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IOdictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * library plumbing * * * * * * * * * * * * * * * //

aldkin_params langmuirStickingFvPatchScalarField::params() const
{
    aldkin_params p;
    aldkin_params_default(&p);
    p.s = s_; p.r = r_; p.kdes = kdes_; p.vbar = vbar_; p.N0 = N0_; p.c0 = c0_;
    p.flatOnly = 1; p.scheme = ALDKIN_SCHEME_EULER;
    return p;
}


void langmuirStickingFvPatchScalarField::resetKinetics()
{
    kin_ = aldkinHandle(params());
    kin_.check(aldkin_set_faces(kin_.get(), theta_.size(), nullptr), "aldkin_set_faces");
    kin_.check(aldkin_set_theta(kin_.get(), ALDKIN_THETA_FLAT, theta_.cdata()), "aldkin_set_theta");
}


// * * * * * * * * * * * * * * * * constructors * * * * * * * * * * * * * * * * //

langmuirStickingFvPatchScalarField::langmuirStickingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    s_(1), r_(0), kdes_(0), vbar_(410), N0_(4e18), c0_(1), D_(-1),
    kin_(),
    theta_(p.size(), 0.0),
    curTimeIndex_(-1)
{
    resetKinetics();
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


langmuirStickingFvPatchScalarField::langmuirStickingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    s_(dict.get<scalar>("s")),
    r_(dict.getOrDefault<scalar>("r", 0.0)),
    kdes_(dict.getOrDefault<scalar>("kdes", 0.0)),
    vbar_(dict.get<scalar>("vbar")),
    N0_(dict.get<scalar>("N0")),
    c0_(dict.get<scalar>("c0")),
    D_(dict.getOrDefault<scalar>("D", -1.0)),
    kin_(),
    theta_(p.size(), 0.0),
    curTimeIndex_(-1)
{
    if (dict.found("theta"))
    {
        theta_ = scalarField("theta", dict, p.size());
    }
    resetKinetics();
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(0.0);
    }
}


langmuirStickingFvPatchScalarField::langmuirStickingFvPatchScalarField
(
    const langmuirStickingFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    s_(ptf.s_), r_(ptf.r_), kdes_(ptf.kdes_), vbar_(ptf.vbar_), N0_(ptf.N0_), c0_(ptf.c0_), D_(ptf.D_),
    kin_(),
    theta_(ptf.theta_, mapper),
    curTimeIndex_(-1)
{
    resetKinetics();
}


langmuirStickingFvPatchScalarField::langmuirStickingFvPatchScalarField
(
    const langmuirStickingFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    s_(ptf.s_), r_(ptf.r_), kdes_(ptf.kdes_), vbar_(ptf.vbar_), N0_(ptf.N0_), c0_(ptf.c0_), D_(ptf.D_),
    kin_(ptf.kin_),
    theta_(ptf.theta_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * * mapping * * * * * * * * * * * * * * * * * //

void langmuirStickingFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
{
    mixedFvPatchScalarField::autoMap(m);
    theta_.autoMap(m);
    resetKinetics();
}


void langmuirStickingFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
    const langmuirStickingFvPatchScalarField& tiptf =
        refCast<const langmuirStickingFvPatchScalarField>(ptf);
    theta_.rmap(tiptf.theta_, addr);
    kin_.check(aldkin_set_theta(kin_.get(), ALDKIN_THETA_FLAT, theta_.cdata()), "aldkin_set_theta");
}


// * * * * * * * * * * * * * * * * * update * * * * * * * * * * * * * * * * * //

void langmuirStickingFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar D = D_;
    if (D < 0)
    {
        const IOdictionary& tp =
            db().lookupObject<IOdictionary>("transportProperties");
        D = tp.get<dimensionedScalar>("DT").value();
    }

    const label ti = db().time().timeIndex();
    const scalarField& Tface = *this;
    if (ti != curTimeIndex_)
    {
        // advance coverage once per time step using the face concentration of the previous step
        // (adsorption on open sites minus first-order desorption; kdes = 0 recovers irreversible Langmuir)
        const scalar dt = db().time().deltaTValue();
        kin_.check(aldkin_step(kin_.get(), Tface.cdata(), dt, nullptr, nullptr), "aldkin_step");
        kin_.check(aldkin_theta(kin_.get(), ALDKIN_THETA_FLAT, theta_.data()), "aldkin_theta");
        curTimeIndex_ = ti;
    }
    else
    {
        kin_.check(aldkin_evaluate(kin_.get(), Tface.cdata(), nullptr, nullptr), "aldkin_evaluate");
    }

    // Wall balance with desorption: K c_face - J_des = D (c_cell - c_face)/delta, J_des = kdes theta N0 / c0
    // -> mixed BC with refValue = J_des/K and valueFraction = K/(K + D/delta)
    const scalarField Dd(D*patch().deltaCoeffs());   // D/delta
    kin_.check
    (
        aldkin_mixed_bc(kin_.get(), Dd.cdata(), refValue().data(), refGrad().data(), valueFraction().data()),
        "aldkin_mixed_bc"
    );

    mixedFvPatchScalarField::updateCoeffs();
}


void langmuirStickingFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("s", s_);
    os.writeEntry("r", r_);
    os.writeEntry("kdes", kdes_);
    os.writeEntry("vbar", vbar_);
    os.writeEntry("N0", N0_);
    os.writeEntry("c0", c0_);
    if (D_ > 0) os.writeEntry("D", D_);
    theta_.writeEntry("theta", os);
    writeEntry("value", os);
}


makePatchTypeField
(
    fvPatchScalarField,
    langmuirStickingFvPatchScalarField
);

} // End namespace Foam
