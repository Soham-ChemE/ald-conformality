#include "langmuirStickingFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IOdictionary.H"

namespace Foam
{

langmuirStickingFvPatchScalarField::langmuirStickingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    s_(1), r_(0), vbar_(410), N0_(4e18), c0_(1), D_(-1),
    theta_(p.size(), 0.0),
    curTimeIndex_(-1)
{
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
    vbar_(dict.get<scalar>("vbar")),
    N0_(dict.get<scalar>("N0")),
    c0_(dict.get<scalar>("c0")),
    D_(dict.getOrDefault<scalar>("D", -1.0)),
    theta_(p.size(), 0.0),
    curTimeIndex_(-1)
{
    if (dict.found("theta"))
    {
        theta_ = scalarField("theta", dict, p.size());
    }
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
    s_(ptf.s_), r_(ptf.r_), vbar_(ptf.vbar_), N0_(ptf.N0_), c0_(ptf.c0_), D_(ptf.D_),
    theta_(ptf.theta_, mapper),
    curTimeIndex_(-1)
{}


langmuirStickingFvPatchScalarField::langmuirStickingFvPatchScalarField
(
    const langmuirStickingFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    s_(ptf.s_), r_(ptf.r_), vbar_(ptf.vbar_), N0_(ptf.N0_), c0_(ptf.c0_), D_(ptf.D_),
    theta_(ptf.theta_),
    curTimeIndex_(-1)
{}


void langmuirStickingFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
{
    mixedFvPatchScalarField::autoMap(m);
    theta_.autoMap(m);
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
}


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
    if (ti != curTimeIndex_)
    {
        // advance coverage once per time step using the face concentration of the previous step
        const scalar dt = db().time().deltaTValue();
        const scalarField& Tface = *this;
        theta_ = min
        (
            scalar(1),
            theta_ + dt*(vbar_/4.0)*s_*(1.0 - theta_)*c0_*max(Tface, scalar(0))/N0_
        );
        curTimeIndex_ = ti;
    }

    const scalarField K((vbar_/4.0)*(s_*(1.0 - theta_) + r_));
    const scalarField Dd(D*patch().deltaCoeffs());   // D/delta

    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = K/(K + Dd);

    mixedFvPatchScalarField::updateCoeffs();
}


void langmuirStickingFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("s", s_);
    os.writeEntry("r", r_);
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
