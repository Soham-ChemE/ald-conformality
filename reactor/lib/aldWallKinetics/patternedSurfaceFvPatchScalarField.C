#include "patternedSurfaceFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IOdictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * profile helpers * * * * * * * * * * * * * * * //

patternedSurfaceFvPatchScalarField::profileTable
patternedSurfaceFvPatchScalarField::readProfile(const dictionary& dict, const word& name)
{
    profileTable tab;
    if (!dict.found(name)) return tab;

    const entry& e = dict.lookupEntry(name, keyType::LITERAL);
    ITstream& is = e.stream();
    const word form(is);
    if (form != "table")
    {
        FatalIOErrorInFunction(dict)
            << name << ": expected 'table ((radius value) ...)', found '" << form << "'"
            << exit(FatalIOError);
    }
    is >> tab;
    if (tab.empty())
    {
        FatalIOErrorInFunction(dict) << name << ": empty table" << exit(FatalIOError);
    }
    for (label i = 1; i < tab.size(); ++i)
    {
        if (tab[i].first() < tab[i-1].first())
        {
            FatalIOErrorInFunction(dict)
                << name << ": radii must be in ascending order" << exit(FatalIOError);
        }
    }
    return tab;
}


// * * * * * * * * * * * * * * * library plumbing * * * * * * * * * * * * * * * //

aldkin_params patternedSurfaceFvPatchScalarField::params() const
{
    aldkin_params p;
    aldkin_params_default(&p);
    p.s = s_; p.r = r_; p.kdes = 0.0; p.vbar = vbar_; p.N0 = N0_; p.c0 = c0_;
    p.d = d_; p.AR = AR_; p.holeDensity = holeDensity_; p.nCells = nCells_;
    p.flatOnly = 0; p.scheme = ALDKIN_SCHEME_EULER;
    return p;
}


scalarField patternedSurfaceFvPatchScalarField::faceRadii() const
{
    const vectorField& Cf = patch().Cf();
    scalarField r(Cf.size());
    forAll(Cf, facei)
    {
        r[facei] = Foam::sqrt(sqr(Cf[facei].y()) + sqr(Cf[facei].z()));   // wedge axis is x
    }
    return r;
}


void patternedSurfaceFvPatchScalarField::setProfile(int which, const profileTable& tab)
{
    scalarField r(tab.size()), v(tab.size());
    forAll(tab, i) { r[i] = tab[i].first(); v[i] = tab[i].second(); }
    kin_.check(aldkin_set_profile(kin_.get(), which, tab.size(), r.cdata(), v.cdata()), "aldkin_set_profile");
}


void patternedSurfaceFvPatchScalarField::resetKinetics()
{
    kin_ = aldkinHandle(params());
    const scalarField r(faceRadii());
    kin_.check(aldkin_set_faces(kin_.get(), r.size(), r.cdata()), "aldkin_set_faces");
    if (!ARProfile_.empty()) setProfile(ALDKIN_PROFILE_AR, ARProfile_);
    if (!holeDensityProfile_.empty()) setProfile(ALDKIN_PROFILE_DENSITY, holeDensityProfile_);
}


// * * * * * * * * * * * * * * * * constructors * * * * * * * * * * * * * * * * //

patternedSurfaceFvPatchScalarField::patternedSurfaceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    s_(1), r_(0), vbar_(410), N0_(4e18), c0_(1), D_(-1),
    d_(1e-7), AR_(50), holeDensity_(0), nCells_(100),
    ARProfile_(), holeDensityProfile_(),
    kin_(),
    curTimeIndex_(-1)
{
    resetKinetics();
    refValue() = 0.0; refGrad() = 0.0; valueFraction() = 1.0;
}


patternedSurfaceFvPatchScalarField::patternedSurfaceFvPatchScalarField
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
    d_(dict.get<scalar>("d")),
    AR_(dict.get<scalar>("AR")),
    holeDensity_(dict.get<scalar>("holeDensity")),
    nCells_(dict.getOrDefault<label>("nCells", 100)),
    ARProfile_(readProfile(dict, "ARProfile")),
    holeDensityProfile_(readProfile(dict, "holeDensityProfile")),
    kin_(),
    curTimeIndex_(-1)
{
    resetKinetics();
    if (dict.found("thetaFlat"))
    {
        const scalarField thetaFlat("thetaFlat", dict, p.size());
        kin_.check(aldkin_set_theta(kin_.get(), ALDKIN_THETA_FLAT, thetaFlat.cdata()), "aldkin_set_theta(thetaFlat)");
    }
    if (dict.found("thetaHole"))
    {
        const scalarField thetaHole("thetaHole", dict, p.size()*nCells_);
        kin_.check(aldkin_set_theta(kin_.get(), ALDKIN_THETA_HOLE, thetaHole.cdata()), "aldkin_set_theta(thetaHole)");
    }
    refValue() = 0.0; refGrad() = 0.0; valueFraction() = 1.0;
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(0.0);
    }
}


patternedSurfaceFvPatchScalarField::patternedSurfaceFvPatchScalarField
(
    const patternedSurfaceFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    s_(ptf.s_), r_(ptf.r_), vbar_(ptf.vbar_), N0_(ptf.N0_), c0_(ptf.c0_), D_(ptf.D_),
    d_(ptf.d_), AR_(ptf.AR_), holeDensity_(ptf.holeDensity_), nCells_(ptf.nCells_),
    ARProfile_(ptf.ARProfile_), holeDensityProfile_(ptf.holeDensityProfile_),
    kin_(),
    curTimeIndex_(-1)
{
    // new faces at their own radii (pattern from the profiles); the flat coverage is mapped, the hole coverages
    // start from 0 as before (mapping is not used in these cases)
    resetKinetics();
    const scalarField thetaFlat(ptf.kin_.theta(ALDKIN_THETA_FLAT, ptf.size()), mapper);
    kin_.check(aldkin_set_theta(kin_.get(), ALDKIN_THETA_FLAT, thetaFlat.cdata()), "aldkin_set_theta(thetaFlat)");
}


patternedSurfaceFvPatchScalarField::patternedSurfaceFvPatchScalarField
(
    const patternedSurfaceFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    s_(ptf.s_), r_(ptf.r_), vbar_(ptf.vbar_), N0_(ptf.N0_), c0_(ptf.c0_), D_(ptf.D_),
    d_(ptf.d_), AR_(ptf.AR_), holeDensity_(ptf.holeDensity_), nCells_(ptf.nCells_),
    ARProfile_(ptf.ARProfile_), holeDensityProfile_(ptf.holeDensityProfile_),
    kin_(ptf.kin_),                 // clones the whole kinetics state
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * * mapping * * * * * * * * * * * * * * * * * //

void patternedSurfaceFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
{
    mixedFvPatchScalarField::autoMap(m);
    scalarField thetaFlat(kin_.theta(ALDKIN_THETA_FLAT, kin_.valid() ? aldkin_n_faces(kin_.get()) : 0));
    thetaFlat.autoMap(m);
    resetKinetics();     // recompute the pattern from the profiles at the new face centres; hole coverages 0
    kin_.check(aldkin_set_theta(kin_.get(), ALDKIN_THETA_FLAT, thetaFlat.cdata()), "aldkin_set_theta(thetaFlat)");
}


void patternedSurfaceFvPatchScalarField::rmap(const fvPatchScalarField& ptf, const labelList& addr)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
    const patternedSurfaceFvPatchScalarField& tiptf = refCast<const patternedSurfaceFvPatchScalarField>(ptf);
    scalarField thetaFlat(kin_.theta(ALDKIN_THETA_FLAT, size()));
    thetaFlat.rmap(tiptf.kin_.theta(ALDKIN_THETA_FLAT, tiptf.size()), addr);
    kin_.check(aldkin_set_theta(kin_.get(), ALDKIN_THETA_FLAT, thetaFlat.cdata()), "aldkin_set_theta(thetaFlat)");
}


// * * * * * * * * * * * * * * * * * update * * * * * * * * * * * * * * * * * //

void patternedSurfaceFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    scalar D = D_;
    if (D < 0)
    {
        const IOdictionary& tp = db().lookupObject<IOdictionary>("transportProperties");
        D = tp.get<dimensionedScalar>("DT").value();
    }

    const label ti = db().time().timeIndex();
    const bool advance = (ti != curTimeIndex_);
    const scalar dt = db().time().deltaTValue();
    const scalarField& Tface = *this;

    // coverages advance once per time step with the face value of the previous step; the further calls within
    // the step only re-evaluate K_eff with the current face value
    if (advance)
    {
        kin_.check(aldkin_step(kin_.get(), Tface.cdata(), dt, nullptr, nullptr), "aldkin_step");
        curTimeIndex_ = ti;
    }
    else
    {
        kin_.check(aldkin_evaluate(kin_.get(), Tface.cdata(), nullptr, nullptr), "aldkin_evaluate");
    }

    const scalarField Dd(D*patch().deltaCoeffs());
    kin_.check
    (
        aldkin_mixed_bc(kin_.get(), Dd.cdata(), refValue().data(), refGrad().data(), valueFraction().data()),
        "aldkin_mixed_bc"
    );

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * output * * * * * * * * * * * * * * * * * //

void patternedSurfaceFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("s", s_); os.writeEntry("r", r_); os.writeEntry("vbar", vbar_);
    os.writeEntry("N0", N0_); os.writeEntry("c0", c0_);
    if (D_ > 0) os.writeEntry("D", D_);
    os.writeEntry("d", d_); os.writeEntry("AR", AR_); os.writeEntry("holeDensity", holeDensity_);
    os.writeEntry("nCells", nCells_);
    if (!ARProfile_.empty())
    {
        os.writeKeyword("ARProfile") << "table " << ARProfile_ << token::END_STATEMENT << nl;
    }
    if (!holeDensityProfile_.empty())
    {
        os.writeKeyword("holeDensityProfile") << "table " << holeDensityProfile_ << token::END_STATEMENT << nl;
    }
    const label n = size();
    scalarField ARField(n, 0.0), holeDensityField(n, 0.0), Keff(n, 0.0);
    kin_.check(aldkin_face_pattern(kin_.get(), ARField.data(), holeDensityField.data()), "aldkin_face_pattern");
    kin_.check(aldkin_keff(kin_.get(), Keff.data()), "aldkin_keff");
    ARField.writeEntry("ARField", os);
    holeDensityField.writeEntry("holeDensityField", os);
    kin_.theta(ALDKIN_THETA_FLAT, n).writeEntry("thetaFlat", os);
    kin_.theta(ALDKIN_THETA_BOTTOM, n).writeEntry("thetaBottom", os);
    kin_.theta(ALDKIN_THETA_MID, n).writeEntry("thetaMid", os);
    Keff.writeEntry("Keff", os);
    kin_.theta(ALDKIN_THETA_HOLE, n*nCells_).writeEntry("thetaHole", os);
    writeEntry("value", os);
}


makePatchTypeField(fvPatchScalarField, patternedSurfaceFvPatchScalarField);

} // End namespace Foam
