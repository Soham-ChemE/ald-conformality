#include "patternedSurfaceFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "mathematicalConstants.H"

namespace Foam
{

patternedSurfaceFvPatchScalarField::patternedSurfaceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    s_(1), r_(0), vbar_(410), N0_(4e18), c0_(1), D_(-1),
    d_(1e-7), AR_(50), holeDensity_(0), nCells_(100),
    thetaFlat_(p.size(), 0.0), thetaHole_(p.size()*nCells_, 0.0), Keff_(p.size(), 0.0),
    curTimeIndex_(-1)
{
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
    thetaFlat_(p.size(), 0.0),
    thetaHole_(p.size()*nCells_, 0.0),
    Keff_(p.size(), 0.0),
    curTimeIndex_(-1)
{
    if (dict.found("thetaFlat")) thetaFlat_ = scalarField("thetaFlat", dict, p.size());
    if (dict.found("thetaHole")) thetaHole_ = scalarField("thetaHole", dict, p.size()*nCells_);
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
    thetaFlat_(ptf.thetaFlat_, mapper),
    thetaHole_(p.size()*ptf.nCells_, 0.0),     // not mapped (mapping is not used in these cases)
    Keff_(ptf.Keff_, mapper),
    curTimeIndex_(-1)
{}


patternedSurfaceFvPatchScalarField::patternedSurfaceFvPatchScalarField
(
    const patternedSurfaceFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    s_(ptf.s_), r_(ptf.r_), vbar_(ptf.vbar_), N0_(ptf.N0_), c0_(ptf.c0_), D_(ptf.D_),
    d_(ptf.d_), AR_(ptf.AR_), holeDensity_(ptf.holeDensity_), nCells_(ptf.nCells_),
    thetaFlat_(ptf.thetaFlat_), thetaHole_(ptf.thetaHole_), Keff_(ptf.Keff_),
    curTimeIndex_(-1)
{}


void patternedSurfaceFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
{
    mixedFvPatchScalarField::autoMap(m);
    thetaFlat_.autoMap(m);
    Keff_.autoMap(m);
    thetaHole_.setSize(size()*nCells_, 0.0);
}


void patternedSurfaceFvPatchScalarField::rmap(const fvPatchScalarField& ptf, const labelList& addr)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
    const patternedSurfaceFvPatchScalarField& tiptf = refCast<const patternedSurfaceFvPatchScalarField>(ptf);
    thetaFlat_.rmap(tiptf.thetaFlat_, addr);
    Keff_.rmap(tiptf.Keff_, addr);
}


// Quasi-steady Knudsen diffusion with Langmuir sink down one hole; mouth density = 1 (relative), closed bottom.
// Returns the mouth flux per unit mouth area per unit mouth density [m/s]; nrel receives the relative density profile.
scalar patternedSurfaceFvPatchScalarField::solveHole(const label facei, scalarField& nrel) const
{
    const label N = nCells_;
    const scalar Lh = AR_*d_;
    const scalar dx = Lh/N;
    const scalar DK = d_*vbar_/3.0;                 // Knudsen diffusivity, cylinder
    const scalar a = DK/(dx*dx);
    List<scalar> lower(N, -a), diag(N, 0.0), upper(N, -a), rhs(N, 0.0);
    for (label i = 0; i < N; ++i)
    {
        const scalar th = thetaHole_[facei*N + i];
        const scalar k = (vbar_/d_)*(s_*(1.0 - th) + r_);   // (4/d)(vbar/4)[...]
        diag[i] = 2.0*a + k;
    }
    diag[N-1] = a + (diag[N-1] - 2.0*a);            // Neumann at the bottom: only the west face conducts
    diag[0] += a;                                   // Dirichlet ghost at the mouth: n_ghost = 2 - n_0
    rhs[0] = 2.0*a;
    // Thomas algorithm
    List<scalar> cp(N, 0.0), dp(N, 0.0);
    cp[0] = upper[0]/diag[0]; dp[0] = rhs[0]/diag[0];
    for (label i = 1; i < N; ++i)
    {
        const scalar m = diag[i] - lower[i]*cp[i-1];
        cp[i] = (i < N-1) ? upper[i]/m : 0.0;
        dp[i] = (rhs[i] - lower[i]*dp[i-1])/m;
    }
    nrel[N-1] = dp[N-1];
    for (label i = N-2; i >= 0; --i) nrel[i] = dp[i] - cp[i]*nrel[i+1];
    // flux into the mouth: D (n_mouth - n_0)/(dx/2) with n_mouth = 1
    return DK*(1.0 - nrel[0])/(0.5*dx);
}


void patternedSurfaceFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    scalar D = D_;
    if (D < 0)
    {
        const IOdictionary& tp = db().lookupObject<IOdictionary>("transportProperties");
        D = tp.get<dimensionedScalar>("DT").value();
    }

    const scalar phi = holeDensity_*constant::mathematical::pi*d_*d_/4.0;
    const label N = nCells_;
    const label ti = db().time().timeIndex();
    const bool advance = (ti != curTimeIndex_);
    const scalar dt = db().time().deltaTValue();
    const scalarField& Tface = *this;

    scalarField nrel(N, 0.0);
    forAll(Keff_, facei)
    {
        const scalar Tf = max(Tface[facei], scalar(0));
        const scalar jhole = solveHole(facei, nrel);          // per unit mouth density
        if (advance)
        {
            // flat part
            thetaFlat_[facei] = min(scalar(1), thetaFlat_[facei]
                + dt*(vbar_/4.0)*s_*(1.0 - thetaFlat_[facei])*c0_*Tf/N0_);
            // hole cells with their local relative density
            for (label i = 0; i < N; ++i)
            {
                scalar& th = thetaHole_[facei*N + i];
                th = min(scalar(1), th + dt*(vbar_/4.0)*s_*(1.0 - th)*c0_*Tf*nrel[i]/N0_);
            }
        }
        Keff_[facei] = (1.0 - phi)*(vbar_/4.0)*(s_*(1.0 - thetaFlat_[facei]) + r_) + phi*jhole;
    }
    if (advance) curTimeIndex_ = ti;

    const scalarField Dd(D*patch().deltaCoeffs());
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = Keff_/(Keff_ + Dd);

    mixedFvPatchScalarField::updateCoeffs();
}


void patternedSurfaceFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("s", s_); os.writeEntry("r", r_); os.writeEntry("vbar", vbar_);
    os.writeEntry("N0", N0_); os.writeEntry("c0", c0_);
    if (D_ > 0) os.writeEntry("D", D_);
    os.writeEntry("d", d_); os.writeEntry("AR", AR_); os.writeEntry("holeDensity", holeDensity_);
    os.writeEntry("nCells", nCells_);
    thetaFlat_.writeEntry("thetaFlat", os);
    scalarField bottom(size(), 0.0), mid(size(), 0.0);
    forAll(bottom, facei)
    {
        bottom[facei] = thetaHole_[facei*nCells_ + nCells_ - 1];
        mid[facei] = thetaHole_[facei*nCells_ + nCells_/2];
    }
    bottom.writeEntry("thetaBottom", os);
    mid.writeEntry("thetaMid", os);
    Keff_.writeEntry("Keff", os);
    thetaHole_.writeEntry("thetaHole", os);
    writeEntry("value", os);
}


makePatchTypeField(fvPatchScalarField, patternedSurfaceFvPatchScalarField);

} // End namespace Foam
