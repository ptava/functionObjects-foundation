/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "sasRefineIndicator.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "IOobjectList.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sasRefineIndicator, 0);
    addToRunTimeSelectionTable(functionObject, sasRefineIndicator, dictionary);
namespace
{

// --- Default values
const word defaultResultName = "sasRefineIndicator";
const scalar sigmaDefault = 0.05;
const scalar weight1Default = 0.0;
const scalar weight2Default = 1.0;
const scalar alphaDefault = 0.6;


// --- Mathematical helpers
inline scalar expm1Safe(const scalar x)
{
    return (mag(x) < 1e-6) ? x + 0.5*sqr(x) : expm1(x);
}

inline scalar oneMinusExpNeg(const scalar x)
{
    return (mag(x) < 1e-6) ? x - 0.5*sqr(x) : -expm1(-x);
}

// --- Auxiliary functions for markCoreSafeScaler
inline scalar qOf(const scalar t)
{
    if (t > 700.0) return 0.0;
    
    const scalar denom = expm1Safe(t);
    if (mag(denom) <= VSMALL) return 2.0;
    
    return 2.0*t/denom;
}

scalar gOf(const scalar t, const scalar alphaEff, const scalar alphaTarget)
{
    if (t <= 0.0) return GREAT;

    const scalar a2  = sqr(alphaEff);
    const scalar a2t = a2*t;
    const scalar q   = qOf(a2t);

    const scalar num  = 2.0*alphaEff*oneMinusExpNeg(t);
    const scalar denA = oneMinusExpNeg(a2t);
    const scalar denB = 1.0 + a2 + (1.0 - a2)*q;

    if (mag(denA) <= VSMALL || mag(denB) <= VSMALL) return GREAT;

    return (num/(denA*denB)) - alphaTarget;
}

// Solves for tStar where beta(tStar) = alpha
scalar solveTStar(const scalar alphaEff, const scalar alphaTarget)
{
    const scalar tMin = 1e-12;
    const scalar tMax = 1e6;
    const label nGrid = 900;

    // 1. Grid Search
    scalar bestT = tMin;
    scalar bestErr = GREAT;
    scalar tLo = tMin;
    scalar gLo = gOf(tLo, alphaEff, alphaTarget);
    scalar tHi = tLo;
    bool bracketFound = false;

    for (label i = 1; i < nGrid; ++i)
    {
        const scalar s = scalar(i)/scalar(nGrid - 1);
        const scalar t = tMin*pow(tMax/tMin, s);
        const scalar gv = gOf(t, alphaEff, alphaTarget);

        if (mag(gv) < bestErr)
        {
            bestErr = mag(gv);
            bestT = t;
        }

        if (gLo*gv <= 0.0)
        {
            tHi = t;
            bracketFound = true;
            break;
        }

        tLo = t;
        gLo = gv;
    }

    if (!bracketFound)
    {
        return bestT;
    }

    // 2. Bisection Refinement
    scalar a = tLo;
    scalar b = tHi;
    scalar fa = gLo;
    
    const label maxIter = 120;
    const scalar tol = 1e-12;

    for (label it = 0; it < maxIter; ++it)
    {
        const scalar m = 0.5*(a + b);
        const scalar fm = gOf(m, alphaEff, alphaTarget);

        if (mag(fm) < tol || (b - a) < tol) return m;

        if (fa*fm <= 0.0)
        {
            b = m;
        }
        else
        {
            a = m;
            fa = fm;
        }
    }

    return 0.5*(a + b);
}

scalar calcKhat(const scalar tStar, const scalar alphaEff)
{
    const scalar a2tStar = sqr(alphaEff)*tStar;
    const scalar qStar = qOf(a2tStar);
    const scalar denomKhat = sqr(alphaEff)*(1.0 - qStar);

    if (mag(denomKhat) > VSMALL)
    {
        return (1.0 + qStar)/denomKhat;
    }
    return 0.0;
}

} // End anonymous namespace

const NamedEnum<sasRefineIndicator::focusRegion, 2>
sasRefineIndicator::focusRegionNames_
{
    "core",
    "periphery"
};

const NamedEnum<sasRefineIndicator::transferFunction, 5>
sasRefineIndicator::transferFunctionNames_
{
    "basic",
    "constant",
    "oddScaler",
    "safeScaler",
    "gaussSink"
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sasRefineIndicator::sasRefineIndicator
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    zone_(dict.lookup("cellZone")),
    focusRegion_(focusRegion::core),
    transferFunction_(transferFunction::constant),
    sigma_(sigmaDefault),
    weight1_(weight1Default),
    weight2_(weight2Default),
    alpha_(alphaDefault),
    controlRunTime_(false),
    resultName_(defaultResultName)
{
    read(dict);

    auto* fldPtr = new volScalarField
    (
        IOobject
        (
            resultName_,
            mesh_.time().name(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, -GREAT)
    );
    mesh_.objectRegistry::store(fldPtr);
}

// * * * * * * * * * * * * * * * *  Helpers  * * * * * * * * * * * * * * * * //

tmp<volScalarField::Internal> sasRefineIndicator::makeTmpInternal
(
    const word& name,
    const dimensionedScalar& init
) const
{
    return tmp<volScalarField::Internal>
    (
        new volScalarField::Internal
        (
            IOobject
            (
                name,
                mesh_.time().name(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            init
        )
    );
}

sasRefineIndicator::optimizationCoeffs sasRefineIndicator::getOptimizationCoeffs
(
    const scalar weight,
    const scalar alpha
) const
{
    // Static cache to avoid recomputing tStar/Khat if weight/alpha don't change
    static HashTable<optimizationCoeffs, word> cache(128);

    // Create a unique key for the current parameters
    const scalar qScale = 1e9;
    const label wKey = label(round(weight*qScale));
    const label aKey = label(round(alpha*qScale));
    const word key(Foam::name(wKey) + "_" + Foam::name(aKey));

    if (cache.found(key))
    {
        return cache[key];
    }

    // Calculate new coefficients
    const scalar alphaEff = (alpha + weight)/(1.0 + weight);
    
    if (!(alphaEff > 0.0 && alphaEff < 1.0))
    {
        FatalErrorInFunction
            << "Derived alphaEff must be in (0,1). Got " << alphaEff 
            << " (alpha=" << alpha << ", weight=" << weight << ")."
            << exit(FatalError);
    }

    optimizationCoeffs coeffs;
    coeffs.tStar = solveTStar(alphaEff, alpha);
    
    if (!(coeffs.tStar > 0.0))
    {
        FatalErrorInFunction
            << "Failed to determine a valid tStar." << exit(FatalError);
    }

    coeffs.Khat = calcKhat(coeffs.tStar, alphaEff);
    
    cache.insert(key, coeffs);
    return coeffs;
}

tmp<volScalarField::Internal> sasRefineIndicator::markCoreBasic
(
    const cellZone& cells,
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2
) const
{
    tmp<volScalarField::Internal> tG =
        makeTmpInternal
        (
            "tmpG",
            dimensionedScalar(dimless, -GREAT)
        );
    volScalarField::Internal& G = tG.ref();

    for(const label i : cells)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply constant value to all positive values.
        scalar d = c2[i] - c1[i];
        G[i] = d > 0.0 ? d : -GREAT;
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markCoreConstant
(
    const cellZone& cells,
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar weight
) const
{
    tmp<volScalarField::Internal> tG =
        makeTmpInternal
        (
            "tmpG",
            dimensionedScalar(dimless, -GREAT)
        );
    volScalarField::Internal& G = tG.ref();

    for(const label i : cells)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply constant value to all positive values.
        scalar d = c2[i] - c1[i];
        G[i] = d > 0.0 ? weight : -GREAT;
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markCoreOddScaler
(
    const cellZone& cells,
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar weight,
    const scalar sigma
) const
{
    const scalar invTwoSigma = 0.5/sqr(sigma);

    tmp<volScalarField::Internal> tG =
        makeTmpInternal
        (
            "tmpG",
            dimensionedScalar(dimless, -GREAT)
        );
    volScalarField::Internal& G = tG.ref();
    const tmp<volScalarField::Internal> td(c2 - c1);
    const volScalarField::Internal& d = td.ref();

    scalar dMaxLocal = -GREAT;
    for (const label i : cells)
    {
        dMaxLocal = max(dMaxLocal, d[i]);
    }
    const scalar dMax = returnReduce(dMaxLocal, maxOp<scalar>());

    if (dMax <= SMALL)
    {
        return tG;
    }

    const scalar shift = weight * dMax;
    for (const label celli : cells)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply an odd, monotonic, sign-preserving function
        const scalar u = d[celli] + shift;
        const scalar x = invTwoSigma * sqr(u);
        G[celli] = u * oneMinusExpNeg(x);
    }
    // Normalise by maximum value of G for d > -shift
    // const scalar maxG = gMax(G);
    // for (const label i : cells)
    // {
    //     const scalar u = d[i] + shift;
    //     const scalar Gi = G[i];
    //     G[i] = (u > 0.0) ? Gi/maxG : Gi;
    // }
    // Normalise by maximum value of G
    const scalar maxG = gMax(G);
    if (maxG > SMALL)
    {
        G /= maxG;
    }

    return tG;
}

tmp<volScalarField::Internal>
sasRefineIndicator::markCoreSafeScaler
(
    const cellZone& cells,
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar weight,
    const scalar sigma,
    const scalar alpha
) const
{
    if (alpha <= 0.0 || alpha > 1.0)
    {
        FatalErrorInFunction
            << "'alpha' must be in the range (0, 1] (got " << alpha << ")."
            << exit(FatalError);
    }

    if (alpha == 1.0)
    {
        return markCoreOddScaler(cells, c1, c2, weight, sigma);
    }

    tmp<volScalarField::Internal> tG =
        makeTmpInternal
        (
            "tmpG",
            dimensionedScalar(dimless, -GREAT)
        );
    volScalarField::Internal& G = tG.ref();
    const tmp<volScalarField::Internal> td(c2 - c1);
    const volScalarField::Internal& d = td.ref();

    scalar dMaxLocal = -GREAT;
    for (const label i : cells)
    {
        dMaxLocal = max(dMaxLocal, d[i]);
    }
    const scalar dMax = returnReduce(dMaxLocal, maxOp<scalar>());

    if (dMax <= SMALL)
    {
        return tG;
    }

    // Retrieve mathematical constants (cached)
    const optimizationCoeffs coeffs = getOptimizationCoeffs(weight, alpha);

    // Calculate scaling factors
    const scalar shift = weight*dMax;
    const scalar uMax  = max((1.0 + weight)*dMax, SMALL);
    const scalar sigmaStar = uMax/max(sqrt(2.0*coeffs.tStar), SMALL);
    const scalar kStar     = coeffs.Khat/sqr(uMax);
    const scalar invTwoSigma2 = 0.5/sqr(sigmaStar);

    // 5. Apply scaler field-wise
    for (const label celli : cells)
    {
        const scalar u = d[celli] + shift;

        // Base gaussian-like scaler
        scalar y = u*oneMinusExpNeg(sqr(u)*invTwoSigma2);

        // Apply damping denominator if u > 0
        if (u > 0.0)
        {
            const scalar denom = 1.0 + kStar*sqr(u);
            y = (denom > SMALL) ? (y / denom) : 0.0;
        }

        G[celli] = y;
    }

    // Normalise by maximum value of G for d > -shift
    // const scalar maxG = gMax(G);
    // for (const label i : cells)
    // {
    //     const scalar u = d[i] + shift;
    //     const scalar Gi = G[i];
    //     G[i] = (u > 0.0) ? Gi/maxG : Gi;
    // }
    // Normalise by maximum value of G
    const scalar maxG = gMax(G);
    if (maxG > SMALL)
    {
        G /= maxG;
    }

    return tG;
}


tmp<volScalarField::Internal> sasRefineIndicator::markCoreGaussSink
(
    const cellZone& cells,
    const volScalarField::Internal& Lvk,
    const volScalarField::Internal& c2,
    const scalar weight,
    const scalar sigma
) const
{
    const scalar invTwoSigma = 0.5/(sqr(sigma));

    tmp<volScalarField::Internal> tG =
        makeTmpInternal
        (
            "tmpG",
            dimensionedScalar(dimless, -GREAT)
        );
    tmp<volScalarField::Internal> tnLvk = Lvk / 
        ( c2 + dimensionedScalar(dimLength, VSMALL));

    auto& G = tG.ref();
    auto& nLvk = tnLvk.ref();

    for(const label i : cells)
    {
        // Calculate normalised von Karman length scale for each cell
        // and apply the Gaussian function to get the indicator value.
        const scalar& nLvki = nLvk[i];
        G[i] = exp(-invTwoSigma * sqr(nLvki - 1))
             - weight * sqr(nLvki - 1);
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markPeripheryGaussSink
(
    const cellZone& cells,
    const volScalarField::Internal& Lvk,
    const scalar LvkRef,
    const scalar weight,
    const scalar sigma
) const
{
    const scalar invTwoSigma = 0.5/(sqr(sigma));

    tmp<volScalarField::Internal> tG =
        makeTmpInternal
        (
            "tmpG",
            dimensionedScalar(dimless, -GREAT)
        );
    tmp<volScalarField::Internal> tnLvk = Lvk / LvkRef;

    auto& G = tG.ref();
    auto& nLvk = tnLvk.ref();

    for(const label i : cells)
    {
        // Calculate normalised von Karman length scale for each cell
        // and apply the Gaussian function to get the indicator value.
        const scalar& nLvki = nLvk[i];
        G[i] = exp(-invTwoSigma * sqr(nLvki - 1))
             - weight * sqr(nLvki - 1);
    }

    return tG;
}

// * * * * * * * * * * * * * *  Calculation  * * * * * * * * * * * * * * * * //

void sasRefineIndicator::calcIndicator()
{
    const cellZone& cells = mesh().cellZones()[zone_];

    if (debug)
    {
        Pout<< "Cells in zone: " << cells.size() << endl;
    }

    auto& fld = mesh_.lookupObjectRef<volScalarField>(resultName_);

    auto& fldI = fld.internalFieldRef();

    switch (focusRegion_)
    {
        case focusRegion::core:
        {
            const volScalarField& C1  = lookupObject<volScalarField>("C1");
            const volScalarField& C2  = lookupObject<volScalarField>("C2");
            const auto& C1I  = C1.internalField();
            const auto& C2I  = C2.internalField();
            switch (transferFunction_)
            {
                case transferFunction::basic:
                {
                    fldI = markCoreBasic(cells, C1I, C2I);
                    break;
                }
                case transferFunction::constant:
                {
                    fldI = markCoreConstant(cells, C1I, C2I, weight2_);
                    break;
                }
                case transferFunction::oddScaler:
                {
                    fldI = markCoreOddScaler(cells, C1I, C2I, weight1_, sigma_);
                    break;
                }
                case transferFunction::safeScaler:
                {
                    // TO DO: read dynamically from dict
                    // if (controlRunTime_)
                    // {
                    //     alpha_ = readScalar(dict().lookup("alpha"));
                    // }
                    fldI = markCoreSafeScaler(cells, C1I, C2I, weight1_, sigma_, alpha_);
                    break;
                }
                case transferFunction::gaussSink:
                {
                    const volScalarField& Lvk = lookupObject<volScalarField>("Lvk");
                    const auto& LvkI = Lvk.internalField();
                    fldI = markCoreGaussSink
                    (
                        cells,
                        LvkI,
                        C2I,
                        weight2_,
                        sigma_
                    );
                    break;
                }
            }
            break;
        }
        case focusRegion::periphery:
        {
            const volScalarField& Lvk = lookupObject<volScalarField>("Lvk");
            const auto& LvkI = Lvk.internalField();
            fldI = markPeripheryGaussSink
            (
                cells,
                LvkI,
                LvkRef_,
                weight2_,
                sigma_
            );
            break;
        }
    }

    // Debug info - range of indicator field
    DebugInfo
        << name() << " (" << resultName_ << "): range = ["
        << gMin(fld.internalField()) << ", "
        << gMax(fld.internalField()) << "]" << nl;

    // Debug info - number of cells with indicator >= 0
    if (debug)
    {
        label nPos = 0;
        forAll(fldI, i)
        {
            if (fldI[i] >= 0)
            {
                nPos++;
            }
        }
        DebugInfo
            << name() << " (" << resultName_ << "): found "
            << returnReduce(nPos, sumOp<label>()) << "/"
            << cells.size()
            << " cells with indicator >= 0: " << endl;
    }
}

// * * * * * * * * * * * * * *  Read/Execute/Write  * * * * * * * * * * * * * //

bool sasRefineIndicator::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    requiredFields_.clear();
    resultName_ = dict.lookupOrDefault<word>("result", defaultResultName);
    sigma_ = dict.lookupOrDefault<scalar>("sigma", sigmaDefault);
    focusRegion_ = focusRegionNames_.read(dict.lookup("focusRegion"));

    switch (focusRegion_)
    {
        case focusRegion::core:
        {
            requiredFields_.append("C1");
            requiredFields_.append("C2");
            requiredFields_.append("Lvk");
            weight1_ = dict.lookupOrDefault<scalar>("weight1", weight1Default);
            weight2_ = dict.lookupOrDefault<scalar>("weight2", weight2Default);
            alpha_ = dict.lookupOrDefault<scalar>("alpha", alphaDefault);
            transferFunction_ = transferFunctionNames_.read(
                dict.lookup("transferFunction")
            );

            break;
        }
        case focusRegion::periphery:
        {
            requiredFields_.append("Lvk");
            weight1_ = dict.lookupOrDefault<scalar>("weight1", weight1Default);
            weight2_ = dict.lookupOrDefault<scalar>("weight2", weight2Default);
            LvkRef_ = dict.lookup<scalar>("LvkRef");
            break;
        }
        default:
        {
            FatalIOErrorInFunction(dict)
                << "Unknown 'focusRegion' " << focusRegionNames_[focusRegion_]
                << exit(FatalIOError);
        }
    }

    if (sigma_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "'sigma' must be > 0 (got " << sigma_ << ")." 
            << exit(FatalIOError);
    }

    if (weight1_ < 0.0 || weight2_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "'weight1' and 'weight2' must be non-negative (got "
            << weight1_ << " and " << weight2_ << ")."
            << exit(FatalIOError);
    }

    if (log)
    {
        Info<< type() << ' ' << name() << ':' << nl
            << "  focusRegion      : " << focusRegionNames_[focusRegion_] << nl
            << "  transferFunction : " << transferFunctionNames_[transferFunction_] << nl
            << "  sigma            : " << sigma_ << nl
            << "  weight1          : " << weight1_ << nl
            << "  weight2          : " << weight2_ << nl
            << "  alpha            : " << alpha_ << nl
            << "  result           : " << resultName_ << nl
            << "  nCells           : " << mesh().cellZones()[zone_].size() << nl
            << endl;
    }

    return true;
}

wordList sasRefineIndicator::fields() const
{
    return requiredFields_;
}

bool sasRefineIndicator::execute()
{
    for (const auto& fieldName : requiredFields_)
    {
        if (!foundObject<volScalarField>(fieldName))
        {
            FatalErrorInFunction
                << "Required field '" << fieldName << "' not found."
                << "Available volScalarField objects in object registry: "
                << obr_.lookupClass<volScalarField>().toc()
                << exit(FatalError);
        }
    }

    calcIndicator();
    return true;
}

bool sasRefineIndicator::write()
{
    Info
        << "FunctionObjects::" << name() 
        << " writing field: " << resultName_ << endl;

    const auto& fld = mesh_.lookupObject<volScalarField>(resultName_);
    fld.write();
    return true;
}

} // End namespace functionObjects
} // End namespace Foam
