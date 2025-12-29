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
const word defaultResultName = "sasRefineIndicator";
const scalar sigmaDefault = 0.05;
const scalar weight1Default = 1.0;
const scalar weight2Default = 1.0;
}

const NamedEnum<sasRefineIndicator::focusRegion, 2>
sasRefineIndicator::focusRegionNames_
{
    "core",
    "periphery"
};

const NamedEnum<sasRefineIndicator::transferFunction, 4>
sasRefineIndicator::transferFunctionNames_
{
    "basic",
    "constant",
    "oddScaler",
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
    focusRegion_(focusRegion::core),
    transferFunction_(transferFunction::constant),
    sigma_(sigmaDefault),
    weight1_(weight1Default),
    weight2_(weight2Default),
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

tmp<volScalarField::Internal> sasRefineIndicator::markCoreBasic
(
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2
) const
{
    DebugInfo
        << "sasRefineIndicator::markCoreBasic called" << nl;

    tmp<volScalarField::Internal> tG
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tmpG",
                mesh_.time().name(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, -GREAT)
        )
    );

    auto& G = tG.ref();
    forAll(G, i)
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
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar weight1
) const
{
    DebugInfo
        << "sasRefineIndicator::markCoreConstant called" << nl;

    tmp<volScalarField::Internal> tG
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tmpG",
                mesh_.time().name(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, -GREAT)
        )
    );

    auto& G = tG.ref();
    forAll(G, i)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply constant value to all positive values.
        scalar d = c2[i] - c1[i];
        G[i] = d > 0.0 ? weight1 : -GREAT;
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markCoreOddScaler
(
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar weight1,
    const scalar sigma
) const
{
    DebugInfo
        << "sasRefineIndicator::markCoreOddScaler called" << nl;

    const scalar invTwoSigma = 0.5/sqr(sigma);

    tmp<volScalarField::Internal> tG
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tmpG",
                mesh_.time().name(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, -GREAT)
        )
    );
    tmp<volScalarField::Internal> td(c2 - c1);

    auto& G = tG.ref();
    auto& d = td.ref();
    
    forAll(G, i)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply an odd, monotonic, sign-preserving function
        const scalar& di = d[i];
        G[i] = di * (weight1 * exp(-invTwoSigma * sqr(di)));
        
    }
    G /= gMax(G);

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markCoreGaussSink
(
    const volScalarField::Internal& Lvk,
    const volScalarField::Internal& c2,
    const scalar weight1,
    const scalar weight2,
    const scalar sigma
) const
{
    DebugInfo
        << "sasRefineIndicator::markCoreGaussSink called" << nl;

    const scalar invTwoSigma = 0.5/(sqr(sigma));

    tmp<volScalarField::Internal> tG
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tmpG",
                mesh_.time().name(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, -GREAT)
        )
    );
    tmp<volScalarField::Internal> tnLvk = Lvk / 
        ( c2 + dimensionedScalar(dimLength, VSMALL));

    auto& G = tG.ref();
    auto& nLvk = tnLvk.ref();

    forAll(G, i)
    {
        // Calculate normalised von Karman length scale for each cell
        // and apply the Gaussian function to get the indicator value.
        const scalar& nLvki = nLvk[i];
        G[i] = weight1 * exp(-invTwoSigma * sqr(nLvki - 1))
             - weight2 * sqr(nLvki - 1);
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markPeripheryGaussSink
(
    const volScalarField::Internal& Lvk,
    const scalar LvkRef,
    const scalar weight1,
    const scalar weight2,
    const scalar sigma
) const
{
    DebugInfo
        << "sasRefineIndicator::markPeripheryGaussSink called" << nl;

    const scalar invTwoSigma = 0.5/(sqr(sigma));

    tmp<volScalarField::Internal> tG
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "tmpG",
                mesh_.time().name(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, -GREAT)
        )
    );
    tmp<volScalarField::Internal> tnLvk = Lvk / LvkRef;

    auto& G = tG.ref();
    auto& nLvk = tnLvk.ref();

    forAll(G, i)
    {
        // Calculate normalised von Karman length scale for each cell
        // and apply the Gaussian function to get the indicator value.
        const scalar& nLvki = nLvk[i];
        G[i] = weight1 * exp(-invTwoSigma * sqr(nLvki - 1))
             - weight2 * sqr(nLvki - 1);
    }

    return tG;
}

// * * * * * * * * * * * * * *  Calculation  * * * * * * * * * * * * * * * * //

void sasRefineIndicator::calcIndicator()
{
    const volScalarField& Lvk = lookupObject<volScalarField>("Lvk");

    auto& fld = mesh_.lookupObjectRef<volScalarField>(resultName_);
    auto& fldI = fld.internalFieldRef();

    const auto& LvkI = Lvk.internalField();

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
                    fldI = markCoreBasic(C1I, C2I);
                    break;
                }
                case transferFunction::constant:
                {
                    fldI = markCoreConstant(C1I, C2I, weight1_);
                    break;
                }
                case transferFunction::oddScaler:
                {
                    fldI = markCoreOddScaler(C1I, C2I, weight2_, sigma_);
                    break;
                }
                case transferFunction::gaussSink:
                {
                    fldI = markCoreGaussSink
                    (
                        LvkI,
                        C2I,
                        weight1_,
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
            fldI = markPeripheryGaussSink
            (
                LvkI,
                LvkRef_,
                weight1_,
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
            << nPos << "/" << fldI.size()
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

    // Testing
    switch (focusRegion_)
    {
        case focusRegion::core:
        {
            requiredFields_.append("C1");
            requiredFields_.append("C2");
            requiredFields_.append("Lvk");
            weight1_ = dict.lookupOrDefault<scalar>("weight1", weight1Default);
            weight2_ = dict.lookupOrDefault<scalar>("weight2", weight2Default);
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
            << "'sigma' must be > 0 (got " << sigma_ << ")." << exit(FatalIOError);
    }

    if (log)
    {
        Info<< type() << ' ' << name() << ':' << nl
            << "  focusRegion      : " << focusRegionNames_[focusRegion_] << nl
            << "  transferFunction : " << transferFunctionNames_[transferFunction_] << nl
            << "  sigma            : " << sigma_ << nl
            << "  weight1          : " << weight1_ << nl
            << "  weight2          : " << weight2_ << nl
            << "  result           : " << resultName_ << nl
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
