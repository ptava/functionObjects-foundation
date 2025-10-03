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

const NamedEnum<sasRefineIndicator::focusRegion, 3>
sasRefineIndicator::focusRegionNames_
{
    "periphery",
    "core",
    "combined"
};

const NamedEnum<sasRefineIndicator::transferFunction, 2>
sasRefineIndicator::transferFunctionNames_
{
    "constant",
    "oddScaler"
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sasRefineIndicator::sasRefineIndicator
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    // Read configuration first (sets zoneSubSetPtr_, resultName_, etc.)
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

tmp<volScalarField::Internal> sasRefineIndicator::markCoreConstant
(
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar coreWeight
) const
{
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
    for (const label i : G)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply constant value to all positive values.
        scalar d = c2[i] - c1[i];
        G[i] = d > 0.0 ? coreWeight : -GREAT;
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markCoreOddScaler
(
    const volScalarField::Internal& c1,
    const volScalarField::Internal& c2,
    const scalar coreWeight,
    const scalar sigma
) const
{
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

    auto& G = tG.ref();
    for (const label i : G)
    {
        // Calculate the difference between c2 and c1 for each cell
        // and apply an odd, monotonic, sign-preserving function
        scalar d = c2[i] - c1[i];
        G[i] = d * (1 + coreWeight * (1 - exp(-invTwoSigma * sqr(d))));
    }

    return tG;
}

tmp<volScalarField::Internal> sasRefineIndicator::markPeripheryGaussSink
(
    const volScalarField::Internal& Lvk,
    const volScalarField::Internal& c2,
    const scalar peripheryWeight1,
    const scalar peripheryWeight2,
    const scalar sigma
) const
{
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

    auto& G = tG.ref();
    for (label i = 0; i < G.size(); ++i)
    {
        // Calculate normalised von Karman length scale for each cell
        // and apply the Gaussian function to get the indicator value.
        scalar nLvk = Lvk[i] / (c2[i] + VSMALL);
        G[i] = peripheryWeight1 * exp(-invTwoSigma * sqr(nLvk - 1))
             - peripheryWeight2 * sqr(nLvk - 1);
    }

    return tG;
}

// * * * * * * * * * * * * * *  Calculation  * * * * * * * * * * * * * * * * //

void sasRefineIndicator::calcIndicator()
{
    // Mandatory fields on main mesh
    if
    (
        !obr_.foundObject<volScalarField>("Lvk")
     || !obr_.foundObject<volScalarField>("C1")
     || !obr_.foundObject<volScalarField>("C2")
    )
    {
        FatalErrorInFunction
            << "Required fields 'Lvk', 'C1', 'C2' must exist in objectRegistry '"
            << obr_.name() << "'. None are optional." << nl
            << "Available volScalarField objects: "
            << obr_.lookupClass<volScalarField>().toc()
            << exit(FatalError);
    }

    const volScalarField& Lvk = obr_.lookupObject<volScalarField>("Lvk");
    const volScalarField& C1  = obr_.lookupObject<volScalarField>("C1");
    const volScalarField& C2  = obr_.lookupObject<volScalarField>("C2");

    auto& fld = mesh_.lookupObjectRef<volScalarField>(resultName_);
    auto& fldI = fld.internalFieldRef();

    const auto& LvkI = Lvk.internalField();
    const auto& C1I  = C1.internalField();
    const auto& C2I  = C2.internalField();

    switch (focusRegion_)
    {
        case focusRegion::core:
            // compute fldI based on the function name selected
            fldI = transferFunction_ == transferFunction::constant
                ? markCoreConstant(C1I, C2I, coreWeight_)
                : markCoreOddScaler(C1I, C2I, coreWeight_, sigma_);
            break;
        case focusRegion::periphery:
            fldI = markPeripheryGaussSink
            (
                LvkI,
                C2I,
                peripheryWeight1_,
                peripheryWeight2_,
                sigma_
            );
            break;
        case focusRegion::combined:
            break;
    }

    DebugInfo
        << type() << " '" << name() << "': indicator range = ["
        << gMin(fld.internalField()) << ", "
        << gMax(fld.internalField()) << "]" << nl;
}

// * * * * * * * * * * * * * *  Read/Execute/Write  * * * * * * * * * * * * * //

bool sasRefineIndicator::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    resultName_ = dict.lookupOrDefault<word>("result", "sasRefineIndicator");
    sigma_ = dict.lookupOrDefault<scalar>("sigma", 0.05);
    focusRegion_ = focusRegionNames_.read(dict.lookup("focusRegion"));

    // Testing
    switch (focusRegion_)
    {
        case focusRegion::core:
        {
            requiredFields_.append("C1");
            requiredFields_.append("C2");
            coreWeight_ = dict.lookupOrDefault<scalar>("coreWeight", 10.0);
            transferFunction_ = transferFunctionNames_.read(
                dict.lookup("transferFunction")
            );
            break;
        }
        case focusRegion::periphery:
        {
            requiredFields_.append("Lvk");
            requiredFields_.append("C2");
            peripheryWeight1_ = dict.lookupOrDefault<scalar>("peripheryWeight1", 1000.0);
            peripheryWeight2_ = dict.lookupOrDefault<scalar>("peripheryWeight2", 10.0);
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
            << "  sigma            : " << sigma_ << nl
            << "  coreWeight       : " << coreWeight_ << nl
            << "  peripheryWeight1 : " << peripheryWeight1_ << nl
            << "  peripheryWeight2 : " << peripheryWeight2_ << nl
            << "  result           : " << resultName_ << nl
            << endl;
    }

    DebugInfo
        << type() << " '" << name() << "': read configuration." << nl;

    return true;
}

wordList sasRefineIndicator::fields() const
{
    return requiredFields_;
}

bool sasRefineIndicator::execute()
{
    calcIndicator();
    return true;
}


bool sasRefineIndicator::write()
{
    Info << " functionObjects::" << type() << ' ' << name();
    Info << " writing field: " << resultName_ << endl;
    const auto& fld = mesh_.lookupObject<volScalarField>(resultName_);
    fld.write();

    return true;
}

} // End namespace functionObjects
} // End namespace Foam
