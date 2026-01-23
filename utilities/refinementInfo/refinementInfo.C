/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "refinementInfo.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(refinementInfo, 0);
    addToRunTimeSelectionTable(functionObject, refinementInfo, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::refinementInfo::writeFileHeader(const label i)
{
    if (headerDone_)
    {
        writeHeader(file(), "===");
    }
    else
    {
        writeHeader(file(), "Refinement Info");
    }

    writeCommented(file(), "Time");

    for (const word& name : parameters_)
    {
        writeTabbed(file(), name);
    }

    file() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::refinementInfo::refinementInfo(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    headerDone_(false),
    parameters_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::refinementInfo::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("refinementParameters") >> parameters_;

    resetName(typeName);

    return true;
}


bool Foam::functionObjects::refinementInfo::execute()
{
    return true;
}

bool Foam::functionObjects::refinementInfo::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());

        for (const word& name : parameters_)
        {
            if (obr_.foundObject<uniformDimensionedScalarField>(name))
            {
                const auto& parameter =
                    obr_.lookupObject<uniformDimensionedScalarField>(name);
                const scalar& value = parameter.value();
                if (std::isnan(value))
                {
                    writeTabbed(file(), "N/A");
                }
                else
                {
                    writeTabbed(file(), Foam::name(value));
                }
            }
            else
            {
                writeTabbed(file(), "N/A");
            }
        }

        file() << endl;
    }

    return true;
}

// ************************************************************************* //
