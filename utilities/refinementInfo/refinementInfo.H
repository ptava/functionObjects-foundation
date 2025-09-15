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

Class
    Foam::functionObjects::refinementInfo

Group
    grpUtilitiesFunctionObjects

Description
    Writes info associated with mesh refinement process.

    Operands:
    \table
      Operand      | Type                | Location
      input        | -                   | -
      output file  | dat                 <!--
               --> | postProcessing/\<FO\>/\<time\>/\<file\>(s)
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    solverInfo1
    {
        // Mandatory entries (unmodifiable)
        type            refinementInfo;
        libs            (myutilityFunctionObjects);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                          | Type | Reqd | Deflt
      type         | Type name: refinementInfo                | word | yes  | -
      libs         | Library name: utilityFunctionObjects | word | yes  | -
                   -->                                    | bool | no   | false
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link logFiles.H \endlink

See also
  - Foam::functionObject
  - Foam::functionObjects::fvMeshFunctionObject
  - Foam::functionObjects::logFiles

SourceFiles
    refinementInfo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_refinementInfo_H
#define functionObjects_refinementInfo_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class refinementInfo Declaration
\*---------------------------------------------------------------------------*/

class refinementInfo 
: 
    public fvMeshFunctionObject,
    public logFiles
{
protected:

    //- Whether the header has been written
    bool headerDone_;

    //- List of refinement process parameters
    wordList parameters_;

    //- Output file header informationS
    void writeFileHeader(const label i);

public:

    //- Runtime type information
    TypeName("refinementInfo");

    //- Constructor from Time and dictionary
    refinementInfo(
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~refinementInfo() = default;

    //- Read refinementInfo settings
    virtual bool read(const dictionary&);

    //- Return the list of fields required
    virtual wordList fields() const
    {
        return wordList::null();
    }

    //- Execute refinementInfo
    virtual bool execute();

    //- Write
    virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
