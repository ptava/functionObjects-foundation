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

Class
    Foam::functionObjects::sasRefineIndicator

Description
    Function object for generating a refinement indicator field for dynamic mesh
    refinement

    Indicator derived from the von Kármán length scale (Lvk) in the kOmegaSSTSAS
    turbulence model formulation.

    Options are provided to enphasise different regions of the flow:

    1. focusRegion == core
        - Procecess Lvk in order to obtain the refinement field that leads to
        refine where the grid filter is active, ideally in the core of detected
        eddies.
        - Transfer function to mark cells: using a constant value function or an
        odd, monotonic, sign-preserving function that behaves linear near 0
        (slope 1) and approaches a scaled line (slope 1+coreWeight) for large
        ∣d|, it uses a Gaussian factor to modulate the slope (markWithConstant
        or markWithOddScaler)

    2. FocusRegion == periphery
        - Procecess nLvk in order to obtain the refinement field that leads to
        refine in the periphery of detected eddies.
        - Transfer function to mark cells: is gaussian-like function centered in
        nLvk=1 minus a quadratic term (markPeripheryGaussSink)

    3. focusRegion == combined
        - Combine the two above

    where:
        - $L_{vk} = max(c_1, c_1)$ where:
        - $c_1$ is the actual length scale;
        - $c_2$ is the high wave number number damper linked to the grid size.
        - $nLvk = Lvk / c2$

    Designed for use with dynamicRefineFvMesh to adapt the mesh around turbulent
    structures.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    sasRefineIndicator.C

\*---------------------------------------------------------------------------*/
#ifndef functionObjects_sasRefineIndicator_H
#define functionObjects_sasRefineIndicator_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "Switch.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class sasRefineIndicator
\*---------------------------------------------------------------------------*/

class sasRefineIndicator
:
    public fvMeshFunctionObject
{
public:

    // Focus region selection
    enum focusRegion
    {
        core,
        periphery,
        combined
    };
    static const NamedEnum<focusRegion, 3> focusRegionNames_;

    enum transferFunction
    {
        constant,
        oddScaler
    };
    static const NamedEnum<transferFunction, 2> transferFunctionNames_;

protected:

    // ---- User options
    focusRegion      focusRegion_;
    scalar           sigma_;
    scalar           coreWeight_;
    scalar           peripheryWeight1_;
    scalar           peripheryWeight2_;
    transferFunction transferFunction_;
    word             resultName_;
    wordList         requiredFields_;

    // ---- Region definition 
    // TODO: address 'momentumError' procedure beacause isn't working
    word regionName_;

    // ---- Created field names (scoped)
    word meshFieldName_;   // field held on the main mesh
    word subFieldName_;    // field on the sub-mesh (if any)

    // ---- Helpers

    //- Transform field via constant value
    tmp<volScalarField::Internal> markCoreConstant
    (
        const volScalarField::Internal& c1,
        const volScalarField::Internal& c2,
        const scalar coreWeight
    ) const;

    //- Transform field via gaussian function
    tmp<volScalarField::Internal> markCoreOddScaler
    (
        const volScalarField::Internal& c1,
        const volScalarField::Internal& c2,
        const scalar coreWeight,
        const scalar sigma
    ) const;

    //- Transform field via gaussian function
    tmp<volScalarField::Internal> markPeripheryGaussSink
    (
        const volScalarField::Internal& Lvk,
        const volScalarField::Internal& c2,
        const scalar peripheryWeight1,
        const scalar peripheryWeight2,
        const scalar sigma
    ) const;

    //- Compute the indicator (writes into created fields)
    void calcIndicator();

public:

    TypeName("sasRefineIndicator");

    // Constructors

    sasRefineIndicator
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Read dictionary, with deprecations handled
    virtual bool read(const dictionary&);

    //- Return the list of fields required
    virtual wordList fields() const;

    //- Execute: compute
    virtual bool execute();

    //- Write
    virtual bool write();

    //- No-op end
    virtual bool end() { return true; }
};

} // End namespace functionObjects
} // End namespace Foam

#endif
