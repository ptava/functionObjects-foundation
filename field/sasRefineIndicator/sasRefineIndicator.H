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

    Refinement indicator based on the von Kármán length scale (Lvk) and is
    designed to serve as an adaptive refinement criterion for a Scale Adaptive
    (SA) turbulente model.

    The porpose is to adaptively refine mesh in regions where vortices are being
    detected by a scale-resolving simulation, refining in vortices-associated
    regions and coarsening outside them.

    Options are provided to enphasise different regions of the flow/vortices:

    1. focusRegion == "core"

        Focuses refinement where the grid filter is most active, ideally in the
        core of detected eddies.

        Available transfer functions:

        a) markCoreConstant:
        apply constant value to all cells where d = c2 - c1 > 0
        \f[
            f(d) = \text{const}
        \f]

        b) markCoreOddScaler:
        apply an odd, monotonic, sign-preserving function, symmetric around
        origin, normalised by its maximum value:
        \f[
          f(d) = d \cdot \left( c - \exp\left( -\frac{d^2}{2 \sigma^2} \right) \right), \quad c \geq 1.0
        \f]
        Notes:
        - As \f$ |d| \to \infty \f$, the exponential vanishes, and
        \f$ f(d) \sim c \cdot d \f$
        - If \f$ c = 1.0 \f$, function behaves cubically near 0 — useful
        for smoothly marking the interface region
        - If \f$ c > 1.0 \f$, the function becomes more linear near 0
        - \f$ \sigma \in [0.0, 1.0] \f$ controls the width of the transition region

        c) markCoreGaussSink:
        gaussian-like function peaking where nLvk = Lvk / c2 = 1, and decays for
        increasing values of nLvk:
        \f[
            f(nLvk) = w_1 \cdot \exp\left( -\frac{1}{2} \left( \frac{nLvk - 1}{\sigma} \right)^2 \right) + w_2 \cdot (nLvk - 1)^2
        \f]
        Parameters:
        - \f$ w_1, w_2 \geq 0 \f$
        - \f$ \sigma \f$ controls the width of the peak
        - \f$ w_1 \f$ controls peak height
        - \f$ w_2 \f$ controls decay rate outside the focus region

2. FocusRegion == "periphery"

        Focuses refinement in the periphery of vortex structures rather than
        their core.

        Available transfer functions:

        d) markPeripheryGaussSink:
        Gaussian-like function peaking at a user-defined reference value
        of Lvk (e.g., a free-stream value). Same as 'c' but different nLvk
        definition:
        \f[
            f(nLvk) = w_1 \cdot \exp\left( -\frac{1}{2} \left( \frac{nLvk - 1}{\sigma} \right)^2 \right) - w_2 \cdot (nLvk - 1)^2
        \f]

    Where:
    - \f$ L_{vk} = \max(c_1, c_2) \f$
    - \f$ c_1 \f$ is the physical turbulence length scale
    - \f$ c_2 \f$ is the grid-related high wave number damper (linked to grid size)

Usage
    Example of function object specification:
    \verbatim
    sasRefineIndicator1
    {
        libs            ( "libmyfieldFunctionObjects.so" );
        type            sasRefineIndicator;
        result          refineIndicatorField;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
        executeAtStart  false;
        regionName      refine_cells;
        focusRegion     core;
        transferFunction oddScaler;
        sigma           1.0;
    }
    \endverbatim

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
#include "fvCellZone.H"

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
        periphery
    };
    static const NamedEnum<focusRegion, 2> focusRegionNames_;

    enum transferFunction
    {
        basic,
        constant,
        oddScaler,
        safeScaler,
        gaussSink
    };
    static const NamedEnum<transferFunction, 5> transferFunctionNames_;

protected:

    // ---- Container for mathematically derived coefficients
    struct optimizationCoeffs
    {
        scalar tStar;
        scalar Khat;
    };

    // ---- User options
    fvCellZone       zone_;
    focusRegion      focusRegion_;
    transferFunction transferFunction_;
    scalar           sigma_;
    scalar           weight1_;
    scalar           weight2_;
    scalar           alpha_;
    scalar           LvkRef_;
    bool             controlRunTime_;
    word             resultName_;

    wordList         requiredFields_;

    // ---- Created field names (scoped)
    word meshFieldName_;   // field held on the main mesh
    word subFieldName_;    // field on the sub-mesh (if any)

    // ---- Helpers

    //- Create an Internal volScalarField (tmp)
    tmp<volScalarField::Internal> makeTmpInternal
    (
        const word& name,
        const dimensionedScalar& init
    ) const;

    //- Calculate or retrieve cached optimization coefficients (tStar, Khat)
    optimizationCoeffs getOptimizationCoeffs
    (
        const scalar weight,
        const scalar alpha
    ) const;

    //- Transform field without modification
    tmp<volScalarField::Internal> markCoreBasic
    (
        const labelList& cells,
        const volScalarField::Internal& c1,
        const volScalarField::Internal& c2
    ) const;

    //- Transform field via constant value
    tmp<volScalarField::Internal> markCoreConstant
    (
        const labelList& cells,
        const volScalarField::Internal& c1,
        const volScalarField::Internal& c2,
        const scalar weight
    ) const;

    //- Transform field via gaussian function
    tmp<volScalarField::Internal> markCoreOddScaler
    (
        const labelList& cells,
        const volScalarField::Internal& c1,
        const volScalarField::Internal& c2,
        const scalar weight,
        const scalar sigma
    ) const;

    //- Transform field via gaussian function with dampening term
    tmp<volScalarField::Internal> markCoreSafeScaler
    (
        const labelList& cells,
        const volScalarField::Internal& c1,
        const volScalarField::Internal& c2,
        const scalar weight,
        const scalar sigma,
        const scalar alpha
    ) const;

    //- Transform field via gaussian function
    tmp<volScalarField::Internal> markCoreGaussSink
    (
        const labelList& cells,
        const volScalarField::Internal& Lvk,
        const volScalarField::Internal& c2,
        const scalar weight,
        const scalar sigma
    ) const;

    tmp<volScalarField::Internal> markPeripheryGaussSink
    (
        const labelList& cells,
        const volScalarField::Internal& Lvk,
        const scalar LvkRef,
        const scalar weight,
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
