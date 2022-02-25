/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::diameterModels::coalescenceModels::Bhole

Description
    Modified Model of Bhole. The coalescence rate is calculated by


    References:
    \verbatim
        "Bubble Coalescence and Break-Up in Air-Sparged Bubble-Columns"
        Prince, M.J., Blanch, H.W.
        AIChE Journal, Vol. 36, Issue 10, 1990, pp. 1485-1499
    \endverbatim

	
	\verbatim
        Bhole, M.R.; Joshi, J.B.; Ramkrishna, D.
		"CFD simulation of bubble columns incorporating population balance modeling"
        Chem. Eng. Sci., Vol. 63, 2008, pp. 2267-2282
    \endverbatim

	
	
	\verbatim
        Shang, X.P.; Ng, B.F.; Wan, M.P.; Ding S.R. 
        Investigation of CFD-PBM simulations based on fixed pivot method:
		Influence of the moment closure,
        The Chemical Engineering Journal
		DOI: 10.1016/j.cej.2019.122882
    \endverbatim

	
Usage
    \table
        Property     | Description             | Required    | Default value
        C1           | Coefficient C1          | no          | 0.089
    \endtable

SourceFiles
    Bhole.C

\*---------------------------------------------------------------------------*/

#ifndef Bhole_H
#define Bhole_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                        Class Bhole Declaration
\*---------------------------------------------------------------------------*/

class Bhole
:
    public coalescenceModel
{
    // Private data

        //- Optional coefficient C1, defaults to 0.089
        dimensionedScalar C1_;

        //- Initial film thickness, defaults to 1e-4m
        dimensionedScalar h0_;

        //- Critical film thickness, defaults to 1e-8m
        dimensionedScalar hf_;

        //- Switch for considering turbulent collisions
        Switch turbulentCollisions_;

        //- Switch for considering buoyancy-induced collisions
        Switch buoyantCollisions_;

        //- Switch for considering buoyancy-induced collisions
        Switch laminarShearCollisions_;


public:

    //- Runtime type information
    TypeName("Bhole");

    // Constructor

        Bhole
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Bhole()
    {}


    // Member Functions

        //- Add to coalescenceRate
        virtual void addToCoalescenceRate
        (
            volScalarField& coalescenceRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
