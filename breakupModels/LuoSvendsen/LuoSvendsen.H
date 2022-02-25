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
    Foam::diameterModels::breakupModels::LuoSvendsen

Description
    Model of Luo-Svendsen (1996). 

    References:
    \verbatim
		Luo, H.; Svendsen, H.F.
        "Theoretical model for drop and bubble breakup in turbulent dispersions"
        AIChE. J., Vol. 42, 1996, pp. 1225-1233
    \endverbatim

Usage
   
SourceFiles
    LuoSvendsen.C

\*---------------------------------------------------------------------------*/

#ifndef LuoSvendsen_H
#define LuoSvendsen_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                  Class LuoSvendsen Declaration
\*---------------------------------------------------------------------------*/

class LuoSvendsen
:
    public breakupModel
{
    // Private data

	const scalar c1_;
	const scalar c2_;
	volScalarField tbreakupRate_;
	volScalarField tfuncs_;
	const volScalarField rhoc_;
	const volScalarField epsilon_;
	const dimensionedScalar sigma_;
	const dimensionedScalar multiplier_;
	scalar nodes_[5];
	scalar A_[5];
	


public:

    //- Runtime type information
    TypeName("LuoSvendsen");

    // Constructor

        LuoSvendsen
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LuoSvendsen()
    {}


    // Member Functions

        //- Set total breakupRate
        virtual void setBreakupRate
        (
            volScalarField& breakupRate,
            const label i
        );
		
		const volScalarField& funcs(const label i, const scalar f, const scalar rxi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
