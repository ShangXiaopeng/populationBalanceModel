/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::diameterModels::daughterSizeDistributionModels::LuoSvendsenDsd

Description
    Uniform binary daughter size distribution. Used for verification and
    validation of the breakup formulation implemented in the
    populationBalanceModel class.

SourceFiles
    LuoSvendsenDsd.C

\*---------------------------------------------------------------------------*/

#ifndef LuoSvendsenDsd_H
#define LuoSvendsenDsd_H

#include "daughterSizeDistributionModel.H"
#include "breakupModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionModels
{

/*---------------------------------------------------------------------------*\
                      Class LuoSvendsenDsd Declaration
\*---------------------------------------------------------------------------*/

class LuoSvendsenDsd
:
    public daughterSizeDistributionModel
{
	
	const scalar c1_;
	const scalar c2_;
	const dimensionedScalar delta_;
	const fvMesh& mesh_;
	const volScalarField rhoc_;
	const volScalarField epsilon_;
	const dimensionedScalar rhocAverage_;
	const dimensionedScalar epsilonAverage_;
	const dimensionedScalar sigma_;
	
public:

    //- Runtime type information
    TypeName("LuoSvendsenDsd");


    // Constructor

        LuoSvendsenDsd
        (
            const breakupModel& breakup,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LuoSvendsenDsd();


    // Member Functions

        //- Return total number of particles assigned to class i when a particle
        //  of class k breaks
        virtual dimensionedScalar calcNik(const label i, const label k) const;
		virtual dimensionedScalar calcNik2(const label i, const label k) const;		
		const dimensionedScalar Bik(const label i, const label k, scalar moment) const;
		const dimensionedScalar Funcs(const label i, const label k, const label jj, const scalar moment) const;
		const dimensionedScalar funcs1(const label k,  const scalar f, const scalar rxi, const scalar moment) const;
		const dimensionedScalar funcs2(const label k,  const scalar f, const scalar rxi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterSizeDistributionModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
