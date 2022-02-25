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

\*---------------------------------------------------------------------------*/

#include "LuoSvendsen.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{
    defineTypeNameAndDebug(LuoSvendsen, 0);
    addToRunTimeSelectionTable
    (
        breakupModel,
        LuoSvendsen,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::breakupModels::LuoSvendsen::
LuoSvendsen
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    breakupModel(popBal, dict),
    c1_(0.9238),
	c2_(2.047),
	tbreakupRate_
	(
		IOobject
		(
			"tbreakupRate",
			popBal.time().timeName(),
			popBal.mesh()
		),
		popBal.mesh(),
		dimensionedScalar("tbreakupRate", dimless, scalar(0))
	),
	tfuncs_
	(
		IOobject
		(
			"tfuncs",
			popBal.time().timeName(),
			popBal.mesh()
		),
		popBal.mesh(),
		dimensionedScalar("tfuncs", dimless, scalar(0))
	),
	rhoc_(popBal.continuousPhase().rho()),
	epsilon_(popBal.continuousTurbulence().epsilon()),
	sigma_("sigma",dimensionSet(1,0,-2,0,0,0,0),readScalar(dict.lookup("sigma"))),
	multiplier_("multiplier",dimensionSet(0,0,0,0,0,0,0),readScalar(dict.lookup("multiplier")))
//	sigma_(popBal.sigmaWithContinuousPhase((*popBal.sizeGroups()[0]).phase()))
{

	nodes_[0] = 0.9061798459;
	nodes_[1] = -0.9061798459;
	nodes_[2] = 0.5384693101;
	nodes_[3] = -0.5384693101;
	nodes_[4] = 0.;
	
	A_[0] = 0.2369268851;
	A_[1] = 0.2369268851;
	A_[2] = 0.4786286705;
	A_[3] = 0.4786286705;
	A_[4] = 0.5688888889;

}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::breakupModels::LuoSvendsen::setBreakupRate
(
    volScalarField& breakupRate,
    const label i
)
{

//    dimensionedScalar epsilonAverage(epsilon_.weightedAverage(popBal_.mesh().V()));

    const sizeGroup& fi = *popBal_.sizeGroups()[i];
	
	const scalar a = 0.;
	const scalar b = 1.;
	const scalar c = 0.;
	const scalar d = 1.;
	
	tbreakupRate_ *= 0.;
	
	for (int ii=0; ii<5; ii++)
	{
		for (int jj=0; jj<5; jj++)
		{
			tbreakupRate_ += A_[ii]*A_[jj]*funcs(i,((b-a)*nodes_[ii]+a+b)/2,((d-c)*nodes_[jj]+c+d)/2);
		}
			
	}
	
	tbreakupRate_ *= (b-a)*(d-c)/4;
	
	breakupRate = 0.5*c1_*(1-fi.phase())*pow(epsilon_/fi.d()/fi.d(), 1./3)*tbreakupRate_*multiplier_;

//////////////////////////////////////////////////////////////////////////////////////////////////


forAll(breakupRate,cellI)
{
breakupRate[cellI] = breakupRate[cellI] * neg(breakupRate.mesh().C()[cellI].z()-0.44);
}


///////////////////////////////////////////////////////////////////////////////////////////////

}



const Foam::volScalarField&
Foam::diameterModels::breakupModels::LuoSvendsen::funcs
(
	const label i,
    const scalar f,
	const scalar rxi
)
{
//    dimensionedScalar epsilonAverage(epsilon_.weightedAverage(popBal_.mesh().V()));

    const dimensionedScalar& di = popBal_.sizeGroups()[i]->d();
	const scalar cf = pow(f, 2./3) + pow(1-f, 2./3) - 1;
	
	tfuncs_ = pow(1+rxi, 2)*pow(rxi, -11./3)*
	exp
	(
		-12*cf*sigma_/c2_/rhoc_/pow(
//						min(epsilon_,10*epsilonAverage), 2./3
						epsilon_, 2./3
					   )/pow(di, 5./3)/pow(rxi, 11./3)
	);
	
	return tfuncs_;

}

// ************************************************************************* //


