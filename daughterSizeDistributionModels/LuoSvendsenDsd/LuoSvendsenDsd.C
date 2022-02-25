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

\*---------------------------------------------------------------------------*/


#include "LuoSvendsenDsd.H"
#include "addToRunTimeSelectionTable.H"
#include "breakupModel.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionModels
{
    defineTypeNameAndDebug(LuoSvendsenDsd, 0);
    addToRunTimeSelectionTable
    (
        daughterSizeDistributionModel,
        LuoSvendsenDsd,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::LuoSvendsenDsd::
LuoSvendsenDsd
(
    const breakupModel& breakup,
    const dictionary& dict
)
:
    daughterSizeDistributionModel(breakup, dict),
	c1_(0.9238),
	c2_(2.047),
	delta_("delta",dimensionSet(0,0,0,0,0,0,0), scalar(1e-15)),
	mesh_(breakup.popBal().mesh()),
	rhoc_(breakup.popBal().continuousPhase().rho()),
	epsilon_(breakup.popBal().continuousTurbulence().epsilon()),
	rhocAverage_
	(
		rhoc_.weightedAverage(mesh_.V())
	),
	epsilonAverage_
	(
		epsilon_.weightedAverage(mesh_.V())
	),
	sigma_("sigma",dimensionSet(1,0,-2,0,0,0,0),readScalar(dict.lookup("sigma")))
//	sigma_("sigma",dimensionSet(1,0,-2,0,0,0,0), scalar(0.07))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::LuoSvendsenDsd::
~LuoSvendsenDsd()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::LuoSvendsenDsd::calcNik
(
    const label i,
    const label k
) const
{
    const dimensionedScalar& xi = breakup_.popBal().sizeGroups()[i]->x();
    const dimensionedScalar& xk = breakup_.popBal().sizeGroups()[k]->x();
    const List<sizeGroup*>& sizeGroups = breakup_.popBal().sizeGroups();

    if (i == 0)
    {
        return (sizeGroups[i+1]->x() - xi)/xk;
    }
    else if (i == k)
    {
        return (xi - sizeGroups[i-1]->x())/xk;
    }
    else
    {
        return (sizeGroups[i+1]->x() - xi)/xk + (xi - sizeGroups[i-1]->x())/xk;
    }
}




Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::LuoSvendsenDsd::calcNik2
(
    const label i,
    const label k
) const
{
    const dimensionedScalar& xi = breakup_.popBal().sizeGroups()[i]->x().value();
    const List<sizeGroup*>& sizeGroups = breakup_.popBal().sizeGroups();
	
	const scalar moment1 = breakup_.popBal().Moment1();
	const scalar moment2 = breakup_.popBal().Moment2();
	
	if (i == 0 && k == 0)
	{
		return 0;
	}
	
    else if (i == 0 && k !=0)
    {
		const dimensionedScalar& xi1 = sizeGroups[i+1]->x().value();

        return (Bik(i,k,moment2)*pow(xi1, moment1)-Bik(i,k,moment1)*pow(xi1, moment2))
		/(pow(xi, moment2)*pow(xi1, moment1)-pow(xi, moment1)*pow(xi1, moment2));
    }
    else if (i == k)
    {
		const dimensionedScalar& xi0 = sizeGroups[i-1]->x().value();
		
        return (Bik(i-1,k,moment2)*pow(xi0, moment1)-Bik(i-1,k,moment1)*pow(xi0, moment2))
		/(pow(xi, moment2)*pow(xi0, moment1)-pow(xi, moment1)*pow(xi0, moment2));
    }
    else
    {
		const dimensionedScalar& xi1 = sizeGroups[i+1]->x().value();
		const dimensionedScalar& xi0 = sizeGroups[i-1]->x().value();
		
        return (Bik(i,k,moment2)*pow(xi1, moment1)-Bik(i,k,moment1)*pow(xi1, moment2))
		/(pow(xi, moment2)*pow(xi1, moment1)-pow(xi, moment1)*pow(xi1, moment2))
		+
		(Bik(i-1,k,moment2)*pow(xi0, moment1)-Bik(i-1,k,moment1)*pow(xi0, moment2))
		/(pow(xi, moment2)*pow(xi0, moment1)-pow(xi, moment1)*pow(xi0, moment2));
    }
}



const Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::LuoSvendsenDsd::Bik
(
    const label i,
    const label k,
	const scalar moment
) const
{
    const scalar& xi = breakup_.popBal().sizeGroups()[i]->x().value();
	const scalar& xi1 = breakup_.popBal().sizeGroups()[i+1]->x().value();
    const scalar& xk = breakup_.popBal().sizeGroups()[k]->x().value();

	const scalar a1 = xi/xk;
	const scalar b1 = xi1/xk;
	const scalar c1 = 0.;
	const scalar d1 = 1.;
	
	const scalar a2 = 0.;
	const scalar b2 = 1.;
	const scalar c2 = 0.;
	const scalar d2 = 1.;
	

	dimensionedScalar tBik1 = 0.;
	dimensionedScalar tBik2 = 0.;
	
	for (int ii=0; ii<5; ii++)
	{
		for (int jj=0; jj<5; jj++)
		{

			tBik1 += A_[ii]*A_[jj]*funcs1(k,((b1-a1)*nodes_[ii]+a1+b1)/2,((d1-c1)*nodes_[jj]+c1+d1)/2,moment);

			tBik2 += A_[ii]*A_[jj]*funcs2(k,((b2-a2)*nodes_[ii]+a2+b2)/2,((d2-c2)*nodes_[jj]+c2+d2)/2);

		}
			
	}
	
	tBik1 *= (b1-a1)*(d1-c1)/4;
//	Info << "tBik2 = " << tBik2*(b2-a2)*(d2-c2)/4 <<endl;
//	tBik2 = max(tBik2*(b2-a2)*(d2-c2)/4, delta_);
	tBik2 *= (b2-a2)*(d2-c2)/4;
//	Info << "tBik1 = " << tBik1 <<endl;

	if (tBik2.value()==0)	return 0;
	else	return 2*pow(xk, moment)*tBik1/tBik2;
	
}

const Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::LuoSvendsenDsd::funcs1
(
	const label k,
    const scalar f,
	const scalar rxi,
	const scalar moment
) const
{
    const dimensionedScalar& dk = breakup_.popBal().sizeGroups()[k]->d();
	const scalar cf = pow(f, 2./3) + pow(1-f, 2./3) - 1;
	const dimensionedScalar Chi = 12*cf*sigma_/c2_/rhocAverage_/pow(epsilonAverage_, 2./3)/pow(dk, 5./3)/pow(rxi, 11./3);

	return pow(f, moment)*pow(1+rxi, 2)*pow(rxi, -11./3)*exp(-Chi);
	
}


const Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::LuoSvendsenDsd::funcs2
(
    const label k,
    const scalar f,
    const scalar rxi
) const
{
    const dimensionedScalar& dk = breakup_.popBal().sizeGroups()[k]->d();
	const scalar cf = pow(f, 2./3) + pow(1-f, 2./3) - 1;
	const dimensionedScalar Chi = 12*cf*sigma_/c2_/rhocAverage_/pow(epsilonAverage_, 2./3)/pow(dk, 5./3)/pow(rxi, 11./3);
	
	return pow(1+rxi, 2)*pow(rxi, -11./3)*exp(-Chi);
	
}




// ************************************************************************* //
