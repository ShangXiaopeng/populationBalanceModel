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

#include "uniformBinary.H"
#include "addToRunTimeSelectionTable.H"
#include "breakupModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionModels
{
    defineTypeNameAndDebug(uniformBinary, 0);
    addToRunTimeSelectionTable
    (
        daughterSizeDistributionModel,
        uniformBinary,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::uniformBinary::
uniformBinary
(
    const breakupModel& breakup,
    const dictionary& dict
)
:
    daughterSizeDistributionModel(breakup, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::uniformBinary::
~uniformBinary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::uniformBinary::calcNik
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
Foam::diameterModels::daughterSizeDistributionModels::uniformBinary::calcNik2
(
    const label i,
    const label k
) const
{
    const dimensionedScalar& xi = breakup_.popBal().sizeGroups()[i]->x().value();
//    const dimensionedScalar& xk = breakup_.popBal().sizeGroups()[k]->x();
    const List<sizeGroup*>& sizeGroups = breakup_.popBal().sizeGroups();

	
	const scalar moment1 = breakup_.popBal().Moment1();
	const scalar moment2 = breakup_.popBal().Moment2();
	
    if (i == 0)
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
Foam::diameterModels::daughterSizeDistributionModels::uniformBinary::Bik
(
    const label i,
    const label k,
	const scalar moment
) const
{
    const dimensionedScalar& xi = breakup_.popBal().sizeGroups()[i]->x();
	const dimensionedScalar& xi1 = breakup_.popBal().sizeGroups()[i+1]->x();
//    const dimensionedScalar& xk = breakup_.popBal().sizeGroups()[k]->x();
	
	dimensionedScalar tBik = 0.;
	
	for (int jj=0; jj<5; jj++)
	{
		tBik += A_[jj]*Funcs(i,k,jj,moment);
			
	}
	
	return (xi1.value()-xi.value())/2.*tBik;
	
}
	
const Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::uniformBinary::Funcs
(
    const label i,
    const label k,
	const label jj,
	const scalar moment
) const
{
	const dimensionedScalar& a = breakup_.popBal().sizeGroups()[i]->x();
	const dimensionedScalar& b = breakup_.popBal().sizeGroups()[i+1]->x();
    const dimensionedScalar& xk = breakup_.popBal().sizeGroups()[k]->x();
	
	return 2./xk.value()*pow((a.value()+b.value())/2.+(b.value()-a.value())/2.*nodes_[jj], moment);
}






// ************************************************************************* //
