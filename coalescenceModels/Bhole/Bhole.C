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

#include "Bhole.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "phaseCompressibleTurbulenceModel.H"

#include "orderedPhasePair.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(Bhole, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        Bhole,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::Bhole::
Bhole
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    C1_("C1", dimless, dict.lookupOrDefault<scalar>("C1", 0.356)),
    h0_("h0", dimLength, dict.lookupOrDefault<scalar>("h0", 1e-4)),
    hf_("hf", dimLength, dict.lookupOrDefault<scalar>("h0", 1e-8)),
    turbulentCollisions_(dict.lookup("turbulentCollisions")),
    buoyantCollisions_(dict.lookup("buoyantCollisions")),
    laminarShearCollisions_(dict.lookup("laminarShearCollisions"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceModels::Bhole::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const label i,
    const label j
)
{

    const sizeGroup& fi = *popBal_.sizeGroups()[i];
    const sizeGroup& fj = *popBal_.sizeGroups()[j];
    const uniformDimensionedVectorField& g =
        popBal_.mesh().lookupObject<uniformDimensionedVectorField>("g");


    const phaseModel& dispersedPhase = fi.phase();
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const scalar Cmu(0.27);
    const scalar Cvm(0.5);
    const dimensionedScalar SMALL("SMALL",dimensionSet(0,1,-1,0,0,0,0),scalar(1e-10));
    const orderedPhasePair pair(dispersedPhase,continuousPhase);
    volScalarField epsilon = popBal_.continuousTurbulence().epsilon();
    volScalarField k = popBal_.continuousTurbulence().k();

//    dimensionedScalar epsilonAverage(epsilon.weightedAverage(popBal_.mesh().V()));


    volScalarField UB = mag(pair.Ur());
    volScalarField Re(pair.Re());
    volScalarField Eo(pair.Eo());


    volScalarField CdRe
    (

        max
        (
           24.0
           *min
            (
                (1.0 + 0.15*pow(Re, 0.687)),
                scalar(3.0)
            ),
            8.0*Eo/(3.0*Eo + 12.0)*Re
        )
    );



    volScalarField StUBi
    (

        4.0/3.0*Cvm/CdRe*fi.d()*Re
        /(Cmu*k/epsilon)

    );

    volScalarField StUBj
    (

        4.0/3.0*Cvm/CdRe*fj.d()*Re
        /(Cmu*k/epsilon)
    );




    dimensionedScalar rij(1.0/(1.0/fi.d() + 1.0/fj.d()));

    const volScalarField collisionEfficiency
    (
        exp
        (
          - sqrt
            (
                pow3(rij)*continuousPhase.rho()
               /(16.0*popBal_.sigmaWithContinuousPhase(fi.phase()))
            )
           *log(h0_/hf_)
           *cbrt(epsilon)/pow(rij, 2.0/3.0)
        )
    );

    if (turbulentCollisions_)
    {
        coalescenceRate +=
            (
                C1_*pi*sqr(fi.d() + fj.d())
               *cbrt(epsilon)
               *sqrt(pow(fi.d(), 2.0/3.0)*(UB/max(SMALL,(UB+StUBi))) + pow(fj.d(), 2.0/3.0))*(UB/max(SMALL,(UB+StUBj)))
            )
           *collisionEfficiency;

/////////////////////////////////////////////////////////////////////////////////////////////////////////


forAll(coalescenceRate,cellI)
{
coalescenceRate[cellI] = coalescenceRate[cellI] * neg(coalescenceRate.mesh().C()[cellI].z()-0.44);
}


///////////////////////////////////////////////////////////////////////////

    }

    if (buoyantCollisions_)
    {
        dimensionedScalar Sij(pi/4.0*sqr(fi.d() + fj.d()));

        coalescenceRate +=
            (
                Sij
               *mag
                (
                    sqrt
                    (
                        2.14*popBal_.sigmaWithContinuousPhase(fi.phase())
                       /(continuousPhase.rho()*fi.d()) + 0.505*mag(g)*fi.d()
                    )
                  - sqrt
                    (
                        2.14*popBal_.sigmaWithContinuousPhase(fi.phase())
                       /(continuousPhase.rho()*fj.d()) + 0.505*mag(g)*fj.d()
                    )
                )
            )
           *collisionEfficiency;
    }

    if (laminarShearCollisions_)
    {
        FatalErrorInFunction
            << "Laminar shear collision contribution not implemented for "
            << this->type() << " coalescence model."
            << exit(FatalError);
    }
}


// ************************************************************************* //
