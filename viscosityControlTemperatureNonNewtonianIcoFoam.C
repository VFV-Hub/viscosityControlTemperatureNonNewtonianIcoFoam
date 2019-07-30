/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    nonNewtonianIcoFoam

Description
    Transient solver for incompressible, laminar flow of non-Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshNoClear.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
	volScalarField viscOld = fluid.nu();
	volScalarField viscNew = fluid.nu();
	double max=1, tol;
	int i,iMax;

	tol=1e-5;
	iMax=100;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"	
//inicio da edição
		i=0;
		max=1;
		while(max>tol && i<iMax)
		{
			max=0;
			viscOld = fluid.nu();

        	fluid.correct();

			viscNew = fluid.nu();
			forAll(mesh.cells(),celli)
			{
				if(max<=(mag(viscNew[celli]-viscOld[celli])))
				{
					max = mag(viscNew[celli]-viscOld[celli]);
				}
				//necessário para o primeiro passo de tempo
				if(max==0)
					max=1;		
			}
			Info<<"Variacao maxima da viscosidade = "<<max<<endl;

//final da edição

        	// Momentum predictor

        	fvVectorMatrix UEqn
        	(
            	fvm::ddt(U)
          	  + fvm::div(phi, U)
          	  - fvm::laplacian(fluid.nu(), U)
          	  - (fvc::grad(U) & fvc::grad(fluid.nu()))
        	);

        	if (piso.momentumPredictor())
        	{
            	solve(UEqn == -fvc::grad(p));
        	}


        	// --- PISO loop
        	while (piso.correct())
        	{
            	volScalarField rAU(1.0/UEqn.A());
            	volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            	surfaceScalarField phiHbyA
            	(
                	"phiHbyA",
                	fvc::flux(HbyA)
              	  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            	);

            	adjustPhi(phiHbyA, U, p);

            	// Update the pressure BCs to ensure flux consistency
            	constrainPressure(p, U, phiHbyA, rAU);

            	// Non-orthogonal pressure corrector loop
            	while (piso.correctNonOrthogonal())
            	{
                	// Pressure corrector

                	fvScalarMatrix pEqn
                	(
                    	fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                	);

                	pEqn.setReference(pRefCell, pRefValue);

                	pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                	if (piso.finalNonOrthogonalIter())
                	{
                    	phi = phiHbyA - pEqn.flux();
                	}
            	}

            	#include "continuityErrs.H"

            	U = HbyA - rAU*fvc::grad(p);
            	U.correctBoundaryConditions();	
        	}

//inicio da edição
			fvScalarMatrix TEqn
			(
				fvm::ddt(T)
			   +fvm::div(phi,T)
			   -fvm::laplacian(DT,T)
			);

			TEqn.solve();
			i++;
		}
		Info<<"numero de iteracoes por causa da viscosidade"<<i<<endl;
//final da edição	

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
