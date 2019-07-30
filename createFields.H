//começo da edição


IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading DT in transportProperties\n" << endl;

dimensionedScalar DT
(
    "DT",
    dimViscosity,
    transportProperties.lookup("DT")
);

Info<< "Reading iMax in transportProperties\n" << endl;

dimensionedScalar iMax
(
    "iMax", 
	dimensionSet(0,0,0,0,0,0,0), 
	transportProperties.lookup("iMax")
);

Info<< "Reading tol in transportProperties\n" << endl;

dimensionedScalar tol
(
    "tol",
	dimensionSet(0,0,0,0,0,0,0),
    transportProperties.lookup("tol")
);

//final da edição

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//começo edição
Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//final edição

#include "createPhi.H"


singlePhaseTransportModel fluid(U, phi);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
