Info<< "Reading field C1\n" << endl;

volScalarField C1
(
    IOobject
    (
        "C1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field C2\n" << endl;
volScalarField C2
(
    IOobject
    (
        "C2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field Fi\n" << endl;

volScalarField Fi
(
    IOobject
    (
        "Fi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField vgradFi
(
    IOobject
    (
        "vgradFi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(Fi)
);

Info<< "Reading properties\n" << endl;

IOdictionary properties
(
    IOobject
    (
        "properties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar DC1("DC1", dimensionSet(0, 2, -1, 0, 0, 0, 0), properties);

dimensionedScalar DC2("DC2", dimensionSet(0, 2, -1, 0, 0, 0, 0), properties);

dimensionedScalar z1("z1", dimensionSet(0, 0, 0, 0, 0, 0, 0), properties);

dimensionedScalar z2("z2", dimensionSet(0, 0, 0, 0, 0, 0, 0), properties);

dimensionedScalar F("F", dimensionSet(0, 0, 1, 0, -1, 1, 0), properties);

dimensionedScalar T("T", dimensionSet(0, 0, 0, 1, 0, 0, 0), properties);

dimensionedScalar R("R", dimensionSet(1, 2, -2, -1, -1, 0, 0), properties);

dimensionedScalar mu1("mu1", dimensionSet(-1, 0, 1, 0, 1, 0, 0), properties);

dimensionedScalar mu2("mu2", dimensionSet(-1, 0, 1, 0, 1, 0, 0), properties);

dimensionedScalar kap("kap", dimensionSet(-1, -3, 3, 0, 0, 2, 0), properties); 

dimensionedScalar Gamma("Gamma", dimensionSet(0, 0, 0, 0, 0, 0, 0), properties);

dimensionedScalar p("p", dimensionSet(0, 0, 0, 0, 0, 0, 0), properties);


volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    kap+F*F*(z1*z1*Gamma*mu1*C1+z2*z2*Gamma*mu2*C2)
);

#include "createFvOptions.H"