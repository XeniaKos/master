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

dimensionedScalar DC3("DC3", dimensionSet(0, 2, -1, 0, 0, 0, 0), properties);

dimensionedScalar F("F", dimensionSet(0, 0, 1, 0, -1, 1, 0), properties);

dimensionedScalar T("T", dimensionSet(0, 0, 0, 1, 0, 0, 0), properties);

dimensionedScalar R("R", dimensionSet(1, 2, -2, -1, -1, 0, 0), properties);

dimensionedScalar mu1("mu1", dimensionSet(-1, 0, 1, 0, 1, 0, 0), properties);

dimensionedScalar mu2("mu2", dimensionSet(-1, 0, 1, 0, 1, 0, 0), properties);

dimensionedScalar mu3("mu3", dimensionSet(-1, 0, 1, 0, 1, 0, 0), properties);

dimensionedScalar z1("z1", dimensionSet(0, 0, 0, 0, 0, 0, 0), properties);

dimensionedScalar z2("z2", dimensionSet(0, 0, 0, 0, 0, 0, 0), properties);

dimensionedScalar z3("z3", dimensionSet(0, 0, 0, 0, 0, 0, 0), properties);

Info<< "Reading field C3\n" << endl;

volScalarField C3
(
    IOobject
    (
        "C3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (-z1*C1-z2*C2)/z3
);

#include "createFvOptions.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Gamma
(
    IOobject
    (
        "Gamma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    p/(tau*tau)
);

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
    F*F*(z1*z1*mu1*Gamma*C1+z2*z2*Gamma*mu2*C2+z3*z3*mu3*Gamma*C3)
);

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    z1*C1+z2*C2+z3*C3
);

volVectorField idif
(
    IOobject
    (
        "idif",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -F*(z1*DC1*Gamma*fvc::grad(C1)+z2*DC2*Gamma*fvc::grad(C2)+z3*DC3*Gamma*fvc::grad(C3))
);

volVectorField itot
(
    IOobject
    (
        "itot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    idif-F*F*fvc::grad(Fi)*(z1*z1*mu1*Gamma*C1+z2*z2*mu2*Gamma*C2+z3*z3*mu3*Gamma*C3)
);