/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    scalarTransportFoam

Group
    grpBasicSolvers

Description
    Passive scalar transport equation solver.

    \heading Solver details
    The equation is given by:

    \f[
        \ddt{T} + \div \left(\vec{U} T\right) - \div \left(D_T \grad T \right)
        = S_{T}
    \f]

    Where:
    \vartable
        T       | Passive scalar
        D_T     | Diffusion coefficient
        S_T     | Source
    \endvartable

    \heading Required fields
    \plaintable
        T       | Passive scalar
        U       | Velocity [m/s]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Passive scalar transport equation solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating scalar transport\n" << endl;

    //#include "CourantNo.H"

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            Gamma=p/(tau*tau);
            surfaceScalarField gradFi=fvc::snGrad(Fi)*mesh.magSf();
            fvScalarMatrix C1Eqn
            (
                p*fvm::ddt(C1)
              - fvm::laplacian(Gamma*DC1, C1)
              - z1*(F/R)*(1/T)*fvm::div(fvc::interpolate(Gamma)*DC1*gradFi,C1)
             ==
                fvOptions(C1)
            );

            C1Eqn.relax();
            fvOptions.constrain(C1Eqn);
            C1Eqn.solve();
            fvOptions.correct(C1);

            fvScalarMatrix C2Eqn
            (
                p*fvm::ddt(C2)
              - fvm::laplacian(Gamma*DC2, C2)
              - z2*(F/R)*(1/T)*fvm::div(fvc::interpolate(Gamma)*DC2*gradFi,C2)
             ==
                fvOptions(C2)
            );

            C2Eqn.relax();
            fvOptions.constrain(C2Eqn);
            C2Eqn.solve();
            fvOptions.correct(C2);

            C3 = (-z1*C1-z2*C2)/z3;

            kappa=F*F*(z1*z1*mu1*Gamma*C1+z2*z2*mu2*Gamma*C2+z3*z3*mu3*Gamma*C3);
            
            k=z1*C1+z2*C2+z3*C3;

            fvScalarMatrix FiEqn
            (
                fvm::laplacian(fvc::interpolate(kappa),Fi)
                + F*fvc::laplacian(z1*Gamma*DC1,C1) 
                + F*fvc::laplacian(z2*Gamma*DC2,C2)
                + F*fvc::laplacian(z3*Gamma*DC3,C3)
            );

            FiEqn.solve();
            vgradFi=fvc::grad(Fi);

            idif=-F*(z1*DC1*Gamma*fvc::grad(C1)+z2*DC2*Gamma*fvc::grad(C2)+z3*DC3*Gamma*fvc::grad(C3));

            itot = idif-F*F*fvc::grad(Fi)*(z1*z1*mu1*Gamma*C1+z2*z2*mu2*Gamma*C2+z3*z3*mu3*Gamma*C3);

        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


