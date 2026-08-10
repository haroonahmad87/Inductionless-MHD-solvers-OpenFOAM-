Info<< "Reading ElectromagneticProperties\n" << endl;
//Additional dictionary for magnetic field definition
IOdictionary electromagneticProperties
(
    IOobject
    (
        "electromagneticProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Magnetic field
dimensionedVector B0( "B0", dimensionSet(1, 0, -2, 0, 0, -1, 0), electromagneticProperties); //electromagneticProperties.lookup("B0") );


