singlePhaseTransportModel laminarTransport(U, phi); //This is a known class in OpenFoam, hence no need of declaration

// Thermal expansion coefficient [1/K]
dimensionedScalar beta("beta", dimless/dimTemperature, laminarTransport);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

//Electrical conductivity
dimensionedScalar sigma("sigma", dimensionSet(-1, -3, 3, 0, 0, 2, 0), laminarTransport);

// Reference density  // *** added by H.A.
dimensionedScalar rho_ref("rho_ref", dimDensity, laminarTransport);
