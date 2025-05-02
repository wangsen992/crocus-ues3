#include "objectRegistryFuncs.H"

Foam::volScalarField& Foam::lookupOrConstruct
(
    const fvMesh& mesh,
    const char* name, 
    IOobject::readOption readOpt=IOobject::MUST_READ,
    IOobject::writeOption writeOpt=IOobject::AUTO_WRITE
)
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
    {
        volScalarField* fPtr
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    readOpt,
                    writeOpt
                ),
                mesh
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volScalarField>(name);
}
