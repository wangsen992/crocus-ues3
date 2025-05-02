
#include "addToRunTimeSelectionTable.H"
#include "fvPatchField.H"
#include "mixedFvPatchField.H"
#include "wrfInletOutletVelocityFvPatchField.H"
#include "surfaceFields.H"
#include "volFields.H"
// #include "nesting_utils.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wrfInletOutletVelocityFvPatchField::wrfInletOutletVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF)
{}


// template<class Type>
// Foam::wrfInletOutletVelocityFvPatchField<Type>::wrfInletOutletVelocityFvPatchField
// (
//     const fvPatch& p,
//     const DimensionedField<Type, volMesh>& iF,
//     const Field<Type>& fld
// )
// :
//     mixedFvPatchField<Type>(p, iF, fld)
// {}


Foam::wrfInletOutletVelocityFvPatchField::wrfInletOutletVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF, dict),
    fieldName_(dict.lookup<word>("field")),
    wrf_case_root_(dict.lookup<string>("wrf_case_root")),
    wrf_case_name_(dict.lookup<string>("wrf_case_name")),
    pwrfTime_
    (
      autoPtr<Time>
      (
        new Time
        (
          Time::controlDictName, 
          dict.lookup<string>("wrf_case_root"),
          dict.lookup<string>("wrf_case_name")
        )
      )
    ),
    pwrfMesh_(),
    testValue_(dict.lookup<vector>("testValue")),
    wrfi_(-1),
    psiOld_(p.size()),
    psiNew_(p.size())

{
    this->refValue()=(Field<vector>(this->patch().size(), this->testValue_));
    operator==(this->refValue());
    this->refGrad() = Zero;
    this->valueFraction() = Zero;
    Info << "Patch " << this->patch().name();
    Info << "[INIT] value = " << average(*this) << endl;
    Info << "[INIT] refValue = " << average(this->refValue()) << endl;
    Info << "[INIT] refGrad = " << average(this->refGrad()) << endl;
    Info << "[INIT] valueFrac = " << average(this->valueFraction()) << endl;
}


Foam::wrfInletOutletVelocityFvPatchField::wrfInletOutletVelocityFvPatchField
(
    const wrfInletOutletVelocityFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper) // Don't map
{
    // Evaluate since value not mapped
    // this->evaluate();
}


Foam::wrfInletOutletVelocityFvPatchField::wrfInletOutletVelocityFvPatchField
(
    const wrfInletOutletVelocityFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wrfInletOutletVelocityFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vector Utarget(0, 0, 0);
    scalar t = this->db().time().value();
	  if (t > 100)
	  {
	  	// Utarget = vector(5 * sin(6.28*t/1000) + 1 * sin(6.28*t/300), 5 * cos(6.28*t/2000) + 0.5 * cos(6.28*t/450),0);
	  	Utarget = 100/100 * vector(5 * sin(6.28*t/1000) + 1 * sin(6.28*t/3000), 5 * cos(6.28*t/2000) + 0.5 * cos(6.28*t/4500),0);
	  }
	  else 
	  {
	  	// Utarget = 100/100 * vector(5 * sin(6.28*t/1000) + 1 * sin(6.28*t/300), 5 * cos(6.28*t/2000) + 0.5 * cos(6.28*t/450),0);
	  	Utarget = t/100 * vector(5 * sin(6.28*t/1000) + 1 * sin(6.28*t/3000), 5 * cos(6.28*t/2000) + 0.5 * cos(6.28*t/4500),0);
	  }
    Info << "t: " << t << "; Utarget = " << Utarget << endl;
    this->refValue()=(Field<vector>(this->patch().size(), Utarget));
    // this->refValue()=(Field<vector>(this->patch().size(), vector(2, 2, 0)));

    const Field<scalar>& phip =
        this->patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "phi"
        );

    // When flux is outward, set to zeroGradient
    this->valueFraction() = 1.0 - pos0(Utarget & this->patch().Sf());
    // this->valueFraction() = 1.0 - pos0(phip + 0.4 * mag(this->patch().Sf()));
    // this->valueFraction() = 1.0 - pos0(phip);
    // When flux absolute value is small, set to zeroGradient

    // Info << "Patch " << this->patch().name();
    // Info << "[INIT] value = " << average(*this) << endl;
    // Info << "[INIT] refValue = " << average(this->refValue()) << endl;
    // Info << "[INIT] refGrad = " << average(this->refGrad()) << endl;
    // Info << "[INIT] valueFrac = " << average(this->valueFraction()) << endl;
    // Info << "[INIT] phi = " << average(phip) << endl;
    // operator==(this->refValue());
    Field<vector>::operator=
    (
        this->valueFraction()*this->refValue()
      +
        (1.0 - this->valueFraction())*
        (
            this->patchInternalField()
          + this->refGrad()/this->patch().deltaCoeffs()
        )
    );
    mixedFvPatchField<vector>::updateCoeffs();
    // Info << "Patch " << this->patch().name();
    // Info << "[Mix] value = " << average(*this) << endl;
    // Info << "[Mix] refValue = " << average(this->refValue()) << endl;
    // Info << "[Mix] refGrad = " << average(this->refGrad()) << endl;
    // Info << "[Mix] valueFrac = " << average(this->valueFraction()) << endl;
    // Info << "[Mix] phi = " << average(phip) << endl;

}


void Foam::wrfInletOutletVelocityFvPatchField::write(Ostream& os) const
{
    // fvPatchField<Type>::write(os);
    mixedFvPatchField<vector>::write(os);
    writeEntry(os, "field", fieldName_);
    writeEntry(os, "wrf_case_root", wrf_case_root_);
    writeEntry(os, "wrf_case_name", wrf_case_name_);
    writeEntry(os, "testValue", this->testValue_);
}

void Foam::wrfInletOutletVelocityFvPatchField::operator=
(
  const fvPatchField<vector>& ptf
)
{
  fvPatchField<vector>::operator=
  (
    this->valueFraction()*this->refValue()
    + (1 - this->valueFraction())*ptf
  );
}


// ************************************************************************* //

namespace Foam
{
    makePatchTypeField
    (
      fvPatchVectorField, 
      wrfInletOutletVelocityFvPatchField
    );
}

// ********************a**************************************************** //
