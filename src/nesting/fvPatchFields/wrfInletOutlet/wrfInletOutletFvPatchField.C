
#include "wrfInletOutletFvPatchField.H"
#include "nesting_utils.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::wrfInletOutletFvPatchField<Type>::wrfInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF)
{}


// template<class Type>
// Foam::wrfInletOutletFvPatchField<Type>::wrfInletOutletFvPatchField
// (
//     const fvPatch& p,
//     const DimensionedField<Type, volMesh>& iF,
//     const Field<Type>& fld
// )
// :
//     mixedFvPatchField<Type>(p, iF, fld)
// {}


template<class Type>
Foam::wrfInletOutletFvPatchField<Type>::wrfInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict),
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
    wrfi_(-1),
    psiOld_(p.size()),
    psiNew_(p.size())

{
}


template<class Type>
Foam::wrfInletOutletFvPatchField<Type>::wrfInletOutletFvPatchField
(
    const wrfInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper) // Don't map
{
    // Evaluate since value not mapped
    // this->evaluate();
}


template<class Type>
Foam::wrfInletOutletFvPatchField<Type>::wrfInletOutletFvPatchField
(
    const wrfInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::wrfInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    label curi(floor(this->db().time().value()/3600.0)+1);
    Time& wrfTime_(pwrfTime_());

    if (curi > wrfi_)
    {
      wrfi_ = curi;

      pwrfMesh_.set
      (
        new fvMesh
        (
          IOobject
          (
            fvMesh::defaultRegion,
            pwrfTime_->timeName(),
            pwrfTime_(),
            IOobject::MUST_READ
          )
        )
      );


      {
        wrfTime_.setTime
        (
          wrfTime_.times()[wrfi_].value(),
          wrfi_
        );
        psiType psi
        (
          IOobject
          (
            fieldName_,
            wrfTime_.timeName(),
            wrfTime_,
            IOobject::MUST_READ
          ),
          pwrfMesh_()
        );

        psiOld_ = interpolate(this->patch().Cf(), psi, "cell");
      }


      {
        wrfTime_.setTime
        (
          wrfTime_.times()[wrfi_+1].value(),
          wrfi_+1
        );
        psiType psi
        (
          IOobject
          (
            fieldName_,
            wrfTime_.timeName(),
            wrfTime_,
            IOobject::MUST_READ
          ),
          pwrfMesh_()
        );

        psiNew_ = interpolate(this->patch().Cf(), psi, "cell");
      }

      pwrfMesh_.clear();

    }

    const fvMesh& mesh(this->patch().boundaryMesh().mesh());
    auto t = mesh.time().value();
    auto told = wrfTime_.times()[wrfi_].value();
    auto tnew = wrfTime_.times()[wrfi_+1].value();


    this->refValue()=((tnew-t)/(tnew-told)*psiOld_ + (t-told)/(tnew-told)*psiNew_);
    const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            "phi"
        );

    this->valueFraction() = 1.0 - pos0(phip);
    Field<Type>::operator=
    (
        this->valueFraction()*this->refValue()
      +
        (1.0 - this->valueFraction())*
        (
            this->patchInternalField()
          + this->refGrad()/this->patch().deltaCoeffs()
        )
    );
    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::wrfInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    // fvPatchField<Type>::write(os);
    mixedFvPatchField<Type>::write(os);
    writeEntry(os, "field", fieldName_);
    writeEntry(os, "wrf_case_root", wrf_case_root_);
    writeEntry(os, "wrf_case_name", wrf_case_name_);
}


// ************************************************************************* //
