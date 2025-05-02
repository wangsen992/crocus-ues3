
#include "wrfFixedValueFvPatchField.H"
#include "nesting_utils.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFvPatchField<Type>(p, iF, fld)
{}


template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
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
    if (dict.found("value"))
    {
        Field<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Essential entry 'value' missing"
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const wrfFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper, false) // Don't map
{
    // Evaluate since value not mapped
    this->evaluate();
}


template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const wrfFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::wrfFixedValueFvPatchField<Type>::updateCoeffs()
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
    auto dt = mesh.time().deltaTValue();
    auto told = wrfTime_.times()[wrfi_].value();
    auto tnew = wrfTime_.times()[wrfi_+1].value();

    // Temporary fix, doesn't support restarting simulations

      this->operator==((tnew-t)/(tnew-told)*psiOld_ + (t-told)/(tnew-told)*psiNew_);
      Info << average(*this) << endl;
      fixedValueFvPatchField<Type>::updateCoeffs();

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::wrfFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "field", fieldName_);
    writeEntry(os, "wrf_case_root", wrf_case_root_);
    writeEntry(os, "wrf_case_name", wrf_case_name_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
