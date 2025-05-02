// Base class for other patchFields, mostly handling the communication
// with WRF class object

#ifndef wrfInletOutletFvPatchField_H
#define wrfInletOutletFvPatchField_H

#include "mixedFvPatchField.H"

namespace Foam
{

template<typename Type>
class wrfInletOutletFvPatchField
  : public mixedFvPatchField<Type>
{
    typedef GeometricField<Type, fvPatchField, volMesh> psiType;

    // autoPtr<WRF>  wrf_;
    word  fieldName_;
    string wrf_case_root_;
    string wrf_case_name_;
    autoPtr<Time> pwrfTime_;
    autoPtr<fvMesh> pwrfMesh_;
    label wrfi_;
    Field<Type> psiOld_;
    Field<Type> psiNew_;

public:
    //- Runtime Type Information
    TypeName("wrfInletOutlet");

    // Constructors

        //- Construct from patch and internal field
        wrfInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wrfInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wrfInletOutletFvPatchField
        //  onto a new patch
        wrfInletOutletFvPatchField
        (
            const wrfInletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        wrfInletOutletFvPatchField
        (
            const wrfInletOutletFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        wrfInletOutletFvPatchField
        (
            const wrfInletOutletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new wrfInletOutletFvPatchField<Type>(*this, iF)
            );
        }
        
    // Member Functions

        // Attributes

            // //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return true;
            }
    
    // Evaluation functions

      //- Update the coefficients associated with the patch field
      virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;

};

}

#ifdef NoRepository
    #include "wrfInletOutletFvPatchField.C"
#endif

#endif
