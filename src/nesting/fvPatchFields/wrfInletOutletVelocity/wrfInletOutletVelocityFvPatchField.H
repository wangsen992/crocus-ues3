// Base class for other patchFields, mostly handling the communication
// with WRF class object

#ifndef wrfInletOutletVelocityFvPatchField_H
#define wrfInletOutletVelocityFvPatchField_H

#include "mixedFvPatchField.H"

namespace Foam
{

class wrfInletOutletVelocityFvPatchField
  : public mixedFvPatchField<vector>
{
    typedef GeometricField<vector, fvPatchField, volMesh> psiType;

    // autoPtr<WRF>  wrf_;
    word  fieldName_;
    string wrf_case_root_;
    string wrf_case_name_;
    autoPtr<Time> pwrfTime_;
    autoPtr<fvMesh> pwrfMesh_;
    vector  testValue_;
    label wrfi_;
    Field<vector> psiOld_;
    Field<vector> psiNew_;

public:
    //- Runtime Type Information
    TypeName("wrfInletOutletVelocity");

    // Constructors

        //- Construct from patch and internal field
        wrfInletOutletVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wrfInletOutletVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wrfInletOutletVelocityFvPatchField
        //  onto a new patch
        wrfInletOutletVelocityFvPatchField
        (
            const wrfInletOutletVelocityFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        wrfInletOutletVelocityFvPatchField
        (
            const wrfInletOutletVelocityFvPatchField&
        ) = delete;

        //- Copy constructor setting internal field reference
        wrfInletOutletVelocityFvPatchField
        (
            const wrfInletOutletVelocityFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new wrfInletOutletVelocityFvPatchField(*this, iF)
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

    //- Member Operators
    virtual void operator=(const fvPatchField<vector>& pvf);

};

}


#endif
