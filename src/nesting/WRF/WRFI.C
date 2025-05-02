#include "WRF.H"
#include "interpolation.H"

using namespace Foam;

template<typename Type>
Type WRF::interpolate(const point& pt, const GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod)
{
  autoPtr<interpolation<Type>> interp
  (
    interpolation<Type>::New(interpMethod, psi)
  );
  point wrf_pt = transform(pt);

  Type interpVal = interp->interpolate
  (
    wrf_pt, 
    tsearcher_->findCell(wrf_pt)
  );
  return interpVal;
  
}

template<typename Type>
Field<Type> WRF::interpolate(const Field<point>& pts, const GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod)
{
  autoPtr<interpolation<Type>> interp
  (
    interpolation<Type>::New(interpMethod, psi)
  );
  pointField wrf_pts = pts;
  Field<Type> interpVals(pts.size());
  for(size_t i=0; i < pts.size(); i++)
  {
    label cellInd = tsearcher_->findCell(wrf_pts[i]);
    if(cellInd < 0)
    {
      Info << "[Debug] point not in domain: " << wrf_pts[i];
      cellInd = tsearcher_->findNearestCell(wrf_pts[i]);
      Info << " distance = " << tsearcher_->mesh().points()[cellInd].z() - wrf_pts[i].z() << endl;
    }
    interpVals[i] = interp->interpolate
    (
      wrf_pts[i],
      tsearcher_->findCell(wrf_pts[i])
    );
  }
  return interpVals;
}
