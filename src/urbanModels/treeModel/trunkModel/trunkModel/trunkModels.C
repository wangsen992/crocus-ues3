#include "trunkModel.H"
#include "trunkSurfaceModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

typedef trunkSurfaceModel<trunkModel> surfaceTrunkModel;

addNamedToRunTimeSelectionTable
(
  trunkModel,
  surfaceTrunkModel,
  treeModel,
  surfaceTrunkModel
);

}
