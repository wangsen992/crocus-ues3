
#include "trunkSurfaceModel.H"

namespace Foam
{
    
template<class BaseTrunkModel>
trunkSurfaceModel<BaseTrunkModel>::trunkSurfaceModel
(
    const treeModel& tree
)
:
    BaseTrunkModel(tree)
{}

}
