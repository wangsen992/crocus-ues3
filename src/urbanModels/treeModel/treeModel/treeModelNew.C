#include "treeModel.H"

namespace Foam
{

autoPtr<treeModel> treeModel::New
(
    const fvMesh& mesh,
    const dictionary& treeModelDict
)
{
    autoPtr<treeModel> treePtr;
    treePtr.set
    (
      new treeModel(mesh, treeModelDict)
    );
    return treePtr;
};

}
