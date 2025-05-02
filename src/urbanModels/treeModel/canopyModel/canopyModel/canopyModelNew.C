
#include "canopyModel.H"
#include "treeModel.H"

namespace Foam
{

autoPtr<canopyModel> canopyModel::New
(
    const treeModel& tree
)
{
    word canopyModelType(tree.dict().subDict("canopy").lookup("type"));

    treeModelConstructorTable::iterator cstrIter = 
        treeModelConstructorTablePtr_->find(canopyModelType);

    if (cstrIter == treeModelConstructorTablePtr_->end())
    {
        FatalErrorInFunction
          << "Unknown canopyModelType type "
          << canopyModelType << endl << endl
          << "Valid canopyModel types are : " << endl
          << treeModelConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return cstrIter()(tree);
}

}
