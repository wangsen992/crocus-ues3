
#include "trunkModel.H"
#include "treeModel.H"

namespace Foam
{

autoPtr<trunkModel> trunkModel::New
(
    const treeModel& tree
)
{
    word trunkModelType(tree.dict().subDict("trunk").lookup("type"));

    treeModelConstructorTable::iterator cstrIter = 
        treeModelConstructorTablePtr_->find(trunkModelType);

    if (cstrIter == treeModelConstructorTablePtr_->end())
    {
        FatalErrorInFunction
          << "Unknown trunkModelType type "
          << trunkModelType << endl << endl
          << "Valid trunkModel types are : " << endl
          << treeModelConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return cstrIter()(tree);
}

}
