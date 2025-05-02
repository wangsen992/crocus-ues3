#include "nesting_utils.H"
#include "polyMesh.H"
// compute cells close to a patch by directly calculating the minimun distance
// of each cell to the patchFaces
// This is a slow method, scales with O(N x M) with N faces and M cells
Foam::HashTable<Foam::scalar, Foam::label> getPatchCloseCells(const Foam::polyMesh& mesh, Foam::string patchName, double dis_lim)
{
  Foam::HashTable<Foam::scalar, Foam::label> tbl;
  // certain distance
  Foam::polyPatch patch(mesh.boundaryMesh()[patchName]);
  Foam::pointField patchFaceCenters(patch.faceCentres());
  Foam::vectorField patchFaceNormals(patch.faceNormals());

  // Find cells cut by rays from face centers
  // by limite the eucleadian distance to set boundarys

  Foam::pointField meshCellCenters(mesh.cellCentres());
  Foam::DynamicList<Foam::label, 10> layerCellsInd;
  Foam::DynamicList<Foam::scalar, 10> layerCellsDis;
  for(int i = 0; i != meshCellCenters.size(); i++)
  {
    Foam::tmp<Foam::vectorField> tDisField = (meshCellCenters[i] - patchFaceCenters);
    Foam::scalar dis = min(mag(tDisField));
    if (dis < dis_lim)
    {
      tbl.set(i, dis);
    }
  }
  return tbl;
}

void combineCloseCellTables(Foam::HashTable<Foam::scalar, Foam::label>& hostTbl, const Foam::HashTable<Foam::scalar, Foam::label>&  srcTbl)
{
    // Note, if already present, update the value only
    for(auto i : srcTbl.sortedToc())
    {
      if(hostTbl.find(i) == hostTbl.end())
      {
          hostTbl.set(i, srcTbl[i]);
      }
      else
      {
          hostTbl.set(i, Foam::max(hostTbl[i], srcTbl[i]));
      }
    }
}
