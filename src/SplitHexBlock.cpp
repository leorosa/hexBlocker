/*
Copyright 2016
Author Leonardo Rosa,
user "leorosa" at github.com

License
    This file is part of hexBlocker.

    hexBlocker is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    hexBlocker is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with hexBlocker.  If not, see <http://www.gnu.org/licenses/>.

    The license is included in the file COPYING.
*/

#include "HexBlocker.h"
#include "HexBlock.h"
#include "HexEdge.h"
#include "HexPatch.h"
//#include "HexReader.h"

//#include <vtkObjectFactory.h>
#include <vtkCollection.h>
#include <vtkIdList.h>
#include <vtkRenderer.h>

void HexBlocker::splitHexBlock(vtkIdType edgeId)
{
    HexEdge *edge;
    vtkSmartPointer<vtkIdList> parallelEdges = vtkSmartPointer<vtkIdList>::New();
    addParallelEdges(parallelEdges, edgeId);

    vtkSmartPointer<vtkIdList> parallelBlocks = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i=0; i<hexBlocks->GetNumberOfItems(); i++)  // get parallel blocks
    {
        HexBlock * b = HexBlock::SafeDownCast(hexBlocks->GetItemAsObject(i));
        for(vtkIdType j=0; j<parallelEdges->GetNumberOfIds(); j++)
        {
            if (b->getParallelEdges(parallelEdges->GetId(j))->GetNumberOfIds() > 0)
            {
                parallelBlocks->InsertUniqueId(i);
                break;
            }
        }
    }

    vtkIdType nv, parBlockId, parEdgeId;
    for(vtkIdType i=0; i<parallelBlocks->GetNumberOfIds(); i++)          // insert new blocks
    {
        parBlockId = parallelBlocks->GetId(i);
        vtkSmartPointer<vtkIdList> baseVertices   = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> middleVertices = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> topVertices    = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> blockVertices0 = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> blockVertices1 = vtkSmartPointer<vtkIdList>::New();
        HexBlock * b = HexBlock::SafeDownCast(hexBlocks->GetItemAsObject(parBlockId));
        for(vtkIdType j=0;j<parallelEdges->GetNumberOfIds();j++)
        {
            parEdgeId = parallelEdges->GetId(j);
            if (b->getParallelEdges(parEdgeId)->GetNumberOfIds() > 0)
            {
                double posm[3], tmpCoords[3];
                edge = HexEdge::SafeDownCast(edges->GetItemAsObject(parEdgeId));
                edge->calcParametricPointOnLine(0.5, posm);

                nv = vertices->GetNumberOfPoints(); // create a new vertice if needed
                vtkIdType ptId;                     // this code is copied from HexBlock::init()
                for(ptId=0; ptId<nv; ptId++)
                {
                    vertices->GetPoint(ptId, tmpCoords);
                    if(    tmpCoords[0] == posm[0]
                        && tmpCoords[1] == posm[1]
                        && tmpCoords[2] == posm[2] )
                    {
                        break;
                    }
                }
                if(ptId==nv)
                {
                    vertices->InsertNextPoint(posm);
                }
// how to guarantee the correct order of the vertices ??
// -> it still doesn't correct the middle point calculated using edges with inverted ends
                middleVertices->InsertUniqueId(ptId);
                baseVertices->InsertUniqueId(edge->getVertIds(0));
                topVertices->InsertUniqueId(edge->getVertIds(1));
            }
        }
        vertices->Modified();

        nv = middleVertices->GetNumberOfIds();
        for(vtkIdType j=0; j<nv; j++)   // construct the vertices list for two new blocks
        {
            blockVertices0->InsertUniqueId(baseVertices->GetId(j));
            blockVertices1->InsertUniqueId(middleVertices->GetId(j));
        }
        for(vtkIdType j=0; j<nv; j++)
        {
            blockVertices0->InsertUniqueId(middleVertices->GetId(j));
            blockVertices1->InsertUniqueId(topVertices->GetId(j));
        }

//      printf("before: ");
//      for(vtkIdType j=0; j<blockVertices0->GetNumberOfIds(); j++) {
//          printf("%d ", blockVertices0->GetId(j)); }
//      printf("\n");
        orderVertices(blockVertices0);
//      printf("after:  ");
//      for(vtkIdType j=0; j<blockVertices0->GetNumberOfIds(); j++) {
//          printf("%d ", blockVertices0->GetId(j)); }
//      printf("\n");

//      printf("before: ");
//      for(vtkIdType j=0; j<blockVertices1->GetNumberOfIds(); j++) {
//          printf("%d ", blockVertices1->GetId(j)); }
//      printf("\n");
        orderVertices(blockVertices1);
//      printf("after:  ");
//      for(vtkIdType j=0; j<blockVertices1->GetNumberOfIds(); j++) {
//          printf("%d ", blockVertices1->GetId(j)); }
//      printf("\n");

        createHexBlock(blockVertices0);
        createHexBlock(blockVertices1);
    }

    for(vtkIdType i=parallelBlocks->GetNumberOfIds()-1; i>=0; i--)      // remove old blocks
    {
        removeHexBlock(parallelBlocks->GetId(i));
    }
}

void HexBlocker::orderVertices(vtkIdList *selectedVertices)
{
    double p[8][3];
    int idx[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    int tmpId;
    vtkSmartPointer<vtkIdList> tmpVertices = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i=0; i<selectedVertices->GetNumberOfIds(); i++)
        tmpVertices->InsertUniqueId(selectedVertices->GetId(i));
    vertices->GetPoint(selectedVertices->GetId(0), p[0]);
    vertices->GetPoint(selectedVertices->GetId(1), p[1]);
    vertices->GetPoint(selectedVertices->GetId(2), p[2]);
    vertices->GetPoint(selectedVertices->GetId(3), p[3]);
    vertices->GetPoint(selectedVertices->GetId(4), p[4]);
    vertices->GetPoint(selectedVertices->GetId(5), p[5]);
    vertices->GetPoint(selectedVertices->GetId(6), p[6]);
    vertices->GetPoint(selectedVertices->GetId(7), p[7]);
// selecting the lower ones in Z
    for(int i=0; i<7; i++)
        for(int j=i+1; j<8; j++)
            if (p[idx[i]][2] > p[idx[j]][2])
            {
                tmpId = idx[j];
                idx[j] = idx[i];
                idx[i] = tmpId;
            }
// selecting the lower ones in Y
    for(int i=0; i<3; i++)
        for(int j=i+1; j<4; j++)
            if (p[idx[i]][1] > p[idx[j]][1])
            {
                tmpId = idx[j];
                idx[j] = idx[i];
                idx[i] = tmpId;
            }
    for(int i=4; i<7; i++)
        for(int j=i+1; j<8; j++)
            if (p[idx[i]][1] > p[idx[j]][1])
            {
                tmpId = idx[j];
                idx[j] = idx[i];
                idx[i] = tmpId;
            }
// selecting the lower ones in X
    if (p[idx[0]][0] > p[idx[1]][0])
    {
        tmpId = idx[0];
        idx[0] = idx[1];
        idx[1] = tmpId;
    }
    if (p[idx[3]][0] > p[idx[2]][0])
    {
        tmpId = idx[2];
        idx[2] = idx[3];
        idx[3] = tmpId;
    }
    if (p[idx[4]][0] > p[idx[5]][0])
    {
        tmpId = idx[4];
        idx[4] = idx[5];
        idx[5] = tmpId;
    }
    if (p[idx[7]][0] > p[idx[6]][0])
    {
        tmpId = idx[6];
        idx[6] = idx[7];
        idx[7] = tmpId;
    }
    for(int i=0; i<8; i++)
        selectedVertices->SetId(i, tmpVertices->GetId(idx[i]));
}

