/*
Copyright 2012, 2013
Author Nicolas Edh,
Nicolas.Edh@gmail.com,
or user "nsf" at cfd-online.com

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

#include "HexBlock.h"
#include <vtkObjectFactory.h>

#include "HexPatch.h"
#include "HexEdge.h"


#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCollection.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>

vtkStandardNewMacro(HexBlock);

//Default constructor
HexBlock::HexBlock()
{
    vertIds = vtkSmartPointer<vtkIdList>::New();
    localEdges = vtkSmartPointer<vtkCollection>::New();
    localPatches = vtkSmartPointer<vtkCollection>::New();

    hexData = vtkSmartPointer<vtkPolyData>::New();
    axesTubes = vtkSmartPointer<vtkTubeFilter>::New();
    globalEdges = vtkSmartPointer<vtkCollection>::New();
    globalPatches = vtkSmartPointer<vtkCollection>::New();

    hexAxisActor = vtkSmartPointer<vtkActor>::New();
    hexBlockActor = vtkSmartPointer<vtkActor>::New();
}

HexBlock::~HexBlock()
{

}

void HexBlock::init(double corner0[3], double corner1[3],
                    vtkSmartPointer<vtkPoints> verts,
                    vtkSmartPointer<vtkCollection> edges,
                    vtkSmartPointer<vtkCollection> patches)
{


    double x0=corner0[0];
    double y0=corner0[1];
    double z0=corner0[2];
    double x1=corner1[0];
    double y1=corner1[1];
    double z1=corner1[2];

    globalVertices = verts;
    globalEdges = edges;
    globalPatches = patches;

    vtkIdType j = globalVertices->GetNumberOfPoints();

    globalVertices->InsertNextPoint(x0,y0,z0);
    globalVertices->InsertNextPoint(x1,y0,z0);
    globalVertices->InsertNextPoint(x1,y1,z0);
    globalVertices->InsertNextPoint(x0,y1,z0);
    globalVertices->InsertNextPoint(x0,y0,z1);
    globalVertices->InsertNextPoint(x1,y0,z1);
    globalVertices->InsertNextPoint(x1,y1,z1);
    globalVertices->InsertNextPoint(x0,y1,z1);


    //set global Ids
    vertIds->SetNumberOfIds(8);
    for(vtkIdType i=0;i<8;i++)
        vertIds->SetId(i,i+j);

    initAll();
}

void HexBlock::init(vtkSmartPointer<vtkIdList> myVertIds,
                    vtkSmartPointer<vtkPoints> verts,
                    vtkSmartPointer<vtkCollection> edges,
                    vtkSmartPointer<vtkCollection> patches)
{
    globalVertices = verts;
    globalEdges = edges;
    globalPatches = patches;

    vertIds = myVertIds;

    initAll();

}

//Extrude from patch and distance
void HexBlock::init(vtkSmartPointer<HexPatch> p,
                    double dist,
                    vtkSmartPointer<vtkPoints> verts,
                    vtkSmartPointer<vtkCollection> edges,
                    vtkSmartPointer<vtkCollection> patches)
{
    HexBlock * fromHex = p->getPrimaryHexBlock();
    globalVertices=verts;
    globalEdges = edges;
    globalPatches = patches;

    vtkIdList * oldIds =  p->vertIds;
    vtkIdList * newIds;
    newIds = vtkSmartPointer<vtkIdList>::New();
    newIds->SetNumberOfIds(4);

    double n_dist[3];
    p->getNormal(n_dist);
    vtkMath::MultiplyScalar(n_dist,dist);
    //Insert the new 4 points into globalVertices
    for(vtkIdType i=0;i<4;i++)
    {
        double oldCoords[3];
        double newCoords[3];
        vtkIdType ptId;
        globalVertices->GetPoint(p->vertIds->GetId(i),oldCoords);
        vtkMath::Add(oldCoords,n_dist,newCoords);
// I would make sense to create 4 new vertices for the extruded block, only if
// the extrusion of a list of patches were permited. However, it is not the
// case, and thus it is more useful to check if new vertices are needed.
        ptId = insertUniqueVertice(globalVertices, newCoords);
        newIds->SetId(i, ptId);
    }

    //We want to keep extruded blocks aligned
    //so the extrude operation is diffierent depending
    //on which patch was selected
    switch (fromHex->getPatchInternalId(p))
    {

    case 0:
        vertIds->InsertNextId(newIds->GetId(0));
        vertIds->InsertNextId(newIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(2));
        vertIds->InsertNextId(newIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(0));
        vertIds->InsertNextId(oldIds->GetId(3));
        vertIds->InsertNextId(oldIds->GetId(2));
        vertIds->InsertNextId(oldIds->GetId(1));
        break;
    case 1:
        vertIds->InsertNextId(newIds->GetId(0));
        vertIds->InsertNextId(newIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(0));
        vertIds->InsertNextId(newIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(2));
        vertIds->InsertNextId(oldIds->GetId(2));
        vertIds->InsertNextId(oldIds->GetId(3));
        break;
    case 2:
        vertIds->InsertNextId(newIds->GetId(0));
        vertIds->InsertNextId(oldIds->GetId(0));
        vertIds->InsertNextId(oldIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(2));
        vertIds->InsertNextId(newIds->GetId(2));
        break;
    case 3:
        vertIds->InsertNextId(oldIds->GetId(0));
        vertIds->InsertNextId(newIds->GetId(0));
        vertIds->InsertNextId(newIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(2));
        vertIds->InsertNextId(oldIds->GetId(2));
        break;
    case 4:
        vertIds->InsertNextId(oldIds->GetId(0));
        vertIds->InsertNextId(oldIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(0));
        vertIds->InsertNextId(oldIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(2));
        vertIds->InsertNextId(newIds->GetId(2));
        vertIds->InsertNextId(newIds->GetId(1));
        break;
    case 5:
        vertIds->InsertNextId(oldIds->GetId(0));
        vertIds->InsertNextId(oldIds->GetId(1));
        vertIds->InsertNextId(oldIds->GetId(2));
        vertIds->InsertNextId(oldIds->GetId(3));
        vertIds->InsertNextId(newIds->GetId(0));
        vertIds->InsertNextId(newIds->GetId(1));
        vertIds->InsertNextId(newIds->GetId(2));
        vertIds->InsertNextId(newIds->GetId(3));
        break;
    default:
        std::cout << "unexpected patch number" << std::endl;
        return;
    }

    initAll();

}

vtkIdType HexBlock::getPatchInternalId(vtkSmartPointer<HexPatch> otherP)
{
    return localPatches->IsItemPresent(otherP)-1;
}

void HexBlock::PrintSelf(ostream &os, vtkIndent indent)
{

}

void HexBlock::exportDict(QTextStream &os)
{
    os << "\t hex (";
    for(vtkIdType j=0; j<vertIds->GetNumberOfIds();j++)
    {
        os << vertIds->GetId(j);
        if(j < vertIds->GetNumberOfIds()-1)
            os << " ";
        else
            os << ") ";
    }

    int nCells[3];
    getNumberOfCells(nCells);
    os << "("
       << nCells[0] <<" "<< nCells[1]<<" " << nCells[2]
       << ") ";
    double gradings[12];

    if(getGradings(gradings))
    {
        os << "simpleGrading (" <<gradings[0] << " "
           << gradings[4] << " " << gradings[8] << ")";
    }
    else
    {
        os << "edgeGrading ( ";
        for(int i =0;i<12;i++)
                os << gradings[i] << " ";
        os << ")" << endl;
    }
}

bool HexBlock::getGradings(double gradings[12])
{
    for(vtkIdType i=0;i<localEdges->GetNumberOfItems();i++)
    {
        HexEdge *e = HexEdge::SafeDownCast(localEdges->GetItemAsObject(i));
        gradings[i] = e->grading;
    }
    return (gradings[0] == gradings[1]
            &&gradings[2] == gradings[3]
            &&gradings[1] == gradings[2] //first four edges
            &&gradings[4] == gradings[5]
            &&gradings[6] == gradings[7]
            &&gradings[4] == gradings[5]//second four edges
            &&gradings[8] == gradings[9]
            &&gradings[10] == gradings[11]
            &&gradings[9] == gradings[10] //last edges
            );
}

void HexBlock::drawLocalaxes()
{
    unsigned char red[3] = {255, 0, 0};
    unsigned char green[3] = {0, 255, 0};
    unsigned char blue[3] = {0, 0, 255};
    vtkSmartPointer<vtkUnsignedCharArray> colors =
            vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");
    colors->InsertNextTupleValue(red);
    colors->InsertNextTupleValue(green);
    colors->InsertNextTupleValue(blue);

    vtkSmartPointer<vtkCellArray> axes =
            vtkSmartPointer<vtkCellArray>::New();
//    axes->SetNumberOfCells(3); //don't do this, this doubles the number of lines why?
//    axes->InitTraversal();


    vtkSmartPointer<vtkLine> xaxes =
            vtkSmartPointer<vtkLine>::New();
    xaxes->GetPointIds()->SetId(0,vertIds->GetId(0));
    xaxes->GetPointIds()->SetId(1,vertIds->GetId(1));
    axes->InsertNextCell(xaxes);

    vtkSmartPointer<vtkLine> yaxes =
            vtkSmartPointer<vtkLine>::New();
    yaxes->GetPointIds()->SetId(0,vertIds->GetId(0));
    yaxes->GetPointIds()->SetId(1,vertIds->GetId(3));
    axes->InsertNextCell(yaxes);

    vtkSmartPointer<vtkLine> zaxes =
            vtkSmartPointer<vtkLine>::New();
    zaxes->GetPointIds()->SetId(0,vertIds->GetId(0));
    zaxes->GetPointIds()->SetId(1,vertIds->GetId(4));
    axes->InsertNextCell(zaxes);

    axesData = vtkSmartPointer<vtkPolyData>::New();
    axesData->SetPoints(globalVertices);
    axesData->SetLines(axes);

    axesData->GetCellData()->SetScalars(colors);

#if VTK_MAJOR_VERSION >= 6
    axesTubes->SetInputData(axesData);
#else
    axesTubes->SetInput(axesData);
#endif    
    axesTubes->SetNumberOfSides(24);
    axesTubes->SetCapping(1);

    // Visualize
    vtkSmartPointer<vtkPolyDataMapper> axesMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    axesMapper->SetInputConnection(axesTubes->GetOutputPort());


    hexAxisActor->SetMapper(axesMapper);
    double x[3];
    globalVertices->GetPoint(vertIds->GetId(0),x);
    hexAxisActor->SetOrigin(x);
    hexAxisActor->SetScale(0.4);
    //hexActor->GetProperty()->SetLineWidth(3);
    hexAxisActor->SetPickable(false);
}

void HexBlock::drawBlock()
{
    vtkSmartPointer<vtkHexahedron> hexahedron =
            vtkSmartPointer<vtkHexahedron>::New();
    for(vtkIdType i=0;i<8;i++)
    {
        hexahedron->GetPointIds()->SetId(i,vertIds->GetId(i));
    }

    vtkSmartPointer<vtkCellArray> hexes =
            vtkSmartPointer<vtkCellArray>::New();
    hexes->InsertNextCell(hexahedron);

    vtkSmartPointer<vtkUnstructuredGrid> uGrid =
             vtkSmartPointer<vtkUnstructuredGrid>::New();
    uGrid->SetPoints(globalVertices);
    uGrid->InsertNextCell(hexahedron->GetCellType(), hexahedron->GetPointIds());

    vtkSmartPointer<vtkDataSetMapper> hexaMapper =
             vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    hexaMapper->SetInputConnection(uGrid->GetProducerPort());
#else
    hexaMapper->SetInputData(uGrid);
#endif

    hexBlockActor->SetMapper(hexaMapper);
    double c[3];
    getCenter(c);
    hexBlockActor->SetOrigin(c);
    hexBlockActor->SetScale(0.4);
    hexBlockActor->GetProperty()->SetColor(0,100,110);
    hexBlockActor->GetProperty()->SetEdgeColor(0,0,0);
    hexBlockActor->GetProperty()->EdgeVisibilityOn();
    hexBlockActor->SetVisibility(true);
    hexBlockActor->GetProperty()->SetOpacity(6.0);

}

void HexBlock::initAll()
{
    initEdges();
    initPatches();
    drawLocalaxes();
    drawBlock();
}

void HexBlock::initEdges()
{
    //clear local collection of old edges
    localEdges->RemoveAllItems();
    //Keep the same order of edges as on
    //docs on blockMesh
    initEdge(vertIds->GetId(1),vertIds->GetId(0)); //0
    initEdge(vertIds->GetId(2),vertIds->GetId(3)); //1
    initEdge(vertIds->GetId(6),vertIds->GetId(7)); //2
    initEdge(vertIds->GetId(5),vertIds->GetId(4)); //3

    initEdge(vertIds->GetId(3),vertIds->GetId(0)); //4
    initEdge(vertIds->GetId(2),vertIds->GetId(1)); //5
    initEdge(vertIds->GetId(6),vertIds->GetId(5)); //6
    initEdge(vertIds->GetId(7),vertIds->GetId(4)); //7

    initEdge(vertIds->GetId(4),vertIds->GetId(0)); //8
    initEdge(vertIds->GetId(5),vertIds->GetId(1)); //9
    initEdge(vertIds->GetId(6),vertIds->GetId(2)); //10
    initEdge(vertIds->GetId(7),vertIds->GetId(3)); //11
}

void HexBlock::initEdge(vtkIdType p0, vtkIdType p1)
{
    vtkSmartPointer<HexEdge> newEdge =
            vtkSmartPointer<HexEdge>::New();
    newEdge->init(p0,p1,globalVertices);
    bool existsInGlobal = false;
    vtkIdType posInGlobal=-1;

    for(vtkIdType i=0;i<globalEdges->GetNumberOfItems();i++)
    {
        HexEdge *e = HexEdge::SafeDownCast(globalEdges->GetItemAsObject(i));
        if(newEdge->equals(e))
        {
//            std::cout<< " edge already exists!" << std::endl;
            existsInGlobal = true;
            posInGlobal=i;
            break;
        }
    }

    if(!existsInGlobal)
    {
//        std::cout << "adding edge ("
//                  << newEdge->vertIds->GetId(0) << ","
//                  << newEdge->vertIds->GetId(1) << ")"
//                  << std::endl;
        globalEdges->AddItem(newEdge);
        localEdges->AddItem(newEdge);
    }
    else
    {
        HexEdge *e = HexEdge::SafeDownCast(globalEdges->GetItemAsObject(posInGlobal));
        localEdges->AddItem(e);
    }

}

void HexBlock::initPatches()
{
    //insert patches like a dice
    //first and last patches define the hexblock.
    //order of points according to OF:
    //Normal is out of domain and according to right hand rule
    //when cycling vertices in patch
    //(dice #1)
    initPatch(0,3,2,1);
    //(dice #2)
    initPatch(0,1,5,4);
    //(dice #3)
    initPatch(0,4,7,3);
    //(dice #4)
    initPatch(1,2,6,5);
    //(dice #5)
    initPatch(3,7,6,2);
    //(dice #6)
    initPatch(4,5,6,7);

    globalPatches->Modified();
}

void HexBlock::initPatch(int id0,int id1,int id2,int id3)
{
    vtkSmartPointer<HexPatch> patch = vtkSmartPointer<HexPatch>::New();

    vtkSmartPointer<vtkIdList> vlist = //patch list of vertIds
            vtkSmartPointer<vtkIdList>::New();

    vlist->InsertId(0,vertIds->GetId(id0));
    vlist->InsertId(1,vertIds->GetId(id1));
    vlist->InsertId(2,vertIds->GetId(id2));
    vlist->InsertId(3,vertIds->GetId(id3));

    patch->init(vlist,globalVertices,this);

    vtkIdType pId=patchIdInGlobalList(patch);
    if(pId >= 0)
    {
        //patch already exist, only add reference to this
        //in the existing one, needed by extrude.
        //std::cout << "item is present" << std::endl;
        HexPatch * existingPatch =
                HexPatch::SafeDownCast(globalPatches->GetItemAsObject(pId));
        existingPatch->setHex(this);
        localPatches->AddItem(existingPatch);
    }
    else //i.e. it doesn't exist in global list
    {
        //add to global and local collection
        globalPatches->AddItem(patch);
        localPatches->AddItem(patch);
    }
}



vtkSmartPointer<vtkIdList> HexBlock::getParallelEdges(vtkIdType edgeId)
{
    vtkSmartPointer<vtkIdList> parallelEdges =
            vtkSmartPointer<vtkIdList>::New();
    // check if we've got a positive edgeId as arg
    if(edgeId<0)
        return parallelEdges;


    HexEdge *edge = HexEdge::SafeDownCast(globalEdges->GetItemAsObject(edgeId));

    //do we have the edge?
    vtkIdType internalEdgeId = localEdges->IsItemPresent(edge) -1;
    bool haveEdge=internalEdgeId > -1;

    if(haveEdge)
    {
        // edges 0-3 are parallel to each other
        // and so are 4-7 and 8-11.
        // the list of parrallel are always edgeId -mod(edgeId) + 0 1 2 3
        // Let's say edge 6 is selected, mainEdgeId=6-6%4=6-2=4 and
        // the parallel edges are mainEdgeId + 0,1,2,3 i.e
        // 4,5,6,7. Find them in the global list and return
        vtkIdType mainEdgeId = internalEdgeId - internalEdgeId%4;

        vtkIdType idInGlobal=-1;
        for(vtkIdType i=0;i<4;i++)
        {
            HexEdge *e = HexEdge::SafeDownCast(localEdges->GetItemAsObject(mainEdgeId+i));
            parallelEdges->InsertNextId(globalEdges->IsItemPresent(e)-1);
        }
/*
        std::cout << "pEdges: (" ;
        for(vtkIdType i=0;i<parallelEdges->GetNumberOfIds();i++)
        {
            std::cout << parallelEdges->GetId(i) << " ";
        }
        std::cout << ")" << endl;
*/
    }

    return parallelEdges;
}

void HexBlock::getNumberOfCells(int nCells[3])
{
    nCells[0] = HexEdge::SafeDownCast(localEdges->GetItemAsObject(0))->nCells;
    nCells[1] = HexEdge::SafeDownCast(localEdges->GetItemAsObject(4))->nCells;
    nCells[2] = HexEdge::SafeDownCast(localEdges->GetItemAsObject(8))->nCells;
}

void HexBlock::setAxesRadius(double rad)
{
    axesTubes->SetRadius(rad);
}



vtkIdType HexBlock::patchIdInGlobalList(vtkSmartPointer<HexPatch> p)
{
    for(vtkIdType i=0;i<globalPatches->GetNumberOfItems();i++)
    {
//        std::cout << "checking with patch: " << i << std::endl;
        HexPatch * op =
                HexPatch::SafeDownCast(globalPatches->GetItemAsObject(i) );
        if(op->equals(p))
        {
//            std::cout << "Patch" << i <<" already exist" << std::endl;
            return i;
        }
    }

    return -1;
}

void HexBlock::rescaleActor()
{
    double x[3];
    globalVertices->GetPoint(vertIds->GetId(0),x);
    hexAxisActor->SetOrigin(x);
    hexAxisActor->SetScale(0.4);

    double c[3];
    getCenter(c);
    hexBlockActor->SetOrigin(c);
    hexBlockActor->SetScale(0.4);
}

void HexBlock::resetColor()
{
    hexBlockActor->GetProperty()->SetColor(0,100,110);
}

void HexBlock::getCenter(double center[])
{
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;

    double pos[3];
    for(vtkIdType i=0;i<8;i++)
    {
        globalVertices->GetPoint(vertIds->GetId(i),pos);
        vtkMath::Add(center,pos,center);
//        std::cout << "boxcenter0: " << center[0] << " "<<  center[1] << " " << center[2] << std::endl;
    }
    vtkMath::MultiplyScalar(center,1.0/8.0);
//    std::cout << "boxcenter1: " << center[0] << " "<<  center[1] << " " << center[2] << std::endl;
}

void HexBlock::changeVertId(vtkIdType from, vtkIdType to)
{
    vtkIdType pos = vertIds->IsId(from);
    if(pos >= 0)
    {
        vertIds->SetId(pos,to);
    }
    drawLocalaxes();
    drawBlock();
}

void HexBlock::reduceVertId(vtkIdType vId)
{
    for(vtkIdType i=0;i<vertIds->GetNumberOfIds();i++)
    {
        vtkIdType oldId = vertIds->GetId(i);
        if(oldId > vId)
            vertIds->SetId(i,oldId-1);
    }
    drawLocalaxes();
    drawBlock();
}

void HexBlock::replacePatch(vtkSmartPointer<HexPatch> oldPatch,
                            vtkSmartPointer<HexPatch> newPatch)
{
    int oldpos = localPatches->IsItemPresent(oldPatch) -1;
    if (oldpos < 0) //not found
        return;
    localPatches->InsertItem(oldpos,newPatch);
}

vtkIdList * HexBlock::commonVertices(HexBlock *hb)
{
    vtkIdList * comVerts;
    for(vtkIdType i=0;i<8;i++)
    {
        //isId returns the id if it's present
        if(this->vertIds->IsId(hb->vertIds->GetId(i)>-1))
            comVerts->InsertNextId(vertIds->GetId(i));
    }
    return comVerts;
}

vtkIdType HexBlock::insertUniqueVertice(vtkSmartPointer<vtkPoints> verts, double pos[3])
{
    double tmpCoords[3];
    vtkIdType ptId;
    vtkIdType nGv=verts->GetNumberOfPoints();
    for(ptId=0; ptId<nGv; ptId++)
    {
        verts->GetPoint(ptId, tmpCoords);
        if(tmpCoords[0]==pos[0] && tmpCoords[1]==pos[1] && tmpCoords[2]==pos[2])
            break;
    }
    if(ptId==nGv)
        verts->InsertNextPoint(pos);
    return ptId;
}

bool HexBlock::hasVertice(vtkIdType vId)
{
    return vertIds->IsId(vId) > -1;
}

bool HexBlock::hasEdge(vtkIdType eId)
{
    HexEdge *ge = HexEdge::SafeDownCast(globalEdges->GetItemAsObject(eId));
    return this->hasEdge(ge);
}

bool HexBlock::hasEdge(HexEdge *e)
{
    return localEdges->IsItemPresent(e)>0 ;
}

bool HexBlock::equals(HexBlock *other)
{
//    std::cout << "\t comparing: this, other" <<std::endl;
    for(vtkIdType i=0;i<vertIds->GetNumberOfIds();i++)
    {
//        std::cout << "\t " << vertIds->GetId(i) << " " << other->vertIds->GetId(i) <<std::endl;
        if( this->vertIds->GetId(i) != other->vertIds->GetId(i))
            return false;
    }

    return true;
}
