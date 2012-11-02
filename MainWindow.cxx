/*
 * Nicos hexBlocker, gpl?
 */



#include "ui_MainWindow.h"
#include "MainWindow.h"
#include "HexBlocker.h"
#include "HexBlock.h" //apparently needed inorder to acces it's actors
#include "InteractorStyleVertPick.h"
#include "InteractorStylePatchPick.h"
#include "ToolBoxWidget.h"
#include "CreateBlockWidget.h"
#include "MoveVerticesWidget.h"
#include "SetBCsWidget.h"
#include "HexBC.h"
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>


#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>

#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkIdList.h>
#include <vtkAreaPicker.h>

#define VTK_CREATE(type, name) \
    vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


// Constructor
MainWindow::MainWindow() 
{
    this->ui = new Ui_MainWindow;
    this->ui->setupUi(this);


    hexBlocker = new HexBlocker();
    hexBlocker->renderer->SetBackground(.2, .3, .4);

    // VTK/Qt wedded
    renwin = this->ui->qvtkWidget->GetRenderWindow();
    renwin->AddRenderer(hexBlocker->renderer);

    // Axes interactor and widget
    axes = vtkSmartPointer<vtkAxesActor>::New();
    widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
    widget->SetOrientationMarker( axes );
    widget->SetInteractor( renwin->GetInteractor() );
    widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
    widget->SetEnabled( 1 );
    widget->InteractiveOff();

    //Area Picker and InteractorStyles
    areaPicker = vtkSmartPointer<vtkAreaPicker>::New();
    styleVertPick = vtkSmartPointer<InteractorStyleVertPick>::New();
    styleVertPick->SetPoints(hexBlocker->vertData);
    styleVertPick->SelectedSphere=hexBlocker->vertSphere;
    stylePatchPick = vtkSmartPointer<InteractorStylePatchPick>::New();
    stylePatchPick->SetPatches(hexBlocker->patches);

    defStyle = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();

    //Qt widgets
    toolbox = new ToolBoxWidget();
    toolbox->setBCsW->hexBCs = hexBlocker->hexBCs;
    this->addDockWidget(Qt::LeftDockWidgetArea,toolbox);


    renwin->Render();

    // Set up action signals and slots
    connect(this->ui->actionOpenFile, SIGNAL(triggered()), this, SLOT(slotOpenFile()));
    connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
    connect(this->ui->actionCreate_HexBlock,SIGNAL(triggered()),this,SLOT(slotOpenCreateHexBlockDialog()));
    connect(toolbox->createBlockW,SIGNAL(apply()),this,SLOT(slotCreateHexBlock()));
    connect(this->ui->actionExtrude_patch,SIGNAL(triggered()),this,SLOT(slotExtrudePatch()));
    connect(this->ui->actionHelp,SIGNAL(triggered()),this,SLOT(slotPrintHexBlocks()));
    connect(this->ui->actionSelectVertices,SIGNAL(triggered()),this,SLOT(slotOpenMoveVerticesDialog()));
    connect(toolbox->moveVerticesW,SIGNAL(moveDone()),this,SLOT(slotResetInteractor()));
    connect(toolbox->moveVerticesW,SIGNAL(moveVertices()),this,SLOT(slotMoveVertices()));
    connect(toolbox,SIGNAL(cancel()),this,SLOT(slotResetInteractor()));
    connect(this->ui->actionSetBCs,SIGNAL(triggered()),this,SLOT(slotOpenSetBCsDialog()));

    connect(toolbox->setBCsW,SIGNAL(startSelectPatches(vtkIdType)),this,SLOT(slotStartSelectPatches(vtkIdType)));
    connect(toolbox->setBCsW,SIGNAL(resetInteractor()), this, SLOT(slotResetInteractor()));
}

MainWindow::~MainWindow()
{
    // The smart pointers should clean up

}

// Action to be taken upon file open 
void MainWindow::slotOpenFile()
{
    hexBlocker->resetBounds();
    renwin->Render();


    //renwin->GetInteractor()->SetInteractorStyle(stylePatchPick);
}

void MainWindow::slotOpenCreateHexBlockDialog()
{
    //set index i toolbox
    toolbox->show();
    toolbox->setCurrentIndex(1);
}

void MainWindow::slotCreateHexBlock()
{
    hexBlocker->createHexBlock(toolbox->createBlockW->c0,toolbox->createBlockW->c1);
    hexBlocker->resetBounds();
    renwin->Render();
}

void MainWindow::slotPrintHexBlocks()
{
    hexBlocker->PrintHexBlocks();
}

void MainWindow::slotExtrudePatch()
{
    hexBlocker->extrudePatch();
    renwin->Render();
}

void MainWindow::slotOpenMoveVerticesDialog()
{
    renwin->GetInteractor()->SetPicker(areaPicker);
    renwin->GetInteractor()->SetInteractorStyle(styleVertPick);

    styleVertPick->StartSelect();

    toolbox->setCurrentIndex(2);

    ui->statusbar->showMessage(tr("Left button to select, press cntrl to rotate"),10000);


}

void MainWindow::slotMoveVertices()
{
    std::cout << "slot is connected, selected IDs: ";
    for(vtkIdType i=0;i<styleVertPick->SelectedList->GetNumberOfIds();i++)
        std::cout << styleVertPick->SelectedList->GetId(i)<< " ";
    std::cout << "."<< std::endl;

    if(toolbox->moveVerticesW->delta)
        hexBlocker->moveVertices(styleVertPick->SelectedList,toolbox->moveVerticesW->dist);

    renwin->Render();
}

void MainWindow::slotResetInteractor()
{
    styleVertPick->SelectedActor->SetVisibility(0);
    styleVertPick->clearSelection();
    renwin->GetInteractor()->SetInteractorStyle(defStyle);

}

void MainWindow::slotOpenSetBCsDialog()
{
    toolbox->setCurrentIndex(3);


}

void MainWindow::slotStartSelectPatches(vtkIdType bcID)
{
    renwin->GetInteractor()->SetInteractorStyle(stylePatchPick);
    connect(stylePatchPick,SIGNAL(selectionDone(vtkIdList *)),
            toolbox->setBCsW,SLOT(slotSelectionDone(vtkIdList*)));

}

void MainWindow::slotExit()
{
    qApp->exit();

}


