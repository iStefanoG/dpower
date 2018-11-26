static char help[] = "This example demonstrates the use of DMNetwork interface for solving a nonlinear electric power grid problem.\n\
                      The available solver options are in the poweroptions file and the data files are in the datafiles directory.\n\
                      The data file format used is from the MatPower package (http://www.pserc.cornell.edu//matpower/).\n\
                      Run this program: mpiexec -n <n> ./pf\n\
                      mpiexec -n <n> ./pfc \n";

/* T
   Concepts: DMNetwork
   Concepts: PETSc SNES solver
*/

#include "power.h"
#include <petscdmnetwork.h>
#include <petsc/private/dmnetworkimpl.h>  /*I  "petscdmnetwork.h"  I*/
#include <petscdmplex.h>
#include <petscsf.h>

PetscErrorCode FormFunction(SNES snes,Vec X, Vec F,void *appctx)
{
  PetscErrorCode ierr;
  DM             networkdm;
  UserCtx_Power  *User=(UserCtx_Power*)appctx;
  Vec            localX,localF;
  PetscInt       nv,ne;
  const PetscInt *vtx,*edges;

  PetscFunctionBegin;
  ierr = SNESGetDM(snes,&networkdm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localF);CHKERRQ(ierr);
  ierr = VecSet(F,0.0);CHKERRQ(ierr);
  ierr = VecSet(localF,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMNetworkGetSubnetworkInfo(networkdm,0,&nv,&ne,&vtx,&edges);CHKERRQ(ierr);
  ierr = FormFunction_Power(networkdm,localX,localF,nv,ne,vtx,edges,User);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SetInitialValues(DM networkdm,Vec X,void* appctx)
{
  PetscErrorCode ierr;
  PetscInt       vStart,vEnd,nv,ne;
  const PetscInt *vtx,*edges;
  Vec            localX;
  UserCtx_Power  *user_power=(UserCtx_Power*)appctx;

  PetscFunctionBegin;
  ierr = DMNetworkGetVertexRange(networkdm,&vStart, &vEnd);CHKERRQ(ierr);

  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = VecSet(X,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = DMNetworkGetSubnetworkInfo(networkdm,0,&nv,&ne,&vtx,&edges);CHKERRQ(ierr);
  ierr = SetInitialGuess_Power(networkdm,localX,nv,ne,vtx,edges,user_power);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMNetworkSetSizes_Parallel(DM dm,PetscInt Nsubnet,PetscInt NsubnetCouple,PetscInt nV[], PetscInt nE[], PetscInt NV[], PetscInt NE[])
{
  PetscErrorCode ierr;
  DM_Network     *network = (DM_Network*) dm->data;
  PetscInt       a[2],b[2],i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  if (Nsubnet <= 0) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Number of subnetworks %D cannot be less than 1",Nsubnet);
  if (NsubnetCouple < 0) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Number of coupling subnetworks %D cannot be less than 0",NsubnetCouple);

  if (Nsubnet > 0) PetscValidLogicalCollectiveInt(dm,Nsubnet,2);
  if (network->nsubnet != 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Network sizes alread set, cannot resize the network");

  network->nsubnet  = Nsubnet + NsubnetCouple;
  Nsubnet          += NsubnetCouple;
  network->ncsubnet = NsubnetCouple;
  ierr = PetscCalloc1(Nsubnet,&network->subnet);CHKERRQ(ierr);
  for(i=0; i < Nsubnet; i++) {
    if (NV) {
      if (NV[i] > 0) PetscValidLogicalCollectiveInt(dm,NV[i],6);
      if (NV[i] > 0 && nV[i] > NV[i]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Subnetwork %D: Local vertex size %D cannot be larger than global vertex size %D",i,nV[i],NV[i]);
    }
    if (NE) {
      if (NE[i] > 0) PetscValidLogicalCollectiveInt(dm,NE[i],7);
      if (NE[i] > 0 && nE[i] > NE[i]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Subnetwork %D: Local edge size %D cannot be larger than global edge size %D",i,nE[i],NE[i]);
    }

    a[0] = nV[i]; a[1] = nE[i];
    ierr = MPIU_Allreduce(a,b,2,MPIU_INT,MPI_SUM,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
    network->subnet[i].Nvtx = b[0]; network->subnet[i].Nedge = b[1];

    network->subnet[i].id = i;

    network->subnet[i].nvtx = nV[i];
    network->subnet[i].vStart = network->nVertices;
    network->subnet[i].vEnd   = network->subnet[i].vStart + network->subnet[i].nvtx;
    network->nVertices += network->subnet[i].nvtx;
    network->NVertices += network->subnet[i].Nvtx;

    network->subnet[i].nedge = nE[i];
    network->subnet[i].eStart = network->nEdges;
    network->subnet[i].eEnd = network->subnet[i].eStart + network->subnet[i].nedge;
    network->nEdges += network->subnet[i].nedge;
    network->NEdges += network->subnet[i].Nedge;
  }
  PetscFunctionReturn(0);
}



PetscErrorCode DMNetworkLayoutSetUp_Parallel(DM dm)
{
  PetscErrorCode ierr;
  DM_Network     *network = (DM_Network*)dm->data;
  PetscInt       dim = 2; /* One dimensional network */
  PetscInt       numCorners=2,spacedim=2;
  double         *vertexcoords=NULL;
  PetscInt       i,j,ndata,ctr=0,nsubnet;
  PetscInt       k,netid,vid;
  PetscInt       *edgelist_couple=NULL;
  PetscSF       sf;
  PetscFunctionBegin;
  MPI_Comm          comm;
  PetscMPIInt       size,rank;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

  if (network->nVertices) {
    ierr = PetscCalloc1(numCorners*network->nVertices,&vertexcoords);CHKERRQ(ierr);
  }

  /* Create the edgelist for the network by concatenating edgelists of the subnetworks */
  nsubnet = network->nsubnet - network->ncsubnet;
  ierr = PetscCalloc1(2*network->nEdges,&network->edges);CHKERRQ(ierr);
  for (i=0; i < nsubnet; i++) {
    for (j = 0; j < network->subnet[i].nedge; j++) {
      network->edges[2*ctr] = network->subnet[i].vStart + network->subnet[i].edgelist[2*j];
      network->edges[2*ctr+1] = network->subnet[i].vStart + network->subnet[i].edgelist[2*j+1];
      ctr++;
    }
  }

  i       = nsubnet; /* netid of coupling subnet */
  nsubnet = network->nsubnet;
  while (i < nsubnet) {
    edgelist_couple = network->subnet[i].edgelist;
    k = 0;
    for (j = 0; j < network->subnet[i].nedge; j++) {
      netid = edgelist_couple[k]; vid = edgelist_couple[k+1];
      network->edges[2*ctr] = network->subnet[netid].vStart + vid; k += 2;
      netid = edgelist_couple[k]; vid = edgelist_couple[k+1];
      network->edges[2*ctr+1] = network->subnet[netid].vStart + vid; k+=2;
      ctr++;
    }
    i++;
  }

#if defined(PETSC_USE_64BIT_INDICES)
  {
    int      *edges;
    PetscInt ii;
    ierr = PetscMalloc1(network->nEdges*numCorners,&edges);CHKERRQ(ierr);
    for (ii=0; ii<network->nEdges*numCorners; ii++) {
      edges[ii] = (int) network->edges[ii];
    }
    ierr = DMPlexCreateFromCellList(PetscObjectComm((PetscObject)dm),dim,network->nEdges,network->nVertices,numCorners,PETSC_FALSE,edges,spacedim,vertexcoords,&network->plex);CHKERRQ(ierr);
    ierr = PetscFree(edges);
  }
#else

  if (size == 1) {
   ierr = DMPlexCreateFromCellList(PetscObjectComm((PetscObject)dm),dim,network->nEdges,network->nVertices,numCorners,PETSC_FALSE,network->edges,spacedim,vertexcoords,&network->plex);CHKERRQ(ierr);
  } else {
    if (rank==0){
   ierr = DMPlexCreateFromCellListParallel(PetscObjectComm((PetscObject)dm),dim,network->nEdges,network->nVertices,numCorners,PETSC_FALSE,network->edges,spacedim,vertexcoords,&sf,&network->plex);CHKERRQ(ierr);
 }
 if (rank==1){
 ierr = DMPlexCreateFromCellListParallel(PetscObjectComm((PetscObject)dm),dim,network->nEdges,network->nVertices-1,numCorners,PETSC_FALSE,network->edges,spacedim,vertexcoords,&sf,&network->plex);CHKERRQ(ierr);
 }
}

 #endif

  if (network->nVertices) {
    ierr = PetscFree(vertexcoords);CHKERRQ(ierr);
  }
  ierr = PetscFree(network->edges);CHKERRQ(ierr);
  ierr = DMPlexGetChart(network->plex,&network->pStart,&network->pEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(network->plex,0,&network->eStart,&network->eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(network->plex,1,&network->vStart,&network->vEnd);CHKERRQ(ierr);
  ierr = PetscSectionCreate(PetscObjectComm((PetscObject)dm),&network->DataSection);CHKERRQ(ierr);
  ierr = PetscSectionCreate(PetscObjectComm((PetscObject)dm),&network->DofSection);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(network->DataSection,network->pStart,network->pEnd);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(network->DofSection,network->pStart,network->pEnd);CHKERRQ(ierr);

  /* Create vertices and edges array for the subnetworks */
  for (j=0; j < network->nsubnet; j++) {
    ierr = PetscCalloc1(network->subnet[j].nedge,&network->subnet[j].edges);CHKERRQ(ierr);
    ierr = PetscCalloc1(network->subnet[j].nvtx,&network->subnet[j].vertices);CHKERRQ(ierr);
    /* Temporarily setting nvtx and nedge to 0 so we can use them as counters in the below for loop.
       These get updated when the vertices and edges are added. */
    network->subnet[j].nvtx = network->subnet[j].nedge = 0;
  }
  network->dataheadersize = sizeof(struct _p_DMNetworkComponentHeader)/sizeof(DMNetworkComponentGenericDataType);
  ierr = PetscCalloc1(network->pEnd-network->pStart,&network->header);CHKERRQ(ierr);
  for (i=network->eStart; i < network->eEnd; i++) {
    network->header[i].index = i;   /* Global edge number */
    for (j=0; j < network->nsubnet; j++) {
      if((network->subnet[j].eStart <= i) && (i < network->subnet[j].eEnd)) {
 network->header[i].subnetid = j; /* Subnetwork id */
 network->subnet[j].edges[network->subnet[j].nedge++] = i;
 break;
      }
    }
    network->header[i].ndata = 0;
    ndata = network->header[i].ndata;
    ierr = PetscSectionAddDof(network->DataSection,i,network->dataheadersize);CHKERRQ(ierr);
    network->header[i].offset[ndata] = 0;
  }

  for(i=network->vStart; i < network->vEnd; i++) {
    network->header[i].index = i - network->vStart;
    for (j=0; j < network->nsubnet; j++) {
      if ((network->subnet[j].vStart <= i-network->vStart) && (i-network->vStart < network->subnet[j].vEnd)) {
 network->header[i].subnetid = j;
 network->subnet[j].vertices[network->subnet[j].nvtx++] = i;
 break;
      }
    }
    network->header[i].ndata = 0;
    ndata = network->header[i].ndata;
    ierr = PetscSectionAddDof(network->DataSection,i,network->dataheadersize);CHKERRQ(ierr);
    network->header[i].offset[ndata] = 0;
  }
  ierr = PetscMalloc1(network->pEnd-network->pStart,&network->cvalue);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMView_Network(DM dm,PetscViewer viewer)
{
  PetscErrorCode ierr;
  DM_Network     *network = (DM_Network*)dm->data;
  PetscBool      iascii;
  PetscMPIInt    rank;
  PetscInt       p,nsubnet,i;
  PetscInt       *edgelist;

  PetscFunctionBegin;

  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank);CHKERRQ(ierr);
  PetscValidHeaderSpecific(dm,DM_CLASSID, 1);
  PetscValidHeaderSpecific(viewer, PETSC_VIEWER_CLASSID, 2);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    nsubnet = network->nsubnet - network->ncsubnet; /* num of subnetworks */
    ierr = PetscViewerASCIIPushSynchronized(viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIISynchronizedPrintf(viewer, "  [%d] nEdges: %D, nVertices: %D, nsubnet: %D\n",rank,network->nEdges,network->nVertices,nsubnet);CHKERRQ(ierr);
    for(i=0; i < nsubnet; i++) {
      edgelist = network->subnet[i].edgelist;
      if (edgelist) {
        ierr = PetscViewerASCIISynchronizedPrintf(viewer, "    [%d]: %D-subnet:\n",rank,i);CHKERRQ(ierr);
        for (p = network->subnet[i].eStart; p < network->subnet[i].eEnd; p++) {
          ierr = PetscViewerASCIISynchronizedPrintf(viewer, "     [%d]: %D ----> %D\n",rank,edgelist[2*p],edgelist[2*p+1]);CHKERRQ(ierr);
        }
      }
    }

    ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopSynchronized(viewer);CHKERRQ(ierr);

  } else SETERRQ1(PetscObjectComm((PetscObject) dm), PETSC_ERR_SUP, "Viewer type %s not yet supported for DMNetwork writing", ((PetscObject)viewer)->type_name);
  PetscFunctionReturn(0);
}



int main(int argc,char ** argv)
{
  PetscErrorCode   ierr;
  char             pfdata_file0[PETSC_MAX_PATH_LEN]="case9.m";
  char             pfdata_file1[PETSC_MAX_PATH_LEN]="case5.m";
  PFDATA           *pfdata0;
  PFDATA           *pfdata1;
  PetscInt         numEdges=0,numVertices=0,NumEdges=PETSC_DETERMINE,NumVertices=PETSC_DETERMINE;
  PetscInt         *edges = NULL;
  PetscInt         i;
  DM               networkdm;
  UserCtx_Power    User;
  PetscLogStage    stage1,stage2;
  PetscMPIInt      rank,size;
  PetscInt         eStart, eEnd, vStart, vEnd,j;
  PetscInt         genj,loadj;
  Vec              X,F;
  Mat              J;
  SNES             snes;

  ierr = PetscInitialize(&argc,&argv,"poweroptions",help);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  {
    /* introduce the const crank so the clang static analyzer realizes that if it enters any of the if (crank) then it must have entered the first */
    /* this is an experiment to see how the analyzer reacts */
    const PetscMPIInt crank = rank;
    /* Create an empty network object */
    ierr = DMNetworkCreate(PETSC_COMM_WORLD,&networkdm);CHKERRQ(ierr);
    /* Register the components in the network */
    ierr = DMNetworkRegisterComponent(networkdm,"branchstruct",sizeof(struct _p_EDGE_Power),&User.compkey_branch);CHKERRQ(ierr);
    ierr = DMNetworkRegisterComponent(networkdm,"busstruct",sizeof(struct _p_VERTEX_Power),&User.compkey_bus);CHKERRQ(ierr);
    ierr = DMNetworkRegisterComponent(networkdm,"genstruct",sizeof(struct _p_GEN),&User.compkey_gen);CHKERRQ(ierr);
    ierr = DMNetworkRegisterComponent(networkdm,"loadstruct",sizeof(struct _p_LOAD),&User.compkey_load);CHKERRQ(ierr);

    ierr = PetscLogStageRegister("Read Data",&stage1);CHKERRQ(ierr);
    PetscLogStagePush(stage1);
    /* READ THE DATA */
    if (crank==0) {
      /*    READ DATA */
      /* Only rank 0 reads the data */
      /*ierr = PetscOptionsGetString(NULL,NULL,"-pfdata",pfdata_file,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);*/
      ierr = PetscOptionsGetString(NULL,NULL,"-pfdata0",pfdata_file0,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);
      ierr=PetscPrintf(PETSC_COMM_SELF,"Hi, I am pfdata_file: %s\n",pfdata_file0);CHKERRQ(ierr);
      ierr = PetscNew(&pfdata0);CHKERRQ(ierr);
      ierr = PFReadMatPowerData(pfdata0,pfdata_file0);CHKERRQ(ierr);
      User.Sbase = pfdata0->sbase;
      numEdges = pfdata0->nbranch;
      numVertices = pfdata0->nbus;
      ierr = PetscMalloc1(2*numEdges,&edges);CHKERRQ(ierr);
      ierr = GetListofEdges_Power(networkdm,pfdata0,edges);CHKERRQ(ierr);

    }
    if (crank==1) {
      ierr = PetscOptionsGetString(NULL,NULL,"-pfdata1",pfdata_file1,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);
      ierr=PetscPrintf(PETSC_COMM_SELF,"Hi, I am pfdata_file: %s\n",pfdata_file1);CHKERRQ(ierr);
      ierr = PetscNew(&pfdata1);CHKERRQ(ierr);
      ierr = PFReadMatPowerData(pfdata1,pfdata_file1);CHKERRQ(ierr);
      User.Sbase = pfdata1->sbase;
      numEdges = pfdata1->nbranch;
      numVertices = pfdata1->nbus;
      ierr = PetscMalloc1(2*numEdges,&edges);CHKERRQ(ierr);
      ierr = GetListofEdges_Power(networkdm,pfdata1,edges);CHKERRQ(ierr);
    }

    /* If external option activated. Introduce error in jacobian */
    ierr = PetscOptionsHasName(NULL,NULL, "-jac_error", &User.jac_error);CHKERRQ(ierr);
    PetscLogStagePop();
    ierr = PetscLogStageRegister("Create network",&stage2);CHKERRQ(ierr);
    PetscLogStagePush(stage2);
    /* Set number of nodes/edges */
    ierr = DMNetworkSetSizes_Parallel(networkdm,1,0,&numVertices,&numEdges,&NumVertices,&NumEdges);CHKERRQ(ierr);
    /* Add edge connectivity */
    ierr = DMNetworkSetEdgeList(networkdm,&edges,NULL);CHKERRQ(ierr);
    /* Set up the network layout */
    ierr = DMNetworkLayoutSetUp_Parallel(networkdm);CHKERRQ(ierr);
    ierr = PetscFree(edges);CHKERRQ(ierr);

    /* Add network components only process 0 has any data to add*/
    if (crank==0) {
      genj=0; loadj=0;
      ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
      for (i = eStart; i < eEnd; i++) {
        ierr = DMNetworkAddComponent(networkdm,i,User.compkey_branch,&pfdata0->branch[i-eStart]);CHKERRQ(ierr);
      }
      ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
      for (i = vStart; i < vEnd; i++) {
        ierr = DMNetworkAddComponent(networkdm,i,User.compkey_bus,&pfdata0->bus[i-vStart]);CHKERRQ(ierr);
        if (pfdata0->bus[i-vStart].ngen) {
          for (j = 0; j < pfdata0->bus[i-vStart].ngen; j++) {
            ierr = DMNetworkAddComponent(networkdm,i,User.compkey_gen,&pfdata0->gen[genj++]);CHKERRQ(ierr);
          }
        }
        if (pfdata0->bus[i-vStart].nload) {
          for (j=0; j < pfdata0->bus[i-vStart].nload; j++) {
            ierr = DMNetworkAddComponent(networkdm,i,User.compkey_load,&pfdata0->load[loadj++]);CHKERRQ(ierr);
          }
        }
        /* Add number of variables */
        ierr = DMNetworkAddNumVariables(networkdm,i,2);CHKERRQ(ierr);
      }
    }
    if (crank==1) {
      genj=0; loadj=0;
      ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
      for (i = eStart; i < eEnd; i++) {
        ierr = DMNetworkAddComponent(networkdm,i,User.compkey_branch,&pfdata1->branch[i-eStart]);CHKERRQ(ierr);
      }
      ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
      for (i = vStart; i < vEnd; i++) {
        ierr = DMNetworkAddComponent(networkdm,i,User.compkey_bus,&pfdata1->bus[i-vStart]);CHKERRQ(ierr);
        if (pfdata1->bus[i-vStart].ngen) {
          for (j = 0; j < pfdata1->bus[i-vStart].ngen; j++) {
            ierr = DMNetworkAddComponent(networkdm,i,User.compkey_gen,&pfdata1->gen[genj++]);CHKERRQ(ierr);
          }
        }
        if (pfdata1->bus[i-vStart].nload) {
          for (j=0; j < pfdata1->bus[i-vStart].nload; j++) {
            ierr = DMNetworkAddComponent(networkdm,i,User.compkey_load,&pfdata1->load[loadj++]);CHKERRQ(ierr);
          }
        }
        /* Add number of variables */
        ierr = DMNetworkAddNumVariables(networkdm,i,2);CHKERRQ(ierr);
      }
    }

    /* Set up DM for use */
    ierr = DMSetUp(networkdm);CHKERRQ(ierr);
    // ierr = DMView_Network(networkdm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    if (crank==0) {
      ierr = PetscFree(pfdata0->bus);CHKERRQ(ierr);
      ierr = PetscFree(pfdata0->gen);CHKERRQ(ierr);
      ierr = PetscFree(pfdata0->branch);CHKERRQ(ierr);
      ierr = PetscFree(pfdata0->load);CHKERRQ(ierr);
      ierr = PetscFree(pfdata0);CHKERRQ(ierr);
    }
    if (crank==1){
      ierr = PetscFree(pfdata1->bus);CHKERRQ(ierr);
      ierr = PetscFree(pfdata1->gen);CHKERRQ(ierr);
      ierr = PetscFree(pfdata1->branch);CHKERRQ(ierr);
      ierr = PetscFree(pfdata1->load);CHKERRQ(ierr);
      ierr = PetscFree(pfdata1);CHKERRQ(ierr);

    }
    PetscLogStagePop();

#if 0
    EDGE_Power     edge;
    PetscInt       key,kk,numComponents;
    VERTEX_Power   bus;
    GEN            gen;
    LOAD           load;

    for (i = eStart; i < eEnd; i++) {
      ierr = DMNetworkGetComponent(networkdm,i,0,&key,(void**)&edge);CHKERRQ(ierr);
      ierr = DMNetworkGetNumComponents(networkdm,i,&numComponents);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d ncomps = %d Line %d ---- %d\n",crank,numComponents,edge->internal_i,edge->internal_j);CHKERRQ(ierr);
    }

    for (i = vStart; i < vEnd; i++) {
      ierr = DMNetworkGetNumComponents(networkdm,i,&numComponents);CHKERRQ(ierr);
      for (kk=0; kk < numComponents; kk++) {
        ierr = DMNetworkGetComponent(networkdm,i,kk,&key,&component);CHKERRQ(ierr);
        if (key == 1) {
          bus = (VERTEX_Power)(component);
          ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d ncomps = %d Bus %d\n",crank,numComponents,bus->internal_i);CHKERRQ(ierr);
        } else if (key == 2) {
          gen = (GEN)(component);
          ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Gen pg = %f qg = %f\n",crank,gen->pg,gen->qg);CHKERRQ(ierr);
        } else if (key == 3) {
          load = (LOAD)(component);
          ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Load pl = %f ql = %f\n",crank,load->pl,load->ql);CHKERRQ(ierr);
        }
      }
    }
#endif
    /* Broadcast Sbase to all processors */
    ierr = MPI_Bcast(&User.Sbase,1,MPIU_SCALAR,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(networkdm,&X);CHKERRQ(ierr);
    ierr = VecDuplicate(X,&F);CHKERRQ(ierr);
    ierr = DMCreateMatrix(networkdm,&J);CHKERRQ(ierr);
    ierr = MatSetOption(J,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
    ierr = SetInitialValues(networkdm,X,&User);CHKERRQ(ierr);

    /* HOOK UP SOLVER */
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
    ierr = SNESSetDM(snes,networkdm);CHKERRQ(ierr);
    ierr = SNESSetFunction(snes,F,FormFunction,&User);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,J,J,FormJacobian_Power,&User);CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
    ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
  /*  ierr = VecView(X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = MatView(J,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);*/
    ierr = VecDestroy(&X);CHKERRQ(ierr);
    ierr = VecDestroy(&F);CHKERRQ(ierr);
    ierr = MatDestroy(&J);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);
    ierr = DMDestroy(&networkdm);CHKERRQ(ierr);
  }
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

   build:
     depends: PFReadData.c pffunctions.c
     requires: !complex double define(PETSC_HAVE_ATTRIBUTEALIGNED)


   test:
     args: -snes_rtol 1.e-3
     localrunfiles: poweroptions case9.m
     output_file: output/power_1.out
     requires: double !complex define(PETSC_HAVE_ATTRIBUTEALIGNED)

   test:
     suffix: 2
     args: -snes_rtol 1.e-3 -petscpartitioner_type simple
     nsize: 4
     localrunfiles: poweroptions case9.m
     output_file: output/power_1.out
     requires: double !complex define(PETSC_HAVE_ATTRIBUTEALIGNED)

TEST*/
