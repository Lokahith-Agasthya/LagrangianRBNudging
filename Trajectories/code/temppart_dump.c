#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//Modified file to use for reading Field temperatures and nudging the second population
#define AMIROOT (!me)
#include "hdf5.h"
#define REAL double
#define npop 9
#define cs2 (1./3.)
#define cs1 sqrt((cs2))
#define cs22 (2.*cs2)
#define cssq (2./9.)
#define pi 3.14159265359

/////////INPUT VALUES/////////
#define nx         1680
#define ny         840
#define scratch    0
#define nstart     0
#define nsteps     5000000
#define nout       20000
#define noutconfig 250
#define rhom1      1.0
#define relax      0.502
#define relaxg     0.502
#define Tup        -0.01
#define Tdown      0.01
#define gravity    (-0.00004)

#define Particle   1
#define Particle_scratch 1
#include "mympi2.h"
///////////////////////////////////////////////////////////////
//chi trova errori contatti Lokahith Agasthya Tel. 3895055105//
///////////////mails: lnagasthya@gmail.com////////////////////
//////////////////////////////////////////////////////////////

typedef struct {
  double p[npop];
} pop_type;

REAL cx[npop];
REAL cy[npop];

int idx0,idx1,idx2,idx3,idx4;
int icount=0,icountconfig=0, icountprof=0;
REAL ww[npop];
int indx,indy;
int istep;
int i,j,k,kk,s,ip,an;

void dumpf(pop_type *f) {
  char fname[128];
  herr_t status;
  double no_halo[NX][NY][npop];
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      int idx1 = j+(NY+2)*i;
      for(k=0;k<npop;k++){
	no_halo[i-1][j-1][k] = f[idx1].p[k];
      }
    }
  }
  H5Eset_current_stack(H5E_DEFAULT);
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hid_t hdf5_status = H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
  hsize_t edimens_3d[3] = {NX,NY,npop};
  hid_t ememspace = H5Screate_simple(3, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  edimens_3d[2] = npop;
  hid_t efilespace = H5Screate_simple( 3, edimens_3d, NULL );
  hsize_t estart_3d[3] = {0,0,0};
  hsize_t estride_3d[3] = {1,1,1};
  hsize_t eblock_3d[3] = {NX,NY,npop};
  hsize_t ecount_3d[3] = {1,1,1};
  hsize_t dstart_3d[3] = {mex * NX,mey * NY,0};
  hsize_t dstride_3d[3] = {1,1,1};
  hsize_t dblock_3d[3] = {NX,NY,npop};
  hsize_t dcount_3d[3] = {1,1,1};
  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
  sprintf(fname,"confnew%d.h5",istep);
  hid_t  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  hid_t dataset_id = H5Dcreate(file_id, "/f", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t  xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  hid_t ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, no_halo);
  status = H5Pclose( xfer_plist );
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  H5Sclose (efilespace);
  H5Sclose (ememspace);
  status = H5Pclose( plist_id );
}

void dumpg(pop_type *f) {
  char fname[128];
  herr_t status;
  double no_halo[NX][NY][npop];
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      int idx1 = j+(NY+2)*i;
      for(k=0;k<npop;k++){
	no_halo[i-1][j-1][k] = f[idx1].p[k];
      }
    }
  }
  H5Eset_current_stack(H5E_DEFAULT);
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hid_t hdf5_status = H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
  hsize_t edimens_3d[3] = {NX,NY,npop};
  hid_t ememspace = H5Screate_simple(3, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  edimens_3d[2] = npop;
  hid_t efilespace = H5Screate_simple( 3, edimens_3d, NULL );
  hsize_t estart_3d[3] = {0,0,0};
  hsize_t estride_3d[3] = {1,1,1};
  hsize_t eblock_3d[3] = {NX,NY,npop};
  hsize_t ecount_3d[3] = {1,1,1};
  hsize_t dstart_3d[3] = {mex * NX,mey * NY,0};
  hsize_t dstride_3d[3] = {1,1,1};
  hsize_t dblock_3d[3] = {NX,NY,npop};
  hsize_t dcount_3d[3] = {1,1,1};
  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
  sprintf(fname,"congnew%d.h5",istep);
  hid_t  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  hid_t dataset_id = H5Dcreate(file_id, "/g", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t  xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  hid_t ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, no_halo);

  status = H5Pclose( xfer_plist );
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  H5Sclose (efilespace);
  H5Sclose (ememspace);
  status = H5Pclose( plist_id );
}

void dumpfpre(pop_type *f) {
  char fname[128];
  herr_t status;
  double no_halo[NX][NY][npop];
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      int idx1 = j+(NY+2)*i;
      for(k=0;k<npop;k++){
	no_halo[i-1][j-1][k] = f[idx1].p[k];
      }
    }
  }
  H5Eset_current_stack(H5E_DEFAULT);
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hid_t hdf5_status = H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
  hsize_t edimens_3d[3] = {NX,NY,npop};
  hid_t ememspace = H5Screate_simple(3, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  edimens_3d[2] = npop;
  hid_t efilespace = H5Screate_simple( 3, edimens_3d, NULL );
  hsize_t estart_3d[3] = {0,0,0};
  hsize_t estride_3d[3] = {1,1,1};
  hsize_t eblock_3d[3] = {NX,NY,npop};
  hsize_t ecount_3d[3] = {1,1,1};
  hsize_t dstart_3d[3] = {mex * NX,mey * NY,0};
  hsize_t dstride_3d[3] = {1,1,1};
  hsize_t dblock_3d[3] = {NX,NY,npop};
  hsize_t dcount_3d[3] = {1,1,1};
  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
  sprintf(fname,"confpre%d.h5",istep);
  hid_t  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  hid_t dataset_id = H5Dcreate(file_id, "/f", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t  xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  hid_t ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, no_halo);
  status = H5Pclose( xfer_plist );
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  H5Sclose (efilespace);
  H5Sclose (ememspace);
  status = H5Pclose( plist_id );
}

void dumpgpre(pop_type *f) {
  char fname[128];
  herr_t status;
  double no_halo[NX][NY][npop];
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      int idx1 = j+(NY+2)*i;
      for(k=0;k<npop;k++){
	no_halo[i-1][j-1][k] = f[idx1].p[k];
      }
    }
  }
  H5Eset_current_stack(H5E_DEFAULT);
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hid_t hdf5_status = H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
  hsize_t edimens_3d[3] = {NX,NY,npop};
  hid_t ememspace = H5Screate_simple(3, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  edimens_3d[2] = npop;
  hid_t efilespace = H5Screate_simple( 3, edimens_3d, NULL );
  hsize_t estart_3d[3] = {0,0,0};
  hsize_t estride_3d[3] = {1,1,1};
  hsize_t eblock_3d[3] = {NX,NY,npop};
  hsize_t ecount_3d[3] = {1,1,1};
  hsize_t dstart_3d[3] = {mex * NX,mey * NY,0};
  hsize_t dstride_3d[3] = {1,1,1};
  hsize_t dblock_3d[3] = {NX,NY,npop};
  hsize_t dcount_3d[3] = {1,1,1};
  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
  sprintf(fname,"congpre%d.h5",istep);
  hid_t  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  hid_t dataset_id = H5Dcreate(file_id, "/g", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t  xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  hid_t ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, no_halo);

  status = H5Pclose( xfer_plist );
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  H5Sclose (efilespace);
  H5Sclose (ememspace);
  status = H5Pclose( plist_id );
}

void pbc (pop_type *f) {

  int procpx, procmx, procpy, procmy;
  MPI_Status status;

  procpx=(mex+1+NPX)%NPX; procmx=(mex-1+NPX)%NPX;
  procpy=(mey+1+NPY)%NPY; procmy=(mey-1+NPY)%NPY;

  MPI_Sendrecv(f+(2+NY)*NX,1,MPI_Poptypex,procpx,10,f,1,MPI_Poptypex,procmx,10,MPI_COMM_ALONG_X,&status);

  MPI_Sendrecv(f+(2+NY),1,MPI_Poptypex,procmx,11,f+(2+NY)*(NX+1),1,MPI_Poptypex,procpx,11,MPI_COMM_ALONG_X,&status);

  MPI_Sendrecv(f+NY,1,MPI_Poptypey,procpy,12,f,1,MPI_Poptypey,procmy,12,MPI_COMM_ALONG_Y,&status);

  MPI_Sendrecv(f+1,1,MPI_Poptypey,procmy,13,f+(1+NY),1,MPI_Poptypey,procpy,13,MPI_COMM_ALONG_Y,&status);

}

void restore_prtcl(prtcl_type *prtcl1,int *prtcl_ids1){
  prtcl_type *prtcl2;
  int prtcl_ids2[Nprtcl];
  prtcl2 = (prtcl_type *) malloc((Nprtcl)*sizeof(prtcl_type));
  //  prtcl_ids2 = (int *) malloc(Nprtcl*sizeof(int));
  int np_loc = 0;
  int send_add = Nprtcl;
  if (me==0){
    np_loc = Nprtcl;
    send_add= 0;
  }
  char fname[128];
  H5Eset_current_stack(H5E_DEFAULT);  
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hid_t hdf5_status = H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
  hsize_t edimens_3d[1] = {np_loc};
  hsize_t fdimens_3d[1] = {Nprtcl};

  hsize_t estart_3d[1] = {0};
  hsize_t estride_3d[1] = {1};
  hsize_t eblock_3d[1] = {np_loc};
  hsize_t ecount_3d[1] = {1};

  hsize_t dstart_3d[1] = {send_add};//offset
  hsize_t dstride_3d[1] = {1};
  hsize_t dblock_3d[1] = {np_loc};//local
  hsize_t dcount_3d[1] = {1};
  hid_t dataset_id,xfer_plist,ret;

  sprintf(fname,"Prtcl_Veldata1/Particles%d.h5",nstart);
  if (nstart>10000000) sprintf(fname,"Prtcl_Veldata2/Particles%d.h5",nstart);
  hid_t  file_id = H5Fopen(fname, H5F_ACC_RDWR,H5P_DEFAULT);
  if (file_id<0){
    printf("Particle file not opened \n");
    exit(8);
  }
  hid_t ememspace = H5Screate_simple(1, edimens_3d, NULL );
  hid_t efilespace = H5Screate_simple( 1, fdimens_3d, NULL );
  hid_t status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
    
  dataset_id = H5Dopen (file_id, "/id", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, ememspace, efilespace, H5P_DEFAULT, prtcl_ids2);
  status = H5Dclose(dataset_id);

  H5Sclose (efilespace);
  H5Sclose (ememspace);
    
  status = H5Fclose(file_id);
  status = H5Pclose( plist_id );

  int loc_prtcl[Nprtcl];
  if (me==0){
    FILE *fin;
    char filename[128];
    sprintf(filename,"prtcl_dump%d.in",nstart);
    fin = fopen(filename,"r");
    fread(prtcl2,sizeof(prtcl_type),Nprtcl,fin);
    fclose(fin);
  }
  MPI_Bcast(prtcl2,Nprtcl,MPI_Prtcltype,0,MPI_COMM_WORLD);
  MPI_Bcast(prtcl_ids2,Nprtcl,MPI_INT,0,MPI_COMM_WORLD);
  nprtcl_local = 0;
  REAL xx,yy;
  
  for (i=0;i<Nprtcl;i++){
    loc_prtcl[i] = 0;
    xx = prtcl2[i].xprtcl;
    yy = prtcl2[i].yprtcl;
    if ((xx > (REAL) (mex*NX))&& (yy > (REAL) (mey*NY))){
      if ((xx < (REAL) (mex*NX + NX))&& (yy < (REAL) (mey*NY + NY))){
	nprtcl_local+=1;
	loc_prtcl[i] = 1;
      }
    }
  }
  j = 0;
  for (i=0;i<Nprtcl;i++){
    if (loc_prtcl[i]==1){
      prtcl1[j].xprtcl = prtcl2[i].xprtcl - (REAL) (mex*NX);
      prtcl1[j].yprtcl = prtcl2[i].yprtcl - (REAL) (mey*NY);
      prtcl1[j].uprtcl = prtcl2[i].uprtcl;
      prtcl1[j].vprtcl = prtcl2[i].vprtcl;
      prtcl1[j].xoprtcl = prtcl2[i].xprtcl;
      prtcl1[j].yoprtcl = prtcl2[i].yprtcl;
      prtcl1[j].uoprtcl = prtcl2[i].uoprtcl;
      prtcl1[j].voprtcl = prtcl2[i].voprtcl;
      //prtcl[j].prtcl_id = prtcl2[i].prtcl_id;
      prtcl_ids1[j] = prtcl_ids2[i];
      j = j+1;
    }
  }
  //  free(prtcl_ids2);
  free(prtcl2);
}

void restoref(pop_type *f){
  double no_halo[NX][NY][npop];
  char fname[128];

  hid_t efilespace;
  hid_t plist_id;               /* property list identifier */
  hid_t file_id, dataset_id;  // identifiers 
  hid_t xfer_plist;
  herr_t ret;
  herr_t status;
  hid_t ememspace;
  hid_t hdf5_status;
  hsize_t  estart_3d[3] = {0,0,0};
  hsize_t  estride_3d[3] = {1,1,1};
  hsize_t  eblock_3d[3] = {NX,NY,npop};
  hsize_t  ecount_3d[3] = {1,1,1};
  hsize_t  dstart_3d[3] = {mex * NX,mey * NY,0}; 
  hsize_t  dstride_3d[3] = {1,1,1};
  hsize_t  dblock_3d[3] = {NX,NY,npop};
  hsize_t  dcount_3d[3] = {1,1,1};
  hsize_t  edimens_3d[3] = {NX,NY,npop};

  H5Eset_current_stack (H5E_DEFAULT);

  plist_id = H5Pcreate (H5P_FILE_ACCESS);
  hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/High_Ra/conf.h5",nstart);
  //sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/benchmark/confnew%d.h5",nstart);
  //sprintf(fname,"confnew%d.h5",nstart);
  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id<0){
    printf("conf not found \n");
    exit(8);
  }
  ememspace = H5Screate_simple( 3, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  edimens_3d[2] = npop;
  
  efilespace = H5Screate_simple( 3, edimens_3d, NULL );

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);

  dataset_id = H5Dopen (file_id, "/f", H5P_DEFAULT);

  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, no_halo);

  status = H5Dclose(dataset_id);
  H5Sclose (efilespace);
  H5Sclose (ememspace);
  status = H5Fclose(file_id);
  status = H5Pclose( plist_id );

  int idx = 0;
  int idy = 0;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      int idx1 = j+(NY+2)*i;

      for(k=0;k<npop;k++){
	f[idx1].p[k] = no_halo[i-1][j-1][k];
      }
    }
  }
  pbc(f);
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      for(k=0;k<npop;k++){
	if ((f[idx1].p[k])==0){
	  printf("Damaar f in the initial function only \n");
	  printf("%d %d %d \n",i,j,k);
	  exit(7);
	}
      }
    }
  }
}


void restoreg(pop_type *f){
  double no_halo[NX][NY][npop];
  char fname[128];
  hid_t efilespace;
  hid_t plist_id;               /* property list identifier */
  hid_t file_id, dataset_id;  // identifiers 
  hid_t xfer_plist;
  herr_t ret;
  herr_t status;
  hid_t ememspace;
  hid_t hdf5_status;
  hsize_t  estart_3d[3] = {0,0,0};
  hsize_t  estride_3d[3] = {1,1,1};
  hsize_t  eblock_3d[3] = {NX,NY,npop};
  hsize_t  ecount_3d[3] = {1,1,1};
  hsize_t  dstart_3d[3] = {mex * NX,mey * NY,0}; 
  hsize_t  dstride_3d[3] = {1,1,1};
  hsize_t  dblock_3d[3] = {NX,NY,npop};
  hsize_t  dcount_3d[3] = {1,1,1};
  hsize_t  edimens_3d[3] = {NX,NY,npop};
  
  H5Eset_current_stack (H5E_DEFAULT);

  plist_id = H5Pcreate (H5P_FILE_ACCESS);
  hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/High_Ra/cong.h5",nstart);
  //sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/benchmark/congnew%d.h5",nstart);
  //sprintf(fname,"congnew%d.h5",nstart);

  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id<0){
    printf("cong not found \n");
    exit(8);
  }
  ememspace = H5Screate_simple( 3, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  edimens_3d[2] = npop;
  
  efilespace = H5Screate_simple( 3, edimens_3d, NULL );

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);

  dataset_id = H5Dopen (file_id, "/g", H5P_DEFAULT);

  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, no_halo);

  status = H5Dclose(dataset_id);
  H5Sclose (efilespace);
  H5Sclose (ememspace);
  status = H5Fclose(file_id);
  status = H5Pclose( plist_id );

  int idx = 0;
  int idy = 0;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      int idx1 = j+(NY+2)*i;

      for(k=0;k<npop;k++){
	f[idx1].p[k] = no_halo[i-1][j-1][k];
      }
    }
  }
  pbc(f);
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      for(k=0;k<npop;k++){
	if ((f[idx1].p[k])==0){
	  printf("Damaar g in the initial function only \n");
	  printf("%d %d %d \n",i,j,k);
	  exit(7);
	}
      }
    }
  }
}

void restorefpre(pop_type *f){
  double no_halo[NX][NY][npop];
  char fname[128];

  hid_t efilespace;
  hid_t plist_id;               /* property list identifier */
  hid_t file_id, dataset_id;  // identifiers 
  hid_t xfer_plist;
  herr_t ret;
  herr_t status;
  hid_t ememspace;
  hid_t hdf5_status;
  hsize_t  estart_3d[3] = {0,0,0};
  hsize_t  estride_3d[3] = {1,1,1};
  hsize_t  eblock_3d[3] = {NX,NY,npop};
  hsize_t  ecount_3d[3] = {1,1,1};
  hsize_t  dstart_3d[3] = {mex * NX,mey * NY,0}; 
  hsize_t  dstride_3d[3] = {1,1,1};
  hsize_t  dblock_3d[3] = {NX,NY,npop};
  hsize_t  dcount_3d[3] = {1,1,1};
  hsize_t  edimens_3d[3] = {NX,NY,npop};

  H5Eset_current_stack (H5E_DEFAULT);

  plist_id = H5Pcreate (H5P_FILE_ACCESS);
  hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/High_Ra/confpre.h5",nstart);
  //sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/benchmark/confpre%d.h5",nstart);
  //sprintf(fname,"confpre%d.h5",nstart);
  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id<0){
    printf("conf not found \n");
    exit(8);
  }
  ememspace = H5Screate_simple( 3, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  edimens_3d[2] = npop;
  
  efilespace = H5Screate_simple( 3, edimens_3d, NULL );

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);

  dataset_id = H5Dopen (file_id, "/f", H5P_DEFAULT);

  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, no_halo);

  status = H5Dclose(dataset_id);
  H5Sclose (efilespace);
  H5Sclose (ememspace);
  status = H5Fclose(file_id);
  status = H5Pclose( plist_id );

  int idx = 0;
  int idy = 0;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      int idx1 = j+(NY+2)*i;

      for(k=0;k<npop;k++){
	f[idx1].p[k] = no_halo[i-1][j-1][k];
      }
    }
  }
  pbc(f);
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      for(k=0;k<npop;k++){
	if ((f[idx1].p[k])==0){
	  printf("Damaar f in the initial function only \n");
	  printf("%d %d %d \n",i,j,k);
	  exit(7);
	}
      }
    }
  }
}


void restoregpre(pop_type *f){
  double no_halo[NX][NY][npop];
  char fname[128];
  hid_t efilespace;
  hid_t plist_id;               /* property list identifier */
  hid_t file_id, dataset_id;  // identifiers 
  hid_t xfer_plist;
  herr_t ret;
  herr_t status;
  hid_t ememspace;
  hid_t hdf5_status;
  hsize_t  estart_3d[3] = {0,0,0};
  hsize_t  estride_3d[3] = {1,1,1};
  hsize_t  eblock_3d[3] = {NX,NY,npop};
  hsize_t  ecount_3d[3] = {1,1,1};
  hsize_t  dstart_3d[3] = {mex * NX,mey * NY,0}; 
  hsize_t  dstride_3d[3] = {1,1,1};
  hsize_t  dblock_3d[3] = {NX,NY,npop};
  hsize_t  dcount_3d[3] = {1,1,1};
  hsize_t  edimens_3d[3] = {NX,NY,npop};
  
  H5Eset_current_stack (H5E_DEFAULT);

  plist_id = H5Pcreate (H5P_FILE_ACCESS);
  hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/High_Ra/congpre.h5",nstart);
  //sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/benchmark/congpre%d.h5",nstart);
  //sprintf(fname,"congpre%d.h5",nstart);

  file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id<0){
    printf("cong not found \n");
    exit(8);
  }
  ememspace = H5Screate_simple( 3, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  edimens_3d[2] = npop;
  
  efilespace = H5Screate_simple( 3, edimens_3d, NULL );

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);

  dataset_id = H5Dopen (file_id, "/g", H5P_DEFAULT);

  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, no_halo);

  status = H5Dclose(dataset_id);
  H5Sclose (efilespace);
  H5Sclose (ememspace);
  status = H5Fclose(file_id);
  status = H5Pclose( plist_id );

  int idx = 0;
  int idy = 0;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      int idx1 = j+(NY+2)*i;

      for(k=0;k<npop;k++){
	f[idx1].p[k] = no_halo[i-1][j-1][k];
      }
    }
  }
  pbc(f);
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      for(k=0;k<npop;k++){
	if ((f[idx1].p[k])==0){
	  printf("Damaar g in the initial function only \n");
	  printf("%d %d %d \n",i,j,k);
	  exit(7);
	}
      }
    }
  }
}

void mirror(REAL *upre, REAL *ubefore){ //copy from upre to ubefore
  for (i=0;i<NX+2;i++){
    for (j=0;j<NY+2;j++){
      idx1=j+(NY+2)*i;
      ubefore[idx1] = upre[idx1];
    }
  }
}


void inithydro(REAL *u,REAL *v,REAL *rho,REAL *temp,pop_type *f,pop_type *g){//INIZIO ROUTINE INITHYDRO
  FILE *fout,*fout1;
  REAL yyytemp,y_c,fl[npop],rhoi;
  int indx,indy;
  char filename[128], filename1[128];
  int c,cont1;
  REAL randnorm = 1.0/RAND_MAX;
#define epstan (10.)
  srand(me);
  if (scratch>0){
    for(i=0;i<NX+2;i++){
      for(j=0;j<NY+2;j++){
	idx1=j+(NY+2)*i;
	rho[idx1]=rhom1;
	yyytemp=((double)(j-1+(NY)*mey));
	y_c= (double)ny*(1)/2.+(randnorm*rand()-0.5)*2.*epstan;
	temp[idx1]=Tdown+(Tup-Tdown)*(tanh((yyytemp-y_c)/epstan)+1.0)/2.0;
	// temp[idx1]=Tdown+(Tup-Tdown)*((double)j-1+NY*mey)/((double)(ny-1));
	v[idx1]=0;
	u[idx1]=0; 
      }
    }
  }else{
    for(i=0;i<NX+2;i++){
      for(j=0;j<NY+2;j++){
	idx1=j+(NY+2)*i;
	
	rho[idx1]=0.0;
	temp[idx1]=0.0;
	for(k=0;k<npop;k++){
	  rho[idx1]=rho[idx1]+f[idx1].p[k];
	}
	for(k=0;k<npop;k++){
	  temp[idx1]=temp[idx1]+g[idx1].p[k];
	}

	rhoi=1./rho[idx1];

	for(k=0;k<npop;k++){
	  fl[k]=f[idx1].p[k]*rhoi;
	}

	u[idx1]=fl[1]-fl[3]+fl[5]-fl[6]-fl[7]+fl[8];
	v[idx1]=fl[5]+fl[2]+fl[6]-fl[7]-fl[4]-fl[8] + (0.5*(-gravity)*temp[idx1]);

      }
    }
  }
  communication(u);
  communication(v);
}//FINE ROUTINE INITHYDRO

void initpop(pop_type *f,pop_type *feq){////ROUTINE INITPOP
  FILE *fout;
  char filename[128];

  for(i=0;i<NX+2;i++){
    for(j=0;j<NY+2;j++){
      idx1=j+(NY+2)*i;
      f[idx1]=feq[idx1];
    }
  }
}/////FINE Routine INITPOP

void forcingconstructWW(REAL *rho,REAL *temp,REAL *forcex,REAL *forcey){
  REAL ran,average;
  char filename[128];
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx0=j+(NY+2)*i;
      average=0.0;//Tdown+(Tup-Tdown)*((double)j-1+mey*NY)/((double)(ny-1));
      forcex[idx0]=0.0;
      forcey[idx0]=-gravity*rho[idx0]*(temp[idx0]-average);
    }
  }
}

void moveplusforcingconstructWW(pop_type *f,pop_type *g,REAL *rho,REAL *temp,REAL *forcex,REAL *forcey,pop_type *memoryup1,pop_type *memoryup2,pop_type *memorydown1,pop_type *memorydown2){//INIZIO ROUTINE MOVEADJUST

  FILE *fout;
  REAL massapost1,massapost2;
  REAL Tnext,densitynow;
  char filename[128];
  if(mey==0){
    for(i=0;i<NX+2;i++){
      idx2=1+(NY+2)*i;  
      memorydown1[i]=f[idx2];
      memorydown2[i]=g[idx2];
    }
  }

  if(mey==NPY-1){
    for(i=0;i<NX+2;i++){
      idx1=NY+(NY+2)*i;  
      memoryup1[i]=f[idx1];
      memoryup2[i]=g[idx1];
    }
  }

  ///////////////////MOVEMENT ///////////////////////

  for(i=NX;i>=1;i--){
    for(j=NY;j>=1;j--){
      idx1=j+(NY+2)*i; 
      idx2=j+(NY+2)*(i-1);
      f[idx1].p[1]=f[idx2].p[1];
      g[idx1].p[1]=g[idx2].p[1];
    }
  }
   

  for(i=1;i<NX+1;i++){
    for(j=NY;j>=1;j--){
      idx1=j+(NY+2)*i; 
      idx2=j-1+(NY+2)*i; 
      f[idx1].p[2]=f[idx2].p[2];
      g[idx1].p[2]=g[idx2].p[2];
    }
  }

  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i; 
      idx2=j+(NY+2)*(i+1); 
      f[idx1].p[3]=f[idx2].p[3];
      g[idx1].p[3]=g[idx2].p[3];
    }
  }

  for(i=NX;i>=1;i--){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i; 
      idx2=j+1+(NY+2)*i; 
      f[idx1].p[4]=f[idx2].p[4];
      g[idx1].p[4]=g[idx2].p[4];
    }
  }
   

  for(i=NX;i>=1;i--){
    for(j=NY;j>=1;j--){
      idx1=j+(NY+2)*i; 
      idx3=j-1+(NY+2)*(i-1); 
      f[idx1].p[5]=f[idx3].p[5];
      g[idx1].p[5]=g[idx3].p[5];
    }
  }
   

  for(i=1;i<NX+1;i++){
    for(j=NY;j>=1;j--){
      idx1=j+(NY+2)*i;  
      idx3=j-1+(NY+2)*(i+1); 
      f[idx1].p[6]=f[idx3].p[6];
      g[idx1].p[6]=g[idx3].p[6];
    }
  }
  
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i; 
      idx3=j+1+(NY+2)*(i+1); 
      f[idx1].p[7]=f[idx3].p[7];
      g[idx1].p[7]=g[idx3].p[7];
    }
  }

  for(i=NX;i>=1;i--){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;  
      idx3=j+1+(NY+2)*(i-1); 
      f[idx1].p[8]=f[idx3].p[8];
      g[idx1].p[8]=g[idx3].p[8];
    }
  }
  ///////////////////////////////////////////////////////////////
  forcingconstructWW(rho,temp,forcex,forcey);  
  ////////////////////////////////////////////////////////////////

  if(mey==NPY-1){

    for(i=1;i<NX+1;i++){ //North wall (Thermalization)

      idx1=NY+(NY+2)*i; 

      idx2=NY-1+(NY+2)*i; 

      massapost1=f[idx1].p[0]+(memoryup1[i].p[2]+memoryup1[i].p[5]+memoryup1[i].p[6])+f[idx1].p[1]+f[idx1].p[3]+f[idx1].p[5]+f[idx1].p[6]+f[idx1].p[2];
 
      f[idx1].p[4]=memoryup1[i].p[2];
  
      f[idx1].p[7]=memoryup1[i].p[5];

      f[idx1].p[8]=memoryup1[i].p[6];

    }
  }

  if(mey==0){

    for(i=1;i<NX+1;i++){ //South wall (Thermalization)

      idx1=1+(NY+2)*i; 

      idx2=2+(NY+2)*i; 

      massapost1=f[idx1].p[0]+(memorydown1[i].p[4]+memorydown1[i].p[8]+memorydown1[i].p[7])+f[idx1].p[1]+f[idx1].p[4]+f[idx1].p[8]+f[idx1].p[7]+f[idx1].p[3];

      f[idx1].p[2]=memorydown1[i].p[4];

      f[idx1].p[5]=memorydown1[i].p[7];

      f[idx1].p[6]=memorydown1[i].p[8];
    }
  }

  if(mey==NPY-1){
    for(i=1;i<NX+1;i++){ 
      idx1=NY+(NY+2)*i;

      g[idx1].p[4]=-g[idx1].p[2] + 2.0*ww[2]*Tup;

      g[idx1].p[7]=-g[idx1].p[5] + 2.0*ww[5]*Tup;

      g[idx1].p[8]=-g[idx1].p[6] + 2.0*ww[6]*Tup;

    }
  }

  if(mey==0){

    for(i=1;i<NX+1;i++){ 
      idx1=1+(NY+2)*i; 

      g[idx1].p[2]=-g[idx1].p[4] + 2.0*ww[4]*Tdown;

      g[idx1].p[5]=-g[idx1].p[7] + 2.0*ww[7]*Tdown;

      g[idx1].p[6]=-g[idx1].p[8] + 2.0*ww[8]*Tdown;

    }
  }
}


void rhocomp(REAL *rho,pop_type *f) {///INIZIO ROUTINE HYDROVAR
  for(i=0;i<NX+2;i++){
    for(j=0;j<NY+2;j++){
      idx1=j+(NY)*i;
      rho[idx1]=0.0;
      for(k=0;k<npop;k++){rho[idx1]=rho[idx1]+f[idx1].p[k];
      }
    }
  }
}

void AVERAGE(REAL *rhopre,REAL *rhopost,REAL *upre,REAL *upost,REAL *vpre,REAL *vpost,REAL *utot,REAL *vtot) {

  REAL mutotpre,mutotpost,muav;
  for(i=0;i<NX+2;i++){
    for(j=0;j<NY+2;j++){

      idx0=j+(NY+2)*i;
      mutotpre=(rhopre[idx0]*upre[idx0]);
      mutotpost=(rhopost[idx0]*upost[idx0]);
      muav=0.5*(mutotpre+mutotpost);
      utot[idx0]=muav/(rhopre[idx0]);
      mutotpre=(rhopre[idx0]*vpre[idx0]);
      mutotpost=(rhopost[idx0]*vpost[idx0]);
      muav=0.5*(mutotpre+mutotpost);
      vtot[idx0]=muav/(rhopre[idx0]);

    }
  }
}


void hydrovar(REAL *u,REAL *v,REAL *rho,REAL *temp,pop_type *f,pop_type *g,REAL *forcex,REAL *forcey) {///INIZIO ROUTINE HYDROVAR

  REAL rhoi;
  REAL fl[npop];
  char filename[128];

  for(i=0;i<NX+2;i++){
    for(j=0;j<NY+2;j++){
      idx1=j+(NY+2)*i;

      rho[idx1]=0.0;
      temp[idx1]=0.0;
      for(k=0;k<npop;k++){
        rho[idx1]=rho[idx1]+f[idx1].p[k];
      }
      for(k=0;k<npop;k++){
        temp[idx1]=temp[idx1]+g[idx1].p[k];
      }

      rhoi=1./rho[idx1];

      for(k=0;k<npop;k++){
        fl[k]=f[idx1].p[k]*rhoi;
      }

      u[idx1]=fl[1]-fl[3]+fl[5]-fl[6]-fl[7]+fl[8] + 0.5*rhoi*forcex[idx1];
      v[idx1]=fl[5]+fl[2]+fl[6]-fl[7]-fl[4]-fl[8] + 0.5*rhoi*forcey[idx1];

    }
  }
}//FINE ROUTINE HYDROVAR


void profile(REAL *rhopre,REAL *temppre,REAL *utot,REAL *vtot) {
  double *upr,*vpr,*usqpr,*vsqpr,*rhopr,*temppr,*vTpr,*dzTpr,*u_prof,*v_prof,*usq_prof,*vsq_prof,*rho_prof,*temp_prof,*vT_prof,*dzT_prof;
  FILE *fout1,*fout2,*fout3,*fout4;
  double norma,norma2,dvbydx,dubydy;
  char filename[128];
  int idy1,idy2,idy0,idxup,idxdn,idxlt,idxrt;

  upr = (double *) malloc((ny+2)*sizeof(double));
  vpr = (double *) malloc((ny+2)*sizeof(double));
  usqpr = (double *) malloc((ny+2)*sizeof(double));
  vsqpr = (double *) malloc((ny+2)*sizeof(double));
  rhopr = (double *) malloc((ny+2)*sizeof(double));
  temppr = (double *) malloc((ny+2)*sizeof(double));
  vTpr = (double *) malloc((ny+2)*sizeof(double));
  dzTpr = (double *) malloc((ny+2)*sizeof(double));

  u_prof = (double *) malloc((ny+2)*sizeof(double));
  v_prof = (double *) malloc((ny+2)*sizeof(double));
  usq_prof = (double *) malloc((ny+2)*sizeof(double));
  vsq_prof = (double *) malloc((ny+2)*sizeof(double));
  rho_prof = (double *) malloc((ny+2)*sizeof(double));
  temp_prof = (double *) malloc((ny+2)*sizeof(double));
  vT_prof = (double *) malloc((ny+2)*sizeof(double));
  dzT_prof = (double *) malloc((ny+2)*sizeof(double));

  for(j=1;j<=ny;j++){
    upr[j]=0;
    vpr[j]=0;
    usqpr[j]=0;
    vsqpr[j]=0;
    rhopr[j]=0;
    temppr[j]=0;
    vTpr[j]=0;
    dzTpr[j]=0;
  }

  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      idx2=j+1+(NY+2)*i;
      indy=j+mey*NY;
      dzTpr[indy]+=(temppre[idx2]-temppre[idx1]);
    }
  }  

  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      idx2=j+1+(NY+2)*i;
      indx=(i)+mex*NX;
      indy=j+mey*NY;
      upr[indy]+=(utot[idx1]);
      vpr[indy]+=(vtot[idx1]);
      usqpr[indy]+=(utot[idx1]*utot[idx1]);
      vsqpr[indy]+=(vtot[idx1]*vtot[idx1]);
      temppr[indy]+=temppre[idx1];
      rhopr[indy]+=rhopre[idx1];
      vTpr[indy]+=(vtot[idx1]*temppre[idx1]);
    }
  }
  //////////


  MPI_Allreduce(upr,u_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(vpr,v_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(usqpr,usq_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(vsqpr,vsq_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(temppr,temp_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(rhopr,rho_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(vTpr,vT_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(dzTpr,dzT_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
 

  if(!me){
    norma=1/((double)(nx));
    sprintf(filename,"profiles.%d.dat",istep);
    fout1=fopen(filename,"w");

    for(j=1;j<ny+1;j++){
      fprintf(fout1,"%d %e %e %e %e %e %e %e %e \n",j,rho_prof[j]*norma,u_prof[j]*norma,v_prof[j]*norma,usq_prof[j]*norma,vsq_prof[j]*norma,temp_prof[j]*norma,vT_prof[j]*norma,dzT_prof[j]*norma);
    }
    fclose(fout1);

  }

  icountprof=icountprof+1;

  free(upr);
  free(vpr);
  free(usqpr);
  free(vsqpr);
  free(rhopr);
  free(temppr);
  free(vTpr);
  free(dzTpr);

  free(u_prof);
  free(v_prof);
  free(usq_prof);
  free(vsq_prof);
  free(rho_prof);
  free(temp_prof);
  free(vT_prof);
  free(dzT_prof);
}


void equilif(REAL *u,REAL *v,REAL *rho,pop_type *feq) {//INIZIO ROUTINE EQUILI

  FILE *fout;
  REAL r1,usq,vsq,sumsq,sumsq2,u2,v2,ui,vi,uv,d0;
  REAL LAPrhoA,LAPrhoB,gradxA,gradxB;
  REAL pres;
  char filename[128];
  
  for(i=0;i<NX+2;i++){
    for(j=0;j<NY+2;j++){

      idx1=j+(NY+2)*i;

      r1=0.0;
      r1=rho[idx1];

      usq=u[idx1]*u[idx1];
      vsq=v[idx1]*v[idx1];

      sumsq=(usq+vsq)/cs22;
      sumsq2=sumsq*(1.-cs2)/cs2;

      u2=usq/cssq;
      v2=vsq/cssq;

      ui=u[idx1]/cs2;
      vi=v[idx1]/cs2;

      uv=ui*vi;

      feq[idx1].p[0]=r1*ww[0]*(1.-sumsq);
      feq[idx1].p[1]=r1*ww[1]*(1.-sumsq+u2+ui);
      feq[idx1].p[2]=r1*ww[2]*(1.-sumsq+v2+vi);
      feq[idx1].p[3]=r1*ww[3]*(1.-sumsq+u2-ui);
      feq[idx1].p[4]=r1*ww[4]*(1.-sumsq+v2-vi);
      feq[idx1].p[5]=r1*ww[5]*(1.+sumsq2+ui+vi+uv);
      feq[idx1].p[6]=r1*ww[6]*(1.+sumsq2-ui+vi-uv);
      feq[idx1].p[7]=r1*ww[7]*(1.+sumsq2-ui-vi+uv);
      feq[idx1].p[8]=r1*ww[8]*(1.+sumsq2+ui-vi-uv);
    }
  }

}//FINE ROUTINE EQUILI

void equilig(REAL *u,REAL *v,REAL *temp,pop_type *feq) {//INIZIO ROUTINE EQUILI

  FILE *fout;
  REAL r1,usq,vsq,sumsq,sumsq2,u2,v2,ui,vi,uv,d0;
  REAL LAPrhoA,LAPrhoB,gradxA,gradxB;
  REAL pres;
  char filename[128];

  for(i=0;i<NX+2;i++){
    for(j=0;j<NY+2;j++){
      idx1=j+(NY+2)*i;

      r1=temp[idx1];

      usq=u[idx1]*u[idx1];
      vsq=v[idx1]*v[idx1];

      sumsq=(usq+vsq)/cs22;
      sumsq2=sumsq*(1.-cs2)/cs2;

      u2=usq/cssq;
      v2=vsq/cssq;

      ui=u[idx1]/cs2;
      vi=v[idx1]/cs2;

      uv=ui*vi;

      feq[idx1].p[0]=r1*ww[0]*(1.-sumsq);
      feq[idx1].p[1]=r1*ww[1]*(1.-sumsq+u2+ui);
      feq[idx1].p[2]=r1*ww[2]*(1.-sumsq+v2+vi);
      feq[idx1].p[3]=r1*ww[3]*(1.-sumsq+u2-ui);
      feq[idx1].p[4]=r1*ww[4]*(1.-sumsq+v2-vi);
      feq[idx1].p[5]=r1*ww[5]*(1.+sumsq2+ui+vi+uv);
      feq[idx1].p[6]=r1*ww[6]*(1.+sumsq2-ui+vi-uv);
      feq[idx1].p[7]=r1*ww[7]*(1.+sumsq2-ui-vi+uv);
      feq[idx1].p[8]=r1*ww[8]*(1.+sumsq2+ui-vi-uv);
    }
  }
}//FINE ROUTINE EQUILIG


void snapshottemp(REAL *rhopre,REAL *temppre,REAL *utot, REAL *vtot) {
  char fname[128];
  char gname[128];
  herr_t status;

  H5Eset_current_stack(H5E_DEFAULT);
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hid_t hdf5_status = H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
  hsize_t edimens_3d[2] = {NX+2,NY+2};
  hid_t ememspace = H5Screate_simple(2, edimens_3d, NULL );
  edimens_3d[0] = nx;
  edimens_3d[1] = ny;
  hid_t efilespace = H5Screate_simple( 2, edimens_3d, NULL );
  hsize_t estart_3d[2] = {1,1};
  hsize_t estride_3d[2] = {1,1};
  hsize_t eblock_3d[2] = {NX,NY};
  hsize_t ecount_3d[2] = {1,1};
  hsize_t dstart_3d[2] = {mex * NX,mey * NY};
  hsize_t dstride_3d[2] = {1,1};
  hsize_t dblock_3d[2] = {NX,NY};
  hsize_t dcount_3d[2] = {1,1};
  hid_t dataset_id,xfer_plist,ret;
  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
  sprintf(fname,"Tempfolder/snaptemp%d.h5",istep);
  hid_t  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id<0){
    printf("snapshottemp not found \n");
    exit(8);
  }

  /* dataset_id = H5Dcreate(file_id, "/vel", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); */
  /* xfer_plist = H5Pcreate (H5P_DATASET_XFER); */
  /* ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE); */
  /* status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, velsq); */
  /* status = H5Pclose( xfer_plist ); */
  /* status = H5Dclose(dataset_id); */

  /* dataset_id = H5Dcreate(file_id, "/vort", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); */
  /* xfer_plist = H5Pcreate (H5P_DATASET_XFER); */
  /* ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE); */
  /* status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, vortsq); */
  /* status = H5Pclose( xfer_plist ); */
  /* status = H5Dclose(dataset_id); */

  dataset_id = H5Dcreate(file_id, "/T", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, temppre);
   
  status = H5Pclose( xfer_plist );
  status = H5Dclose(dataset_id);

  dataset_id = H5Dcreate(file_id, "/u", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, utot);
   
  status = H5Pclose( xfer_plist );
  status = H5Dclose(dataset_id);


  dataset_id = H5Dcreate(file_id, "/v", H5T_NATIVE_DOUBLE, efilespace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, vtot);
   
  status = H5Pclose( xfer_plist );
  status = H5Dclose(dataset_id);


  H5Sclose (ememspace);
  H5Sclose (efilespace);


  status = H5Fclose(file_id);
  icountconfig=icountconfig+1;
  //  free(vortsq);
  // free(velsq);
}


void snapshotold(REAL *rhopre,REAL *temppre,REAL *utot,REAL *vtot) {

  int plusx1,plusy1;
  int minusx1,minusy1;

  int plusx2,plusy2;
  int minusx2,minusy2;


  REAL flusso,muav,mutotpre,mutotpost,mforcing,A,B,g1,g2,PPP,INTERA,INTERB;
  REAL LAPrhoA,LAPrhoB,gradxA,gradxB;

  REAL s1,s2,aaa,bbb,uav,div,div2,term;


  FILE *fout1;
  char filename[128];

  sprintf(filename,"Snapfolder/all_xynew.%d.%d.%d.dat",istep,mex,mey);
  fout1=fopen(filename,"w");

  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx0=j+(NY+2)*i;
      indx=(i)+mex*NX;
      indy=j+mey*NY;
      fprintf(fout1,"%d %d %e %e %e %e \n",indx,indy,rhopre[idx0],utot[idx0],vtot[idx0],temppre[idx0]);
    }
    fprintf(fout1,"\n");
  }
  fclose(fout1);

  icountconfig=icountconfig+1;

}

void computetau(REAL *rho,REAL *temp,REAL *tau,REAL *taug){//INIZIO ROUTINE COLLIS
  REAL omeganow;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      tau[idx1]=relax; 
      taug[idx1]=relaxg; 
    }
  }

}//FINE ROUTINE COMPUTETAU


void collisf(pop_type *f,pop_type *feq,REAL *tau,REAL *u,REAL *v,REAL *forcex,REAL *forcey){//INIZIO ROUTINE COLLIS
  REAL omeganow,source,factor;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      omeganow=1./tau[idx1];
      for(k=0;k<npop;k++){
	source = (1.0- 0.5*omeganow)*ww[k];
	factor = ((cx[k]-u[idx1])*forcex[idx1] + (cy[k]-v[idx1])*forcey[idx1])/cs2;
	factor = factor + ((cx[k]*u[idx1] + cy[k]*v[idx1])*(cx[k]*forcex[idx1]+cy[k]*forcey[idx1]))/(cs2*cs2);
        f[idx1].p[k]=f[idx1].p[k]*(1.-omeganow)+omeganow*feq[idx1].p[k] + source*factor;
      }
    }
  }
}//FINE ROUTINE COLLIS

void collisg(pop_type *f,pop_type *feq,REAL *tau){//INIZIO ROUTINE COLLIS
  REAL omeganow;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      omeganow=1./tau[idx1];
      for(k=0;k<npop;k++){
        f[idx1].p[k]=f[idx1].p[k]*(1.-omeganow)+omeganow*feq[idx1].p[k];
      }
    }
  }
}//FINE ROUTINE COLLIS

void computeUeq(REAL *forcex,REAL *forcey,REAL *upre,REAL *vpre,REAL *rhopre,REAL *eqfieldxf,REAL *eqfieldyf,REAL *eqfieldxg,REAL *eqfieldyg,REAL *tau){
  REAL momx, momy, invrho;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){

      idx1=j+(NY+2)*i;

      eqfieldxf[idx1]=upre[idx1]+tau[idx1]*forcex[idx1]/(rhopre[idx1]);

      eqfieldyf[idx1]=vpre[idx1]+tau[idx1]*forcey[idx1]/(rhopre[idx1]);

      eqfieldxg[idx1]=upre[idx1]+0.5*forcex[idx1]/(rhopre[idx1]);

      eqfieldyg[idx1]=vpre[idx1]+0.5*forcey[idx1]/(rhopre[idx1]);
    }
  }

}//COMPUTEUeq


void dellar(){
  /////////////////////WEIGHTS////////////////

  ww[0]=16./36.;
  ww[1]=4./36.;
  ww[2]=4./36.;
  ww[3]=4./36.;
  ww[4]=4./36.;
  ww[5]=1./36.;
  ww[6]=1./36.;
  ww[7]=1./36.;
  ww[8]=1./36.;
  /////////////////////LATTICE SPEEDS/////////////////////
  cx[0]=0.;
  cy[0]=0.;

  cx[1]=1.;
  cy[1]=0.;

  cx[2]=0.;
  cy[2]=1.;

  for(i=3;i<=4;i++){
    cx[i]=-cx[i-2];
    cy[i]=-cy[i-2];
  }


  cx[5]=1.;
  cy[5]=1.;

  cx[6]=-1.0;
  cy[6]=1.0;

  for(i=7;i<=8;i++){
    cx[i]=-cx[i-2];
    cy[i]=-cy[i-2];
  } 
}//FINE ROUTINE DELLAR

void prtcl_dump(prtcl_type *prtcl1){
  int Nprtcl_test;
  int i,j,send_add,*send_addarr,*send_addarrt; //send_add is number of particles stored in process with rank < me so that we know where to add the prtcl information in the global array
  char filename[128];
  int *nprtcl_array,*nprtcl_arraytemp;
  FILE *fout1;
  REAL xx,yy,uu,vv;
  send_add = 0;
  send_addarr = (int *) malloc((NP)*sizeof(int));
  send_addarrt = (int *) malloc((NP)*sizeof(int));
  nprtcl_array = (int *) malloc((NP)*sizeof(int));
  nprtcl_arraytemp = (int *) malloc((NP)*sizeof(int));
  
  prtcl_type *prtcl_global,*prtcl2;
  prtcl_global = (prtcl_type *) malloc((Nprtcl)*sizeof(prtcl_type));
  prtcl2 = (prtcl_type *) malloc((nprtcl_local)*sizeof(prtcl_type));
  for (i=0;i<nprtcl_local;i++){
    prtcl2[i].xprtcl = ((REAL) (mex*NX)) + prtcl1[i].xprtcl;
    prtcl2[i].yprtcl = ((REAL) (mey*NY)) + prtcl1[i].yprtcl;
    prtcl2[i].uprtcl = prtcl1[i].uprtcl;
    prtcl2[i].vprtcl = prtcl1[i].vprtcl;
  }
   
  for (i=0;i<NP;i++){
    nprtcl_arraytemp[i] = 0;
    send_addarrt[i] = 0;
  }
  nprtcl_arraytemp[me] = nprtcl_local;

  MPI_Reduce(&nprtcl_local,&Nprtcl_test,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  
  MPI_Allreduce(nprtcl_arraytemp,nprtcl_array,NP,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  for (i=0;i<me;i++){
    send_add+=nprtcl_array[i];
  }
  send_addarrt[me]=send_add;
  MPI_Allreduce(send_addarrt,send_addarr,NP,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Gatherv(prtcl2,nprtcl_local,MPI_Prtcltype,prtcl_global,nprtcl_array,send_addarr,MPI_Prtcltype,0,MPI_COMM_WORLD);

  if (me==0){
    sprintf(filename,"prtcl_dump%d.in",istep,mex,mey);
    fout1=fopen(filename,"w");
    //  printf("Entering particle loop %d \n",nprtcl_local);

    fwrite(prtcl_global,sizeof(prtcl_type),Nprtcl,fout1);
    fclose(fout1);
  }
  free(prtcl_global);
  free(nprtcl_array);
  free(nprtcl_arraytemp);
  free(prtcl2);  
}

void prtcl_write_array(hid_t file_id, const char *datasetname, hid_t datatype, void *array, 
		       hid_t status, hid_t ememspace, hid_t efilespace, hid_t xfer_plist, hid_t dataset_id, hid_t ret) {
  // Create dataset
  dataset_id = H5Dcreate(file_id, datasetname, datatype, efilespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
  status = H5Dwrite(dataset_id, datatype, ememspace, efilespace, xfer_plist, array);
  status = H5Pclose(xfer_plist);
  status = H5Dclose(dataset_id);
}

void prtcl_write(prtcl_type *prtcl1,REAL *temp,int *prtcl_ids1){
  int Nprtcl_test;
  int i,j,send_add; //send_add is number of particles stored in process with rank < me so that we know where to add the prtcl information in the global array
  char fname[128];
  int nprtcl_array[NP],nprtcl_arraytemp[NP];
  FILE *fout11;
  REAL xx,yy,uu,vv,temp_prtcl[nprtcl_local],temp_prtclnow;
  REAL xxo,yyo,uuo,vvo;
  //  send_addarr = (int *) malloc((NP)*sizeof(int));
  //  send_addarrt = (int *) malloc((NP)*sizeof(int));
  //  nprtcl_array = (int *) malloc((NP)*sizeof(int));
  //nprtcl_arraytemp = (int *) malloc((NP)*sizeof(int));
  
  REAL prtcl2_x[nprtcl_local],prtcl2_y[nprtcl_local];
  REAL prtcl2_u[nprtcl_local],prtcl2_v[nprtcl_local];
  REAL prtcl2_xo[nprtcl_local],prtcl2_yo[nprtcl_local];
  REAL prtcl2_uo[nprtcl_local],prtcl2_vo[nprtcl_local];
  int prtcl_ids2[nprtcl_local];
  
  /* prtcl2_x = (REAL *) malloc((nprtcl_local)*sizeof(REAL)); */
  /* prtcl2_y = (REAL *) malloc((nprtcl_local)*sizeof(REAL)); */
  /* temp_prtcl = (REAL *) malloc((nprtcl_local)*sizeof(REAL)); */
  /* prtcl_ids2 = (int *) malloc(nprtcl_local*sizeof(int)); */
  communication(temp);
  for (i=0;i<nprtcl_local;i++){
    xx = prtcl1[i].xprtcl;
    yy = prtcl1[i].yprtcl;
    prtcl2_x[i] = ((REAL) (mex*NX)) + xx;
    prtcl2_y[i] = ((REAL) (mey*NY)) + yy;
    bilinear(xx,yy,temp,&temp_prtclnow);
    temp_prtcl[i] = temp_prtclnow;
    prtcl2_u[i] = prtcl1[i].uprtcl;
    prtcl2_v[i] = prtcl1[i].vprtcl;
    prtcl_ids2[i] = prtcl_ids1[i];

    xxo = prtcl1[i].xoprtcl;
    yyo = prtcl1[i].yoprtcl;
    prtcl2_xo[i] = ((REAL) (mex*NX)) + xxo;
    prtcl2_yo[i] = ((REAL) (mey*NY)) + yyo;
    prtcl2_uo[i] = prtcl1[i].uoprtcl;
    prtcl2_vo[i] = prtcl1[i].voprtcl;

  }

  int offs_arrt[NP],offs_arr[NP];
  for (i=0;i<NP;i++){
    nprtcl_arraytemp[i] = 0;
    offs_arrt[i] = 0;
  }
  nprtcl_arraytemp[me] = nprtcl_local;

  MPI_Allreduce(nprtcl_arraytemp,nprtcl_array,NP,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  send_add = 0;
  for (i=0;i<me;i++){
    send_add+=nprtcl_array[i];
  }

  H5Eset_current_stack(H5E_DEFAULT);
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hid_t hdf5_status = H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
  hsize_t edimens_3d[1] = {nprtcl_local};
  hsize_t fdimens_3d[1] = {Nprtcl};
  hid_t ememspace = H5Screate_simple(1, edimens_3d, NULL );
  hid_t efilespace = H5Screate_simple( 1, fdimens_3d, NULL );

  hsize_t estart_3d[1] = {0};
  hsize_t estride_3d[1] = {1};
  hsize_t eblock_3d[1] = {nprtcl_local};
  hsize_t ecount_3d[1] = {1};

  hsize_t dstart_3d[1] = {send_add};//offset
  hsize_t dstride_3d[1] = {1};
  hsize_t dblock_3d[1] = {nprtcl_local};//local
  hsize_t dcount_3d[1] = {1};
  hid_t dataset_id,xfer_plist,ret;
  hid_t status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
  status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);

  sprintf(fname,"Prtcl_Veldata1/Particles%d.h5",istep);
  if (istep>10000000) sprintf(fname,"Prtcl_Veldata2/Particles%d.h5",istep);
  hid_t  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id<0){
    printf("Particle file not created \n");
    exit(8);
  }

  prtcl_write_array(file_id, "/id", H5T_NATIVE_INT, prtcl_ids2, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/x", H5T_NATIVE_DOUBLE, prtcl2_x, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/y", H5T_NATIVE_DOUBLE, prtcl2_y, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/u", H5T_NATIVE_DOUBLE, prtcl2_u, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/v", H5T_NATIVE_DOUBLE, prtcl2_v, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/xo", H5T_NATIVE_DOUBLE, prtcl2_xo, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/yo", H5T_NATIVE_DOUBLE, prtcl2_yo, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/uo", H5T_NATIVE_DOUBLE, prtcl2_uo, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/vo", H5T_NATIVE_DOUBLE, prtcl2_vo, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  prtcl_write_array(file_id, "/temp", H5T_NATIVE_DOUBLE, temp_prtcl, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
  
  
  status = H5Fclose(file_id);

  H5Sclose (efilespace);
  H5Sclose (ememspace);
  
  status = H5Pclose(plist_id);
}

int main(int argc, char *argv[]) {
  
  FILE *fout4,*fout,*fout6,*fout7,*fin,*fout1,*fout2,*fout3,*fout5,*fout8,*fout9;

  REAL utotpre,utotpost,forcing;
  pop_type *f,*feq,*g,*geq,*memoryup1,*memoryup2,*memorydown1,*memorydown2;
  REAL *upre,*vpre,*upost,*vpost,*rhopre,*temppre,*rhopost,*temppost,*forcex,*forcey,*utot,*vtot,*tempread;
  REAL *uread,*vread;
  REAL *eqfieldxf,*eqfieldxg,*eqfieldyf,*eqfieldyg;
  REAL *tau,*taug;

  REAL r1,r2,r3,A,B,g1,g2,totalenergy,totalenergy1,totaltemp,totaltemp1;;

  int i;


  myInit(argc,argv);

#if Particle
  prtcl_type *prtcl;
  int *prtcl_ids;
  REAL *ubefore, *vbefore, *dtu, *dtv, *dxu, *dxv, *dyu, *dyv;
#endif


  f = (pop_type *) malloc((NX+2)*(NY+2)*sizeof(pop_type));
  feq = (pop_type *) malloc((NX+2)*(NY+2)*sizeof(pop_type));
  g = (pop_type *) malloc((NX+2)*(NY+2)*sizeof(pop_type));
  geq = (pop_type *) malloc((NX+2)*(NY+2)*sizeof(pop_type));
  memoryup1 = (pop_type *) malloc((NX+2)*sizeof(pop_type));
  memoryup2 = (pop_type *) malloc((NX+2)*sizeof(pop_type));
  memorydown1 = (pop_type *) malloc((NX+2)*sizeof(pop_type));
  memorydown2 = (pop_type *) malloc((NX+2)*sizeof(pop_type));
  upre = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));   //upre means before the u velocity before the collision step -- It's a part of Lattice Boltzmann evolution
  vpre = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  rhopre = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  upost = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  vpost = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  rhopost = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  temppre = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  temppost = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  forcex = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  forcey = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  eqfieldxf = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  eqfieldyf = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  eqfieldxg = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  eqfieldyg = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  utot = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  vtot = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  tau = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  taug = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  tempread = (REAL *) malloc(NX * NY * sizeof(REAL));
  uread = (REAL *) malloc(NX * NY * sizeof(REAL));
  vread = (REAL *) malloc(NX * NY * sizeof(REAL));
  dellar();

#if Particle
    ubefore = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));   //ubefore means the upre from the previous time-step. Storing this because we need it for Maxey-Riley Equation
    vbefore = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
    dtu = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));   //ubefore means the upre from the previous time-step. Storing this because we need it for Maxey-Riley Equation
    dtv = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
    dxu = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));   //ubefore means the upre from the previous time-step. Storing this because we need it for Maxey-Riley Equation
    dxv = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
    dyu = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));   //ubefore means the upre from the previous time-step. Storing this because we need it for Maxey-Riley Equation
    dyv = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
    
    prtcl = (prtcl_type *) malloc(nprtcl_max*sizeof(prtcl_type));
    prtcl_ids = (int *) malloc(nprtcl_max*sizeof(int));
#endif

  if(scratch) {
    inithydro(upre,vpre,rhopre,temppre,f,g);
    equilif(upre,vpre,rhopre,feq);
    initpop(f,feq);
    equilig(upre,vpre,temppre,geq);
    initpop(g,geq);
  } else {
    restoref(f);
    restoreg(g);
    restorefpre(feq);
    restoregpre(geq);
    inithydro(upre,vpre,rhopre,temppre,feq,geq);
  }
  //////////////////////////MAIN LOOP//////////////////////////
#if Particle

#if Particle_scratch
  init_prtcl(prtcl,upre,vpre,prtcl_ids,temppre);
  //printf("Particle Initialised \n");
#endif

  if (Particle_scratch==0) restore_prtcl(prtcl,prtcl_ids);
#endif

  icount=0;
  icountconfig=0;
  icountprof=0;
  if (me==0) {
    fout6=fopen("evolution.dat","w");
  }

  for(istep=nstart+1;istep<=nsteps+nstart;istep++){ //nsteps
    pbc(f);
    pbc(g);
    moveplusforcingconstructWW(f,g,rhopre,temppre,forcex,forcey,memoryup1,memoryup2,memorydown1,memorydown2);
    ////////////////////////////////////////
    /////and the move-machinery is over/////
    ////////////////////////////////////////
    pbc(f);
    pbc(g);
#if Particle
    mirror(upre,ubefore);
    mirror(vpre,vbefore);
#endif
    hydrovar(upre,vpre,rhopre,temppre,f,g,forcex,forcey);

    /* if (me==0){ */
    /*   printf("Loop %d hydro stuff done \n",istep); */
    /* } */

#if Particle
    propagate_prtcl(prtcl,upre,vpre,ubefore,vbefore,dtu,dtv,dxu,dxv,dyu,dyv,temppre,prtcl_ids);
    toggle_prtcl(&prtcl,&prtcl_ids);
#endif
    if (me==0){
      printf("Loop %d particles propagated \n",istep);
    }


    if (istep==(nstart+nsteps)){
      dumpfpre(f);
      dumpgpre(g);

    }
    if (istep%2000000==(0)){
      dumpfpre(f);
      dumpgpre(g);
    }

    computetau(rhopre,temppre,tau,taug);

    //    computeUeq(forcex,forcey,upre,vpre,rhopre,eqfieldxf,eqfieldyf,eqfieldxg,eqfieldyg,tau);
    equilif(upre,vpre,rhopre,feq);
    collisf(f,feq,tau,upre,vpre,forcex,forcey);

    equilig(upre,vpre,temppre,geq);
    collisg(g,geq,taug);

    pbc(f);
    pbc(g);
    hydrovar(upost,vpost,rhopost,temppost,f,g,forcex,forcey);

    AVERAGE(rhopre,rhopost,upre,upost,vpre,vpost,utot,vtot);
    if (istep%2000000==(0)){
      dumpf(f);
      dumpg(g);
#if Particle
      prtcl_dump(prtcl);
#endif
    }
    
    if((istep%nout)==0){  ////TO GET THE PROFILES
      MPI_Barrier(MPI_COMM_WORLD);
      //profile(rhopre,temppre,utot,vtot);
      //snapshottemp(rhopre,temppre,upre,vpre);

      totalenergy=0.0;
      totalenergy1=0.0;
      totaltemp=0.0;
      totaltemp1=0.0;   
      for(i=1;i<NX+1;i++){for(j=1;j<NY+1;j++){
	  idx1=j+(NY+2)*i;
	  totalenergy+=0.5*rhopre[idx1]*(upre[idx1]*upre[idx1]+vpre[idx1]*vpre[idx1]);
	  totaltemp+=temppre[idx1];
	}}
      MPI_Reduce(&totalenergy,&totalenergy1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&totaltemp,&totaltemp1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      if (me==0) {
	totalenergy1/=(((double)nx)*((double)(ny)));
	totaltemp1/=(((double)nx)*((double)(ny)));
	fprintf(fout6,"%d %e %e\n",istep,totalenergy1,totaltemp1);
	fflush(fout6);

      }

    }
  }//END OF MAIN LOOP
  if (me==0){
    fclose(fout6);
  }
  dumpf(f);
  dumpg(g);
#if Particle
  prtcl_dump(prtcl);
#endif
  myFinalize();
}
