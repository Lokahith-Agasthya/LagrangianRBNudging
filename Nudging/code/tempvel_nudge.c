#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define AMIROOT (!me)
#include "hdf5.h"
#define REAL double
#define npop 9          ////NUMERO DELLE POPOLAZIONI
#define cs2 (1./3.)
#define cs1 sqrt((cs2))
#define cs22 (2.*cs2)
#define cssq (2./9.)
#define pi 3.14159265359
#define alpha 0.01
/////////INPUT VALUES/////////
#define nx         1680
#define ny         840
#define scratch    1         //Set to 1 if initialising from scratch
#define contin     0
#define nstart     0
#define nsteps     3000000
#define nout       20000
#define noutconfig 20000
#define ntempread  250
#define rhom1      1.0
#define relax      0.502
#define relaxg     0.502
#define Tup        -0.01
#define Tdown      0.01
#define gravity    (-0.00004)

#define Particle   1

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
int indx,indy,idxread,indx1;
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
  sprintf(fname,"../dump_Guo2/conf.h5",nstart+1);
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
  sprintf(fname,"../dump_Guo2/cong.h5",nstart+1);
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
	if ((f[idx1].p[k])==0.0){
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
  sprintf(fname,"../dump_Guo2/confpre.h5",nstart);
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
  sprintf(fname,"../dump_Guo2/congpre.h5",nstart);
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
	if ((f[idx1].p[k])==0.0){
	  printf("Damaar g in the initial function only \n");
	  printf("%d %d %d \n",i,j,k);
	  exit(7);
	}
      }
    }
  }
}

void inithydro(REAL *u,REAL *v,REAL *rho,REAL *temp,pop_type *f,pop_type *g,int *nudge_vol,REAL *tempreadint){//INIZIO ROUTINE INITHYDRO
  FILE *fout,*fout1;
  REAL yyytemp,y_c,fl[npop],rhoi;
  int indx,indy;
  char filename[128], filename1[128];
  int c,cont1m,idxread;
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
    if (contin==0){
      for(i=0;i<NX+2;i++){
	for(j=0;j<NY+2;j++){
	  idx1=j+(NY+2)*i;
	  idxread = j+NY*i;
	  rho[idx1]=0.0;
	  temp[idx1]=0.0;
	  u[idx1] = 0.0;
	  v[idx1] = 0.0;
	  rho[idx1] = rhom1;
	  if (nudge_vol[idx1]==0){
	    for (k=0;k<npop;k++){
	      g[idx1].p[k] = 0.0;
	    }
	  }else{
	    rho[idx1]=rhom1;
	    for (k=0;k<npop;k++){
	      temp[idx1]+= g[idx1].p[k];
	    }
	  }
	}
      }
    }else{
      REAL inv_fac;
      for(i=0;i<NX+2;i++){
	for(j=0;j<NY+2;j++){
	  idx1=j+(NY+2)*i;
	  idxread = (j-1)+NY*(i-1);
	  rho[idx1]=0.0;
	  temp[idx1]=0.0;
	  for(k=0;k<npop;k++){
	    rho[idx1]=rho[idx1]+f[idx1].p[k];	    
	    temp[idx1]=temp[idx1]+g[idx1].p[k];
	  }
	  if (nudge_vol[idx1]>0){
	    inv_fac = 1.0/(1.0 + 0.5*alpha);
	    temp[idx1] = (temp[idx1] + 0.5*alpha*tempreadint[idxread])*inv_fac;
	  }

	  for(k=0;k<npop;k++){
	    fl[k]=f[idx1].p[k]*rhoi;
	  }
	  rhoi=1./rho[idx1];
	  
	  u[idx1]=fl[1]-fl[3]+fl[5]-fl[6]-fl[7]+fl[8];
	  v[idx1]=fl[5]+fl[2]+fl[6]-fl[7]-fl[4]-fl[8] + (0.5*(-gravity)*temp[idx1]);
	}
      }
    }
  }
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
      idxread = (j-1)+NY*(i-1);
      average=0.0;//Tdown+(Tup-Tdown)*((double)j-1+mey*NY)/((double)(ny-1));
      forcey[idx0]=-gravity*rho[idx0]*(temp[idx0]-average);
      forcex[idx0]=0.0;
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
      idx1=2+(NY+2)*i;
      memorydown1[i]=f[idx2];
      memorydown2[i]=g[idx2];
      g[idx2] = g[idx1];
    }
  }

  if(mey==NPY-1){
    for(i=0;i<NX+2;i++){
      idx1=NY+(NY+2)*i;  
      idx2=(NY-1)+(NY+2)*i;
      memoryup1[i]=f[idx1];
      memoryup2[i]=g[idx1];
      g[idx1] = g[idx2];
    }
  }

  ///////////////////MOVEMENT of bulk ///////////////////////
  if ((mey>0)&&(mey<NPY-1)){
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
  }

  /////Bottom row zero gradient boundary condition//////
  if (mey==0){
    for(i=NX;i>=1;i--){
      for(j=NY;j>=1;j--){
	idx1=j+(NY+2)*i; 
	idx2=j+(NY+2)*(i-1);
	f[idx1].p[1]=f[idx2].p[1];
	if (j>1){
	g[idx1].p[1]=g[idx2].p[1];
	}
      }
    }
    
    for(i=1;i<NX+1;i++){
      for(j=NY;j>=1;j--){
	idx1=j+(NY+2)*i; 
	idx2=j-1+(NY+2)*i; 
	f[idx1].p[2]=f[idx2].p[2];
	if (j>1){
	g[idx1].p[2]=g[idx2].p[2];
	}
      }
    }

    for(i=1;i<NX+1;i++){
      for(j=1;j<NY+1;j++){
	idx1=j+(NY+2)*i; 
	idx2=j+(NY+2)*(i+1); 
	f[idx1].p[3]=f[idx2].p[3];
	if (j>1){
	g[idx1].p[3]=g[idx2].p[3];
	}
      }
    }

    for(i=NX;i>=1;i--){
      for(j=1;j<NY+1;j++){
	idx1=j+(NY+2)*i; 
	idx2=j+1+(NY+2)*i; 
	f[idx1].p[4]=f[idx2].p[4];
	if (j>1){
	g[idx1].p[4]=g[idx2].p[4];
	}
      }
    }
   
    for(i=NX;i>=1;i--){
      for(j=NY;j>=1;j--){
	idx1=j+(NY+2)*i; 
	idx3=j-1+(NY+2)*(i-1); 
	f[idx1].p[5]=f[idx3].p[5];
	if (j>1){
	g[idx1].p[5]=g[idx3].p[5];
	}
      }
    }
   
    for(i=1;i<NX+1;i++){
      for(j=NY;j>=1;j--){
	idx1=j+(NY+2)*i;  
	idx3=j-1+(NY+2)*(i+1); 
	f[idx1].p[6]=f[idx3].p[6];
	if (j>1){
	g[idx1].p[6]=g[idx3].p[6];
	}
      }
    }
  
    for(i=1;i<NX+1;i++){
      for(j=1;j<NY+1;j++){
	idx1=j+(NY+2)*i; 
	idx3=j+1+(NY+2)*(i+1); 
	f[idx1].p[7]=f[idx3].p[7];
	if (j>1){
	g[idx1].p[7]=g[idx3].p[7];
	}
      }
    }

    for(i=NX;i>=1;i--){
      for(j=1;j<NY+1;j++){
	idx1=j+(NY+2)*i;  
	idx3=j+1+(NY+2)*(i-1); 
	f[idx1].p[8]=f[idx3].p[8];
	if (j>1){
	g[idx1].p[8]=g[idx3].p[8];
	}
      }
    }
  }

  /////Top row zero gradient boundary condition//////
  if (mey==NPY-1){
    for(i=NX;i>=1;i--){
      for(j=NY;j>=1;j--){
	idx1=j+(NY+2)*i; 
	idx2=j+(NY+2)*(i-1);
	f[idx1].p[1]=f[idx2].p[1];
	if (j<NY){
	g[idx1].p[1]=g[idx2].p[1];
	}
      }
    }
    
    for(i=1;i<NX+1;i++){
      for(j=NY;j>=1;j--){
	idx1=j+(NY+2)*i; 
	idx2=j-1+(NY+2)*i; 
	f[idx1].p[2]=f[idx2].p[2];
	if (j<NY){
	g[idx1].p[2]=g[idx2].p[2];
	}
      }
    }

    for(i=1;i<NX+1;i++){
      for(j=1;j<NY+1;j++){
	idx1=j+(NY+2)*i; 
	idx2=j+(NY+2)*(i+1); 
	f[idx1].p[3]=f[idx2].p[3];
	if (j<NY){
	g[idx1].p[3]=g[idx2].p[3];
	}
      }
    }

    for(i=NX;i>=1;i--){
      for(j=1;j<NY+1;j++){
	idx1=j+(NY+2)*i; 
	idx2=j+1+(NY+2)*i; 
	f[idx1].p[4]=f[idx2].p[4];
	if (j<NY){
	g[idx1].p[4]=g[idx2].p[4];
	}
      }
    }
   
    for(i=NX;i>=1;i--){
      for(j=NY;j>=1;j--){
	idx1=j+(NY+2)*i; 
	idx3=j-1+(NY+2)*(i-1); 
	f[idx1].p[5]=f[idx3].p[5];
	if (j<NY){
	g[idx1].p[5]=g[idx3].p[5];
	}
      }
    }
   
    for(i=1;i<NX+1;i++){
      for(j=NY;j>=1;j--){
	idx1=j+(NY+2)*i;  
	idx3=j-1+(NY+2)*(i+1); 
	f[idx1].p[6]=f[idx3].p[6];
	if (j<NY){
	g[idx1].p[6]=g[idx3].p[6];
	}
      }
    }
  
    for(i=1;i<NX+1;i++){
      for(j=1;j<NY+1;j++){
	idx1=j+(NY+2)*i; 
	idx3=j+1+(NY+2)*(i+1); 
	f[idx1].p[7]=f[idx3].p[7];
	if (j<NY){
	g[idx1].p[7]=g[idx3].p[7];
	}
      }
    }

    for(i=NX;i>=1;i--){
      for(j=1;j<NY+1;j++){
	idx1=j+(NY+2)*i;  
	idx3=j+1+(NY+2)*(i-1); 
	f[idx1].p[8]=f[idx3].p[8];
	if (j<NY){
	g[idx1].p[8]=g[idx3].p[8];
	}
      }
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

  /* if(mey==NPY-1){ */
  /*   for(i=1;i<NX+1;i++){  */
  /*     idx1=NY+(NY+2)*i; */

  /*     /\* g[idx1].p[4]=-g[idx1].p[2] + 2.0*ww[2]*Tup; *\/ */

  /*     /\* g[idx1].p[7]=-g[idx1].p[5] + 2.0*ww[5]*Tup; *\/ */

  /*     /\* g[idx1].p[8]=-g[idx1].p[6] + 2.0*ww[6]*Tup; *\/ */
      
  /*   } */
  /* } */

  /* if(mey==0){ */

  /*   for(i=1;i<NX+1;i++){  */
  /*     idx1=1+(NY+2)*i;  */

  /*     g[idx1].p[2]=-g[idx1].p[4] + 2.0*ww[4]*Tdown; */

  /*     g[idx1].p[5]=-g[idx1].p[7] + 2.0*ww[7]*Tdown; */

  /*     g[idx1].p[6]=-g[idx1].p[8] + 2.0*ww[8]*Tdown; */

  /*   } */
  /* } */
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


void hydrovar(REAL *u,REAL *v,REAL *rho,REAL *temp,pop_type *f,pop_type *g,int *nudge_vol,REAL *tempreadint,REAL *forcex,REAL *forcey,double *nud_term) {///INIZIO ROUTINE HYDROVAR

  REAL rhoi;
  REAL fl[npop];
  char filename[128];
  REAL Q;
  REAL lamda = 8.0*pi/((REAL) nx);
  for(i=0;i<NX+2;i++){
    for(j=0;j<NY+2;j++){
      idx1=j+(NY+2)*i;
      idxread = (j-1)+NY*(i-1);
      rho[idx1]=0.0;
      temp[idx1]=0.0;
      for(k=0;k<npop;k++){
        rho[idx1]=rho[idx1]+f[idx1].p[k];
      }
      for(k=0;k<npop;k++){
        temp[idx1]=temp[idx1]+g[idx1].p[k];
      }
      if (nudge_vol[idx1]>0){
	temp[idx1] = (temp[idx1] + 0.5*nud_term[idx1]);
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
  double *vorsqpr,*vorsq_prof;
  double *vorlocal;
  FILE *fout1,*fout2,*fout3,*fout4;
  double norma,norma2,dvbydx,dubydy;
  char filename[128];
  int idy1,idy2,idy0,idxup,idxdn,idxlt,idxrt;

  vorlocal = (double *) malloc((NX+2)*(NY+2)*sizeof(REAL));

  for (i=1;i<NX+1;i++){
    for (j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      idxup = j+1 + (NY+2)*i;
      idxdn = j-1 + (NY+2)*i;
      idxlt = j + (NY+2)*(i-1);
      idxrt = j + (NY+2)*(i+1);
      dvbydx = 0.5*(vtot[idxrt] - vtot[idxlt]);
      dubydy = 0.5*(utot[idxup] - utot[idxdn]);
      vorlocal[idx1] = dvbydx - dubydy;
    }
  }
  
  if (mey==0){
    for (i=1;i<NX+1;i++){
      j = 1;
      idx1=j+(NY+2)*i;
      idxup = j+1 + (NY+2)*i;
      idxlt = j + (NY+2)*(i-1);
      idxrt = j + (NY+2)*(i+1);
      dvbydx = 0.5*(vtot[idxrt] - vtot[idxlt]);
      dubydy = (utot[idxup] - utot[idx1]);
      vorlocal[idx1] = dvbydx - dubydy;
    }
  }

  if (mey==NPY-1){
    for (i=1;i<NX+1;i++){
      j = NY;
      idx1=j+(NY+2)*i;
      idxdn = j-1 + (NY+2)*i;
      idxlt = j + (NY+2)*(i-1);
      idxrt = j + (NY+2)*(i+1);
      dvbydx = 0.5*(vtot[idxrt] - vtot[idxlt]);
      dubydy = (utot[idx1] - utot[idxdn]);
      vorlocal[idx1] = dvbydx - dubydy;
    }
  }

  upr = (double *) malloc((ny+2)*sizeof(double));
  vpr = (double *) malloc((ny+2)*sizeof(double));
  usqpr = (double *) malloc((ny+2)*sizeof(double));
  vsqpr = (double *) malloc((ny+2)*sizeof(double));
  rhopr = (double *) malloc((ny+2)*sizeof(double));
  temppr = (double *) malloc((ny+2)*sizeof(double));
  vTpr = (double *) malloc((ny+2)*sizeof(double));
  dzTpr = (double *) malloc((ny+2)*sizeof(double));
  vorsqpr = (double *) malloc((ny+2)*sizeof(double));

  u_prof = (double *) malloc((ny+2)*sizeof(double));
  v_prof = (double *) malloc((ny+2)*sizeof(double));
  usq_prof = (double *) malloc((ny+2)*sizeof(double));
  vsq_prof = (double *) malloc((ny+2)*sizeof(double));
  rho_prof = (double *) malloc((ny+2)*sizeof(double));
  temp_prof = (double *) malloc((ny+2)*sizeof(double));
  vT_prof = (double *) malloc((ny+2)*sizeof(double));
  dzT_prof = (double *) malloc((ny+2)*sizeof(double));
  vorsq_prof = (double *) malloc((ny+2)*sizeof(double));

  for(j=1;j<=ny;j++){
    upr[j]=0;
    vpr[j]=0;
    usqpr[j]=0;
    vsqpr[j]=0;
    rhopr[j]=0;
    temppr[j]=0;
    vTpr[j]=0;
    dzTpr[j]=0;
    vorsqpr[j]=0;
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
      idxread = (j-1)+NY*(i-1);
      upr[indy]+=(utot[idx1]);
      vpr[indy]+=(vtot[idx1]);
      usqpr[indy]+=(utot[idx1]*utot[idx1]);
      vsqpr[indy]+=(vtot[idx1]*vtot[idx1]);
      temppr[indy]+=temppre[idx1];
      rhopr[indy]+=rhopre[idx1];
      vTpr[indy]+=(vtot[idx1]*temppre[idx1]);
      vorsqpr[indy]+=(vorlocal[idx1]*vorlocal[idx1]);
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
  MPI_Allreduce(vorsqpr,vorsq_prof,ny+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
 

  if(!me){
    norma=1/((double)(nx));
    sprintf(filename,"profiles.%d.dat",istep);
    fout1=fopen(filename,"w");

    for(j=1;j<ny+1;j++){
      fprintf(fout1,"%d %e %e %e %e %e %e %e %e %e\n",j,rho_prof[j]*norma,u_prof[j]*norma,v_prof[j]*norma,usq_prof[j]*norma,vsq_prof[j]*norma,temp_prof[j]*norma,vT_prof[j]*norma,dzT_prof[j]*norma);
    }
    fclose(fout1);

  }

  icountprof=icountprof+1;

  free(vorlocal);
  free(upr);
  free(vpr);
  free(usqpr);
  free(vsqpr);
  free(rhopr);
  free(temppr);
  free(vTpr);
  free(dzTpr);
  free(vorsqpr);

  free(u_prof);
  free(v_prof);
  free(usq_prof);
  free(vsq_prof);
  free(rho_prof);
  free(temp_prof);
  free(vT_prof);
  free(dzT_prof);
  free(vorsq_prof);
}

void media(REAL *rho,REAL *temp) {///INIZIO ROUTINE MEDIA

  REAL nu,kappa,Delta,L,Ra;
   FILE *fout;

  if (me==0){
    fout=fopen("Rayleigh.dat","w");

    nu=cs2*(relax-0.5);
    kappa=cs2*(relaxg-0.5);

    Delta=(Tdown-Tup);
    L=(REAL)(ny-1); 

    Ra=-(gravity/rhom1)*L*L*L*Delta/(nu*kappa);

    fprintf(fout,"%e \n",Ra);

  fclose(fout);
  }

}//FINE ROUTINE MEDIA

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



void read_prtcl(REAL *tempreadint,int *nudge_vol,int *chk_id,int *bc_chkx,int *bc_chky,REAL *diff_x,REAL *prtcl_xint,REAL *diff_y,REAL *prtcl_yint,REAL *diff_T,REAL *prtcl_Tint,REAL *prtcl_x,REAL *prtcl_y,REAL *prtcl_T,int *prtcl_ids){

  int prtcl_xidx[Nprtcl],prtcl_yidx[Nprtcl];
  int x_now,y_now,i1,j1,tmp_xnow,xx,yy;
  REAL xposnow,yposnow,Tnow;
  int procmx,procpx,procmy,procpy;
  procpx=(mex+1+NPX)%NPX; 
  procmx=(mex-1+NPX)%NPX;
  procpy=(mey+1+NPY)%NPY; 
  procmy=(mey-1+NPY)%NPY;
  if ((istep%ntempread)==1){
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

    sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/High_Ra/St0.5/beta1/Prtcl_Veldata1/Particles%d.h5",istep+ntempread-1);
    /* if (istep>2000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata2/Particles%d.h5",istep+ntempread-1); */
    /* if (istep>4000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata3/Particles%d.h5",istep+ntempread-1); */
    /* if (istep>6000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata4/Particles%d.h5",istep+ntempread-1); */
    /* if (istep>8000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata5/Particles%d.h5",istep+ntempread-1); */
    /* if (istep>10000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata6/Particles%d.h5",istep+ntempread-1); */
    /* if (istep>12000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata7/Particles%d.h5",istep+ntempread-1); */
    /* if (istep>14000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata8/Particles%d.h5",istep+ntempread-1); */
    hid_t  file_id = H5Fopen(fname, H5F_ACC_RDWR,H5P_DEFAULT);
    if (file_id<0){
      printf("Particle file not opened %d \n",istep);
      exit(8);
    }
    hid_t ememspace = H5Screate_simple(1, edimens_3d, NULL );
    hid_t efilespace = H5Screate_simple( 1, fdimens_3d, NULL );
    hid_t status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
    status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
    REAL prtcl2_x[Nprtcl];
    REAL prtcl2_y[Nprtcl];
    REAL prtcl2_T[Nprtcl];
    /* for (i=0;i<Nprtcl;i++){ */
    /*   prtcl2_x[i] = 0.; */
    /*   prtcl2_y[i] = 0.; */
    /*   prtcl2_T[i] = 0.;       */
    /* } */
    
    dataset_id = H5Dopen (file_id, "/x", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, prtcl2_x);
    status = H5Pclose( plist_id );
    status = H5Dclose(dataset_id);

    plist_id = H5Pcreate (H5P_FILE_ACCESS);
    hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    dataset_id = H5Dopen (file_id, "/y", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, prtcl2_y);
    status = H5Pclose( plist_id );
    status = H5Dclose(dataset_id);

    int prtcl2_id[Nprtcl];
    plist_id = H5Pcreate (H5P_FILE_ACCESS);
    hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    dataset_id = H5Dopen (file_id, "/id", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, ememspace, efilespace, H5P_DEFAULT, prtcl2_id);
    status = H5Pclose( plist_id );
    status = H5Dclose(dataset_id);

    plist_id = H5Pcreate (H5P_FILE_ACCESS);
    hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    dataset_id = H5Dopen (file_id, "/temp", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, prtcl2_T);
    status = H5Pclose( plist_id );
    status = H5Dclose(dataset_id);

    H5Sclose (efilespace);
    H5Sclose (ememspace);
    
    status = H5Fclose(file_id);

    MPI_Bcast(&prtcl2_x,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&prtcl2_y,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&prtcl2_T,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&prtcl2_id,Nprtcl,MPI_INT,0,MPI_COMM_WORLD);
    
    np_loc = 0;
    //prtcl2_id[chk_id[i]] = prtcl_ids[i] is what we need
    for (i=0;i<Nprtcl;i++){
      for (j=0;j<Nprtcl;j++){
	if (prtcl2_id[j]==prtcl_ids[i]){
	  chk_id[i] = j;
	  break;
	}
      }
    }

    REAL xread,yread,Tread,uread,vread;
    for (i=0;i<Nprtcl;i++){
      bc_chky[i] = 0; //if 1, then it will pass through the y-boundary at some point of time
      //if 2, then it has actually passed in current time step. Need to change it to the read temperature for nudging
      xread = prtcl2_x[chk_id[i]];
      yread = prtcl2_y[chk_id[i]];
      diff_x[i] =(xread - prtcl_xint[i])/((double) ntempread);
      if (abs(xread - prtcl_xint[i]) > ((double) nx/2.0)){
	if (xread > prtcl_xint[i]){
	  xread = xread - (REAL) nx;
	  diff_x[i] =(xread - prtcl_xint[i])/((double) ntempread);
	}else{
	  xread = xread + (REAL) nx;
	  diff_x[i] =(xread - prtcl_xint[i])/((double) ntempread);
	}
      }
      xposnow = prtcl_xint[i] + diff_x[i];
      if (xposnow>(REAL)nx){
	xposnow = xposnow - (REAL)nx;
      }
      if (xposnow < 0.0){
	xposnow = xposnow + (REAL)nx;
      }
      prtcl_xint[i] = xposnow;

      diff_y[i] =(yread - prtcl_yint[i])/((double) ntempread);
      if (abs(yread - prtcl_yint[i]) > ((double) ny/2.0)){
	bc_chky[i] = 1;
	if (yread > prtcl_yint[i]){
	  yread = yread - (REAL) ny + 1.0;
	  diff_y[i] =(yread - prtcl_yint[i])/((double) ntempread);
	}else{
	  yread = yread + (REAL) ny - 1.0;
	  diff_y[i] =(yread - prtcl_yint[i])/((double) ntempread);
	}
      }
      yposnow = prtcl_yint[i] + diff_y[i];
      if (yposnow>(REAL)ny){
	yposnow = yposnow - (REAL)ny + 1.0;
	bc_chky[i] = 2;
      }
      if (yposnow < 1.0){
	yposnow = yposnow + (REAL)ny - 1.0;
	bc_chky[i] = 2;
      }
      prtcl_yint[i] = yposnow;

      prtcl_xidx[i] = round(xposnow);
      prtcl_yidx[i] = round(yposnow);
      
      
      Tread = prtcl2_T[chk_id[i]];
      diff_T[i] = (Tread - prtcl_Tint[i])/((double) ntempread);
      Tnow = prtcl_Tint[i] + diff_T[i];
            
      if (bc_chky[i]==1){
	Tnow = prtcl_Tint[i];
      }
      if (bc_chky[i]==2){
	Tnow = prtcl2_T[chk_id[i]];
	bc_chky[i]=1;
      }
      prtcl_Tint[i] = Tnow;
    }

    for (i=0;i<NX+2;i++){
      for(j=0;j<NY+2;j++){
	idx1=j+(NY+2)*i;
	nudge_vol[idx1]=0;
      }
    }

    for (i=1;i<NX+1;i++){
      for (j=1;j<NY+1;j++){
	idxread = (j-1)+NY*(i-1);
	tempreadint[idxread] = 0.0;
      }
    }

    for (k=0;k<Nprtcl;k++){
      x_now = prtcl_xidx[k] - mex*NX;
      y_now = prtcl_yidx[k] - mey*NY;
      for (i1=-3;i1<4;i1++){
	for (j1=-3;j1<4;j1++){
	  xx = abs(i1);
	  yy = abs(j1);
	  i = x_now+i1;
	  j = y_now+j1;
	  if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	    idx1=j+(NY+2)*i;
	    idxread = (j-1)+NY*(i-1);
	    tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	    nudge_vol[idx1] += 1;
	  }
	}
      }
      if ((prtcl_xidx[k]<4)&&(mex==NPX-1)){
	tmp_xnow = NX + prtcl_xidx[k];
	for (i1=-3;i1<1;i1++){
	  for (j1=-3;j1<4;j1++){
	    i = tmp_xnow+i1;
	    j = y_now+j1;
	    if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	      //printf("%d %d %d %d \n",i,j,prtcl_xidx[k],istep);
	      idx1=j+(NY+2)*i;
	      idxread = (j-1)+NY*(i-1);
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }      
	  }
	}
      }
      if ((prtcl_xidx[k]>nx-2)&&(mex==0)){
	tmp_xnow = prtcl_xidx[k]-nx;
	for (i1=0;i1<4;i1++){
	  for (j1=-3;j1<4;j1++){
	    i = tmp_xnow+i1;
	    j = y_now+j1;
	    if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	      //printf("%d %d %d %d After \n",i,j,prtcl_xidx[k],istep);
	      idx1=j+(NY+2)*i;
	      idxread = (j-1)+NY*(i-1);
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }
	  }
	}
      }
      prtcl_x[k] = prtcl2_x[k];
      prtcl_y[k] = prtcl2_y[k];
      prtcl_T[k] = prtcl2_T[k];
    }
  }else if((istep%ntempread)==0){
    for (i=0;i<Nprtcl;i++){
      prtcl_xint[i]=prtcl_x[chk_id[i]];
      prtcl_yint[i]=prtcl_y[chk_id[i]];
      prtcl_Tint[i]=prtcl_T[chk_id[i]];
      prtcl_xidx[i] = round(prtcl_x[chk_id[i]]);
      prtcl_yidx[i] = round(prtcl_y[chk_id[i]]);
    }
    for (i=0;i<NX+2;i++){
      for(j=0;j<NY+2;j++){
	idx1=j+(NY+2)*i;
	nudge_vol[idx1]=0;
      }
    }
    for (i=1;i<NX+1;i++){
      for (j=1;j<NY+1;j++){
	idxread = (j-1)+NY*(i-1);
	tempreadint[idxread] = 0.0;
      }
    }
    for (k=0;k<Nprtcl;k++){
      x_now = prtcl_xidx[k] - mex*NX;
      y_now = prtcl_yidx[k] - mey*NY;
      for (i1=-3;i1<4;i1++){
	for (j1=-3;j1<4;j1++){
	  i = x_now+i1;
	  j = y_now+j1;
	  if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	    idx1=j+(NY+2)*i;
	    idxread = (j-1)+NY*(i-1);
	    tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	    nudge_vol[idx1] += 1;
	  }
	}
      }
      if ((prtcl_xidx[k]<4)&&(mex==NPX-1)){
	int tmp_xnow = NX+x_now;
	for (i1=-3;i1<1;i1++){
	  for (j1=-3;j1<4;j1++){
	    i = tmp_xnow+i1;
	    j = y_now+j1;
	    if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	      idx1=j+(NY+2)*i;
	      idxread = (j-1)+NY*(i-1);
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }      
	  }
	}
      }
      if ((prtcl_xidx[k]>nx-2)&&(mex==0)){
	int tmp_xnow = x_now-NX;
	for (i1=0;i1<4;i1++){
	  for (j1=-3;j1<4;j1++){
	    i = tmp_xnow+i1;
	    j = y_now+j1;
	    if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	      idx1=j+(NY+2)*i;
	      idxread = (j-1)+NY*(i-1);
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }
	  }
	}
      }
    }
  }else{
    for (i=0;i<Nprtcl;i++){
      xposnow = prtcl_xint[i] + diff_x[i];
      if (xposnow>(REAL)nx){
	xposnow = xposnow - (REAL)nx;
      }
      if (xposnow < 0.0){
	xposnow = xposnow + (REAL)nx;
      }
      prtcl_xint[i] = xposnow;

      yposnow = prtcl_yint[i] + diff_y[i];
      if (yposnow>(REAL)ny){
	yposnow = yposnow - (REAL)ny + 1.0;
	bc_chky[i] = 2;
      }
      if (yposnow < 1.0){
	yposnow = yposnow + (REAL)ny - 1.0;
	bc_chky[i] = 2;
      }
      prtcl_yint[i] = yposnow;
      prtcl_xidx[i] = round(xposnow);
      prtcl_yidx[i] = round(yposnow);

      Tnow = prtcl_Tint[i] + diff_T[i];
      if (bc_chky[i]==1){
	Tnow = prtcl_Tint[i];
      }
      if (bc_chky[i]==2){
	Tnow = prtcl_T[chk_id[i]];
	bc_chky[i]=1;
      }      
      prtcl_Tint[i] = Tnow;
    }
    for (i=0;i<NX+2;i++){
      for(j=0;j<NY+2;j++){
	idx1=j+(NY+2)*i;
	nudge_vol[idx1]=0;
      }
    }
    for (i=1;i<NX+1;i++){
      for (j=1;j<NY+1;j++){
	idxread = (j-1)+NY*(i-1);
	tempreadint[idxread] = 0.0;
      }
    }

    for (k=0;k<Nprtcl;k++){
      x_now = prtcl_xidx[k] - mex*NX;
      y_now = prtcl_yidx[k] - mey*NY;
      for (i1=-3;i1<4;i1++){
	for (j1=-3;j1<4;j1++){
	  i = x_now+i1;
	  j = y_now+j1;
	  if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	    idx1=j+(NY+2)*i;
	    idxread = (j-1)+NY*(i-1);
	    tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	    nudge_vol[idx1] += 1;
	  }
	}
      }
      if ((prtcl_xidx[k]<4)&&(mex==NPX-1)){
	int tmp_xnow = NX+x_now;
	for (i1=-3;i1<1;i1++){
	  for (j1=-3;j1<4;j1++){
	    i = tmp_xnow+i1;
	    j = y_now+j1;
	    if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	      idx1=j+(NY+2)*i;
	      idxread = (j-1)+NY*(i-1);
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }      
	  }
	}
      }
      if ((prtcl_xidx[k]>nx-2)&&(mex==0)){
	int tmp_xnow = x_now-NX;
	for (i1=0;i1<4;i1++){
	  for (j1=-3;j1<4;j1++){
	    i = tmp_xnow+i1;
	    j = y_now+j1;
	    if (((i>0)&&(i<=NX))&&((j>0)&&(j<=NY))){
	      idx1=j+(NY+2)*i;
	      idxread = (j-1)+NY*(i-1);
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl_Tint[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }
	  }
	}
      }
    }
  }
}


void snapshottemp(REAL *rhopre,REAL *temppre,REAL *utot, REAL *vtot) {
  char fname[128];
  char gname[128];
  herr_t status;

  //  double *vortsq,*velsq;
  int idxup,idxdn,idxlt,idxrt;
  double dvbydx,dubydy;
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
    printf("Couldn't create snaptemp");
    return(8);
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

void snapshot(REAL *rhopre,REAL *temppre,REAL *utot,REAL *vtot) {

  int plusx1,plusy1;
  int minusx1,minusy1;

  int plusx2,plusy2;
  int minusx2,minusy2;


  REAL flusso,muav,mutotpre,mutotpost,mforcing,A,B,g1,g2,PPP,INTERA,INTERB;
  REAL LAPrhoA,LAPrhoB,gradxA,gradxB;

  REAL s1,s2,aaa,bbb,uav,div,div2,term;


  FILE *fout1;
  char filename[128];

  sprintf(filename,"all_xynew.%d.%d.%d.dat",istep,mex,mey);
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

void collisg(pop_type *f,pop_type *feq,REAL *tau,int *nudge_vol,REAL *tempreadint,REAL *temppre,double *nud_term){//INIZIO ROUTINE COLLIS

  REAL omeganow;
  REAL B;
  REAL Q;
  int count1 = 0;
  Q = 0.0;
  for(i=1;i<NX+1;i++){
    for(j=1;j<NY+1;j++){
      idx1=j+(NY+2)*i;
      idxread = (j-1) + NY*(i-1);
      omeganow=1./tau[idx1];
      if (nudge_vol[idx1]>0){
	count1 = count1 + nudge_vol[idx1];
	B = 1.0 - 0.5*omeganow;
	Q = alpha*(temppre[idx1] - tempreadint[idxread]);
	for(k=0;k<npop;k++){
	  f[idx1].p[k]=f[idx1].p[k]*(1.-omeganow)+omeganow*feq[idx1].p[k] - B*ww[k]*Q;
	}
	nud_term[idx1] = -1.*Q;
      }else{
	B = 0.0;
	for(k=0;k<npop;k++){
	  f[idx1].p[k]=f[idx1].p[k]*(1.-omeganow)+omeganow*feq[idx1].p[k] + B*ww[k]*Q;
	}
      }
    }
  }
}//FINE ROUTINE COLLIS

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


void dellar(){///////////////////MODELLO BIDIMENSIONALE PER LE VELOCITA'

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
}

void set_nudge_vol(int *nudge_vol,REAL *prtcl_x,REAL *prtcl_y,REAL *prtcl_T,int *prtcl_ids,REAL *prtcl_xint,REAL *prtcl_yint,REAL *prtcl_Tint,REAL *tempreadint,double *nud_term){
  if (contin==0){
    int nsquares,nsqx,nsqy,kx,ky,i1,j1,count,countred;
    int lsq = 6;
    double temp1,temp2;
    int xx,yy;
    nsqy = 40;
    nsqx = 40;
    nsquares = nsqx*nsqy;
    if (nsquares!=Nprtcl){
      printf("nsquares ne Nprtcl");
      exit(8);
    }
    int xsqcoord[nsqx],ysqcoord[nsqy];
    count = 0;
    for (i=0;i<NX+2;i++){
      for (j=0;j<NY+2;j++){
	indx1 = j + (NY+2)*i;
	nudge_vol[indx1] = 0;
	nud_term[indx1] = 0;
      }
    }
    xsqcoord[0] = 5;
    for(i=1;i<nsqx;i++){
      xsqcoord[i] = xsqcoord[i-1]+42+(i%2);
    }
    xsqcoord[nsqx-1]+=2;
    xsqcoord[nsqx-1]+=1;

    ysqcoord[0] = 10;
    for(i=1;i<nsqy;i++){
      ysqcoord[i] = ysqcoord[i-1]+21;
    }
    count = 0;
    for (i1=0;i1<nsqx;i1++){
      for (j1=0;j1<nsqy;j1++){
	for (i=1;i<NX+1;i++){
	  for (j=1;j<NY+1;j++){
	    indx1 = j + (NY+2)*i;
	    indx=i+mex*NX;
	    indy=j+mey*NY;
	    if ((abs(indx-xsqcoord[i1])<4)&&(abs(indy-ysqcoord[j1])<4)){
	      nudge_vol[indx1]=1;
	    }
	  }
	}
      }
    }
    double prtcl2_x[Nprtcl],prtcl2_y[Nprtcl],prtcl2_T[Nprtcl];
    int prtcl2_id[Nprtcl];
    if (me==0){
      FILE *fin;
      char fname[128];

      sprintf(fname,"../../St0.5/beta1/init_prtclx0.in",nstart);
      fin = fopen(fname,"r");
      fread(prtcl2_x,sizeof(double),Nprtcl,fin);
      fclose(fin);

      sprintf(fname,"../../St0.5/beta1/init_prtcly0.in",nstart);
      fin = fopen(fname,"r");
      fread(prtcl2_y,sizeof(double),Nprtcl,fin);
      fclose(fin);

      sprintf(fname,"../../St0.5/beta1/init_prtclT0.in",nstart);
      fin = fopen(fname,"r");
      fread(prtcl2_T,sizeof(double),Nprtcl,fin);
      fclose(fin);

      sprintf(fname,"../../St0.5/beta1/init_prtcl_id0.in",nstart);
      fin = fopen(fname,"r");
      fread(prtcl2_id,sizeof(int),Nprtcl,fin);
      fclose(fin);
    
    }
  

    MPI_Bcast(prtcl2_x,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(prtcl2_y,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(prtcl2_T,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(prtcl2_id,Nprtcl,MPI_INT,0,MPI_COMM_WORLD);
 
    for (i=0;i<Nprtcl;i++){
      prtcl_x[i] = prtcl2_x[i];
      prtcl_y[i] = prtcl2_y[i];
      prtcl_T[i] = prtcl2_T[i];
      prtcl_ids[i] = prtcl2_id[i];
      prtcl_xint[i] = prtcl2_x[i];
      prtcl_yint[i] = prtcl2_y[i];
      prtcl_Tint[i] = prtcl2_T[i];
    }
  }else{
    double prtcl2_x[Nprtcl],prtcl2_y[Nprtcl],prtcl2_T[Nprtcl];
    int prtcl2_id[Nprtcl];
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

    //sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata1/Particles%d.h5",nstart);
    sprintf(fname,"/mnt/petaStor/agasthya/LagrangianNudging/High_Ra/St0.5/beta1/Prtcl_Veldata1/Particles%d.h5",istep+ntempread-1);
    /* if (nstart>2000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata2/Particles%d.h5",nstart); */
    /* if (nstart>4000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata3/Particles%d.h5",nstart); */
    /* if (nstart>6000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata4/Particles%d.h5",nstart); */
    /* if (nstart>8000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata5/Particles%d.h5",nstart); */
    /* if (nstart>10000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata6/Particles%d.h5",nstart); */
    /* if (nstart>12000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata7/Particles%d.h5",nstart); */
    /* if (nstart>14000000) sprintf(fname,"../Prtcl_dump1_Guo/Prtcl_Veldata8/Particles%d.h5",nstart); */
    hid_t  file_id = H5Fopen(fname, H5F_ACC_RDWR,H5P_DEFAULT);
    if (file_id<0){
      printf("Particle file not opened %d \n",istep);
      exit(8);
    }
    hid_t ememspace = H5Screate_simple(1, edimens_3d, NULL );
    hid_t efilespace = H5Screate_simple( 1, fdimens_3d, NULL );
    hid_t status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d );
    status = H5Sselect_hyperslab( efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
    
    dataset_id = H5Dopen (file_id, "/x", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, prtcl2_x);
    status = H5Pclose( plist_id );
    status = H5Dclose(dataset_id);

    plist_id = H5Pcreate (H5P_FILE_ACCESS);
    hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    dataset_id = H5Dopen (file_id, "/y", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, prtcl2_y);
    status = H5Pclose( plist_id );
    status = H5Dclose(dataset_id);

    plist_id = H5Pcreate (H5P_FILE_ACCESS);
    hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    dataset_id = H5Dopen (file_id, "/id", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, ememspace, efilespace, H5P_DEFAULT, prtcl2_id);
    status = H5Pclose( plist_id );
    status = H5Dclose(dataset_id);

    plist_id = H5Pcreate (H5P_FILE_ACCESS);
    hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    dataset_id = H5Dopen (file_id, "/temp", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ememspace, efilespace, H5P_DEFAULT, prtcl2_T);

    status = H5Pclose( plist_id );
    status = H5Dclose(dataset_id);

    H5Sclose (efilespace);
    H5Sclose (ememspace);
    
    status = H5Fclose(file_id);

    MPI_Bcast(&prtcl2_x,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&prtcl2_y,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&prtcl2_T,Nprtcl,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&prtcl2_id,Nprtcl,MPI_INT,0,MPI_COMM_WORLD);
    double xposnow,yposnow;
    int prtcl_xidx[Nprtcl],prtcl_yidx[Nprtcl],x_now,x_nowbc;
    for (i=0;i<NX+2;i++){
      for (j=0;j<NY+2;j++){
	indx1 = j + (NY+2)*i;
	nudge_vol[indx1] = 0;
	nud_term[indx1] = 0;
      }
    }
    for (k=0;k<Nprtcl;k++){
      xposnow = prtcl2_x[k];
      yposnow = prtcl2_y[k];
      prtcl_xidx[k] = round(xposnow);
      prtcl_yidx[k] = round(yposnow);
      for (i=1;i<NX+1;i++){
	for (j=1;j<NY+1;j++){
	  idx1=j+(NY+2)*i;
	  indx = mex * NX + i;
	  indy = mey * NY + j;
	  idxread = (j-1)+NY*(i-1);
	  x_now = prtcl_xidx[k];
	  x_nowbc = nx - (x_now % nx);
	  if (abs(indy-prtcl_yidx[k])<4){
	    if  (abs(indx-x_now)<4) {
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl2_T[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }
	    if (abs( (x_now - nx) - indx)<4) {
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl2_T[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }
	    if (abs( (x_now + nx) - indx)<4) {
	      tempreadint[idxread] = (tempreadint[idxread]*(REAL)nudge_vol[idx1] + prtcl2_T[k])/((REAL) nudge_vol[idx1]+1.0);
	      nudge_vol[idx1] += 1;
	    }
	  }
	}
      }
      prtcl_x[k] = prtcl2_x[k];
      prtcl_y[k] = prtcl2_y[k];
      prtcl_T[k] = prtcl2_T[k];
      prtcl_ids[k] = prtcl2_id[k];
      prtcl_xint[k] = prtcl2_x[k];
      prtcl_yint[k] = prtcl2_y[k];
      prtcl_Tint[k] = prtcl2_T[k];
    }
  }
}

int main(int argc, char *argv[]) {

  FILE *fout4,*fout,*fout6,*fout7,*fin,*fout1,*fout2,*fout3,*fout5,*fout8,*fout9;
  REAL forcing;
  pop_type *f,*feq,*g,*geq,*memoryup1,*memoryup2,*memorydown1,*memorydown2;
  REAL *upre,*vpre,*rhopre,*temppre,*forcex,*forcey,*tempread,*tempreadint;
  REAL *tau,*taug;
  int *nudge_vol,*prtcl_ids,*chk_id,*bc_chkx,*bc_chky;
  double *nud_term;
  REAL *diff_x,*diff_y,*diff_T;
  REAL r1,r2,r3,A,B,g1,g2,totalenergy,totalenergy1,totaltemp,totaltemp1;
  REAL *prtcl_x,*prtcl_y,*prtcl_T,*prtcl_xint,*prtcl_yint,*prtcl_Tint;
  int i;
  myInit(argc,argv);
  printf("Rank %d: NX=%d, NY=%d, Nprtcl=%d\n", me, NX, NY, Nprtcl);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  f = (pop_type *) malloc((NX+2)*(NY+2)*sizeof(pop_type));
  feq = (pop_type *) malloc((NX+2)*(NY+2)*sizeof(pop_type));
  g = (pop_type *) malloc((NX+2)*(NY+2)*sizeof(pop_type));
  geq = (pop_type *) malloc((NX+2)*(NY+2)*sizeof(pop_type));
  memoryup1 = (pop_type *) malloc((NX+2)*sizeof(pop_type));
  memoryup2 = (pop_type *) malloc((NX+2)*sizeof(pop_type));
  memorydown1 = (pop_type *) malloc((NX+2)*sizeof(pop_type));
  memorydown2 = (pop_type *) malloc((NX+2)*sizeof(pop_type));
  upre = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  vpre = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  rhopre = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  temppre = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  forcex = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  forcey = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  tau = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  taug = (REAL *) malloc((NX+2)*(NY+2)*sizeof(REAL));
  tempread = (REAL *) malloc(NX * NY * sizeof(REAL));
  tempreadint = (REAL *) malloc(NX * NY * sizeof(REAL));
  nudge_vol = (int *) malloc((NX+2) * (NY+2) * sizeof(int));
  nud_term = (double *) malloc((NX+2)*(NY+2)*sizeof(double));
  //Particle Stuff
  /* prtcl_mex = (int *) malloc(Nprtcl*sizeof(int)); */
  /* prtcl_mey = (int *) malloc(Nprtcl*sizeof(int)); */
  chk_id = (int *) malloc(Nprtcl*sizeof(int));
  prtcl_ids = (int *) malloc(Nprtcl*sizeof(int));
  prtcl_x = (REAL *) malloc(Nprtcl*sizeof(REAL));
  prtcl_y = (REAL *) malloc(Nprtcl*sizeof(REAL));
  prtcl_T = (REAL *) malloc(Nprtcl*sizeof(REAL));
  prtcl_xint = (REAL *) malloc(Nprtcl*sizeof(REAL));
  prtcl_yint = (REAL *) malloc(Nprtcl*sizeof(REAL));
  prtcl_Tint = (REAL *) malloc(Nprtcl*sizeof(REAL));
  diff_x = (REAL *) malloc(Nprtcl*sizeof(REAL));
  diff_y = (REAL *) malloc(Nprtcl*sizeof(REAL));
  diff_T = (REAL *) malloc(Nprtcl*sizeof(REAL));
  bc_chkx = (int *) malloc(Nprtcl*sizeof(int));
  bc_chky = (int *) malloc(Nprtcl*sizeof(int));
  dellar();
  set_nudge_vol(nudge_vol,prtcl_x,prtcl_y,prtcl_T,prtcl_ids,prtcl_xint,prtcl_yint,prtcl_Tint,tempreadint,nud_term);
  printf("Nudge vol set %d \n", me);
  fflush(stdout);

  if(scratch) {
    inithydro(upre,vpre,rhopre,temppre,f,g,nudge_vol,tempreadint);
    equilif(upre,vpre,rhopre,feq);
    initpop(f,feq);
    equilig(upre,vpre,temppre,geq);
    initpop(g,geq);
  } else {
    restoreg(g);
    restoref(f);
    restoregpre(geq);
    restorefpre(feq);
    inithydro(upre,vpre,rhopre,temppre,feq,geq,nudge_vol,tempreadint);
    if (contin==0){
      equilif(upre,vpre,rhopre,feq);
      initpop(f,feq);
    }
    /* equilig(upre,vpre,temppre,geq); */
    /* initpop(g,geq); */
  }
  printf("Initialised %d \n", me);
  fflush(stdout);

  //////////////////////////MAIN LOOP//////////////////////////
  
  icount=0;
  icountconfig=0;
  icountprof=0;



  if (me==0) {
    fout6=fopen("evolution.dat","w");
  }

  for(istep=1+nstart;istep<=nsteps+nstart;istep++){ //nsteps
    //if ((me==0)&&((istep%10)==0)) printf("%d \n",istep);
    //if ((me==0)) printf("%d \n",istep);
    pbc(f);
    pbc(g);
    read_prtcl(tempreadint,nudge_vol,chk_id,bc_chkx,bc_chky,diff_x,prtcl_xint,diff_y,prtcl_yint,diff_T,prtcl_Tint,prtcl_x,prtcl_y,prtcl_T,prtcl_ids);
    moveplusforcingconstructWW(f,g,rhopre,temppre,forcex,forcey,memoryup1,memoryup2,memorydown1,memorydown2);
    ////////////////////////////////////////
    /////and the move-machinery is over/////
    ////////////////////////////////////////
    pbc(f);
    pbc(g);
    hydrovar(upre,vpre,rhopre,temppre,f,g,nudge_vol,tempreadint,forcex,forcey,nud_term);
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
    collisg(g,geq,taug,nudge_vol,tempreadint,temppre,nud_term);

    pbc(f);
    pbc(g);
    
    if((istep%noutconfig)==0){ 
      snapshottemp(rhopre,temppre,upre,vpre);
    }
    if ((istep%nout)==0) profile(rhopre,temppre,upre,vpre);
      
    if((istep%5000)==0){  ////TO GET THE PROFILES
      //     media(rhopre,temppre);
      //
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
    if ((istep%2000000)==0){
      dumpf(f);
      dumpg(g);
    }
  }//END OF MAIN LOOP

  if (me==0){
    fclose(fout6);
  }
  dumpf(f);
  dumpg(g);
  myFinalize();
}
