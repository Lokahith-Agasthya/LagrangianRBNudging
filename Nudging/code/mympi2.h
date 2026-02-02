#include <mpi.h>
typedef struct{
  double xprtcl,yprtcl;
  double uprtcl,vprtcl;
  double xoprtcl,yoprtcl;
  double uoprtcl,voprtcl;
} prtcl_type;

int dims[2], period[2], coords[2];
int NP,NPX,NPY, me, mex, mey,i;
int NX, NY;
int NXP2, NYP2;
int NXNY;
int NXP6, NYP6,*prtcl_ids;
#define Nprtcl 40*40
MPI_Datatype MPI_Poptype, MPI_Poptypey, MPI_Poptypex,MPI_Prtcltype;
MPI_Comm MPI_COMM_CART, MPI_COMM_ALONG_X, MPI_COMM_ALONG_Y;

/* MPI Init here */
void myInit(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  //NP and me respectively store the number of processes and the rank (address) of the current process
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  dims[0]=0;
  dims[1]=0;
  MPI_Dims_create(NP, 2, dims); //Create an appropriate 2D array of NP nodes. dims stores the no. of processors
  //in each dimension

  MPI_Type_contiguous( 9, MPI_DOUBLE, &MPI_Poptype); //Creates a new datatype which takes 9 array elements and makes it one contiguous data. Sending the 9 entries separately and sending the contiguous data once are identical
  MPI_Type_commit(&MPI_Poptype); //Compulsory syntax line to be used before the new data type can be used

  MPI_Type_contiguous( 8, MPI_DOUBLE, &MPI_Prtcltype);
  MPI_Type_commit(&MPI_Prtcltype);

  /*condizioni periodiche lungo x*/
  period[0]=1;
  period[1]=1;

  MPI_Cart_create( MPI_COMM_WORLD, 2, dims, period, 0, &MPI_COMM_CART);

  /*Costruisce i sub-comunicatori lungo X e Y*/
  coords[0]=1; coords[1]=0;
  MPI_Cart_sub(MPI_COMM_CART, coords, &MPI_COMM_ALONG_X);
  coords[0]=0; coords[1]=1;
  MPI_Cart_sub(MPI_COMM_CART, coords, &MPI_COMM_ALONG_Y);
  //  if (MPI_COMM_ALONG_X == MPI_COMM_NULL) fprintf(stderr,"rank %d: ALONG_X is NULL\n", me);
  // if (MPI_COMM_ALONG_Y == MPI_COMM_NULL) fprintf(stderr,"rank %d: ALONG_Y is NULL\n", me);

  /* Rank lungo X ed Y*/
  MPI_Comm_rank(MPI_COMM_ALONG_X, &mex);
  MPI_Comm_rank(MPI_COMM_ALONG_Y, &mey);


  NPX=dims[0];
  NPY=dims[1];

  NX = nx/NPX;
  NY = ny/NPY;

  if ((NX*NPX!=nx)||(NY*NPY!=ny)){
    if (me==0){
      fprintf(stderr,"Error: chosen dimension %d %d not divisible on the processor mesh %d %d\n",nx,ny,NPX,NPY);
    }
    exit(8);
  }

  MPI_Type_vector( (2+NX), 1, 2+NY, MPI_Poptype, &MPI_Poptypey);
  MPI_Type_vector( 1, (NY+2), (2+NY)*(2+NX), MPI_Poptype, &MPI_Poptypex);
  MPI_Type_commit(&MPI_Poptypey);
  MPI_Type_commit(&MPI_Poptypex);


  if (me==0){
    fprintf(stderr,"SIZEX = %d nx = %d on %d processors\n",nx,NX,NPX);
    fprintf(stderr,"SIZEY = %d ny = %d on %d processors\n",ny,NY,NPY);
  }
  for (i=0;i<2;i++){
    fprintf(stderr,"direction %d number of procs %d\n",i,dims[i]);
  }
  
  NXP6 = NX + 6;
  NYP6 = NY + 6;
  NXP2 = NX + 2;
  NYP2 = NY + 2;
  NXNY=NXP2*NYP2;
  // fprintf(stderr, "rank %d: NX=%d NY=%d NXP6=%d NYP6=%d\n", me, NX, NY, NXP6, NYP6);
  // printf("Going to the next line %d \n",me);
  MPI_Barrier(MPI_COMM_WORLD);

}


void communication(REAL *u){
  REAL sendbufposx[NY];
  REAL sendbufnegx[NY];
  REAL sendbufposy[NX+2];
  REAL sendbufnegy[NX+2];
  REAL recvbufposy[NX+2];
  REAL recvbufnegy[NX+2];
  int indx1,procpx,procmx;
  MPI_Status status1,status2;
  for (int i = 1;i<NY+1;i++){
    indx1 = i + (NY+2)*NX;
    sendbufposx[i-1] = u[indx1];
    indx1 = i + (NY+2);
    sendbufnegx[i-1] = u[indx1];
  }
  procpx=(mex+1+NPX)%NPX; procmx=(mex-1+NPX)%NPX;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Sendrecv(&sendbufposx,NY,MPI_DOUBLE,procpx,10,u+1,NY,MPI_DOUBLE,procmx,10,MPI_COMM_ALONG_X,&status1);
  MPI_Sendrecv(&sendbufnegx,NY,MPI_DOUBLE,procmx,11,u+(NX+1)*(NY+2),NY,MPI_DOUBLE,procpx,11,MPI_COMM_ALONG_X,&status2);

  for (int i = 0;i<NX+2;i++){
    indx1 = NY + (NY+2)*i;
    sendbufposy[i] = u[indx1];
    indx1 = 1 + (NY+2)*i;
    sendbufnegy[i] = u[indx1];
  }
  procpx=(mey+1+NPY)%NPY; procmx=(mey-1+NPY)%NPY;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Sendrecv(&sendbufposy,NX+2,MPI_DOUBLE,procpx,10,&recvbufnegy,NX+2,MPI_DOUBLE,procmx,10,MPI_COMM_ALONG_Y,&status1);
  MPI_Sendrecv(&sendbufnegy,NX+2,MPI_DOUBLE,procmx,11,&recvbufposy,NX+2,MPI_DOUBLE,procpx,11,MPI_COMM_ALONG_Y,&status2);

  for (int i=0;i<NX+2;i++){
    indx1 = NY+1 + (NY+2)*i;
    u[indx1] = recvbufposy[i];
    indx1 = (NY+2)*i;
    u[indx1] = recvbufnegy[i];
  }
}


void myFinalize()
{
  MPI_Finalize(); 
}
