#include <mpi.h>
typedef struct{
  double xprtcl,yprtcl;
  double uprtcl,vprtcl;
  double xoprtcl,yoprtcl;
  double uoprtcl,voprtcl;
} prtcl_type;

static prtcl_type *prtcl_sendposx = NULL;
static prtcl_type *prtcl_sendnegx = NULL;
static int *prtcl_id_sendposx = NULL;
static int *prtcl_id_sendnegx = NULL;

#define Nprtcl 1600
#define tracer 0 //Set tracer to 1 and heavyprtcl to 0 for tracer evolution
#define heavyprtcl 1 //Set heavyprtcl to 1 and tracer to 0 to evolve according to Maxey-Riley equation
#define taulag 2150.0
#define onebytaulag (1/taulag)
#define beta 0.1 //relevant only if heavyprtcl is set to 1 and tracer set to 0

#define Prtcl_Writefull 1 //Flag to write the full particle information. Use only for benchmarking

extern int istep;

int dims[2], period[2], coords[2];
int NP,NPX,NPY, me, mex, mey,i;
int NX, NY;
int NXP2, NYP2;
int NXNY;
int NXP6, NYP6;
int proc1,proc2,conta;
int nprtcl_local;
int nprtcl_max;
int nprtcl_max_comm;
int i,j,idx1,idxright,idxleft,idxtop,idxbelow;

MPI_Datatype MPI_Poptype, MPI_Poptypey, MPI_Poptypex,MPI_Prtcltype;
MPI_Comm MPI_COMM_CART, MPI_COMM_ALONG_X, MPI_COMM_ALONG_Y;

/* MPI Init here */
void myInit(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NP);

  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  dims[0]=0;
  dims[1]=0;
  MPI_Dims_create(NP, 2, dims);

  MPI_Type_contiguous( 9, MPI_DOUBLE, &MPI_Poptype); //Creates a new datatype which takes 9 array elements and makes it one contiguous data. Sending the 9 entries separately and sending the contiguous data once are identical
  MPI_Type_commit(&MPI_Poptype);
  
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
    fprintf(stderr,"SIZEY = %d ny = %d on %d processors\n",ny,NY,NPY);}
  for (i=0;i<2;i++){
    fprintf(stderr,"direction %d number of procs %d\n",i,dims[i]);
  }
  NXP6 = NX + 6;
  NYP6 = NY + 6;
  NXP2 = NX + 2;
  NYP2 = NY + 2;
  NXNY=NXP2*NYP2;

  nprtcl_local = (int) Nprtcl/NP;
  nprtcl_max = (int) (nprtcl_local*5.0);
  nprtcl_max_comm = (int) (nprtcl_local * 2);
  if (nprtcl_local==0){
    nprtcl_max = Nprtcl;
    nprtcl_max_comm = (int) (Nprtcl/2);
  }
}

void myFinalize()
{
  MPI_Finalize(); 
}


void init_prtcl(prtcl_type *prtcl,REAL *u, REAL *v,int *prtcl_ids,REAL *temp){
  nprtcl_local = 0;
  int i,j,indx1,indx,indy;
  int nsquares,nsqx,nsqy,kx,ky,i1,j1,count,countred;
  int lsq = 6;
  REAL temp1,temp2,xx,yy,tempup,tempvp;
  int *prtcl_ids2;
  nsqy = 40;
  nsqx = 40;
  int xsqcoord[nsqx],ysqcoord[nsqy];
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
  REAL del = 0.000001;
  for (i1=0;i1<nsqx;i1++){
    for (j1=0;j1<nsqy;j1++){
      for (i=1;i<NX+1;i++){
        for (j=1;j<NY+1;j++){
          indx1 = j + (NY+2)*i;
          indx=i+mex*NX;
          indy=j+mey*NY;
          if ((abs(indx-xsqcoord[i1])<1)&&(abs(indy-ysqcoord[j1])<1)){
            prtcl[nprtcl_local].xprtcl = (REAL) i;
            prtcl[nprtcl_local].yprtcl = (REAL) j;
            if (i==1){
              prtcl[nprtcl_local].xprtcl = (REAL) i + del;    
            }	 
            if (j==1){
              prtcl[nprtcl_local].yprtcl = (REAL) j + del;    
            }	 
            if (i==NX){
              prtcl[nprtcl_local].xprtcl = (REAL) i - del;
            }
            if (j==NY){
              prtcl[nprtcl_local].yprtcl = (REAL) j - del;
            }
            prtcl[nprtcl_local].xoprtcl = prtcl[nprtcl_local].xprtcl;
            prtcl[nprtcl_local].yoprtcl = prtcl[nprtcl_local].yprtcl;
            nprtcl_local+=1;
            count+=1;
          }
        }
      }
    }
  }

  MPI_Reduce(&count,&countred,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
  if (me==0){
    printf("count = %d, nsquares = %d Nprtcl = %d \n",countred,nsquares,Nprtcl);
  }
  printf("local = %d, max = %d, me = %d \n",nprtcl_local,nprtcl_max,me);
  REAL prtcl_T[nprtcl_local];
  REAL prtcl_x[nprtcl_local],prtcl_y[nprtcl_local],prtcl_u[nprtcl_local],prtcl_v[nprtcl_local];
  communication(u);
  communication(v);
  communication(temp);

  if (tracer==1){ 
    for (int i = 0;i<nprtcl_local;i++){
      xx = prtcl[i].xprtcl;
      yy = prtcl[i].yprtcl;
      bilinear(xx,yy,u,&tempup);
      bilinear(xx,yy,v,&tempvp);
      prtcl[i].uprtcl = tempup;
      prtcl[i].vprtcl = tempvp;
      prtcl[i].uoprtcl = prtcl[i].uprtcl;
      prtcl[i].voprtcl = prtcl[i].vprtcl;
      bilinear(xx,yy,temp,&tempvp);    
      prtcl_T[i] = tempvp;
      prtcl_x[i] = ((REAL) (mex*NX)) + xx;
      prtcl_y[i] = ((REAL) (mey*NY)) + yy;
      prtcl_u[i] = prtcl[i].uprtcl;
      prtcl_v[i] = prtcl[i].vprtcl;
    }

  } else if (heavyprtcl==1){ 
    for (int i = 0;i<nprtcl_local;i++){
      xx = prtcl[i].xprtcl;
      yy = prtcl[i].yprtcl;
      bilinear(xx,yy,u,&tempup);
      bilinear(xx,yy,v,&tempvp);
      prtcl[i].uprtcl = 0.98*tempup;  //This can be changed according to the initial condition we want. 
      //For now I am just making it go in the opposite direction to the fluid with a certain factor. 
      prtcl[i].vprtcl = 0.98*tempvp;
      prtcl[i].uoprtcl = prtcl[i].uprtcl;        
      prtcl[i].voprtcl = prtcl[i].vprtcl; 
      bilinear(xx,yy,temp,&tempvp);    
      prtcl_T[i] = tempvp;
      prtcl_x[i] = ((REAL) (mex*NX)) + xx;
      prtcl_y[i] = ((REAL) (mey*NY)) + yy;
      prtcl_u[i] = prtcl[i].uprtcl;
      prtcl_v[i] = prtcl[i].vprtcl;
    }
  }
  

  int pid_global[Nprtcl],np_global[NP],np_globalt[NP];
  int offs_arrt[NP],offs_arr[NP];
  for (i=0;i<NP;i++){
    np_globalt[i] = 0;
    offs_arrt[i] = 0;
  }
  np_globalt[me] = nprtcl_local;
  MPI_Allreduce(np_globalt,np_global,NP,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
  int  offs = 0;
  for (i=0;i<me;i++){
    offs+=np_global[i];
  }
  
  offs_arrt[me] = offs;
  //  printf("%d %d %d \n",nprtcl_local,offs,me);
  MPI_Allreduce(offs_arrt,offs_arr,NP,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);

  if (me == 0){
    for (i=0;i<Nprtcl;i++){
      pid_global[i] = i;
    }
  }

  
  prtcl_ids2 = (int *) malloc(nprtcl_local*sizeof(int));
  MPI_Scatterv(pid_global,np_global,offs_arr,MPI_INTEGER,prtcl_ids2,nprtcl_local,MPI_INTEGER,0,MPI_COMM_WORLD);
  for (i=0;i<nprtcl_local;i++){
    prtcl_ids[i] = prtcl_ids2[i];
  }
  REAL Tpglobal[Nprtcl],xpglob[Nprtcl],ypglob[Nprtcl],upglob[Nprtcl],vpglob[Nprtcl];
  MPI_Gatherv(prtcl_x,nprtcl_local,MPI_DOUBLE,xpglob,np_global,offs_arr,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(prtcl_y,nprtcl_local,MPI_DOUBLE,ypglob,np_global,offs_arr,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(prtcl_T,nprtcl_local,MPI_DOUBLE,Tpglobal,np_global,offs_arr,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(prtcl_u,nprtcl_local,MPI_DOUBLE,upglob,np_global,offs_arr,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(prtcl_v,nprtcl_local,MPI_DOUBLE,vpglob,np_global,offs_arr,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(prtcl_ids,nprtcl_local,MPI_INT,pid_global,np_global,offs_arr,MPI_INT,0,MPI_COMM_WORLD);
  if (me==0){
    FILE *fout11;
    char fname[128];
    sprintf(fname,"init_prtclx%d.in",nstart);
    fout11 = fopen(fname,"w");
    fwrite(xpglob,sizeof(REAL),Nprtcl,fout11);
    fclose(fout11);

    sprintf(fname,"init_prtcly%d.in",nstart);
    fout11 = fopen(fname,"w");
    fwrite(ypglob,sizeof(REAL),Nprtcl,fout11);
    fclose(fout11);

    sprintf(fname,"init_prtclT%d.in",nstart);
    fout11 = fopen(fname,"w");
    fwrite(Tpglobal,sizeof(REAL),Nprtcl,fout11);
    fclose(fout11);

    sprintf(fname,"init_prtclu%d.in",nstart);
    fout11 = fopen(fname,"w");
    fwrite(upglob,sizeof(REAL),Nprtcl,fout11);
    fclose(fout11);

    sprintf(fname,"init_prtclv%d.in",nstart);
    fout11 = fopen(fname,"w");
    fwrite(vpglob,sizeof(REAL),Nprtcl,fout11);
    fclose(fout11);

    sprintf(fname,"init_prtcl_id%d.in",nstart);
    fout11 = fopen(fname,"w");
    fwrite(pid_global,sizeof(int),Nprtcl,fout11);
    fclose(fout11);

  }
}

void bilinear(REAL xx, REAL yy,REAL *u,REAL *tempu){
  int x1 = (int) xx;
  int x2 = x1 + 1;
  int y1 = (int) yy;
  int y2 = y1 + 1;
  int indx11 = y1 + (NY+2)*x1;
  int indx12 = y2 + (NY+2)*x1;
  int indx21 = y1 + (NY+2)*x2;
  int indx22 = y2 + (NY+2)*x2;
  REAL u11 = u[indx11];
  REAL u12 = u[indx12];
  REAL u21 = u[indx21];
  REAL u22 = u[indx22];
  REAL val1 = u11*((REAL)x2 - xx)*((REAL)y2 - yy);
  REAL val2 = u21*(xx - (REAL)x1)*((REAL)y2 - yy);
  REAL val3 = u12*((REAL)x2 - xx)*(yy - (REAL)y1);
  REAL val4 = u22*(xx - (REAL)x1)*(yy - (REAL)y1);
  REAL val = val1 + val2 + val3 + val4;
  *tempu = val1 + val2 + val3 + val4;
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
  MPI_Sendrecv(&sendbufposx,NY,MPI_DOUBLE,procpx,10,u+1,NY,MPI_DOUBLE,procmx,10,MPI_COMM_ALONG_X,&status1);
  MPI_Sendrecv(&sendbufnegx,NY,MPI_DOUBLE,procmx,11,u+1+(NX+1)*(NY+2),NY,MPI_DOUBLE,procpx,11,MPI_COMM_ALONG_X,&status2);

  for (int i = 0;i<NX+2;i++){
    indx1 = NY + (NY+2)*i;
    sendbufposy[i] = u[indx1];
    indx1 = 1 + (NY+2)*i;
    sendbufnegy[i] = u[indx1];
  }
  procpx=(mey+1+NPY)%NPY; procmx=(mey-1+NPY)%NPY;
  MPI_Sendrecv(&sendbufposy,NX+2,MPI_DOUBLE,procpx,10,&recvbufnegy,NX+2,MPI_DOUBLE,procmx,10,MPI_COMM_ALONG_Y,&status1);
  MPI_Sendrecv(&sendbufnegy,NX+2,MPI_DOUBLE,procmx,11,&recvbufposy,NX+2,MPI_DOUBLE,procpx,11,MPI_COMM_ALONG_Y,&status2);

  for (int i=0;i<NX+2;i++){
    indx1 = NY+1 + (NY+2)*i;
    u[indx1] = recvbufposy[i];
    indx1 = (NY+2)*i;
    u[indx1] = recvbufnegy[i];
  }
}

void propagate_prtcl(prtcl_type *prtcl, REAL *u, REAL *v, REAL *ubef, REAL *vbef, REAL *dtu, REAL *dtv, REAL *dxu, REAL *dxv, REAL *dyu, REAL *dyv, REAL *temppre,int *prtcl_ids){
  REAL delta_t = 1.0;
  REAL tempup,tempvp,tempdtu,tempdtv,tempdxu,tempdxv,tempdyu,tempdyv;
  communication(u);
  communication(v);

  /* prtcl[i].uprtcl = tempup; */
  /* prtcl[i].vprtcl = tempvp; */
  if (tracer==1){ 
    for (int i = 0;i<nprtcl_local;i++){
      REAL xx = prtcl[i].xprtcl;
      REAL yy = prtcl[i].yprtcl;
      prtcl[i].uoprtcl = prtcl[i].uprtcl;
      prtcl[i].voprtcl = prtcl[i].vprtcl;
      bilinear(xx,yy,u,&tempup);
      bilinear(xx,yy,v,&tempvp);
      prtcl[i].xoprtcl = xx;
      prtcl[i].yoprtcl = yy;
      prtcl[i].xprtcl = xx + 1.5*delta_t*tempup-0.5*delta_t*prtcl[i].uoprtcl;
      prtcl[i].yprtcl = yy + 1.5*delta_t*tempvp-0.5*delta_t*prtcl[i].voprtcl;
      prtcl[i].uprtcl = tempup;
      prtcl[i].vprtcl = tempvp;
    }
  } else if (heavyprtcl==1){ 
    for (i=1;i<NX+1;i++){
      for (j=1;j<NY+1;j++){
        idx1=j+(NY+2)*i;
        idxright = j+(NY+2)*(i+1);
        idxleft = j+(NY+2)*(i-1);
        idxtop = (j+1)+(NY+2)*i;
        idxbelow = (j-1)+(NY+2)*i;

        dtu[idx1] = u[idx1] - ubef[idx1];
        dtv[idx1] = v[idx1] - vbef[idx1];
        dxu[idx1] = (u[idxright] - u[idxleft])*0.5;
        dxv[idx1] = (v[idxright] - v[idxleft])*0.5;
        dyu[idx1] = (u[idxtop] - u[idxbelow])*0.5;
        dyv[idx1] = (v[idxtop] - v[idxbelow])*0.5;
      }
    }
    if (mey==0){
      for(i=1;i<NX+1;i++){
        idx1 = 1+(NY+2)*i;
        idxtop = 2+(NY+2)*i;
        dyu[idx1] = (u[idxtop] - u[idx1]);
        dyv[idx1] = (v[idxtop] - v[idx1]);
      }
    }
    if (mey==NPY-1){
      for(i=1;i<NX+1;i++){
        idx1 = NY+(NY+2)*i;
        idxbelow = NY-1+(NY+2)*i;
        dyu[idx1] = (u[idx1] - u[idxbelow]);
        dyv[idx1] = (v[idx1] - v[idxbelow]);
      }
    }
    communication(dtu);
    communication(dtv);
    communication(dxu);
    communication(dxv);
    communication(dyu);
    communication(dyv);
    
    for (int i = 0;i<nprtcl_local;i++){
      REAL xx = prtcl[i].xprtcl;
      REAL yy = prtcl[i].yprtcl;
      prtcl[i].uoprtcl = prtcl[i].uprtcl;
      prtcl[i].voprtcl = prtcl[i].vprtcl;
      bilinear(xx,yy,u,&tempup);
      bilinear(xx,yy,v,&tempvp);
      bilinear(xx,yy,dtu,&tempdtu);
      bilinear(xx,yy,dtv,&tempdtv);
      bilinear(xx,yy,dxu,&tempdxu);
      bilinear(xx,yy,dxv,&tempdxv);
      bilinear(xx,yy,dyu,&tempdyu);
      bilinear(xx,yy,dyv,&tempdyv);
      prtcl[i].xoprtcl = xx;
      prtcl[i].yoprtcl = yy;
      prtcl[i].uprtcl = prtcl[i].uprtcl + delta_t*(beta*(tempdtu + tempup*tempdxu + tempvp*tempdyu) - onebytaulag*(prtcl[i].uprtcl - tempup));
      prtcl[i].vprtcl = prtcl[i].vprtcl + delta_t*(beta*(tempdtv + tempup*tempdxv + tempvp*tempdyv) - onebytaulag*(prtcl[i].vprtcl - tempvp));
      prtcl[i].xprtcl = xx + delta_t*prtcl[i].uprtcl;
      prtcl[i].yprtcl = yy + delta_t*prtcl[i].vprtcl;
    }
  }
  if ((istep%noutconfig)==0){
    //#if Particle
    MPI_Barrier(MPI_COMM_WORLD);
    prtcl_write(prtcl,temppre,prtcl_ids);
    //#endif 
    MPI_Barrier(MPI_COMM_WORLD);

    if (Prtcl_Writefull==1){
      REAL prtcl_dtu[nprtcl_local],prtcl_dtv[nprtcl_local],prtcl_dxu[nprtcl_local],prtcl_dxv[nprtcl_local],prtcl_dyu[nprtcl_local],prtcl_dyv[nprtcl_local];
      REAL prtcl_fu[nprtcl_local],prtcl_fv[nprtcl_local];
      if (tracer==1){
        for (i=1;i<NX+1;i++){
          for (j=1;j<NY+1;j++){
            idx1=j+(NY+2)*i;
            idxright = j+(NY+2)*(i+1);
            idxleft = j+(NY+2)*(i-1);
            idxtop = (j+1)+(NY+2)*i;
            idxbelow = (j-1)+(NY+2)*i;
    
            dtu[idx1] = u[idx1] - ubef[idx1];
            dtv[idx1] = v[idx1] - vbef[idx1];
            dxu[idx1] = (u[idxright] - u[idxleft])*0.5;
            dxv[idx1] = (v[idxright] - v[idxleft])*0.5;
            dyu[idx1] = (u[idxtop] - u[idxbelow])*0.5;
            dyv[idx1] = (v[idxtop] - v[idxbelow])*0.5;
          }
        }
        if (mey==0){
          for(i=1;i<NX+1;i++){
            idx1 = 1+(NY+2)*i;
            idxtop = 2+(NY+2)*i;
            dyu[idx1] = (u[idxtop] - u[idx1]);
            dyv[idx1] = (v[idxtop] - v[idx1]);
          }
        }
        if (mey==NPY-1){
          for(i=1;i<NX+1;i++){
            idx1 = NY+(NY+2)*i;
            idxbelow = NY-1+(NY+2)*i;
            dyu[idx1] = (u[idx1] - u[idxbelow]);
            dyv[idx1] = (v[idx1] - v[idxbelow]);
          }
        }
        communication(dtu);
        communication(dtv);
        communication(dxu);
        communication(dxv);
        communication(dyu);
        communication(dyv);
      }
      for (int i = 0;i<nprtcl_local;i++){
	REAL xx = prtcl[i].xoprtcl;
	REAL yy = prtcl[i].yoprtcl;
	bilinear(xx,yy,u,&tempup);
	bilinear(xx,yy,v,&tempvp);
	bilinear(xx,yy,dtu,&tempdtu);
	bilinear(xx,yy,dtv,&tempdtv);
	bilinear(xx,yy,dxu,&tempdxu);
	bilinear(xx,yy,dxv,&tempdxv);
	bilinear(xx,yy,dyu,&tempdyu);
	bilinear(xx,yy,dyv,&tempdyv);          
	prtcl_fu[i] = tempup;
	prtcl_fv[i] = tempvp;
	prtcl_dtu[i] = tempdtu;
	prtcl_dtv[i] = tempdtv;
	prtcl_dxu[i] = tempdxu;
	prtcl_dxv[i] = tempdxv;
	prtcl_dyu[i] = tempdyu;
	prtcl_dyv[i] = tempdyv;
      }
      int i,j,send_add; //send_add is number of particles stored in process with rank < me so that we know where to add the prtcl information in the global array
      char fname[128];
      int nprtcl_array[NP],nprtcl_arraytemp[NP];
      FILE *fout11;
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
      hid_t  file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
      if (file_id<0){
        printf("Particle file does not exist \n");
        exit(8);
      }
      prtcl_write_array(file_id, "/dtu", H5T_NATIVE_DOUBLE, prtcl_dtu, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
      prtcl_write_array(file_id, "/dtv", H5T_NATIVE_DOUBLE, prtcl_dtv, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
      prtcl_write_array(file_id, "/dxu", H5T_NATIVE_DOUBLE, prtcl_dxu, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
      prtcl_write_array(file_id, "/dxv", H5T_NATIVE_DOUBLE, prtcl_dxv, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
      prtcl_write_array(file_id, "/dyu", H5T_NATIVE_DOUBLE, prtcl_dyu, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
      prtcl_write_array(file_id, "/dyv", H5T_NATIVE_DOUBLE, prtcl_dyv, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
      prtcl_write_array(file_id, "/fu", H5T_NATIVE_DOUBLE, prtcl_fu, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
      prtcl_write_array(file_id, "/fv", H5T_NATIVE_DOUBLE, prtcl_fv, status, ememspace, efilespace, xfer_plist, dataset_id, ret);
      status = H5Fclose(file_id);

      H5Sclose (efilespace);
      H5Sclose (ememspace);
      
      status = H5Pclose(plist_id);
    
    
    }
  }
}
  void toggle_prtcl(prtcl_type **prtcl1,int **prtcl_ids1){

    prtcl_type *prtcl = *prtcl1;
    int *prtcl_ids = *prtcl_ids1;

    int sendposx = 0;
    int sendnegx = 0;
    int keephere = 0;
    int recvposx = 0;
    int recvnegx = 0;
    int procpx, procmx;
    MPI_Status status1,status2;

    if ((prtcl_sendposx == NULL)||(prtcl_sendnegx == NULL)){
      prtcl_sendposx = (prtcl_type *) realloc(prtcl_sendposx,sizeof(prtcl_type) * nprtcl_max_comm);
      prtcl_sendnegx = (prtcl_type *) realloc(prtcl_sendnegx,sizeof(prtcl_type) * nprtcl_max_comm);
    }
    if ((prtcl_id_sendposx == NULL)||(prtcl_id_sendnegx == NULL)){
      prtcl_id_sendposx = (int *) realloc(prtcl_id_sendposx,sizeof(int)*nprtcl_max_comm);
      prtcl_id_sendnegx = (int *) realloc(prtcl_id_sendnegx,sizeof(int)*nprtcl_max_comm);
    }
    for (int i = 0;i<nprtcl_local;i++){
      if (prtcl[i].xprtcl >= (REAL)NX){
	if (sendposx >= nprtcl_max_comm){
	  nprtcl_max_comm += 1;
	  prtcl_sendposx = (prtcl_type *) realloc(prtcl_sendposx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_sendnegx = (prtcl_type *) realloc(prtcl_sendnegx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_id_sendposx = (int *) realloc(prtcl_id_sendposx,sizeof(int)*nprtcl_max_comm);
	  prtcl_id_sendnegx = (int *) realloc(prtcl_id_sendnegx,sizeof(int)*nprtcl_max_comm);
	}
	if (mex==(NPX-1)){
	  prtcl[i].xprtcl = prtcl[i].xprtcl - (REAL) NX;
	}else{
	  prtcl[i].xprtcl = prtcl[i].xprtcl - (REAL) NX;
	}
	prtcl_sendposx[sendposx] = prtcl[i];
	prtcl_id_sendposx[sendposx] = prtcl_ids[i];
	sendposx = sendposx+1;
      } else if ((prtcl[i].xprtcl <= 0.0)&&(mex==0)){
	if (sendnegx >= nprtcl_max_comm){
	  nprtcl_max_comm += 1;
	  prtcl_sendposx = (prtcl_type *) realloc(prtcl_sendposx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_sendnegx = (prtcl_type *) realloc(prtcl_sendnegx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_id_sendposx = (int *) realloc(prtcl_id_sendposx,sizeof(int)*nprtcl_max_comm);
	  prtcl_id_sendnegx = (int *) realloc(prtcl_id_sendnegx,sizeof(int)*nprtcl_max_comm);
	}
	prtcl[i].xprtcl = prtcl[i].xprtcl + (REAL) NX;
	prtcl_sendnegx[sendnegx] = prtcl[i];
	prtcl_id_sendnegx[sendnegx] = prtcl_ids[i];
	sendnegx = sendnegx+1;
      } else if ((prtcl[i].xprtcl <= 0.0)&&(mex>0))  {
	if (sendnegx >= nprtcl_max_comm){
	  nprtcl_max_comm += 1;
	  prtcl_sendposx = (prtcl_type *) realloc(prtcl_sendposx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_sendnegx = (prtcl_type *) realloc(prtcl_sendnegx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_id_sendposx = (int *) realloc(prtcl_id_sendposx,sizeof(int)*nprtcl_max_comm);
	  prtcl_id_sendnegx = (int *) realloc(prtcl_id_sendnegx,sizeof(int)*nprtcl_max_comm);
	}
	prtcl[i].xprtcl = prtcl[i].xprtcl + (REAL) NX;
	prtcl_sendnegx[sendnegx] = prtcl[i];
	prtcl_id_sendnegx[sendnegx] = prtcl_ids[i];
	sendnegx = sendnegx+1;
      }else {
	prtcl[keephere] = prtcl[i];
	prtcl_ids[keephere] = prtcl_ids[i];
	keephere = keephere+1;
      }
    }

    procpx=(mex+1+NPX)%NPX; procmx=(mex-1+NPX)%NPX;
    MPI_Sendrecv(&sendposx,1,MPI_INT,procpx,10,&recvnegx,1,MPI_INT,procmx,10,MPI_COMM_ALONG_X,&status1);
    MPI_Sendrecv(&sendnegx,1,MPI_INT,procmx,11,&recvposx,1,MPI_INT,procpx,11,MPI_COMM_ALONG_X,&status2);

    int flag = 0;
    if ((recvposx > nprtcl_max_comm)||(recvnegx > nprtcl_max_comm)){
      flag = 1;
      if (recvposx >= recvnegx){
	nprtcl_max_comm = recvposx + 1;
      } else {
	nprtcl_max_comm = recvnegx + 1;
      }
    }

    if (flag>0){
      prtcl_sendposx = (prtcl_type *) realloc(prtcl_sendposx,sizeof(prtcl_type) * nprtcl_max_comm);
      prtcl_sendnegx = (prtcl_type *) realloc(prtcl_sendnegx,sizeof(prtcl_type) * nprtcl_max_comm);
      prtcl_id_sendposx = (int *) realloc(prtcl_id_sendposx,sizeof(int)*nprtcl_max_comm);
      prtcl_id_sendnegx = (int *) realloc(prtcl_id_sendnegx,sizeof(int)*nprtcl_max_comm);
    }
  
    int np_loc_temp = keephere + recvnegx + recvposx;
    if (np_loc_temp > nprtcl_max){
      nprtcl_max = np_loc_temp + 1;
      prtcl = *prtcl1 = (prtcl_type *) realloc(prtcl,sizeof(prtcl_type) * nprtcl_max);
      prtcl_ids = *prtcl_ids1 = (int *) realloc(prtcl_ids,sizeof(int)*nprtcl_max);
    }
    nprtcl_local = np_loc_temp;
    if ((prtcl==NULL)||(prtcl_ids==NULL)){
      printf("Problem Reallocing");
    }


    MPI_Sendrecv(prtcl_sendposx,sendposx,MPI_Prtcltype,procpx,12,prtcl+keephere,recvnegx,MPI_Prtcltype,procmx,12,MPI_COMM_ALONG_X,&status1);
    MPI_Sendrecv(prtcl_sendnegx,sendnegx,MPI_Prtcltype,procmx,13,prtcl+keephere+recvnegx,recvposx,MPI_Prtcltype,procpx,13,MPI_COMM_ALONG_X,&status2);
    MPI_Sendrecv(prtcl_id_sendposx,sendposx,MPI_INT,procpx,14,prtcl_ids+keephere,recvnegx,MPI_INT,procmx,14,MPI_COMM_ALONG_X,&status1);
    MPI_Sendrecv(prtcl_id_sendnegx,sendnegx,MPI_INT,procmx,15,prtcl_ids+keephere+recvnegx,recvposx,MPI_INT,procpx,15,MPI_COMM_ALONG_X,&status2);

    sendposx = 0;
    sendnegx = 0;
    keephere = 0;
    recvposx = 0;
    recvnegx = 0;

    for (int i = 0;i<nprtcl_local;i++){
      if ((prtcl[i].yprtcl >= (REAL)NY)&&(mey<NPY-1)){
	if (sendposx >= nprtcl_max_comm){
	  nprtcl_max_comm += 1;
	  prtcl_sendposx = (prtcl_type *) realloc(prtcl_sendposx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_sendnegx = (prtcl_type *) realloc(prtcl_sendnegx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_id_sendposx = (int *) realloc(prtcl_id_sendposx,sizeof(int)*nprtcl_max_comm);
	  prtcl_id_sendnegx = (int *) realloc(prtcl_id_sendnegx,sizeof(int)*nprtcl_max_comm);
	}
	prtcl[i].yprtcl = prtcl[i].yprtcl - (REAL) NY;
	prtcl_sendposx[sendposx] = prtcl[i];
	prtcl_id_sendposx[sendposx] = prtcl_ids[i];
	sendposx = sendposx+1;
      }  else if ((prtcl[i].yprtcl <= 0.0)&&(mey>0)){
	if (sendnegx >= nprtcl_max_comm){
	  nprtcl_max_comm += 1;
	  prtcl_sendposx = (prtcl_type *) realloc(prtcl_sendposx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_sendnegx = (prtcl_type *) realloc(prtcl_sendnegx,sizeof(prtcl_type) * nprtcl_max_comm);
	  prtcl_id_sendposx = (int *) realloc(prtcl_id_sendposx,sizeof(int)*nprtcl_max_comm);
	  prtcl_id_sendnegx = (int *) realloc(prtcl_id_sendnegx,sizeof(int)*nprtcl_max_comm);
	}
	prtcl[i].yprtcl = prtcl[i].yprtcl + (REAL) NY;
	prtcl_sendnegx[sendnegx] = prtcl[i];
	prtcl_id_sendnegx[sendnegx] = prtcl_ids[i];
	sendnegx = sendnegx+1;
      }else if ((prtcl[i].yprtcl<1.0)&&(mey==0)){
	prtcl[i].yprtcl = prtcl[i].yprtcl + 2.0*(1.0 - prtcl[i].yprtcl);
  	prtcl[i].vprtcl = - prtcl[i].vprtcl;
	prtcl[keephere] = prtcl[i];
	prtcl_ids[keephere] = prtcl_ids[i];
	keephere = keephere+1;
      } else if ((prtcl[i].yprtcl >= (REAL)NY)&&(mey==NPY-1)) {
	prtcl[i].yprtcl = prtcl[i].yprtcl + 2.0*((REAL)NY - prtcl[i].yprtcl);
  	prtcl[i].vprtcl = - prtcl[i].vprtcl;
	prtcl[keephere] = prtcl[i];
	prtcl_ids[keephere] = prtcl_ids[i];
	keephere = keephere+1;
      } else {
	prtcl[keephere] = prtcl[i];
	prtcl_ids[keephere] = prtcl_ids[i];
	keephere = keephere+1;
      }
    }
    procpx=(mey+1+NPY)%NPY; procmx=(mey-1+NPY)%NPY;
    MPI_Sendrecv(&sendposx,1,MPI_INT,procpx,14,&recvnegx,1,MPI_INT,procmx,14,MPI_COMM_ALONG_Y,&status1);
    MPI_Sendrecv(&sendnegx,1,MPI_INT,procmx,15,&recvposx,1,MPI_INT,procpx,15,MPI_COMM_ALONG_Y,&status2);
    flag = 0;
    if ((recvposx > nprtcl_max_comm)||(recvnegx > nprtcl_max_comm)){
      flag = 1;
      if (recvposx >= recvnegx){
	nprtcl_max_comm = recvposx + 1;
      } else {
	nprtcl_max_comm = recvnegx + 1;
      }
    }

    if (flag>0){
      prtcl_sendposx = (prtcl_type *) realloc(prtcl_sendposx,sizeof(prtcl_type) * nprtcl_max_comm);
      prtcl_sendnegx = (prtcl_type *) realloc(prtcl_sendnegx,sizeof(prtcl_type) * nprtcl_max_comm);
      prtcl_id_sendposx = (int *) realloc(prtcl_id_sendposx,sizeof(int)*nprtcl_max_comm);
      prtcl_id_sendnegx = (int *) realloc(prtcl_id_sendnegx,sizeof(int)*nprtcl_max_comm);
    }
  
  
    np_loc_temp = keephere + recvnegx + recvposx;
    if (np_loc_temp > nprtcl_max){
      nprtcl_max = np_loc_temp + 1;
      prtcl = *prtcl1 = (prtcl_type *) realloc(prtcl,sizeof(prtcl_type) * nprtcl_max);
      prtcl_ids = *prtcl_ids1 = (int *) realloc(prtcl_ids,sizeof(int) * nprtcl_max);
    }
    nprtcl_local = np_loc_temp;
    if (prtcl==NULL){
      printf("Problem Reallocing");
    }

    MPI_Sendrecv(prtcl_sendposx,sendposx,MPI_Prtcltype,procpx,16,prtcl+keephere,recvnegx,MPI_Prtcltype,procmx,16,MPI_COMM_ALONG_Y,&status1);
    MPI_Sendrecv(prtcl_sendnegx,sendnegx,MPI_Prtcltype,procmx,17,prtcl+keephere+recvnegx,recvposx,MPI_Prtcltype,procpx,17,MPI_COMM_ALONG_Y,&status2);
    MPI_Sendrecv(prtcl_id_sendposx,sendposx,MPI_INT,procpx,18,prtcl_ids+keephere,recvnegx,MPI_INT,procmx,18,MPI_COMM_ALONG_Y,&status1);
    MPI_Sendrecv(prtcl_id_sendnegx,sendnegx,MPI_INT,procmx,19,prtcl_ids+keephere+recvnegx,recvposx,MPI_INT,procpx,19,MPI_COMM_ALONG_Y,&status2);
  }
