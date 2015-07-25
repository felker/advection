#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "visit_writer.h"

//select solver options
#define X2_PERIODIC 1 //flag to wrap phi coordinate and automatically mesh entire circle
#define AUTO_TIMESTEP 1 //flag for automatically setting dt such that max{cfl_array} = CFL
#define CFL 0.8 //if set to 1.0, 1D constant advection along grid is exact. 
#define OUTPUT_INTERVAL 1 //how many timesteps to dump simulation data. 0 for only last step, 1 for every step
#undef SECOND_ORDER //flag to turn on van Leer flux limiting

//select problem (mutually exclusive)
#undef RADIAL_INWARD
#undef SEMI_CLOCK
#undef HORIZONTAL
#undef GAUSS_RADIAL
#define GAUSS_CLOCK

#ifdef SEMI_CLOCK
#define X2_PERIODIC 0
#endif

double stream_function(double x, double y);
double X_physical(double, double);
double Y_physical(double, double);
void vector_coordinate_to_physical(double vr, double vphi, double phi, double *vx, double *vy);
void vector_physical_to_coordinate(double vx, double vy, double phi, double *vr, double *vphi);
void velocity_physical(double x_b,double y_b,double *vx,double *vy);
void velocity_coordinate(double r_b,double phi_b,double *vr,double *vphi);
double initial_condition(double x, double y);
double bc_x1i(double x, double y);
double bc_x1f(double x, double y,double t);
double bc_x2i(double x, double y);
double bc_x2f(double x, double y);
double flux_PLM(double ds,double *imu);
float find_max(float a[], int n);
float find_min(float a[], int n); 
float sum(float a[], int n);
double gaussian(double x_0, double y_0,double x,double y);

int main(int argc, char **argv){
  int i,j,k,n; 
  int nsteps=150;
  double dt =0.01;
  
  /* Computational (2D polar) grid coordinates */
  int nx1 = 200;
  int nx2 = 200;

  //number of ghost cells on both sides of each dimension
  //only need 1 for piecewise constant method
  //need 2 for piecewise linear reconstruction
  int num_ghost = 2;
  int nx1_r = nx1;
  int nx2_r = nx2;
  nx1 += 2*num_ghost; 
  nx2 += 2*num_ghost; 

  /*non-ghost indices */
  int is = num_ghost;
  int ie= is+nx1_r; 
  int js = num_ghost;
  int je = js+nx2_r; 

  //convention: these ranges refer to the non-ghost cells
  //however, the ghost cells have real coordinate interpretations
  //this means that we must be sure that the number of ghost cells makes sense with the range of coordinates
  //this is a big problem if the phi polar coordinate runs the whole range [0,2pi) for example

  //further, all mesh structures (zonal and nodal) will be nx1 x nx2, although we may not fill the ghost entries with anything meaningful
  //this is to standardize the loop indexing from is:ie

  /*another convention: when the phi coordinate is periodic, do we repeat the boundary mesh points? yes for now */

  double lx1 = 2.0;//these values are inclusive [x1_i, x1_f]
  double lx2 = M_PI;
  double x1_i = 0.5;
  double x2_i = 0.0;

  double dx2 = lx2/(nx2_r-1);   
  double x2_f = x2_i + lx2;

  if(X2_PERIODIC){
    dx2 = 2*M_PI/(nx2_r);   
    lx2 = dx2*(nx2_r-1); 
    x2_i = 0.0;
  }

  double dx1 = lx1/(nx1_r-1);   
  double x1_f = x1_i + lx1;

  printf("dx1=%lf dx2=%lf \n",dx1,dx2); 
  
  /*Cell centered (zonal) values of computational coordinate position */
  double *x1 = (double *) malloc(sizeof(double)*nx1); 
  double *x2 = (double *) malloc(sizeof(double)*nx2);
  x1[is] = x1_i;
  x2[js] = x2_i;
  for(i=is+1; i<ie; i++){ 
    x1[i] = x1[i-1] + dx1;
    //    printf("%lf\n",x1[i]);
  }
  for(i=js+1; i<je; i++){
    x2[i] = x2[i-1] + dx2;
    //printf("%lf\n",x2[i]);
  } 

  /*Mesh edge (nodal) values of computational coordinate position */
  double *x1_b = (double *) malloc(sizeof(double)*(nx1+1)); 
  double *x2_b = (double *) malloc(sizeof(double)*(nx2+1));
  x1_b[is] = x1_i - dx1/2;
  x2_b[js] = x2_i - dx2/2;
  for(i=is+1; i<=ie; i++){ 
    x1_b[i] = x1_b[i-1] + dx1;
    //    printf("%lf\n",x1_b[i]);
  }
  for(i=js+1; i<=je; i++){
    x2_b[i] = x2_b[i-1] + dx2;
    //printf("%lf\n",x2_b[i]);
  } 

  /*Cell centered (zonal) values of physical coordinate position */
  //These must be 2D arrays since coordinate transformation is not diagonal
  //indexed by x1 rows and x2 columns
  double **x = (double **) malloc(sizeof(double *)*nx1); 
  double **y = (double **) malloc(sizeof(double *)*nx1);
  double *dataX = (double *) malloc(sizeof(double)*nx1*nx2);
  double *dataY = (double *) malloc(sizeof(double)*nx1*nx2);
  for(i=0; i<nx1; i++){
    x[i] = &(dataX[nx2*i]);
    y[i] = &(dataY[nx2*i]);
  }

  /*precompute coordinate mappings */
  for(i=is; i<ie; i++){
    for(j=js; j<je; j++){
      x[i][j] = X_physical(x1[i],x2[j]); 
      y[i][j] = Y_physical(x1[i],x2[j]); 
    }
  }
  /*Mesh edge (nodal) values of physical coordinate position */
  double **x_b = (double **) malloc(sizeof(double *)*(nx1+1)); 
  double **y_b = (double **) malloc(sizeof(double *)*(nx1+1));
  double *dataXb = (double *) malloc(sizeof(double)*(nx1+1)*(nx2+1));
  double *dataYb = (double *) malloc(sizeof(double)*(nx1+1)*(nx2+1));
  for(i=0; i<=nx1; i++){
    x_b[i] = &(dataXb[(nx2+1)*i]);
    y_b[i] = &(dataYb[(nx2+1)*i]);
  }
  for(i=is; i<=ie; i++){
    for(j=js; j<=je; j++){
      x_b[i][j] = X_physical(x1_b[i],x2_b[j]); 
      y_b[i][j] = Y_physical(x1_b[i],x2_b[j]); 
    }
  }
  /*Edge normal vectors in Cartesian coordinates */
  // point radially outward and +\phi
  

  /*Coordinate cell capacity */
  double **kappa;
  double *datakappa;
  kappa = (double **) malloc(sizeof(double *)*nx1);
  datakappa = (double *) malloc(sizeof(double)*nx1*nx2);
  for(i=0; i<nx1; i++){
    kappa[i] = &(datakappa[nx2*i]);
  }
  
  for(i=is; i<ie; i++){
    for(j=js; j<je; j++){
      kappa[i][j] = x1[i]*dx1*dx2/(dx1*dx2); // C_ij/(dx1*dx2)
      //capacity in ghost cells
      kappa[ie][j] = (x1[ie-1]+dx1)*dx1*dx2/(dx1*dx2); 
    }
    kappa[i][je] = (x1[i]+dx1)*dx1*dx2/(dx1*dx2); 
  }
 
  /*Average normal edge velocities */
  //now we move from cell centered quantities to edge quantities 
  //the convention in this code is that index i refers to i-1/2 edge
  double **U, **V;
  double *dataU, *dataV; 
  U = (double **) malloc(sizeof(double *)*nx1);
  V = (double **) malloc(sizeof(double *)*nx1);
  dataU = (double *) malloc(sizeof(double)*nx1*nx2);
  dataV = (double *) malloc(sizeof(double)*nx1*nx2);
  for(i=0; i<nx1; i++){
    U[i] = &(dataU[nx2*i]);
    V[i] = &(dataV[nx2*i]);
  }
  
  /*Option #1 for computing edge velocities: path integral of Cartesian stream function */
  /* Stream function */ //probably dont need this array
  //this variable is cell centered stream 
  /* double **stream;
  double *datastream;
  stream = (double **) malloc(sizeof(double *)*nx1);
  datastream = (double *) malloc(sizeof(double)*nx1*nx2);
  for(i=0; i<nx1; i++){
    stream[i] = &(datastream[nx2*i]);
  }
  for(i=is; i<ie; i++){
    for(j=js; j<je; j++){
      stream[i][j] = stream_function(x[i][j],y[i][j]);
    }
  }
 
   for(i=is; i<ie; i++){
    for(j=js; j<je; j++){ //go an additional step to capture nx2_r+1/2
      //average the corner stream functions
      U[i][j]= (stream_function(X_physical(x1[i]-dx1/2,x2[j]+dx2/2),Y_physical(x1[i]-dx1/2,x2[j]+dx2/2)) - 
		stream_function(X_physical(x1[i]-dx1/2,x2[j]-dx2/2),Y_physical(x1[i]-dx1/2,x2[j]-dx2/2)))/(dx2);
      V[i][j]= -(stream_function(X_physical(x1[i]+dx1/2,x2[j]-dx2/2),Y_physical(x1[i]+dx1/2,x2[j]-dx2/2)) - 
      	 stream_function(X_physical(x1[i]-dx1/2,x2[j]-dx2/2),Y_physical(x1[i]-dx1/2,x2[j]-dx2/2)))/dx1;
      if(j==10)
	printf("i=%d, u,v=%lf,%lf\n",i,U[i][j],V[i][j]); 
    }
    j=je;
    U[i][j]= (stream_function(X_physical(x1[i]-dx1/2,x2[j-1]+3*dx2/2),Y_physical(x1[i]-dx1/2,x2[j-1]+3*dx2/2)) - 
	      stream_function(X_physical(x1[i]-dx1/2,x2[j-1]+dx2/2),Y_physical(x1[i]-dx1/2,x2[j-1]+dx2/2)))/(dx2);
    V[i][j]= -(stream_function(X_physical(x1[i]+dx1/2,x2[j-1]+dx2/2),Y_physical(x1[i]+dx1/2,x2[j-1]+dx2/2)) - 
	       stream_function(X_physical(x1[i]-dx1/2,x2[j-1]+dx2/2),Y_physical(x1[i]-dx1/2,x2[j-1]+dx2/2)))/dx1;
  }
  i=ie;
  for(j=js; j<je; j++){
    U[i][j]= (stream_function(X_physical(x1[i-1]+dx1/2,x2[j]+dx2/2),Y_physical(x1[i-1]+dx1/2,x2[j]+dx2/2)) - 
	      stream_function(X_physical(x1[i-1]+dx1/2,x2[j]-dx2/2),Y_physical(x1[i-1]+dx1/2,x2[j]-dx2/2)))/(dx2);
    V[i][j]= -(stream_function(X_physical(x1[i-1]+3*dx1/2,x2[j]-dx2/2),Y_physical(x1[i-1]+3*dx1/2,x2[j]-dx2/2)) - 
	      stream_function(X_physical(x1[i-1]+dx1/2,x2[j]-dx2/2),Y_physical(x1[i-1]+dx1/2,x2[j]-dx2/2)))/dx1;
  } 
  */
  double ux,vy,temp; 
  for(i=is; i<=ie; i++){
    for(j=js; j<=je; j++){ 
      /*Option #2: directly specify velocity in Cartesian coordinate basis as a function of cartesian position*/   
      
      //radial face i-1/2
      velocity_physical(X_physical(x1_b[i],x2_b[j]+dx2/2),Y_physical(x1_b[i],x2_b[j]+dx2/2),&ux,&vy);
      // Average normal edge velocity: just transform face center velocity to local orthonormal basis? 
      vector_physical_to_coordinate(ux,vy,x2_b[j]+dx2/2,&U[i][j],&temp); 
      
#ifdef HORIZONTAL
      //EXACT SOLUTION FOR EDGE VELOCITY FOR HORIZONTAL FLOW
      //      printf("U_before = %lf\n",U[i][j]);
      U[i][j] = (sin(x2_b[j]+dx2) - sin(x2_b[j]))/(dx2); 
      //      printf("U_after = %lf\n",U[i][j]);
#endif

      //phi face j-1/2
      velocity_physical(X_physical(x1_b[i]+dx1/2,x2_b[j]),Y_physical(x1_b[i]+dx1/2,x2_b[j]),&ux,&vy);
      vector_physical_to_coordinate(ux,vy,x2_b[j],&temp,&V[i][j]); 
#if defined(SEMI_CLOCK) || defined(GAUSS_CLOCK)
      velocity_coordinate(x1_b[i],x2_b[j],&U[i][j],&V[i][j]);
#endif
      //      printf("U,V = %lf,%lf\n",U[i][j],V[i][j]);
    }
  } 

  /* check normalization of velocities  */
  /* these are naturally not normalized because the cell has finite volume so the edge velocities arent the velocity of the same point */
  double norm; 
  double max_deviation =0.0; 
  for(i=is; i<=ie; i++){
    for(j=js; j<=je; j++){ 
      norm = sqrt(U[i][j]*U[i][j] + V[i][j]*V[i][j]); 
      if (fabs(norm-1.0) > max_deviation)
	max_deviation = fabs(norm-1.0); 
      //  printf("%0.12lf\n",norm);
    }
  }
  //  printf("maximum deviation from 1.0 = %0.12lf\n",max_deviation); 



  /*Option #3: specify velocity in polar coordinates */


  /*Check CFL condition, reset timestep */
  float **cfl_array;
  float *dataCFL;
  cfl_array = (float **) malloc(sizeof(float *)*nx1);
  dataCFL = (float *) malloc(sizeof(float)*nx1*nx2);
  for(i=0; i<nx1; i++){
    cfl_array[i] = &(dataCFL[nx2*i]);
  }
  for(i=1; i<nx1; i++){
    for(j=1; j<nx2; j++){ //based on edge velocities or cell centered u,v?
      if (i >=is && i< ie && j >=js && j <je)
	cfl_array[i][j] = fabs(U[i][j])*dt/dx1 + fabs(V[i][j])*dt/(x1_b[i]*dx2); //use boundary radius
      else
	cfl_array[i][j] =0.0; 
    }
  }
  //find maximum CFL value in domain
  float  max_cfl = find_max(dataCFL,nx1*nx2); 
  printf("Largest CFL number = %lf\n",max_cfl); 
  if (max_cfl > CFL || AUTO_TIMESTEP){//reset timestep if needed
    dt = CFL*dt/max_cfl; 
    for(i=1; i<nx1; i++){
      for(j=1; j<nx2; j++){ 
	if (i >=is && i< ie && j >=js && j <je)
	  cfl_array[i][j] = fabs(U[i][j])*dt/dx1 + fabs(V[i][j])*dt/(x1_b[i]*dx2);
	else
	  cfl_array[i][j] =0.0; 
      }
    } 
  }
  max_cfl = find_max(dataCFL,nx1*nx2); 
  printf("Largest CFL number = %lf\n",max_cfl); 
  
#ifdef SEMI_CLOCK
  nsteps = M_PI/dt;
  printf("nsteps = %d, dt = %lf, t_final = %lf\n",nsteps,dt,nsteps*dt);
#endif
  
#ifdef GAUSS_CLOCK
  nsteps = 2*M_PI/dt;
  //turn dt down until mod(2*MPI,nsteps) = 0
  double remainder = 2*M_PI -nsteps*dt;
  double extra = remainder/nsteps; 
  double dt_new = dt +extra; 
  printf("nsteps = %d, dt = %lf, dt_new =%lf, remainder = %lf, t_final = %lf t_final_new =%lf\n",nsteps,dt,dt_new,remainder,nsteps*dt,nsteps*dt_new);
  max_cfl *= dt_new/dt; 
  dt = dt_new; 
  printf("Largest CFL number = %lf\n",max_cfl); 
#endif

  /*Conserved variable on the computational coordinate mesh*/
  double **Q;
  double *dataQ;
  //make contiguous multiarray                                                                                                
  Q = (double **) malloc(sizeof(double *)*nx1);
  dataQ = (double *) malloc(sizeof(double)*nx1*nx2);
  for(i=0; i<nx1; i++){
    Q[i] = &(dataQ[nx2*i]);
  }
  /*Initial condition */
  //specified as a function of cartesian physical cooridnates, as is the stream function
  for(i=is; i<ie; i++){
    for(j=js; j<je; j++){
      Q[i][j] = initial_condition(x[i][j],y[i][j]); 
    }
  }

  //net fluctuations/fluxes
  double U_plus,U_minus,V_plus,V_minus;
  double **net_fluctuation;
  double *dataFlux;
  net_fluctuation = (double **) malloc(sizeof(double *)*nx1);
  dataFlux = (double *) malloc(sizeof(double)*nx1*nx2);
  for(i=0; i<nx1; i++){
    net_fluctuation[i] = &(dataFlux[nx2*i]);
  }

  /* Using Visit VTK writer */
  char filename[20];
  int dims[] = {nx1_r+1, nx2_r+1, 1}; //dont output ghost cells. //nodal variables have extra edge point
  int nvars = 2;
  int vardims[] = {1, 3}; //Q is a scalar, velocity is a 3-vector 
  int centering[] = {0, 1}; // Q is cell centered, velocity is defined at edges
  const char *varnames[] = {"Q", "edge_velocity"};
  /* Curvilinear mesh points stored x0,y0,z0,x1,y1,z1,...*/
  //An array of size nI*nJ*nK*3 . These points are nodal, not zonal... unfortunatley
  float *pts = (float *) malloc(sizeof(float)*(nx1_r+1)*(nx2_r+1)*3); 
  //The array should be layed out as (pt(i=0,j=0,k=0), pt(i=1,j=0,k=0), ...
  //pt(i=nI-1,j=0,k=0), pt(i=0,j=1,k=0), ...).
  int index=0; 
  for(k=0; k<1; k++){
    for(j=js; j<=je; j++){
      for(i=is; i<=ie; i++){
	pts[index] = x_b[i][j];
	pts[++index] = y_b[i][j];
	pts[++index] = 0.0;
	index++;
      }
    }
  }

  /* pack U and V into a vector */
  float *edge_vel = (float *) malloc(sizeof(float)*(nx1_r+1)*(nx2_r+1)*3); //An array of size nI*nJ*nK*3 
  index=0; 
  for(k=0; k<1; k++){
    for(j=js; j<=je; j++){
      for(i=is; i<=ie; i++){
	vector_coordinate_to_physical(U[i][j],V[i][j], x2[j], &ux,&vy);
	//edge_vel[index] = U[i][j];
	//edge_vel[++index] = V[i][j]; 
	edge_vel[index] = ux;
	edge_vel[++index] = vy;
	edge_vel[++index] = 0.0;
	index++;
      }
    }
  } 

  //  vars       An array of variables.  The size of vars should be nvars.
  //                 The size of vars[i] should be npts*vardim[i].
  float *realQ; 
  realQ =(float *) malloc(sizeof(double)*nx1_r*nx2_r);
  float *vars[] = {(float *) realQ, (float *)edge_vel};

  /*-----------------------*/
  /* Main timestepping loop */
  /*-----------------------*/
  for (n=0; n<nsteps; n++){
    /*Boundary conditions */
    //bcs are specified along a computational coord direction, but are a function of the physical coordinate of adjacent "real cells"
    for(k=0;k<num_ghost; k++){
      for (j=js; j<je; j++){
	Q[k][j] = bc_x1i(x[is][j],y[is][j]);
	Q[nx1-1-k][j] = bc_x1f(x[ie-1][j],y[ie-1][j],(n+1)*dt);
      }
      for (i=is; i<ie; i++){
	if(X2_PERIODIC){
	  Q[i][k] = Q[i][je-1-k];
	  Q[i][nx2-1-k] = Q[i][js+k];
	}
	else{
	Q[i][k] = bc_x2i(x[i][js],y[i][js]);
	Q[i][nx2-1-k] = bc_x2f(x[i][je-1],y[i][je-1]);
	}
      }
    }

    double flux_limiter =0.0; 
    double *qmu = (double *) malloc(sizeof(double)*3); //manually copy array for computing slope limiters
    /* Donor cell upwinding */
    for (i=is; i<ie; i++){
      for (j=js; j<je; j++){
	/* First coordinate */
	U_plus = fmax(U[i][j],0.0); // max{U_{i-1/2,j},0.0} LHS boundary
	U_minus = fmin(U[i+1][j],0.0); // min{U_{i+1/2,j},0.0} RHS boundary
	/*Fluctuations: A^+ \Delta Q_{i-1/2,j} + A^- \Delta Q_{i+1/2,j} */
	//	net_fluctuation[i][j] = dt/(kappa[i][j]*dx1)*(U_plus*(Q[i][j] - Q[i-1][j]) + U_minus*(Q[i+1][j] - Q[i][j]));
	/* First order fluxes: F_i+1/2 - F_i-1/2 */
	net_fluctuation[i][j] = dt/(kappa[i][j]*dx1)*(x1_b[i+1]*(fmax(U[i+1][j],0.0)*Q[i][j] + U_minus*Q[i+1][j])-x1_b[i]*(U_plus*Q[i-1][j] + fmin(U[i][j],0.0)*Q[i][j]));
#ifdef SECOND_ORDER
	/* Second order fluxes */
	if (U[i+1][j] > 0.0){ //middle element is always the upwind element
	  qmu[0] = Q[i-1][j];
	  qmu[1] = Q[i][j];  
	  qmu[2] = Q[i+1][j];  
	  flux_limiter= flux_PLM(dx1,qmu);
	}
	else{
	  qmu[0] = Q[i+2][j];
	  qmu[1] = Q[i+1][j];  
	  qmu[2] = Q[i][j];  
	  flux_limiter= flux_PLM(dx1,qmu);
	  /*	  if (flux_limiter != 0.0){
	    printf("i,j: %d,%d F_{i+1/2} flux limiter: %lf\n",i,j,flux_limiter);
	    printf("Q0 = %lf Q1= %lf Q2 = %lf\n",qmu[0],qmu[1],qmu[2]); } */
	}
	//F^H_{i+1/2,j}
	net_fluctuation[i][j] -= dt/(kappa[i][j]*dx1)*(x1_b[i+1]*(1-dt*fabs(U[i+1][j])/(dx1))*fabs(U[i+1][j])*flux_limiter/2);
	//	net_fluctuation[i][j] -= dt/(kappa[i][j]*dx1)*((kappa[i][j]/x1_b[i+1]-dt*fabs(U[i+1][j])/(x1_b[i+1]*dx1))*fabs(U[i+1][j])*flux_limiter/2);
	//	net_fluctuation[i][j] += dt/(kappa[i][j]*dx1)*(x1_b[i+1]*(kappa[i][j]/x1_b[i+1]-dt*fabs(U[i+1][j])/(x1_b[i+1]*dx1))*fabs(U[i+1][j])*flux_limiter/2);
	if (U[i][j] > 0.0){
	  qmu[0] = Q[i-2][j];  //points to the two preceeding bins; 
	  qmu[1] = Q[i-1][j];  
	  qmu[2] = Q[i][j];  
	  flux_limiter= flux_PLM(dx1,qmu);
	}
	else{
	  qmu[0] = Q[i+1][j]; //centered around current bin
	  qmu[1] = Q[i][j];  
	  qmu[2] = Q[i-1][j];  
	  flux_limiter= flux_PLM(dx1,qmu);
	}
	//F^H_{i-1/2,j}
	net_fluctuation[i][j] += dt/(kappa[i][j]*dx1)*(x1_b[i]*(1-dt*fabs(U[i][j])/(dx1))*fabs(U[i][j])*flux_limiter/2);
	//	net_fluctuation[i][j] += dt/(kappa[i][j]*dx1)*((kappa[i-1][j]/x1_b[i]-dt*fabs(U[i][j])/(x1_b[i]*dx1))*fabs(U[i][j])*flux_limiter/2);
	//net_fluctuation[i][j] -= dt/(kappa[i][j]*dx1)*(x1_b[i]*(kappa[i-1][j]/x1_b[i]-dt*fabs(U[i][j])/(x1_b[i]*dx1))*fabs(U[i][j])*flux_limiter/2);
#endif
	/* Second coordinate */
	V_plus = fmax(V[i][j],0.0); // max{V_{i,j-1/2},0.0} LHS boundary
	V_minus = fmin(V[i][j+1],0.0); // min{V_{i,j+1/2},0.0} RHS boundary
	/*Fluctuations: B^+ \Delta Q_{i,j-1/2} + B^- \Delta Q_{i,j+1/2} */
	//net_fluctuation[i][j] += dt/(kappa[i][j]*dx2)*(V_plus*(Q[i][j] - Q[i][j-1]) + V_minus*(Q[i][j+1] - Q[i][j]));	
	/* Fluxes: G_i,j+1/2 - G_i,j-1/2 */
	net_fluctuation[i][j] += dt/(kappa[i][j]*dx2)*((fmax(V[i][j+1],0.0)*Q[i][j] + V_minus*Q[i][j+1])-(V_plus*Q[i][j-1] + fmin(V[i][j],0.0)*Q[i][j]));
#ifdef SECOND_ORDER
	/* Second order fluxes */
	if (V[i][j+1] > 0.0){
	  qmu[0] = Q[i][j-1];  //points to the two preceeding bins; 
	  qmu[1] = Q[i][j];  
	  qmu[2] = Q[i][j+1];  
	  flux_limiter= flux_PLM(dx2,qmu);
	}
	else{
	  qmu[0] = Q[i][j+2]; //centered around current bin
	  qmu[1] = Q[i][j+1];  
	  qmu[2] = Q[i][j];  
	  flux_limiter= flux_PLM(dx2,qmu);
	  /*	  if (flux_limiter != 0.0){
	    printf("i,j: %d,%d F_{i+1/2} flux limiter: %lf\n",i,j,flux_limiter);
	    printf("Q0 = %lf Q1= %lf Q2 = %lf\n",qmu[0],qmu[1],qmu[2]); } */
	}
	//G^H_{i,j+1/2}
	//net_fluctuation[i][j] -= dt/(kappa[i][j]*dx2)*((1-dt*fabs(V[i][j+1])/(kappa[i][j]*dx2))*fabs(V[i][j+1])*flux_limiter/2);
	net_fluctuation[i][j] -= dt/(kappa[i][j]*dx2)*((1-dt*fabs(V[i][j+1])/(dx2))*fabs(V[i][j+1])*flux_limiter/2);
	if (V[i][j] > 0.0){
	  qmu[0] = Q[i][j-2];  //points to the two preceeding bins; 
	  qmu[1] = Q[i][j-1];  
	  qmu[2] = Q[i][j];  
	  flux_limiter = flux_PLM(dx2,qmu);
	}
	else{
	  qmu[0] = Q[i][j+1]; //centered around current bin
	  qmu[1] = Q[i][j];  
	  qmu[2] = Q[i][j-1];  
	  flux_limiter= flux_PLM(dx2,qmu);
	}
	//G^H_{i,j-1/2}
	//net_fluctuation[i][j] -= dt/(kappa[i][j]*dx2)*((1-dt*fabs(V[i][j])/(kappa[i][j]*dx2))*fabs(V[i][j])*flux_limiter/2);
	net_fluctuation[i][j] += dt/(kappa[i][j]*dx2)*((1-dt*fabs(V[i][j])/(dx2))*fabs(V[i][j])*flux_limiter/2);
#endif
      }
    }

    /*Apply fluctuations */
    for (i=is; i<ie; i++)
      for (j=js; j<je; j++){
	Q[i][j] -= net_fluctuation[i][j];
      }
    
    /*Source terms */
    for (i=is; i<ie; i++)
      for (j=js; j<je; j++){
      }

    /*Output */
    //for now, explicitly copy subarray corresponding to real zonal info:
    index=0; 
    for (k=0; k<1; k++){
      for (j=js; j<je; j++){
	for (i=is; i<ie; i++){
	  //index =(j-num_ghost)*nx2_r + (i-num_ghost); 
	  realQ[index] = (float) Q[i][j];//*kappa[i][j]; //\bar{q}=qk density in computational space
	  index++;
	}
      }
    }
    //debug only horizontal flow
    /*    if (find_max(realQ,nx1_r*nx2_r) > 1.1){
      printf("Q greater than 1.0!\n"); 
      for (i=1;i<nx1; i++){
	for (j=1; j<nx2; j++){
	  if (Q[i][j] > 1.0){
	    printf("i=%d j=%d Q[i][j] = %0.10lf\n",i,j,Q[i][j]); 
	    return(0); 
	  }
	}
      }
      }*/

    sprintf(filename,"advect-%.3d.vtk",n); 
    if(!OUTPUT_INTERVAL){
      if (n==nsteps-1) //for only the final result
	write_curvilinear_mesh(filename,3,dims, pts, nvars,vardims, centering, varnames, vars);}
    else{
      if (!(n%OUTPUT_INTERVAL)) //HAVENT CHECKED THIS
	write_curvilinear_mesh(filename,3,dims, pts, nvars,vardims, centering, varnames, vars);}
      
      printf("step: %d time: %lf max{Q} = %0.7lf min{Q} = %0.7lf sum{Q} = %0.7lf \n",
	     n+1,(n+1)*dt,find_max(realQ,nx1_r*nx2_r),find_min(realQ,nx1_r*nx2_r),sum(realQ,nx1_r*nx2_r));
  }
  return(0); 
}

/* Map to physical (cartesian) coordinates */
double X_physical(double x1, double x2){
  double x_cartesian = x1*cos(x2); 
  return(x_cartesian); 
}
double Y_physical(double x1, double x2){
  double y_cartesian = x1*sin(x2); 
  return(y_cartesian); 
}

double initial_condition(double x, double y){
#if defined(GAUSS_CLOCK) || defined(GAUSS_RADIAL)
  //if ((x<-1.0) && (x>=-1.5) && (y>1.0) && (y<1.5)){
  //  return(1.0);
  //}

  //create a smooth gaussian in 2D

  return(gaussian(-1,1.5,x,y));
#endif 

  return(0.0);
}

//for 2D polar coordinates:
//bc at innermost radius
double bc_x1i(double x, double y){
  //  if ((y > 0.2) && (y < 0.3 ))
  //return(1.0); 
  return(0.0);
}
//bc at outermost radius
double bc_x1f(double x, double y, double t){
#ifdef HORIZONTAL
  if ((x<-1.5) && (x>=-2.0) && (y>1.5) && (y<=2.0)){
    return(1.0);
  }
#endif
#ifdef RADIAL_INWARD
    return(1.0);
#endif
  return(0.0);
}
//bc at phi=0.0
double bc_x2i(double x, double y){
  return(0.0);
}
//bc at phi_final
double bc_x2f(double x, double y){
#ifdef SEMI_CLOCK
  //  if (x <= -1.0 && x>= -2.0)
    return(1.0);
#endif
  return(0.0);
}

float find_max(float a[], int n) {
  int i,index;
  float max; 
  max = a[0];
  index = 0;
  for (i = 1; i < n; i++) {
    if (a[i] > max) {
      index = i;
      max = a[i];
    }
  }
  return(max); 
}

float find_min(float a[], int n) {
  int i,index;
  float min; 
  min = a[0];
  index = 0;
  for (i = 1; i < n; i++) {
    if (a[i] < min) {
      index = i;
      min = a[i];
    }
  }
  return(min); 
}

float sum(float a[], int n) {
  int i,index;
  float sum; 
  sum = a[0];
  for (i = 1; i < n; i++) {
    sum+= a[i]; 
  }
  return(sum); 
}

/*Transform polar vector (i.e. edge velocity) from orthonormal local basis to physical basis */
void vector_coordinate_to_physical(double vr, double vphi, double phi, double *vx, double *vy){
  *vx = vr*cos(phi) -sin(phi)*vphi; 
  *vy = vr*sin(phi) +cos(phi)*vphi; 
  return;
}

/*Transform Cartesian vector to orthonormal local basis */
void vector_physical_to_coordinate(double vx, double vy, double phi, double *vr, double *vphi){
  *vr = vx*cos(phi) +sin(phi)*vy; 
  *vphi = -vx*sin(phi) +cos(phi)*vy; 
  return;
}
/* VELOCITY OPTIONS: PICK ONLY ONE */
void velocity_physical(double x_b,double y_b,double *vx,double *vy){
#ifdef HORIZONTAL
  *vx=1.0;
  *vy=0.0;
#endif
 
#ifdef RADIAL_INWARD
  double angle = atan2(y_b,x_b); 
  *vx = -cos(angle);
  *vy = -sin(angle);
#endif

  //these produce nonuniform radial velocities-- why do we get uniform scaled edge velocities with atan2(y,x) in stream function?
  //bc the act of differencing to get proper velocity requires discrete divergence free condition
  /*  *vx = -x_b/(x_b*x_b + y_b*y_b);
   *vy = -y_b/(x_b*x_b + y_b*y_b);*/
  return; 
}

void velocity_coordinate(double r_b,double phi_b,double *vr,double *vphi){
#if defined(SEMI_CLOCK)
  *vr =0.0;
  *vphi = -1.0;
#endif
#if defined(GAUSS_CLOCK)
  *vr =0.0;
  *vphi = -r_b;
#endif
  return;
}

/*Stream function in physical coordinates */
double stream_function(double x, double y){
  double radius = sqrt(x*x + y*y); 
  double phi = atan2(y,x);
  return(y); //stream1: velocity is moving in the y-direction
  //  return(radius);  //stream2: velocity is moving clockwise
  //return(-phi); //stream3: radially inward velocity
  //  return(phi); //stream3: radially outward velocity
  //return(x+y); //stream4 diagonally to the bottom right  
}

/* a duplication of the function in ATHENA FullRT_flux.c */
double flux_PLM(double ds,double *imu){
  //the upwind slope
  double delq1 = (imu[2] - imu[1])/ds;
  double delq2 = (imu[1] - imu[0])/ds;
  double dqi0;
  if(delq1*delq2 >0.0) //minmod function
    dqi0 = 2.0*delq1*delq2/(delq1+delq2);
  else
    dqi0 = 0.0;
  //unknown why ds is in this function. might have to do with nonuniform grids
  return(ds*dqi0); 
}


double gaussian(double x_0, double y_0,double x,double y){
  // set standard deviation to 1.0
  double sigma = 0.25;
  double kernel, s = 2.0 * sigma * sigma;
  // sum is for normalization
  double sum = 0.0; 

  kernel= exp(-1/s*((x_0 -x)*(x_0 -x) + (y_0 -y)*(y_0 -y)));
 
  // normalize the Kernel
  kernel /= sqrt(2*M_PI)*sigma;
  return(kernel);
}
