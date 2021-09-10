#include "auto_f2c.h"
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*   Grazing Managment model				                                      */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
  /* Local variables */
  doublereal B,W,Bx,Wx, yy,zz,L;
  doublereal P,eta,kappa,mu,nu,lambda,gamma,rho,DB,DW;


  P        = par[0];
  eta      = par[1];
  kappa    = par[2];
  mu       = par[3];
  nu       = par[4];
  lambda   = par[5];
  gamma    = par[6];
  rho      = par[7];
  DB       = par[8];
  DW       = par[9];

  yy       = par[18];
  zz 	   = par[19];
  L  	   = 1;

  B      = u[0];
  W      = u[1];
  Bx     = u[2];
  Wx     = u[3];

  f[0] = L * Bx;
  f[1] = L * Wx;

  f[2] = -L * ( lambda*W*B*(1+eta*B)*(1+eta*B)*(1-B/kappa) - mu*B + zz*Bx) / DB;
  f[3] = -L * ( P - nu*W/(1+rho*B/kappa) - gamma*W*B*(1+eta*B)*(1+eta*B)  + (yy+zz)*Wx) / DW;

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par)
{
  /* Initialize the equation parameters */
  doublereal B0,W0;
  doublereal P,eta,kappa,mu,nu,lambda,gamma,rho,DB,DW;

  P      = (doublereal)50.0;
  eta	   = (doublereal)1.5;
  kappa	 = (doublereal)1.0;
  mu     = (doublereal)2.0;
  nu     = (doublereal)5.0;
  lambda = (doublereal)0.1;
  gamma  = (doublereal)4.0;
  rho    = (doublereal)0.1;
  DB     = (doublereal)0.01;
  DW     = (doublereal)10.0;	

  par[0] = P;
  par[1] = eta;
  par[2] = kappa;
  par[3] = mu;
  par[4] = nu;
  par[5] = lambda;
  par[6] = gamma;
  par[7] = rho;
  par[8] = DB;
  par[9] = DW;

  par[18] = (doublereal).00000000000000000009;	/* yy  */
  par[19] = (doublereal).00000000000000000010;	/* zz  */

  B0 = 0.0;
  W0 = P/nu;
/* Initialize the solution */
  u[0] = B0;
  u[1] = W0;
  u[2] = 0.0;
  u[3] = 0.0;

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc)
{
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint)
{
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par)
{
  integer tmp;
  extern doublereal getp();
  extern doublereal getpuwe();

  if(par[10]>0)
  	par[11]=6.2832/par[10];
  else
	par[11]=0;
  par[12]=u[0]*u[0];

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
