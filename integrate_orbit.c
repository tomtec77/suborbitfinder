/**
 * @file integrate_orbit.c
 * @brief Functions to perform the integration of the galaxy orbit
 */
#include "allvars.h"
#include "readparameterfile.h"
#include "proto.h"

/** Minimum allowed subhalo mass (in code units). */
#define MINSUBMASS 1.0e-7

/** Minimum value allowed for code units. */
#define MINUNITS 1.0e-6

/** Minimum value to check for energy/angular momentum conservation. 
    For Target simulation, this can be as small as 1e-6. For DolagClusters 
    simulations, this needs to be 1e-4. */
#define MINCHECK 1.0e-4

/**
 * @brief Integrate orbit of a type 2 satellite.
 * @param p Index of selected galaxy
 * @param central Index of the corresponding central galaxy
 * @param tmax Maximum integration time
 * @param dt Integration timestep
 * @param eps Softening to use in the potential
 * @return 0 on success, -1 if energy increases (try reducing the 
 * softening)
 *
 * Note that this function uses its own system of units, where the unit of 
 * mass is the host halo mass, the unit length is the host halo radius
 * and the unit velocity is set so that G = 1 (i.e. it should be the virial
 * velocity).
 * Important: elements Posrel, Vrel, Orbit(X,Y,Z) in struct GALAXY are 
 * assumed to store information in physical units.
 * Note new galaxy type 4: disrupted subhalo
 */
int integrate_type2_orbit(int p, int central, double tmax, double dt,
    double eps)
{
  int elements, k, nstep;
  int trigger=0;
  double r[3], v[3], a[3], j[3];
  double mhalo, lambda_c, v2, vhalo, j2, j2_0, djdt;
  double j2prev;
  double r2, dr, rscale, rgx, rdyn;
  double ekin, epot, e_fin, e_in;
  double omega;
  double t, t_out, dtout;
  double tdyn;
  double tdyn_orbital;
  double mass, minit, radius, rtidal, rperi, rapo;
  double unitlength, unitmass, unittime, unitvel, unitenergy;
  double m_of_r, mstrip;
  struct GALAXY galaxySave;
  
  EmergencyStop = 0;

  if (GalaxyA[p].Type != 2) {
    fprintf(stderr, "Error (integrate_type2_orbit): "
	    "function called for galaxy not type 2\n"
	    "GalaxyA[%d].Type = %d - Exit\n", p, GalaxyA[p].Type);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  galaxySave = GalaxyA[p];

  /* NEW: satellites of type 1 galaxies can be stripped away, and they
     orbit the central galaxy of the FOF. Update the definition of
     central and the relative position and velocity accordingly */
  if (GalaxyA[p].Relocated == 1)
    central = FirstGalInFOFGroup_A[GalaxyA[p].ParentGroup];

  /* Use the host halo mass as unit mass */
  unitmass = (double)GalaxyA[central].Mvir;
  /*printf("p, central, unitmass: %d %d %g (PartMass %g)\n", p, central, 
    GalaxyA[central].Mvir, PartMass); fflush(stdout);*/
  
  /* Use the host halo virial radius as unit length */
  unitlength = (double)GalaxyA[central].Rvir;
  
  /* Unit velocity in code units (usually km/s) */
  unitvel = sqrt(G*unitmass/unitlength); 

  /* Unit time in main code units */
  unittime = (unitlength/unitvel);

  unitenergy = unitmass*unitvel*unitvel;

  tmax /= unittime;    /* Max time to internal units */
  dtout = tmax/STEPS;  /* Time interval to store the orbit */

#ifdef DEBUG
  printf("Unit velocity:   %g km/s (virial velocity %g)\n", unitvel,
	 GalaxyA[central].Vvir);
  printf("Unit time:       %g Myr\n", unittime*UnitTime_in_s / 
	 SEC_PER_MEGAYEAR / Hubble_h);
  /* Note that, unlike testorbit, we don't need to multiply by the 
     conversion from km to kpc because we use G in code units */

  printf("Maximum time:    %g Myr\n", tmax*unittime*UnitTime_in_s /
	 SEC_PER_MEGAYEAR / Hubble_h);
  printf("Output timestep: %g Myr\n", dtout*unittime*UnitTime_in_s /
	 SEC_PER_MEGAYEAR / Hubble_h);
  fflush(stdout);
#endif  

  if (unittime < MINUNITS || unitvel < MINUNITS || unitmass < MINUNITS
      || unitlength < MINUNITS) {
    fprintf(stderr, "Error (integrate_type2_orbit): "
	    "error while setting up units - Exit\n"); 
    fprintf(stderr, 
	    "  p, central: %d %d\n"
	    "  Unit length:   %g\n"
	    "  Unit time:     %g\n"
	    "  Unit velocity: %g\n"
	    "  Unit mass:     %g\n",
	    p, central, unitlength, unittime, unitvel, unitmass); 
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  /* TOMAS 2013-04-24: with the new units, this is 1 */
  /*rhalo = GalaxyA[central].Rvir*UnitLength_in_kpc/ULENGTH_ORBIT;*/
  
  /* Initial conditions for integration: start with the last recorded 
     satellite position, velocity, mass and outer radius */
  mass   = (double)GalaxyA[p].OrbitMass[STEPS-1]/unitmass;
  minit  = mass;
  radius = (double)GalaxyA[p].Rdm/unitlength;

  /* Coordinates and velocities relative to the central galaxy */
  r2 = 0.0;
  v2 = 0.0;
  for (k=0; k<3; k++) {
    r[k] = (double)GalaxyA[p].Posrel[k]/unitlength;
    v[k] = (double)GalaxyA[p].Vrel[k]*UnitVelocity_in_cm_per_s/unitvel/
      1.0e5;
    r2 += r[k]*r[k];
    v2 += v[k]*v[k];
  }
  rperi = sqrt(r2);
  rapo  = rperi;
  /*if (rperi > 1) {
    printf("WARNING: galaxy %d initial radius r=%g exceeds host "
    "Rvir\n", p, rperi); fflush(stdout);
    }*/

  /* Angular momentum */
  j[0] = r[1]*v[2] - r[2]*v[1];
  j[1] = r[2]*v[0] - r[0]*v[2];
  j[2] = r[0]*v[1] - r[1]*v[0];
  djdt = 0.0;
  j2 = 0.0;
  for (k=0; k<3; k++)
    j2 += j[k]*j[k];

#ifdef DEBUG
  printf("Angular momentum (previous): %g\n", 
	 sqrt(GalaxyA[p].Jorb[0]*GalaxyA[p].Jorb[0] + 
	      GalaxyA[p].Jorb[1]*GalaxyA[p].Jorb[1] + 
	      GalaxyA[p].Jorb[2]*GalaxyA[p].Jorb[2]));
  printf("Angular momentum (initial):  %g\n", sqrt(GalaxyA[p].J2init));
  printf("Angular momentum:            %g (%g)\n",
	 sqrt(j2)*unitlength*unitvel, sqrt(j2));
#endif

  /* Convert the initial angular momentum to local units */
  j2_0 = (double)GalaxyA[p].J2init/(unitlength*unitlength)/
    (unitvel*unitvel);

  /* Time-scale for tidal stripping: we will use the subhalo 
     dynamical time */
  /* TOMAS 2013-03-26: for a SIS host halo, if the outer radius of the
     subhalo is the tidal radius as in Tormen, Diaferio & Syer 1998, the
     dynamical time scales linearly with the subhalo's distance to the
     centre. I don't know if this is appropriate... Maybe as a first order
     treatment we should use the dynamical time at infall */
  omega = sqrt(j2)/r2;
  tdyn_orbital = 2*M_PI/omega;

  if (TSTimescaleSelect == TSTIME_SUBHALO)
    tdyn = sqrt(radius*radius*radius/mass);
  if (TSTimescaleSelect == TSTIME_SUBINFALL) {
    rdyn = (double)GalaxyA[p].Rvir/unitlength;
    tdyn = sqrt(rdyn*rdyn*rdyn*unitmass/((double)GalaxyA[p].Mvir));
  }
  dt = tmax/dt;
#ifdef DEBUG
  printf("Subhalo dynamical time:  %g Myr\n", 
	 tdyn*unittime*UnitTime_in_s/SEC_PER_MEGAYEAR/Hubble_h);
  /*printf("Orbital dynamical time:  %g Myr\n\n", 
    tdyn_orbital*unittime/SEC_PER_MEGAYEAR);*/
  printf("Timestep:                %g Myr\n", dt*unittime*UnitTime_in_s / 
	 SEC_PER_MEGAYEAR/Hubble_h);
  printf("Output timestep:         %g Myr\n", dtout*unittime*UnitTime_in_s/ 
	 SEC_PER_MEGAYEAR/Hubble_h);
  printf("Softening:               %g h^-1 kpc\n", unitlength*eps);
  fflush(stdout);
#endif

  /* Store the initial conditions, in PHYSICAL coordinates */
  GalaxyA[p].OrbitX[0] = r[0]*unitlength;
  GalaxyA[p].OrbitY[0] = r[1]*unitlength;
  GalaxyA[p].OrbitZ[0] = r[2]*unitlength;
  GalaxyA[p].OrbitVx[0] = v[0]*unitvel;
  GalaxyA[p].OrbitVy[0] = v[1]*unitvel;
  GalaxyA[p].OrbitVz[0] = v[2]*unitvel;
  GalaxyA[p].OrbitMass[0] = mass*unitmass;
  GalaxyA[p].OrbitType[0] = GalaxyA[p].Type;
  
  /*
   * Integrate the orbit using a kick-drift-kick scheme 
   */
  t     = 0.0;
  t_out = 0.0;
  nstep = 0;
  ekin  = 0.5*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  /*epot  = log(sqrt(r2)/rhalo)/rhalo;*/
  /* Softened isothermal sphere potential */
  epot = potential_sis(r[0], r[1], r[2], 1, 1, eps);
  /*log(sqrt(r2)) + 0.5*log(1 + r2/(eps*eps)) - 
    log(sqrt(r2)/eps) + eps*atan(sqrt(r2)/eps)/sqrt(r2);*/
  e_in  = ekin + epot;

#ifdef DEBUG
  Nprints = 0;
  output(e_in*unitenergy, ekin*unitenergy, epot*unitenergy);
#endif

  /* Initial acceleration. */
  for (k=0; k<3; k++)
    a[k] = accel(r[k], v[k], r2, v2, mass, (double)GalaxyA[p].Coulomb, eps);

  while (t < tmax) {

    for (k=0; k<3; k++) {          /* First kick, then drift */
      v[k] += 0.5*dt*a[k];
      r[k] += dt*v[k];
    }
    
    r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    
    for (k=0; k<3; k++) {          /* Second kick */
      a[k] = accel(r[k], v[k], r2, v2, mass, (double)GalaxyA[p].Coulomb, 
		   eps);
      v[k] += 0.5*dt*a[k];
    }
    
    t += dt;
    
    /* Find the pericentre of the orbit */
    if (sqrt(r2) <= rperi) rperi = sqrt(r2);

    /* Find also the apocentre */
    if (sqrt(r2) >= rapo) rapo = sqrt(r2);
    
    /* When the galaxy loses all of its initial specific angular momentum,
       a merger has happened */
    j2prev = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    j[0] = r[1]*v[2] - r[2]*v[1];
    j[1] = r[2]*v[0] - r[0]*v[2];
    j[2] = r[0]*v[1] - r[1]*v[0];
    j2 = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];

    /* Check energy */
    ekin  = 0.5*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    epot = potential_sis(r[0], r[1], r[2], 1, 1, eps);

    /* If dynamical friction is enabled, a galaxy should lose energy,
       never gain it. Jumps in energy may be caused by an inadequate 
       choice of softening or integration timestep */
    /* TOMAS 2013-05-14 Apparently, reducing the softening does not solve
       the problem for some cases */
    if ((ekin+epot)-e_in > MINCHECK) {
      printf("\nWARNING (integrate_type2_orbit): energy increasing "
	     "for galaxy %d (%g) (snap %d)\n", p, (ekin+epot)-e_in,
	     Snapshot); 
      /*printf("Retrying with smaller softening... (current %g)\n\n", eps);*/
      printf("Retrying with smaller timestep... (current %g)\n\n", dt);
      fflush(stdout);
      if (dt < 1.0e-13) {
        printf("WARNING: very small timestep - perhaps you should try "
               "with a smaller EpsGrav...\n"); fflush(stdout);
      }
      GalaxyA[p] = galaxySave;
      //return -1;
      return -2;
    }

#ifdef DEBUG
    /*if (Snapshot == 43 && p == 21) {
      printf("Checking angular momentum: j/j0 = %g (previous), %g "
	     "(current), r = %g, E = %g\n", sqrt(j2prev/j2_0), 
	     sqrt(j2/j2_0), sqrt(r2), ekin+epot); fflush(stdout);
	     }*/
#endif
    if (sqrt(j2/j2_0) < SAT_MINJORB) {
      if (GalaxyA[p].Type == 2) {
	GalaxyA[p].Type = 3;
	GalaxyA[p].TypeT[Snapshot] = 3;
      }
      GalaxyA[p].MergTimeActual += (float)(t*unittime);
      GalaxyA[p].MergTimeSat    += (float)(t*unittime);
      GalaxyA[p].MergerType     = 1;
      break;
    }

    /*
     * Control for jumps in the angular momentum. As the galaxy sinks
     * towards the halo centre, accelerations become large and the 
     * timestep may have to be adjusted. Exit and retry the orbit
     */
    djdt = (sqrt(j2) - sqrt(j2prev))/dt;
    if (djdt > MINCHECK) {
      printf("\nWARNING: increasing j (djdt=%g) for galaxy %d (type %d) "
	     "(snap %d)\n", djdt, p, GalaxyA[p].Type, Snapshot); 
      printf("r = %g - v = %g - j/j0 = %g jprev/j0 = %g djdt = %g\n", 
	     sqrt(r2)*unitlength, sqrt(v2)*unitvel, sqrt(j2/j2_0), 
	     sqrt(j2prev/j2_0), djdt);
      printf("E_in = %g, E_curr = %g\n", e_in, ekin+epot);
      printf("Retrying with smaller timestep... (current %g)\n\n", dt);
      fflush(stdout);
      if (dt < 1.0e-13) {
        printf("WARNING: very small timestep - perhaps you should try "
               "with a smaller EpsGrav...\n"); fflush(stdout);
      }
      GalaxyA[p] = galaxySave;
      return -2;
      //exit(EXIT_FAILURE);
    }

    /* TOMAS 2012-08-28 NEW: consider also mergers by proximity */
    if (ProximityMergerOn == 1) {
      rgx = 0.1;  /* Consider the typical size of a central galaxy
		     to be 10 per cent the halo size */
                  /* TOMAS 2013-04-25: this is based on the Shen et al
		     (2003) results for z ~ 0, so its probably not
		     accurate for high redshift */
      if (rperi < rgx) {
	/*printf("rgx = %g kpc\n", rgx*ULENGTH_ORBIT);*/
	if (GalaxyA[p].Type == 2) {
	  GalaxyA[p].Type = 3;
	  GalaxyA[p].TypeT[Snapshot] = 3;
	}
	GalaxyA[p].MergTimeActual += (float)(t*unittime);
	GalaxyA[p].MergTimeSat    += (float)(t*unittime);
	GalaxyA[p].MergerType     = 2;
	break;
      }
    }
    
    /* Mass loss due to tidal stripping */
    if (TidalStrippingOn == 1) {
      /* Mass of host halo within the position of the satellite */
      m_of_r = sqrt(r2);
  
      /* Calculate the tidal radius */
      rtidal = sqrt(r2)*pow(mass/m_of_r, 1.0/3.0);

      if (rtidal < radius) {
	/* Find the mass beyond rtidal in the satellite */
	mstrip = mass*(1.0 - rtidal/radius);
    
	/*printf("p, mass, mstrip, mstrip*dt/tdyn: %d %g %g %g\n", 
	  p, mass, mstrip, mstrip*dt/tdyn); fflush(stdout);*/

	/* Remove a fraction dt/tstrip of the mass beyond rtidal,
	   and update the bounding radius of the satellite */
	if (dt/tdyn < 1) mstrip *= dt/tdyn;
	if (fabsf(mass - mstrip)*unitmass > MINSUBMASS) {
	  radius *= (mass - mstrip)/mass;
	  mass -= mstrip;
	  GalaxyA[p].StrippedMass += (float)(unitmass*mstrip);
  	}
	else {
	  GalaxyA[p].StrippedMass += (float)mass;
	  radius = 0.0;
	  mass = 0.0;
	}
	if (mass < 0.0) {
	  fprintf(stderr, "Error (integrate_type2_orbit): mass = %g < 0"
		  " for galaxy %d - Exit\n", mass, p); fflush(stderr);
	  exit(EXIT_FAILURE);
	}
      }
    }

    /* Flag galaxies as disrupted */
    if (mass/minit <= SAT_DISRUPT) {
      GalaxyA[p].Type = 4;
      GalaxyA[p].TypeT[Snapshot] = 4;
      break; /* OJO no se si esta bien cortar aca */
    }

    /* Store the data, in physical coordinates for pos and vel */
    if (t >= t_out) {
      nstep++;
      if (nstep >= STEPS) break;
      GalaxyA[p].OrbitX[nstep] = (float)(r[0]*unitlength);
      GalaxyA[p].OrbitY[nstep] = (float)(r[1]*unitlength);
      GalaxyA[p].OrbitZ[nstep] = (float)(r[2]*unitlength);
      GalaxyA[p].OrbitVx[nstep] =(float)(v[0]*unitvel*1.0e5/
					 UnitVelocity_in_cm_per_s);
      GalaxyA[p].OrbitVy[nstep] = (float)(v[1]*unitvel*1.0e5/
					  UnitVelocity_in_cm_per_s);
      GalaxyA[p].OrbitVz[nstep] = (float)(v[2]*unitvel*1.0e5/
					  UnitVelocity_in_cm_per_s);
      GalaxyA[p].OrbitMass[nstep] = (float)(mass*unitmass);
      GalaxyA[p].OrbitType[nstep] = GalaxyA[p].Type;

      /*printf("r_DM, r_tidal: %g %g\n", radius, rtidal); fflush(stdout);*/

#ifdef DEBUG
    output(e_in*unitenergy, ekin*unitenergy, epot*unitenergy);
#endif

      t_out += dtout;    /* Advance output clock */
    }
      
    /* Alternative: control the angular momentum loss 
       (but see the observation above) */
    /*if (fabsf(j2/j2_0) > 1) {
      if (EmergencyStop == 0) {
	fprintf(stderr, "Error (integrate_type2_orbit): specific angular "
		"momentum increasing for galaxy %d\n"
		"j/j_0 = %g\n\n", p, sqrt(j2/j2_0)); 
	fflush(stderr);
      }
      EmergencyStop++;
      }*/
  } /* Close while t < tmax */
  
  /* If galaxy has merged, or completely disrupted (zero mass), fill out 
     the data with zeros and the galaxy type along the orbit with the final 
     state, but save the last DM radius and bound mass */
  if ((GalaxyA[p].Type == 3) && nstep <= STEPS-1) {
    for (k=nstep; k<STEPS; k++) {
      GalaxyA[p].OrbitX[k] = 0.0;
      GalaxyA[p].OrbitY[k] = 0.0;
      GalaxyA[p].OrbitZ[k] = 0.0;
      GalaxyA[p].OrbitVx[k] = 0.0;
      GalaxyA[p].OrbitVy[k] = 0.0;
      GalaxyA[p].OrbitVz[k] = 0.0;
      GalaxyA[p].OrbitMass[k] = GalaxyA[p].OrbitMass[nstep-1];
      GalaxyA[p].OrbitType[k] = GalaxyA[p].Type;
    }
  }

  /* Update position, velocity, DM bounding radius of galaxy */
  GalaxyA[p].Posrel[0] = GalaxyA[p].OrbitX[STEPS-1];
  GalaxyA[p].Posrel[1] = GalaxyA[p].OrbitY[STEPS-1];
  GalaxyA[p].Posrel[2] = GalaxyA[p].OrbitZ[STEPS-1];
  GalaxyA[p].Vrel[0] = GalaxyA[p].OrbitVx[STEPS-1];
  GalaxyA[p].Vrel[1] = GalaxyA[p].OrbitVy[STEPS-1];
  GalaxyA[p].Vrel[2] = GalaxyA[p].OrbitVz[STEPS-1];

  GalaxyA[p].Rdm = (float)(radius*unitlength);
  GalaxyA[p].Mdm = GalaxyA[p].OrbitMass[STEPS-1];

  GalaxyA[p].Jorb[0] = 
    GalaxyA[p].Posrel[1]*GalaxyA[p].Vrel[2] - 
    GalaxyA[p].Posrel[2]*GalaxyA[p].Vrel[1];
  GalaxyA[p].Jorb[1] = 
    GalaxyA[p].Posrel[2]*GalaxyA[p].Vrel[0] - 
    GalaxyA[p].Posrel[0]*GalaxyA[p].Vrel[2];
  GalaxyA[p].Jorb[2] = 
    GalaxyA[p].Posrel[0]*GalaxyA[p].Vrel[1] - 
    GalaxyA[p].Posrel[1]*GalaxyA[p].Vrel[0];
  
  /* Advance merging clock for non-merged / non-disrupted galaxies */
  if (GalaxyA[p].Type < 3) {
    GalaxyA[p].MergTimeActual += (float)(tmax*unittime);
    GalaxyA[p].MergTimeSat    += (float)(tmax*unittime);
  }

  GalaxyA[p].EnergyLoss += (float)((e_in - ekin - epot)*unitenergy);
  
  GalaxyA[p].Rpericentre = (float)(rperi*unitlength);
  GalaxyA[p].Rapocentre  = (float)(rapo*unitlength);

#ifdef DEBUG
  if (GalaxyA[p].Type >= 3)
    printf("Final time:              %g Myr\n", t*unittime*UnitTime_in_s /
	   SEC_PER_MEGAYEAR / Hubble_h);
  else
    printf("Final time:              %g Myr\n", tmax*unittime*UnitTime_in_s/
	   SEC_PER_MEGAYEAR / Hubble_h);
  if (DynamicalFrictionOn == 1) {
    printf("Expected merger time:    %g Myr (final time = %g times the "
	   "expected time)\n", 
	   GalaxyA[p].MergTime*UnitTime_in_s/SEC_PER_MEGAYEAR/Hubble_h, 
	   t*unittime/GalaxyA[p].MergTime);
    printf("Angular momentum (final): %g\n", 
	   sqrt(GalaxyA[p].Jorb[0]*GalaxyA[p].Jorb[0] + 
		GalaxyA[p].Jorb[1]*GalaxyA[p].Jorb[1] + 
		GalaxyA[p].Jorb[2]*GalaxyA[p].Jorb[2]));
    printf("Angular momentum loss:   %g\n", 
	   sqrt((GalaxyA[p].Jorb[0]*GalaxyA[p].Jorb[0] +
		 GalaxyA[p].Jorb[1]*GalaxyA[p].Jorb[1] +
		 GalaxyA[p].Jorb[2]*GalaxyA[p].Jorb[2])/GalaxyA[p].J2init));
    if (GalaxyA[p].Type == 3)
      printf("Merger happened:      Yes\n");
    else if (GalaxyA[p].Type <= 2) 
      printf("Merger happened:      No\n");
    else
      printf("Subhalo disrupted!\n");
  }
  printf("Orbital pericentre:   %g h^-1 kpc\n", 
	 rperi*unitlength*UnitLength_in_kpc);
  printf("Orbital apocentre:    %g h^-1 kpc\n", 
         rapo*unitlength*UnitLength_in_kpc);
  printf("Final galaxy mass:    %g h^-1 M_Sun (%.1f per cent of initial "
	 "mass)\n", GalaxyA[p].OrbitMass[STEPS-1]*UnitMass_in_Msun, 
	 GalaxyA[p].OrbitMass[STEPS-1]*100.0/minit/unitmass);
  printf("Total mass lost:      %g h^-1 M_Sun\n", 
	 GalaxyA[p].StrippedMass*UnitMass_in_Msun);
  fflush(stdout);
#endif

  return 0;
}


/**
 * Calculate a single component of the acceleration for the galaxy, 
 * including dynamical friction if enabled.
 * @param x Position coordinate of galaxy (either x, y or z)
 * @param dxdt The corresponding velocity component of galaxy
 * @param r2 Square of the current halocentric radius of the galaxy
 * @param v2 Square of current galaxy velocity
 * @param mass Mass of the galaxy
 * @param lnC Coulomb logarithm at infall
 * @param eps Softening length used in the gravitational potential
 * @return Acceleration in internal code units
 */
double accel(double x, double dxdt, double r2, double v2, double mass, 
	     double lnC, double eps)
{
  double a, rad, vel, vel2;

  /* TOMAS 2013-04-24: with the new choice of system of units, the host 
     halo virial quantities are all 1 */ 
  rad = sqrt(r2);
  vel = sqrt(v2);
  
  /* Gravitational potential term */
  /*a = -x/(r0*r2);*/
  /* NEW 2013-04-18: to avoid excessively large accelerations close to 
     the centre, we will use a softened isothermal sphere instead */
  a = -(1.0/rad)*(1.0 - eps*atan(rad/eps)/rad) * x/rad;
  
  /* Dynamical friction term */
  if (DynamicalFrictionOn == 1) {
    a -= (mass/r2)*lnC * (1.0/v2) *
      (erf(vel) - 0.5*sqrt(M_PI)*vel*exp(-v2)) * dxdt/vel;
  }
  
  return a;
}


#ifdef DEBUG
/**
 * Output information on the orbit integration to the standard output
 * @param ekin Kinetic energy of the galaxy
 * @param epot Potential energy of the galaxy
 */
void output(double e0, double ekin, double epot)
{
  if (Nprints == 0) {
    printf("===================================================================\n");
    printf("|  E_kinetic    | E_potential   |   E_total    |  E var.          |\n");
    printf("===================================================================\n");
  }
  printf("| %le  | %le | %le | %le   |\n", ekin, epot, ekin+epot, 
	 (ekin+epot)-e0);
  Nprints++;

  return;
}
#endif
