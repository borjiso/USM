/**
 * @file example.c
 * @author Borja Jimeno Soto
 * @date 2/06/2017
 */
#include "hoko.h"
#include <math.h>
#include "aux_fun.h"
extern int Yes_return;

#define EXAMPLE_PI 3.1415926535898
#define EARTH_RADIUS 6378.137 /* This is the WGS84 value in km. */
#define MJD_AT_DEC_31_1957 36203

#define ALL_OUTPUT 1
#define XYZ_OUTPUT 2
#define KEPLER_MEAN_OUTPUT 3
#define KEPLER_OSC_OUTPUT 4

int main (void)
{
	printf("Empiezo.\n");
  int L, NTOCH;
  double MIZ[6];
  struct bxprog BX;
  struct bixprog * BIX;
  double * TKNK;
  double epoch_in =0, inclination_deg=0, ra_asc_node_deg=0,
    eccentricity=0, arg_perigee_deg=0, mean_anomaly_deg=0,
    semimajor_axis_earth_radii=0, deg2rad=0, rad2deg=0, nu=0, snu=0, cnu=0, ballistic_coef=0;
  int prmodel=0, pock=0, pzdachi=0;
  int object_in=0, rev_num_in=0, I,J,i,j;
  char input_file[256];
  int NUM_STEPS_TO_PROPAGATE;
  double STEP_SIZE_IN_DAYS;
  int output_type;

  rad2deg = 180 / EXAMPLE_PI;
  deg2rad = EXAMPLE_PI / 180;

  /* Read the data from the input.txt file. */
printf("1.\n");
  snprintf(input_file, 256, "%s", "usm_input.txt");
  if (get_elset_input(input_file,
		      &epoch_in,
		      &object_in,
		      &rev_num_in,
		      &STEP_SIZE_IN_DAYS,
		      &NUM_STEPS_TO_PROPAGATE,
		      &semimajor_axis_earth_radii,
		      &inclination_deg,
		      &ra_asc_node_deg,
		      &eccentricity,
		      &arg_perigee_deg,
		      &mean_anomaly_deg,
		      &ballistic_coef,
          &prmodel,
          &pock,
          &pzdachi,
		      &output_type) == 1)
    {
      printf("Error reading input from %s.\n",input_file);
      return 1;
    }
	printf("epoch_in: %lf\n",epoch_in);
	printf("object_in: %d\n",object_in);
	printf("rev_num_in: %d\n",rev_num_in);
	printf("STEP_SIZE_IN_DAYS: %lf\n",STEP_SIZE_IN_DAYS);
	printf("NUM_STEPS_TO_PROPAGATE: %d\n",NUM_STEPS_TO_PROPAGATE);
	printf("semimajor_axis_earth_radii: %lf\n",semimajor_axis_earth_radii);
	printf("inclination_deg: %lf\n",inclination_deg);
	printf("ra_asc_node_deg: %lf\n",ra_asc_node_deg);
	printf("eccentricity: %lf\n",eccentricity);
	printf("arg_perigee_deg: %lf\n",arg_perigee_deg);
	printf("mean_anomaly_deg: %lf\n",mean_anomaly_deg);
	printf("ballistic_coef: %d\n",ballistic_coef);
	printf("pock: %d\n",pock);
	printf("pzdachi: %d\n",pzdachi);
	printf("output_type: %d\n",output_type);
	printf("prmodel: %d\n",prmodel);

  /* Allocate memory for the TKNK and BIX arrays. */
  if ((TKNK = (double *) calloc(NUM_STEPS_TO_PROPAGATE, sizeof(double))) == NULL)
    {
      printf("Cannot Allocate space for TKNK array.\n");
      return 1;
    }

  if ((BIX = /*(struct bixprog *)*/ calloc(NUM_STEPS_TO_PROPAGATE, sizeof(struct bixprog))) == NULL)
    {
      printf("Cannot Allocate space for BIX array.\n");
      return 1;
    }

  /* Populate the MIZ array with test data. */
  MIZ[0] = semimajor_axis_earth_radii * EARTH_RADIUS;
  MIZ[1] = eccentricity;
  MIZ[2] = inclination_deg * deg2rad;
  MIZ[3] = arg_perigee_deg * deg2rad;
  MIZ[4] = ra_asc_node_deg * deg2rad;
  MIZ[5] = mean_anomaly_deg * deg2rad;

  if (output_type == ALL_OUTPUT)
    {

      printf("These are the original Keplerian Elements in radians.\n");
      for (i = 0; i < 6; i++)
	{
	  printf("%lf\n", MIZ[i]);
	}
      printf("\n");
    }

  /************************************************************
   Prepare for the first call to PROGNOZ.
  *************************************************************/

  /* Populate TKNK which stores epochs. */
  for (I = 0; I < 1; I++)
    {
      TKNK[I] =  epoch_in - MJD_AT_DEC_31_1957;
    }

  /* Initial orbit determination as first measurement. */
  BX.DT       =  TKNK[0];
  BX.NO       =  object_in;
  BX.PZADACHI =  0;
  BX.POCK     =  -2;
  BX.DT       =  TKNK[0];
  BX.U        =  arg_lat(MIZ[1], MIZ[5], MIZ[3]);
  BX.A        =  MIZ[0];
  BX.I        =  MIZ[2];
  BX.DBY      =  MIZ[4];
  BX.L        =  MIZ[1]*cos(MIZ[3]);
  BX.H        =  MIZ[1]*sin(MIZ[3]);
  BX.KB       =  ballistic_coef;

  L           =  PROGNOZ(1,&BX,BIX,TKNK);

  /* Check for errors. */
  if (L > 0)
    {
      printf("Error return %d from PROGNOZ.\n", L);
    }
  else if (output_type == ALL_OUTPUT)
    {
      printf("Here are the contents of the output structure from PROGNOZ.\n");
      printf("\n");
      printf("Output Elements --------------------------\n");
      printf("Satellite %ld, Epoch: %lf\n", BIX[0].bix.NO, BIX[0].bix.DT);
      printf("\nECI Vector:\n");
      for (i = 0; i < 6; i++)
	printf("%lf\n",BIX[0].X[i]);

      printf("\nKeplerian Elements:\n");
      for (i = 0; i < 6; i++)
	printf("%lf\n", BIX[0].Y[i]);

      printf("\nNonsingular Mean Elements:\n");
      for ( i = 0; i < 6; i++)
	printf("%lf\n", BIX[0].Z[i]);

      printf("\nEccentricity: %lf\n", BIX[0].E);
      printf("Argument of Perigee, radians: %lf\n", BIX[0].W);
      printf("Apogee Altitude: %lf\n", BIX[0].HA);
      printf("Perigee Altitude: %lf\n", BIX[0].HP);
      printf("Geocentric range, in km: %lf\n", BIX[0].RAD);
      printf("Latitude: %lf\n", BIX[0].ALTITUDE);
      printf("Longitude: %lf\n", BIX[0].LONGITUDE);
    }

  /**********************************************************
   Begin subsequent calls to PROGNOZ.
  **********************************************************/

  /* Populate TKNK[] with additional epochs. */
  for (I = 0; I < NUM_STEPS_TO_PROPAGATE; I++)
    {
      TKNK[I] = epoch_in - MJD_AT_DEC_31_1957 + STEP_SIZE_IN_DAYS * I;
    }

  /* Populate the input structure. */
  BX.DT       =  TKNK[0];
  BX.NO       =  object_in;
  BX.PZADACHI =  0;
  BX.POCK     =  -2;
  BX.DT       =  TKNK[0];
  BX.U        =  arg_lat(MIZ[1], MIZ[5], MIZ[3]);
  BX.A        =  MIZ[0];
  BX.I        =  MIZ[2];
  BX.DBY      =  MIZ[4];
  BX.L        =  MIZ[1]*cos(MIZ[3]);
  BX.H        =  MIZ[1]*sin(MIZ[3]);
  BX.KB       =  ballistic_coef;

  L           =  PROGNOZ(NUM_STEPS_TO_PROPAGATE,&BX,BIX,TKNK);

  /* Check for errors. */
  if (L > 0)
    {
      printf("Error return %d from PROGNOZ.\n", L);
    }
  else if (output_type == ALL_OUTPUT)
    {
      for (i = 0; i < NUM_STEPS_TO_PROPAGATE; i++)
	{
	  printf("\nHere are the contents of the output structure from PROGNOZ.\n");
	  printf("\n");
	  printf("Output Elements --------------------------\n");
	  printf("Satellite %ld, Epoch: %lf\n", BIX[i].bix.NO, BIX[i].bix.DT);
	  printf("\nECI Vector:\n");
	  for (j = 0; j < 6; j++)
	    printf("%lf\n",BIX[i].X[j]);

	  printf("\nKeplerian Elements:\n");
	  for (j = 0; j < 6; j++)
	    printf("%lf\n", BIX[i].Y[j]);

	  printf("\nNonsingular Mean Elements:\n");
	  for (j = 0; j < 6; j++)
	    printf("%lf\n", BIX[i].Z[j]);

	  printf("\nEccentricity: %lf\n", BIX[i].E);
	  printf("Argument of Perigee, radians: %lf\n", BIX[i].W);
	  printf("Apogee Altitude: %lf\n", BIX[i].HA);
	  printf("Perigee Altitude: %lf\n", BIX[i].HP);
	  printf("Geocentric range, in km: %lf\n", BIX[i].RAD);
	  printf("Latitude: %lf\n", BIX[i].ALTITUDE);
	  printf("Longitude: %lf\n", BIX[i].LONGITUDE);
	}
    }
  else if (output_type == XYZ_OUTPUT)
    {
      for (i = 0; i < NUM_STEPS_TO_PROPAGATE; i++)
	{
	  printf("%09.9lf %09.9lf %09.9lf %09.9lf\n",
		 BIX[i].bix.DT,
		 BIX[i].X[0],
		 BIX[i].X[1],
		 BIX[i].X[2]);
	}
    }
  else if (output_type == KEPLER_MEAN_OUTPUT)
    {
      for (i = 0; i < NUM_STEPS_TO_PROPAGATE; i++)
	{
	  printf("%09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf\n",
		 BIX[i].bix.DT,
		 BIX[i].Y[0],
		 BIX[i].Y[1],
		 BIX[i].Y[2],
		 BIX[i].Y[3],
		 BIX[i].Y[4],
		 BIX[i].Y[5]);
        }
    }
  else if (output_type == KEPLER_OSC_OUTPUT)
    {
      /*       printf("%09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf\n",
	     TKNK[0] + MJD_AT_DEC_31_1957,
	     MIZ[0],
	     MIZ[2] * rad2deg,
	     MIZ[4] * rad2deg,
	     MIZ[1]*cos(MIZ[3]),
	     MIZ[1]*sin(MIZ[3]),
	     arg_lat(MIZ[1], MIZ[5], MIZ[3]) * rad2deg); */

      for (i = 0; i < NUM_STEPS_TO_PROPAGATE; i++)
	{
	  printf("%09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf\n",
		 BIX[i].bix.DT + MJD_AT_DEC_31_1957,
		 BIX[i].bix.A,
		 BIX[i].bix.I * rad2deg,
		 BIX[i].bix.DBY * rad2deg,
		 BIX[i].bix.L,
		 BIX[i].bix.H,
		 BIX[i].bix.U * rad2deg);
	}
    }

  free(BIX);
  free(TKNK);

  return 0;
}

double arg_lat(double e, double M, double g){
	double zeta;
	double Eo = M;
	double En1 = 0;
	double En = Eo;
	do{
		En1 = En-(f_prima(e,En,M)/f(e,En,M));
	}while(abs(En1-En)<10E-12);
	double tgf = sqrt((1+e)/(1-e)*tan(En1/2));
	double f = atan(tgf)*2;
	zeta = f+g;
	return zeta;
}

double f(double e, double E, double M){
	return (E -M -e*sin(E));
}

double f_prima(double e, double E, double M){
	return (1-e*cos(E));
}

int get_elset_input(char * filename_in,
		    double * epoch_in,
		    int * object_in,
		    int * rev_num_in,
		    double * step_size_in_days,
		    int * number_of_steps,
		    double * semimajor_axis_earth_radii,
		    double * inclination_deg,
		    double * ra_asc_node_deg,
		    double * eccentricity,
		    double * arg_perigee_deg,
		    double * mean_anomaly_deg,
		    double * ballistic_coef,
		    int * prmodel,
		    int * pock,
		    int * pzadachi,
		    int * output_type){
	printf("Abriendo fichero...\n");
	FILE *fe = fopen(filename_in, "r");
printf("Fichero abierto.\n");
  if(fe==NULL){
    return 1;
  }else{
	int x;
	printf("Obteniendo datos...\n");
	//fscanf(fe, " epoch_in = %ld", epoch_in);
	//if(x == 0){return 1;}
 	//printf("1 lectura\n");
	x = fscanf(fe, " object_in = %d", object_in);
	if(x == 0){return 1;}
  	printf("2 lectura\n");
	x = fscanf(fe, " rev_num_in = %d", rev_num_in);
	if(x == 0){return 1;}
  	printf("3 lectura\n");
	x = fscanf(fe, " step_size_in_days = %lf", step_size_in_days);
	if(x == 0){return 1;}
  	printf("4 lectura\n");
	x = fscanf(fe, " number_of_steps = %d", number_of_steps);
	if(x == 0){return 1;}
  	printf("5 lectura\n");
	x = fscanf(fe, " semimajor_axis_km = %lf", 		semimajor_axis_earth_radii);
	if(x == 0){return 1;}
  	printf("6 lectura\n");
	x = fscanf(fe, " inclination_deg = %lf", inclination_deg);
	if(x == 0){return 1;}
  	printf("7 lectura\n");
	x = fscanf(fe, " ra_asc_node_deg = %lf", ra_asc_node_deg);
	if(x == 0){return 1;}
 	printf("8 lectura\n");
	x = fscanf(fe, " eccentricity = %lf", eccentricity);
	if(x == 0){return 1;}
  	printf("9 lectura\n");
	x = fscanf(fe, " arg_perigee_deg = %lf", arg_perigee_deg);
	if(x == 0){return 1;}
 	printf("10 lectura\n");
	x = fscanf(fe, " mean_anomaly_deg = %lf", mean_anomaly_deg);
	if(x == 0){return 1;}
  	printf("11 lectura\n");
	x = fscanf(fe, " ballistic_coef = %lf", ballistic_coef);
	if(x == 0){return 1;}
 	printf("12 lectura\n");
	x = fscanf(fe, " pock = %d", pock);
	if(x == 0){return 1;}
  	printf("13 lectura\n");
	x = fscanf(fe, " pzadachi = %d", pzadachi);
	if(x == 0){return 1;}
 	printf("14 lectura\n");
	x = fscanf(fe, " output_type = %d", output_type);
	if(x == 0){return 1;}
	printf("15 lectura\n");
	x = fscanf(fe, " prmodel = %d", prmodel);
	if(x == 0){return 1;}
 	printf("16 lectura\n");
	fscanf(fe, " epoch_in = %lf", epoch_in);
	if(x == 0){return 1;}
 	printf("1 lectura\n");
	printf("antes\n");
	if(fe !=NULL){
		printf("cierra\n");
		fclose(fe);
	}
 return 0;
  }
}
