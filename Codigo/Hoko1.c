/**
 * @file Hoko1.c
 * @author Borja Jimeno Soto
 * @date 2/06/2017
 */
#include "hoko.h"
int LM = 2, LT = 0, LKL = ((LLUN + 1) * (LLUN + 2) - 6),
    LKS = ((LSUN + 1) * (LSUN + 2) - 6), JREZ[MAXREZ], L, M, PSUNDAB = 0,
    PATMOSFERA = 0, PLUNASUN = 0;

double

    A,
    E, I, SI, CI, SI2, CI2, DT0, RO, KA[6][5], KC[5][2], S0, HP, ALTITUDE,
    LONGITUDE, DTREENTRY;
double far *AT1, F135 = 150.0, FT = 150.0, KP = 3.3, KB0 = 0, KOTP = 0,
                 AX = 1.0, EX = 2.0, IX = 10.0, AG = 1.0, EG = 2.0, IG = 10.0;

double QQ[LMAX - 1], PP[LMAX - 1], RR[LMAX - 1], TT[LMAX - 1]

    ;
double GME = 398600.44, RZ = 6378.136, C20 = -1.08262683E-3;
#if defined(GEM_10)
double HK1[300] = {
    -484.16544E-6, 0,          0,          0,          2.43404E-6, -1.3991E-6,
    9.5838E-7,     0,          2.02855E-6, 2.5197E-7,  8.9272E-7,  -6.2346E-7,
    7.0028E-7,     1.41250E-6, 5.4112E-7,  0,          -5.3521E-7, -4.6926E-7,
    3.5208E-7,     6.6404E-7,  9.885E-7,   -2.0179E-7, -1.9531E-7, 2.9883E-7,
    6.862E-8,      0,          -5.117E-8,  -9.379E-8,  6.5146E-7,  -3.2769E-7,
    -4.6712E-7,    -2.0298E-7, -2.8754E-7, 4.99E-8,    1.5617E-7,  -6.5983E-7,
    -1.51E-7,      0,          -7.293E-8,  2.3E-8,     4.935E-8,   -3.5387E-7,
    5.697E-8,      3.32E-9,    -1.0039E-7, -4.6157E-7, -2.5833E-7, -5.373E-7,
    2.71E-9,       -2.4213E-7, 9.312E-8,   0,          2.7044E-7,  1.0196E-7,
    3.2437E-7,     1.0813E-7,  2.3109E-7,  -2.1615E-7, -2.8455E-7, -1.2984E-7,
    1.498E-8,      4.312E-8,   -3.617E-7,  1.3055E-7,  -7.170E-9,  1.688E-9,
    5.100E-8,      0,          5.0E-9,     3.1E-8,     8.8E-8,     4.2E-8,
    -9.00E-9,      -0.830E-7,  -2.37E-7,   6.8E-8,     -1.5E-8,    6.9E-8,
    -5.2E-8,       3.00E-7,    0.72E-7,    0.74E-7,    -0.106E-6,  0.127E-6,
    0.2754E-7,     0,          0.16E-6,    0.8E-8,     0.27E-7,    -0.357E-7,
    -0.162E-6,     -0.9E-7,    -9.8E-9,    1.36E-8,    -0.7E-8,    -5.6E-8,
    3.9E-8,        2.16E-7,    -1.03E-7,   0.748E-7,   0.199E-6,   -1.2E-8,
    -5.5E-8,       0.91E-7,    0.5261E-7,  0,          0.89E-7,    -1.3E-7,
    -8.5E-8,       -1.3E-8,    -1.8E-8,    -1.61E-8,   -9.7E-8,    -7.7E-8,
    -6.5E-8,       -3.1E-8,    -3.8E-8,    -8.4E-8,    0.4E-8,     1.8E-8,
    4.3E-8,        -6.9E-8,    1.2E-8,     -4.8E-8,    0.1E-7,     -2.2E-8,
    -4.8E-8,       0,          5.5E-9,     -2.4E-9,    3.1E-8,     -9.1E-8,
    -5.1E-8,       -1.3E-7,    -4.8E-8,    -7.4E-8,    4.7E-8,     7.1E-8,
    -3.5E-9,       3.1E-8,     1.4E-8,     -8.5E-8,    1.1E-8,     2.5E-8,
    -3.3E-8,       3.9E-8,     -6.5E-8,    -1.2E-9,    4.9E-8,     -7.0E-8,
    3.9E-8,        0,          -6.9E-8,    -5.3E-8,    1.8E-9,     -4.2E-9,
    5.9E-8,        2.6E-8,     -8.0E-8,    -2.2E-8,    4.4E-8,     6.5E-9,
    -2.1E-9,       2.6E-8,     -1.9E-8,    4.6E-8,     -2.6E-8,    2.5E-8,
    3.7E-8,        1.2E-8,     -4.7E-9,    5.1E-8,     1.7E-8,     -3.8E-9,
    -5.6E-9,       -1.3E-8,    4.4E-8,     0,          -3.6E-8,    2.6E-8,
    2.7E-8,        -5.7E-8,    -2.2E-8,    7.0E-8,     -1.6E-8,    -2.3E-9,
    6.0E-8,        5.0E-8,     -3.7E-8,    3.0E-9,     -1.2E-9,    -3.3E-9,
    -1.8E-8,       -3.8E-9,    2.1E-8,     4.5E-8,     3.1E-8,     -3.1E-8,
    -3.9E-8,       -1.1E-8,    -3.2E-8,    9.1E-8,     -6.0E-8,    7.0E-8,
    -2.3E-8,       0,          -5.5E-9,    4.0E-8,     -3.6E-8,    3.3E-8,
    3.3E-8,        -6.8E-9,    -6.4E-9,    -2.8E-10,   2.0E-8,     -1.7E-8,
    -9.5E-9,       -1.8E-9,    2.3E-8,     -2.2E-8,    -4.5E-8,    -4.3E-9,
    3.2E-8,        1.0E-8,     4.7E-8,     1.3E-8,     2.1E-8,     -3.7E-8,
    9.8E-9,        -3.1E-8,    2.8E-8,     4.2E-8,     -5.1E-8,    -5.4E-9,
    1.5E-9,        0,          -2.8E-9,    2.7E-9,     3.3E-9,     -2.0E-8,
    2.2E-8,        2.4E-8,     -4.1E-8,    -3.4E-9,    6.1E-10,    5.8E-10,
    2.4E-8,        -4.5E-8,    6.5E-8,     1.7E-8,     -1.6E-8,    2.9E-8,
    1.1E-8,        3.4E-8,     1.9E-8,     2.9E-9,     8.0E-11,    6.8E-10,
    -3.4E-8,       1.7E-8,     -2.2E-8,    -2.2E-9,    3.9E-9,     -2.5E-8,
    -2.1E-8,       -4.5E-9,    -7.3E-9,    0,          2.0E-8,     3.7E-9,
    -9.4E-9,       2.6E-8,     -1.1E-8,    -2.0E-8,    3.2E-8,     3.3E-8,
    -1.1E-8,       -1.4E-9,    -6.4E-9,    -2.9E-8,    -2.4E-9,    -7.7E-9,
    -1.9E-8,       8.3E-9,     -1.8E-8,    -4.5E-8,    9.0E-9,     -2.1E-9,
    2.3E-8,        6.3E-9,     1.8E-8,     8.8E-9,     1.2E-8,     -6.3E-9,
    -1.9E-8,       -3.8E-8,    -1.3E-8,    -2.6E-8,    -2.6E-8,    7.7E-9};
#else
double HK1[LGEM];
#endif

/*
double const    HK1[300]
={-484.1656973E-6,0,0,0,2.43690E-6,-1.3953E-6,9.5800E-7,0,2.03500E-6,2.4600E-7,8.8900E-7,-
  6.2600E-7,7.1800E-7,1.3800E-6,5.4200E-7,0,-5.2800E-7,-4.5300E-7,3.3300E-7,6.5500E-7,9.890E-7,-2.0500E-7,-
  1.9300E-7,3.0200E-7,6.900E-8,0,-7.600E-8,-7.800E-8,6.3900E-7,-3.0700E-7,-4.5900E-7,-2.17E-7,-3.20E-7,
  6.2E-8,1.91E-7,-6.74E-7,-1.51E-7,0,-5.7E-8,2.2E-8,3E-8,-3.57E-7,5.5E-8,-6E-9,-
  1.0400E-7,-4.4300E-7,-2.5100E-7,-5.240E-7,15.0E-9,-2.3200E-7,9.300E-8,0,2.7100E-7,9.4000E-8,3.3000E-7,
  1.2400E-7,2.2500E-7,-2.0900E-7,-2.6500E-7,-1.1900E-7,-2.00E-9,7.000E-9,-3.38E-7,1.31E-7,15.0E-9,31.0E-9,
  5.100E-8,0,5.0E-9,3.1E-8,8.8E-8,4.2E-8,-9.00E-9,-0.830E-7,-2.37E-7,6.8E-8,-1.5E-8,6.9E-8,-5.2E-8,3.00E-7,
  0.72E-7,0.74E-7,-0.106E-6,0.127E-6,0.2754E-7,0,0.16E-6,0.8E-8,0.27E-7,-0.357E-7,-0.162E-6,-0.9E-7,-9.8E-9,
  1.36E-8,-0.7E-8,-5.6E-8,3.9E-8,2.16E-7,-1.03E-7,0.748E-7,0.199E-6,-1.2E-8,-5.5E-8,0.91E-7,0.5261E-7,0,
  0.89E-7,-1.3E-7,-8.5E-8,-1.3E-8,-1.8E-8,-1.61E-8,-9.7E-8,-7.7E-8,-6.5E-8,-3.1E-8,-3.8E-8,-8.4E-8,0.4E-8,
  1.8E-8,4.3E-8,-6.9E-8,1.2E-8,-4.8E-8,0.1E-7,-2.2E-8,-4.8E-8,0,5.5E-9,-2.4E-9,3.1E-8,-9.1E-8,-5.1E-8,-
  1.3E-7,-4.8E-8,-7.4E-8,4.7E-8,7.1E-8,-3.5E-9,3.1E-8,1.4E-8,-8.5E-8,1.1E-8,2.5E-8,-3.3E-8,3.9E-8,-6.5E-8,-
  1.2E-9,4.9E-8,-7.0E-8,3.9E-8,0,-6.9E-8,-5.3E-8,1.8E-9,-4.2E-9,5.9E-8,2.6E-8,-8.0E-8,-2.2E-8,4.4E-8,6.5E-9,
  -2.1E-9,2.6E-8,-1.9E-8,4.6E-8,-2.6E-8,2.5E-8,3.7E-8,1.2E-8,-4.7E-9,5.1E-8,1.7E-8,-3.8E-9,-5.6E-9,-1.3E-8,
  4.4E-8,0,-3.6E-8,2.6E-8,2.7E-8,-5.7E-8,-2.2E-8,7.0E-8,-1.6E-8,-2.3E-9,6.0E-8,5.0E-8,-3.7E-8,3.0E-9,-1.2E-9
  ,-3.3E-9,-1.8E-8,-3.8E-9,2.1E-8,4.5E-8,3.1E-8,-3.1E-8,-3.9E-8,-1.1E-8,-3.2E-8,9.1E-8,-6.0E-8,7.0E-8,-
  2.3E-8,0,-5.5E-9,4.0E-8,-3.6E-8,3.3E-8,3.3E-8,-6.8E-9,-6.4E-9,-2.8E-10,2.0E-8,-1.7E-8,-9.5E-9,-1.8E-9,
  2.3E-8,-2.2E-8,-4.5E-8,-4.3E-9,3.2E-8,1.0E-8,4.7E-8,1.3E-8,2.1E-8,-3.7E-8,9.8E-9,-3.1E-8,2.8E-8,4.2E-8,-
  5.1E-8,-5.4E-9,1.5E-9,0,-2.8E-9,2.7E-9,3.3E-9,-2.0E-8,2.2E-8,2.4E-8,-4.1E-8,-3.4E-9,6.1E-10,5.8E-10,2.4E-8
  ,-4.5E-8,6.5E-8,1.7E-8,-1.6E-8,2.9E-8,1.1E-8,3.4E-8,1.9E-8,2.9E-9,8.0E-11,6.8E-10,-3.4E-8,1.7E-8,-2.2E-8,-
  2.2E-9,3.9E-9,-2.5E-8,-2.1E-8,-4.5E-9,-7.3E-9,0,2.0E-8,3.7E-9,-9.4E-9,2.6E-8,-1.1E-8,-2.0E-8,3.2E-8,3.3E-8
  ,-1.1E-8,-1.4E-9,-6.4E-9,-2.9E-8,-2.4E-9,-7.7E-9,-1.9E-8,8.3E-9,-1.8E-8,-4.5E-8,9.0E-9,-2.1E-9,2.3E-8,
  6.3E-9,1.8E-8,8.8E-9,1.2E-8,-6.3E-9,-1.9E-8,-3.8E-8,-1.3E-8,-2.6E-8,-2.6E-8,7.7E-9};
*/
/*########################## ������ GEM_input ###########################*/
//   ����������: ���� ������������� ��������������� ����
//   ����� ���������: GEM_input().
//	 �������� ����������:
//		true - ������������ �������.
//      false - ������ ����� �������������
/*=========================================================================
  Destination: Input of the Earth gravity field coefficients
  Call: GEM_input().
  Input data:none.
  Output data:
        true, if inputting was succesful,
    false, if it was detected input error.
=========================================================================*/
;
/**
@brief Input of the Earth gravity field coefficients
@returns bool true, if inputting was succesful, false, if it was detected input
error.

*/
bool JVC_DLL GEM_input(void) {
  FILE *f1;
  int I, l, m;
  double c, s;
  char string[125];
  if ((f1 = fopen("gem.txt", "r")) == NULL) {
    return (false);
  } else {
    if (fscanf(f1, "%s\n", string) != 1)
      return (false);
    if (fscanf(f1, "%d%d\n", &l, &m) != 2)
      return (false);
    if (fscanf(f1, "%lf\n", &GME) != 1)
      return (false);
    if (fscanf(f1, "%lf\n", &RZ) != 1)
      return (false);
    //            l_max=MIN(LGEM,(l+1)*(l+2)-6);
    for (I = 0; I < LGEM; I++)
      HK1[I] = 0;
    while (fscanf(f1, "%d%d%lf%lf\n", &l, &m, &c, &s) == 4) {
      if (l <= LMAX) {
        HK1[l * (l + 1) - 6 + 2 * m] = c;
        HK1[l * (l + 1) - 6 + 2 * m + 1] = s;
      }
    }
    C20 = sqrt(5) * (HK1[0]);
    fclose(f1);
    return (true);
  }
}

/*########################## ������ GEM_input ###########################*/

/*########################## ������ SIGN ###########################*/
;
/**
@brief Calculation of the sing of a number.
@param X double
@returns int 1, if X is greater than 0, -1 if X less than 0, 0, if X is equal 0.

*/
int JVC_DLL sign(double X) {
  if (X == 0)
    return (0);
  else
    return ((X > 0) ? 1 : -1);
}
/*########################## ������ SIGN ###########################*/
/*########################## ������ SIGN ###########################*/
;
/**
@brief
@param i
@returns bool true, if DEMO is true, false, if DEMO is false;

*/
bool JVC_DLL demo(int i) {
  return (DEMO);
}
/*########################## ������ SIGN ###########################*/

/*########################## ������ ������ ###########################*/
/* ======================================================================
   ����������: ���������� �������� ��������������� �������� �� ������
               ������� ��������� �������.
   ����� ���������: KEPLER(M,EKC).
   ������� ����������:
                � - ������� �������� (� ��������),
                E - �������������� (�/�).
         �������� ����������:
                EKC - ��������������� �������� (� ��������).
   ------------------------------------------------------------------- */
/* ======================================================================
   Destination: Calculation of Kepler equation.
   Call: KEPLER(M,EKC).
   Input data:
                � is mean anomaly in radian,
                E is eccentricity.
   Output data:
                E is eccentric anomaly in radian.
   ------------------------------------------------------------------- */

;
/**
@brief Calculation of Kepler equation.
@param M mean anomaly in radian.
@param E X eccentricity.
@returns double eccentric anomaly in radian.

*/
double JVC_DLL KEPLER(double M, double EKC) {
  double E, E0 = M;
  do {
    E = E0;
    E0 = M + EKC * sin(E0);
  } while (fabs(E - E0) > 1.0E-8);
  return (E0);
}
/*########################## ������ ������ ###########################*/
/*########################## ������ ALFDEL ###########################*/
/* ================================================================

         ����������: ������ ������� ����������� � ��������� ������.
         ����� ���������: ALFDEL(T,ALF,DEL).
         ������� ����������:
                T - ����� �� 0-�� ������ 1958�. ( � ������).
         �������� ����������:
                ALF,DEL - �������� ������� ����������� � ��������� ������ (�
   ��������).
   ------------------------------------------------------------------- */
/* ================================================================

         Destination: Calculation of right ascension and declination of the Sun.
         Call: ALFDEL(T,ALF,DEL).
         Input data:
                T is time in days from 0h UTC Dec 31 1957.
         Output data:
                ALF,DEL are right ascension and declination of the Sun in
   radiian.
   ------------------------------------------------------------------- */

;

/**
@brief Calculation of right ascension and declination of the Sun.
@param T is time in days from 0h UTC Dec 31 1957.
@param ALF right ascension of the sun in radian.
@param DEL right declination of the sun in radian.

@returns void
*/
void JVC_DLL ALFDEL(double T, double far *ALF, double far *DEL) {
  double TE, EC, MCC, LC, SLC;
  EC = 0.01675104 - 0.0000418 * (TE = (21183.5 + T + TDT_TAI + DAT) / 36525.0);
  MCC = 6.2565835 + 628.30193 * TE;
  LC = 4.881627921 + 628.3319503 * TE + 2 * EC * sin(MCC) +
       1.25 * EC * EC * sin(2 * MCC);
  *DEL = asin(sin(EC = 0.40931975 - 0.00022711097 * TE) * (SLC = sin(LC)));
  *ALF = asin(SLC * cos(EC) / cos(*DEL));
  if (cos(LC) <= 0)
    *ALF = PI - *ALF;
}
/*########################## ����� ALFDEL ###########################*/
/* ================================================================
         ����������: �������� ���������� ��������� E0[] � ���� EK[]
         ��������� : a(��), e, i(���), w(���), ���(���), �(���).
         ����      : X(��), Y(��), Z(��), Vx(��/c), Vy(��/c), Vz(��/c).
         ����� ���������: GHCK(E0,EK).
         ������������ ��������� � ��������� KEPLER � ���������������
         �������������� ���������� GME.
         ---------------------------------------------------------------- */
/* =====================================================================
         Destination: Transformation of Kepler elements to coordinates and
   velocities Call: GHCK(E0,EK). Input data: array E0[0:5], where E0[0]=a(km),
   E0[1]=e,E0[2]=i(���),E0[3]=w(���),E0[4]=���(���),E0[5]=�(���). Output data:
   array Ek[0:5], where EK[0]=X(km),EK[1]=Y(km),EK[2]=Z(km),
     EK[3]=Vx(km/c),EK[4]=Vy(km/c),EK[5]=Vz(km/c).
         ---------------------------------------------------------------- */

/**
@brief Transformation of Kepler elements to coordinates and velocities.
@param E0[] array E0[0:5], where
E0[0]=a(km),E0[1]=e,E0[2]=i(���),E0[3]=w(���),E0[4]=���(���),E0[5]=�(���)
@param XL X coordinate of the Moon in km.
@returns double
*/
double JVC_DLL GHCK(double far E0[], double far EK[]) {
  double A, B, E, F1, F2, K1, K2, K3, K4, W, U;
  E = KEPLER(E0[5], E0[1]);
  W = 2 * atan2(sqrt((1 + E0[1]) / (1 - E0[1])) * sin(0.5 * E), cos(0.5 * E));
  U = E0[3] + W;
  A = E0[0] * (1 - E0[1] * cos(E));
  EK[0] = A * ((K1 = cos(E0[4])) * (F1 = cos(U)) -
               (K2 = sin(E0[4])) * (F2 = sin(U)) * (K3 = cos(E0[2])));
  EK[1] = A * (K2 * F1 + K1 * F2 * K3);
  EK[2] = A * F2 * (K4 = sin(E0[2]));
  B = (A = sqrt(GME / (E0[0] * (1 - E0[1] * E0[1])))) * (1 + E0[1] * cos(W));
  EK[3] = (A *= E0[1] * sin(W)) * (K1 * F1 - K2 * F2 * K3) -
          B * (K1 * F2 + K2 * F1 * K3);
  EK[4] = A * (K2 * F1 + K1 * F2 * K3) - B * (K2 * F2 - K1 * F1 * K3);
  EK[5] = A * F2 * K4 + B * F1 * K4;
  return (U);
}
/* ------------------------- Function ���� ----------------------------- */

/* ==================================Function
   XBL================================= ����������: �������� �� ���� �[] �
   ���������� �������� ������ E[]. ����       : X(��), Y(��), Z(��), Vx(��/c),
   Vy(��/c), Vz(��/c). ���������� : a(��),���,���, P, Q,lambda(���). �����
   ���������: XBL(X,E). ������������ ��������������� �������������� ����������
   GME.

         ------------------------------------------------------------------- */
;
/* =====================================================================
         Destination: Transformation of coordinates and velocities to
   nonsingular elements Call: XBL(X,E) Input data: array X[0:5], where
   X[0]=X(km),X[1]=Y(km),X[2]=Z(km), X[3]=Vx(km/c),X[4]=Vy(km/c),X[5]=Vz(km/c).
         Output data: array E[0:5],
         ???  where E[0]=a(km), E[1]=,E[2]=,E[3]=,E[4]= ,E[5]= (���).
         ---------------------------------------------------------------- */

/**
 * @brief Transformation of coordinates and velocities to nonsingular elements.
 * @param  X array X[0:5], where X[0]=X(km),X[1]=Y(km),X[2]=Z(km),
 * X[3]=Vx(km/c),X[4]=Vy(km/c),X[5]=Vz(km/c)
 * @param  array E[0:5], ???  where E[0]=a(km), E[1]=,E[2]=,E[3]=,E[4]= ,E[5]=
 * (���).
 * @return double E1
 */

double JVC_DLL XBL(double far X[], double far E[]) {
  double C1, C2, C3, C, R, DR, MO, F1, F2, F3, F, A1, A2, B1, B2, E1, A;

  C1 = X[1] * X[5] - X[2] * X[4];
  C2 = X[2] * X[3] - X[0] * X[5];
  C3 = X[0] * X[4] - X[1] * X[3];
  C = sqrt(C1 * C1 + C2 * C2 + C3 * C3);
  DR = (X[0] * X[3] + X[1] * X[4] + X[2] * X[5]) /
       (R = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]));
  F1 = -(MO = GME / R) * X[0] + C3 * X[4] - C2 * X[5];
  F2 = -MO * X[1] - C3 * X[3] + C1 * X[5];
  F3 = -MO * X[2] + C2 * X[3] - C1 * X[4];
  F = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
  F1 /= F;
  F2 /= F;
  F3 /= F;
  A = sqrt(2 * C * (C + C3));
  E[0] = C * (C / GME) / (1 - (F / GME) * (F / GME));
  E[3] = -C2 / A;
  E[4] = C1 / A;
  A1 = F1 - 2 * E[4] * C * F3 / A;
  B1 = F2 + 2 * E[3] * C * F3 / A;
  B2 = (X[1] - C2 * X[2] / (C + C3)) / R;
  A2 = (X[0] - C1 * X[2] / (C + C3)) / R;
  E[1] = F * A1 / GME;
  E[2] = F * B1 / GME;
  A = 1 + sqrt(1 - E[1] * E[1] - E[2] * E[2]);
  F1 = R * GME * (E[2] + B2 - C * DR * E[1] / GME / A) / (C * C);
  F2 = R * GME * (E[1] + A2 + C * DR * E[2] / GME / A) / (C * C);
  E1 = atan2(F1, F2);
  E[5] = E1 - E[1] * F1 + E[2] * F2;
  return (E1);
}

/*########################## Function XBL ###########################*/
/*########################## Function LBK   ###########################*/
/* ======================================================================
         ����������: �������� ����������� ��������� ������ E[] � ��������� E1[].
         ���������� : a(��),���,���, P, Q, lambda(���).
         ���������  : a(��), e(���), i(���), w(���), ���(���), �(���).
         ����� ���������: LBK(E,E1)
         ����������� �������:
         ������� ��������:
                E={26569.13206,.005,.8660254038e-2,.224143868,.1294095226,1.047197551},
    PI=3.14159265359;
   �������� ��������: E1={26569.13206,.01,.52359878,.52359878,.52359878,.0}.
         ------------------------------------------------------------------- */
/* =====================================================================
         Destination: Transformation of nonsingular elements to Kepler elements
     Call: LBK(E,EK).
         Input data: array E[0:5],
???  where E[0]=a(km), E[1]=,E[2]=,E[3]=,E[4]= ,E[5]= (���).
     Output data: array EK[0:5],
     where EK[0]=a(km),EK[1]=e(km),EK[2]=i(rad),
     EK[3]=w(rad),EK[4]=W(rad),EK[5]=L(rad).
         ---------------------------------------------------------------- */
/**
 * @brief Transformation of nonsingular elements to Kepler elements.
 * @param  E  array E[0:5], ???  where E[0]=a(km), E[1]=,E[2]=,E[3]=,E[4]=
 * ,E[5]= (���).
 * @param  E1 array EK[0:5], where
 * EK[0]=a(km),EK[1]=e(km),EK[2]=i(rad),EK[3]=w(rad),EK[4]=W(rad),EK[5]=L(rad).
 * @return double U
 */
;
double JVC_DLL LBK(double far E[], double far E1[]) {
  double EA, U;
  E1[0] = E[0];
  if ((E1[1] = hypot(E[1], E[2])) < 1.0E-8)
    E1[1] = 1.0E-8;
  E1[3] = (E[2] == 0 && E[1] == 0) ? 0 : atan2(E[2], E[1]);
  E1[2] = hypot(E[3], E[4]);
  if ((E1[2] = 2.0 * asin(E1[2])) < 1.0E-8)
    E1[2] = 1.0E-8;
  E1[4] = (E[4] == 0 && E[3] == 0) ? 0 : atan2(E[4], E[3]);
  E1[5] = E[5] - E1[3];
  E1[3] -= E1[4];
  if (E1[3] < 0)
    E1[3] += PI2;
  if (E1[4] < 0)
    E1[4] += PI2;
  EA = E1[5] = fmod(E1[5], PI2);
  for (U = EA + 0.05; fabs(U - EA) > 1.0E-8; EA = E1[5] + E1[1] * sin(EA))
    U = EA;
  ;
  U = E1[3] +
      2 * atan2(sqrt((1 + E1[1]) / (1 - E1[1])) * sin(EA / 2), cos(EA / 2));
  U = fmod(U, PI2);
  if (U < 0)
    U += PI2;
  return (U);
}
/*########################## Function LBK ###########################*/
//		double ST(int D,double T0)

/*########################## Function ST   ###########################*/
/* ======================================================================
   ����������: ������ cpe����� ��������� ������� �� �������� ������
               �������.
   ����� ���������: ST(D,T0).
   ST(D,T0) - �������� ����������-��������.
         ������� ����������:
                D - ����� ����� �����, ������������� �� 0-�� ������ 1958 �. ��
   ��������� �������, T0 - ����������� ��������� ����� �� �������� ������.
   �������� ����������:
    ST(D,T0) - �������� �����, ��������������� ��������� ������� �������
    (D,T0).
   ------------------------------------------------------------------- */
/* =====================================================================
         Destination: Calculation of sidereal time
     Call: ST(D,T0).
         Input data: D - days from 0h UTC Dec 31, 1957;
                 T0 - current time from 0h UTC in days.
     Output data: sidereal time in radian.
         ---------------------------------------------------------------- */

;
/**
 * @brief Calculation of sideral time.
 * @param  D  days from 0h UTC Dec 31, 1957.
 * @param  T0 current time from 0h UTC in days.
 * @return double sideral time in radian.
 */
double JVC_DLL ST(int D, double T0) {
  double T1, A;
#if defined(EPOCHA_IS_QUASINERTIAL)
  T1 = (D + 2921);
  A = 0.277987616 + 0.00273781191 * T1 + 1.00273781191 * (T0 + DUT1);
  return (modf(A, &T1) * PI2);
#else
  T1 = (D + T0 + 21183.5 + DUT1) / 36525.0;
  A = 0.276919398147 + 100.0021359027 * T1 + T1 * T1 * 0.107523148148E-5;
  return (modf(A + T0, &A) * PI2);
#endif
};
/*########################## function ST ###########################*/

/*########################## function KBL   ###########################*/
/* ======================================================================
   ����������: �������� ���������� ��������� ������ E[] � ���������� E1[].
         ���������  : a(��), e(���), i(���), w(���), ���(���), �(���).
         ���������� : a(��), ���, ���, P, Q, lambda(���).
         ����� ���������: KBL(E,E1)
         ����������� �������:
         ������� ����������:
   E={26569.13206,.01,.52359878,.52359878,.52359878,.0},. PI=3.14159265359;
         �������� ����������:
                E1={26569.13206,.005,.8660254038e-2,.224143868,.1294095226,1.047197551};
         ------------------------------------------------------------------- */
/* =====================================================================
         Destination: Transformation of Kepler elements to nonsingular elements
     Call: LBK(E,EK).
         Input data: array E[0:5],
     where EK[0]=a(km),EK[1]=e(km),EK[2]=i(rad),
     EK[3]=w(rad),EK[4]=W(rad),EK[5]=L(rad).
     Output data: array EK[0:5],
???  where E[0]=a(km), E[1]=,E[2]=,E[3]=,E[4]= ,E[5]= (���).
         ---------------------------------------------------------------- */

;
/**
 * @brief Transformation of Kepler elements to nonsingular elements.
 * @param  E  array E[0:5], where
 * EK[0]=a(km),EK[1]=e(km),EK[2]=i(rad),EK[3]=w(rad),EK[4]=W(rad),EK[5]=L(rad).
 * @param  E1 array EK[0:5], ???  where E[0]=a(km), E[1]=,E[2]=,E[3]=,E[4]=
 * ,E[5]= (���).
 * @return A
 */
double JVC_DLL KBL(double far E[], double far E1[]) {
  double A;
  E1[0] = E[0];
  E1[1] = E[1] * cos(A = E[3] + E[4]);
  E1[2] = E[1] * sin(A);
  E1[5] = A + E[5];
  E1[5] = fmod(E1[5], PI2);
  if (E1[5] < 0)
    E1[5] += PI2;
  E1[3] = (A = sin(0.5 * E[2])) * cos(E[4]);
  E1[4] = A * sin(E[4]);
  return (A);
}
/*########################## Function KBL ###########################*/

/*########################## Function C202 ###########################*/
/* ======================================================================
   ����������: ���������� �������������, ������������ ���������� ��
                                                         ������� ������� ��
   ������ ��������� ���������. ����� ���������: C202(). ������� ����������:
                ������� ������� �(��), �������������� E, ����� � �������
   ���������� SI,CI, ����������� ��� ������ ��������� ��������� C20,
   �������������� ������ ����� RZ. �������� ����������: ������ �������������
   KC[2:6, 1:2]. ����������� �������: ������� ����������: �=26569.13206, E=.01,
   SI=sin(PI/6), CI=cos(PI/6),C20=-.001082626848, RZ=6378.137. ��������
   ����������: KC[2,2]=-.935383e-11, KC[3,2]=.78254139e-13,
   KC[4,1]=.26227595e-7, KC[4,2]=-.93553246e-9,
   KC[5,1]=-.129622723e-7,KC[5,2]=-.2055e-12, KC[6,1]=.5761024e-8,
   KC[6,2]=.93519593e-9.
         ------------------------------------------------------------------- */
/* ======================================================================
   Destination: Calculation of slow-changing functions for taken into
               account perturbations of order two due to the second zonal
  garmonics, (see formula (1.3.4) Topic 5.2, Report 1). Call: C202(); It is used
  variables, declared as global: A(axis in km),E (eccentricity), SI(sin(I)), CI2
  (cos(I)),CI2(cos(I/2)), where I is inclination. Output data: static array
  KC[0:4,0:1].
  ------------------------------------------------------------------- */

;
/**
 * @brief Calculation of slow-changing functions for taken into
account perturbations of order two due to the second zonal
garmonics, (see formula (1.3.4) Topic 5.2, Report 1). Call: C202(); It is used
variables, declared as global: A(axis in km),E (eccentricity), SI(sin(I)), CI2
(cos(I)),CI2(cos(I/2)), where I is inclination. Output data: static array
KC[0:4,0:1].
 * @return  void
 */
void JVC_DLL C202(void) {
  double

      Z,
      A1, H2, K, K1, K2;
  ;
  Z = SI * SI;
  H2 = sqrt(1 - E * E);
  A1 = RZ / A;
  A1 = A1 * A1 * A1 * A1;
  K = 3 * C20 * C20 * A1 / (32 * sqr(sqr(sqr(H2))));
  K1 = K * (14 - 15 * Z);
  KC[0][1] = -K1 * Z * H2 * H2 * E;
  KC[1][1] = 0.5 * K1 * E * E * SI * CI * CI2;
  K1 = 5 * Z - 4;
  K2 = 3 * Z - 2;
  KC[2][0] = 0.25 * K *
             (24 * K1 * K1 + Z * (136 - 170 * Z) +
              (56 - 36 * Z - 45 * Z * Z) * E * E + 24 * K1 * K2 * H2);
  KC[2][1] =
      0.25 * K * (Z * (60 * Z - 56) + (56 - 316 * Z + 270 * Z * Z) * E * E);
  KC[3][0] = K * CI * (40 * Z - 36 - (4 + 5 * Z) * E * E + 12 * H2 * K2);
  KC[3][1] = K * CI * (30 * Z - 14) * E * E;
  KC[4][0] = 0.25 * K * H2 *
             (80 - 200 * Z + 130 * Z * Z + (40 - 40 * Z - 25 * Z * Z) * E * E +
              16 * H2 * K2 * K2);
  KC[4][1] = 0.25 * K * H2 * (56 - 60 * Z) * Z * (1 - 2.5 * E * E)

      ;
}
/*########################## ����� C202 ###########################*/
/*########################## ������ BIN ###########################*/
/* ======================================================================
         ����������: ������ �� ��������� ������� b,a,k �������� �������
               b =b*(a!/(k!*(a-k)!)
   ����� ���������: BIN(B,A,K).
   BIN(B,A,K) - �������� ����������-��������.
   ����������� �������:
   ������� ����������:
   B=5, A= 5, K=3.
   �������� ����������: BIN=50.
   ------------------------------------------------------------------- */
/* ======================================================================
   Destination: Calculation of function b =b*(a!/(k!*(a-k)!).
   Call: BIN(B,A,K).
   Input data: B, A, K.
   Output data: result of function value calculation.
  ------------------------------------------------------------------- */

;
/**
 * @brief Calculation of function b =b*(a!/(k!*(a-k)!).
 * @param  B
 * @param  A
 * @param  K
 * @return   double result of function value calculation.
 */
double JVC_DLL BIN(double B, int A, int K) {
  double S;
  int I, A1;
  S = ((A < K) ? 0 : B);
  I = A - K;
  if (K > I)
    K = I;
  A1 = A + 1;
  for (I = 1; I <= K; I++)
    S *= (A1 - I) / (double)I;
  return (S);
}
/*########################## ����� BIN ###########################*/
/*########################## ������ FINC ###########################*/
/* ======================================================================
   ����������: ������ �� ��������� ������������� ���������� L,M,P �
               ���������� i �������� ������� ���������� F (i) � ��
               �����������.
   ����� ���������: FINC(L,M,P,S,C,F1,F2).
   S=sin(.5*i); C=cos(.5*i);
   ������������ ��������� � ��������� BIN.
   ����������� �������:
   ������� ����������:
   �) L= 6, M= 5, P=0, I=PI/4;
   b) L= 6, M= 6, P=6, I=PI/4.
   �������� ����������:
   �) F1=1.34382208,F2=-1.43932458;
   b) F1=0.2389e-4,F2=0.34605e-3.
   ------------------------------------------------------------------- */
/* ======================================================================
   Destination: Calculation of inclination function and its derivative.
   Call: FINC(L,M,P,S,C,F1,F2),
         where S=sin(.5*i); C=cos(.5*i);
   Input data: L, M, P, C, S.
   Output data: FF1 - value of inclination function,
                FF2 - derivative of inclination function.
  ------------------------------------------------------------------- */

;
/**
 * @brief Calculation of inclination function and its derivative.
 * @param  L
 * @param  M
 * @param  P
 * @param  S   sin(.5*i)
 * @param  C   cos(.5*i)
 * @param  FF1 value of inclination function.
 * @param  FF2 derivative of inclination function.
 * @return     void
 */
void JVC_DLL FINC(int L, int M, int P, double S, double C, double far *FF1,
                  double far *FF2) {
  int I, J, J1, J2, K, R, N, LM;

  double D, F, G, CK2, SK2, CK, SK;
  K = 2 * L;
  LM = L - M;
  *FF1 = 0;
  *FF2 = 0;
  N = L + 1;
  G = BIN(1.0, L, P) / pow(2, L) * sqrt((1 + sign(M)) * (K + 1));
  if (((LM + 1) / 2) % 2)
    G *= -1;
  CK2 = C * C;
  SK2 = S * S;
  for (R = 1; R <= M; R++)
    G *= sqrt((double)(L + R) / (N - R));
  R = 2 * P;
  N = K - R;
  J1 = ((LM > R) ? LM - R : 0);
  J2 = ((LM > N) ? N : LM);
  CK = pow(C, LM + N - 2 * J1 + 2);
  SK = pow(S, L + M - N + 2 * J1 - 2);
  for (J = J1; J <= J2; J++) {
    D = BIN(1.0, N, J);
    if (J % 2)
      D = -D;
    F = BIN(D, R, LM - J);
    I = N - 2 * J;
    K = LM + I;
    I = L + M - I;
    CK /= CK2;
    D = F * CK;
    SK *= SK2;
    *FF1 += D * SK;
    *FF2 -= F * K * CK / C * SK * S;
    if (I > 0)
      *FF2 += D * I * C * SK / S;
  };
  *FF1 *= G;
  *FF2 *= G * 0.5;
}
/*########################## ����� FINC ###########################*/
/*########################## ������ HANSEN   ###########################*/
/* ======================================================================
   ����������: ������ �� ��������� ������������� ���������� N,Q,J �
               ��������������� e �������� ������� ������� H(e) � ��
               �����������.
   ����� ���������: HANSEN(N,P,Q,E,A,A1).
   ����������� �������:
   ������� ����������:
   �) N=-3, P= 0, Q=0, E=.75;
   b) N=-4, P=-1, Q=0, E=.75.
   �������� ����������:
   �) A=3.45567,A1=17.77204 ;
   b) A=5.92401,A1=58.67595.
   ------------------------------------------------------------------- */
/* ======================================================================
   Destination: Calculation of Hansen function and its derivative.
   Call: HANSEN(N,P,Q,E,A,A1).
   Input data: N, P, Q, E (eccentricity).
   Output data: A - value of Hansen function,
                A1 - derivative of Hansen function.
  ------------------------------------------------------------------- */

;
/**
 * @brief Calculation of Hansen function and its derivative.
 * @param  N  eccentricity
 * @param  P  eccentricity
 * @param  Q  eccentricity
 * @param  E  eccentricity
 * @param  A  value of Hansen function
 * @param  A1 derivative of Hansen function
 * @return    void
 */
void JVC_DLL HANSEN(int N, int P, int Q, double E, double far *A,
                    double far *A1) {
  int S, K, M;

  double BT, D, NU, V, V1, V2, V3, V4, V5, U, U1, U2, U3, U4, U5, W, W1, W2, W3,
      W4, W5, B, C;
  if (E < 0.5) {
    BT = E / (1 + sqrt(1 - E * E));
    if (Q >= P) {
      K = Q;
      M = P;
    } else {
      K = -Q;
      M = -P;
    };
    D = BT * BT;
    NU = K / (1 + D);
    V = 0;
    V3 = 0;
    for (S = 0; S <= (K - M); S++) {
      if (S == 0) {
        U = 1.0;
        U3 = 0.0;
      } else {
        U = (V * (NU - N + M + S - 2) - BT * NU * W) * BT / S;
        U3 = ((V * (1 - D) - 2 * BT * W) * NU / (1 + D) + V * (M - N + S - 2) +
              (V3 * (NU - N + M + S - 2) - BT * NU * W3) * BT) /
             S;
      };
      V4 = W3 = V3;
      V3 = U3;
      V1 = W = V;
      V = U;
      ;
    };
    *A = 0.0;
    *A1 = 0.0;
    V2 = 0.0;
    V5 = 0.0;
    S = 0;
    U1 = U;
    U2 = 1.0;
    U4 = U3;
    U5 = 0.0;
    do {
      if (S > 0) {
        U1 = (V1 * (NU - N + K + S - 2) - BT * NU * W1) * BT / (S + K - M);
        U2 = (V2 * (S - 2 - NU - N - M) + BT * NU * W2) * BT / S;
        U4 = ((V1 * (1 - D) - 2 * BT * W1) * NU / (1 + D) +
              V1 * (K + S - 2 - N) +
              (V4 * (NU - N + K + S - 2) - BT * NU * W4) * BT) /
             (S + K - M);
        U5 = (-(V2 * (1 - D) - 2 * BT * W2) * NU / (1 + D) +
              V2 * (S - 2 - N - M) +
              (V5 * (S - 2 - NU - N - M) + BT * NU * W5) * BT) /
             S;
      };
      W4 = V4;
      V4 = U4;
      W5 = V5;
      V5 = U5;
      B = U1 * U2;
      C = U1 * U5 + U2 * U4 - 2 * BT * (N + 1) * B / (1 + D);
      *A += B;
      *A1 += C;
      W1 = V1;
      V1 = U1;
      W2 = V2;
      V2 = U2;
      S++;
    } while ((fabs(B) > fabs(*A) * 1.0E-3) || (fabs(C) > fabs(*A1) * 1.0E-3));
    /*	  (*A)*=U=pow(1+D,-(N+1.0));*/
    U = 1.0;
    if (N + 1) {
      V = (N + 1 > 0) ? 1.0 / (1 + D) : 1 + D;
      for (S = 0; S < abs(N + 1); S++)
        U *= V;
    }
    (*A) *= U;

    (*A1) *= U / (1 + sqrt(1 - E * E)) / sqrt(1 - E * E);
  } else {
    *A = XNEWC1(N, P, Q, E);
    *A1 = (*A - XNEWC1(N, P, Q, E - 0.001)) / 0.001;
  }
}
/*########################## ����� HANSEN ###########################*/
/*########################## ������ XNEWC1 ###########################*/
/**
 * @brief Calculation of Hansen function.
 * @param  n eccentricity
 * @param  m eccentricity
 * @param  k eccentricity
 * @param  E eccentricity
 * @return   double value of Hansen function.
 */
double JVC_DLL XNEWC1(int n, int m, int k, double E) {
  /*      double XNEWC1(N,M,K,E,X) */
  double far *ars, far *y, far *z, far *w, X, hk, e2, pk1, pk2, pn2, pkn, anm,
      anj, anp, pknm, ann, ai, s;
  int mm, kk, km, ns, nt, nta, i, j;
  /*
  C  ���������� �������� ������������ �������  X = X(n,m,k;e)
  C            (�� �������� ������� - �������)
  C  N, M, K - �������� �������� ��������  n, m, k
  C        E - �������� ( �������������� )
  C      ����������� �� ��, ��� � �  XHANSV
  */
  /* ======================================================================
     Destination: Calculation of Hansen function.
     Call: X=XNEWC1(N,M,K,E,X).
     Input data: N, M, K, E (eccentricity).
     Output data: X is value of Hansen function,
    ------------------------------------------------------------------- */

  if ((ars = (double far *)calloc(1, 100 * sizeof(double))) == NULL ||
      (y = (double far *)calloc(1, 300 * sizeof(double))) == NULL ||
      (z = (double far *)calloc(1, 100 * sizeof(double))) == NULL ||
      (w = (double far *)calloc(1, 100 * sizeof(double))) == NULL) {
    perror("\n������ ��� ������� ������������ ������� ������������.\n");
    return (2);
  }

  if (k <= m) {
    mm = m;
    kk = k;
  } else {
    mm = -m;
    kk = -k;
  }
  km = kk - mm;
  if (km > 200 || E > 0.95) {
    /*
     print*,'      ������ ������������� ������� �� ��������,'
     print*,' ��������� k-m ��� E   -   ��� ��������� ���������'
     print*,'  k-m =',km,'   e =',E
    */
    return (9999);
  }
  hk = 0.5 * kk;
  ns = 7.0 + 97.0 * E;
  nt = ns + km;
  e2 = E * E;
  /*C  ���������� ��������������� ����������*/
  ars[1] = -1.5;
  for (i = 2; i <= ns; i++)
    ars[i] = (-2.5 + i) * ars[i - 1] / i;
  /*C ���������� ������������ ��������� H������*/
  pn2 = 0.5 * n;
  pk1 = kk - pn2;
  pk2 = pk1 - pn2;
  pkn = (kk + n + 4) * 0.25;
  y[nt + 1] = 1.0;
  y[nt] = pk1 - 1.0;
  for (i = 2; i <= nt; i++)
    y[nt + 1 - i] =
        ((pk1 - i) * y[nt + 2 - i] + 0.25 * (pk2 - i) * y[nt + 3 - i]) / i;
  /*C ���������� ������ ��������� H������*/
  z[1] = y[1];
  anm = ns - mm;
  for (j = 1; j <= ns; j++) {
    y[1] = y[j + 1];
    anj = anm - j;
    anp = anj - pn2;
    pknm = pkn - 0.75 * anj;
    ann = (anp - pn2) * 0.25;
    y[2] = anp * z[1] - (pknm - 1) * y[1];
    for (i = 2; i <= j; i++) {
      ai = 1.0 / i;
      y[i + 1] = (anp * z[i] + ann * w[i - 1] - (pknm - i) * y[i]) * ai;
      s = 0.0;
      for (nta = 2; nta <= i; nta++)
        s += ars[nta] * y[i + 1 - nta];
      y[i + 1] += hk * ai * s;
    }
    for (i = 1; i <= j; i++) {
      w[i] = z[i];
      z[i] = y[i];
    }
    z[j + 1] = y[j + 1];
  }
  /*C ���������� ������������ �������*/
  X = y[ns + 1];
  for (i = 1; i <= ns; i++)
    X = y[ns + 1 - i] + e2 * X;
  X *= pow(E, km);
  free(ars);
  free(y);
  free(z);
  free(w);
  return (X);
}
/*########################## ����� XNEWC1 ###########################*/
/*########################## ������ QPRT ###########################*/
/* ======================================================================
   ����������: ���������� �������� �������, ������������ �������� �������
          � ��������� ������������������ ���������� �� ��������� ��������.
   ����� ���������: QPRT().
         ������� ����������:
         ������� ������� �(��.), �������������� E, ���������� I(���.),
   ���������- ��� ����� ����������� �������� L�, �������������� ������ �����
   RZ(��.), ������ ������������� ���������� ������������� HK1[1:300,1:2].
         �������� ����������: QQ[0:L�-2],PP[0:L�-2],RR[0:L�-2],TT[0:L�-2].
         ����������� �������:
         ������� ����������:
         A=26569.13206, E=.01, I=.5235987756, L�=7.
   �������� ����������:
   QQ {0., .17992277e-7, .2847690667e-10, -.2105617865e-14, 0., 0.},
   PP {.218608e-11, -.18e-7, .28487e-10, .21057e-14, .14506e-16},
   RR {-.257069e-8, -.16749e-10, .137e-12, .2758e-16, .105e-18, 0.},
   TT {.3037e-9, -.1435e-8, .1456e-11, .4095e-16, 0., 0.}.
   ------------------------------------------------------------------- */
/* ======================================================================
         Destination: Calculation of slow-changing functions for taken into
     account perturbations due to the zonal garmonics,
     (see formula (1.3.10) Topic 5.2, Report 1).
         Call: QPRT().
         It is used variables, declared as global:
         �(semimajor axis, km), E (eccentricity),
     SI2(sin(I/2)),CI2 (cos(I/2)), where I is inclination.
         Otput data: static arrays with results of functioon calculation
                QQ[0:LMAX-2],
        PP[0:LMAX-2],
                RR[0:LMAX-2],
        TT[0:LMAX-2].
   ------------------------------------------------------------------- */

;
/**
 * Calculation of slow-changing functions for taken into account perturbations
 * due to the zonal garmonics.
 * @return  [description]
 */
void JVC_DLL QPRT(void) {

  int I0, K, K1, K2, M, N, N1, A1, A3, J, L;
  double Y = 1.0, E1, F, DF, B, A2, X0 = 1.0, X1, X2, DX0 = 0.0, DX1, DX2, EE1,
         EE2 = 1.0, H0, H1, H2, DH0, DH1, DH2;

  ;
  X1 = (3 * SI * SI - 2) / 4.0;
  DX1 = 1.5 * SI * CI;
  E1 = 1 - E * E;
  EE1 = 0.5 * E / E1;
  for (K = 0; K <= LM - 2; K++) {
    ;
    QQ[K] = PP[K] = RR[K] = TT[K] = 0.0;
    M = K + 1;
    K1 = M;
    K2 = K * K;
    H0 = 0.0;
    DH0 = 0.0;
    if (K > 0)
      EE2 *= EE1;
    H1 = EE2 / sqrt(E1);
    DH1 = ((E == 0) && (K == 0)) ? 0.0 : H1 * (K + (K + 1) * E * E) / (E * E1);
    L = ((K1 < 3) ? 3 : K1) /*INT(MAX(K1,3.0))*/
        ;
    if ((L - K) % 2)
      L++;
    I0 = L;
    while (I0 <= LM) {
      N = I0 - 2;
      N1 = 2 * N;
      A1 = (N1 - 1) * (I0 * I0 - K2);
      A2 =
          (N1 + 1) * (2 * (N * N + K2 + N - 1) - SI * SI * (N1 - 1) * (N1 + 3));
      A3 = (N1 + 3) * (sqr(N - 1) - K2);
      X2 = -(A2 * X1 + A3 * X0) / A1;
      DX2 = -(A2 * DX1 + A3 * DX0 -
              2 * (N1 + 1) * (N1 - 1) * (N1 + 3) * SI * CI * X1) /
            A1;
      B = onedeg((I0 + 1) / 2 + (I0 + K) / 2) * sqrt(2.0 * I0 + 1.0);
      F = B * X2;
      DF = B * DX2;
      X0 = X1;
      X1 = X2;
      DX0 = DX1;
      DX1 = DX2;
      for (J = M; J < I0; J++) {
        H2 = J * ((2 * J - 1) * H1 - (J - 1) * H0) / ((J * J - K2) * E1);
        DH2 = 2 * E * H2 / E1 +
              J * ((2 * J - 1) * DH1 - (J - 1) * DH0) / ((J * J - K2) * E1);
        H0 = H1;
        H1 = H2;
        DH0 = DH1;
        DH1 = DH2;
      };
      M = I0;
      B = pow(RZ / A, I0) * HK1[I0 * (I0 + 1) - 6];
      QQ[K] += onedeg(I0) * B * F * H2;
      PP[K] += B * F * DH2;
      RR[K] += B * DF * H2;
      TT[K] += B * F * H2 * 2.0 * (I0 + 1);
      I0 += 2;
    };
    N = ((K == 0) ? 1 : 2);
    QQ[K] *= 2.0 * K / E;
    PP[K] *= N;
    RR[K] *= N / (2.0 * sqrt(E1) * CI2);
    TT[K] *= N;
    Y *= (2 * K1 - 1.0) / (2.0 * K1);
    X1 = Y * pow(SI, K1);
    DX1 = Y * pow(SI, K) * CI * K1;
    DX0 = 0.0;
    X0 = 0.0;
  }

  ;
}
/*########################## ����� QPRT ###########################*/
#if HOKO

double AT[525] =
    {-1.829910e+01, 7.009000e-01,  1.153430e+02,  -5.101900e+00, 6.258100e-02,
     -1.672100e-04, -6.828000e-01, 5.576200e-03,  9.523800e-07,  -4.384000e+00,
     8.063000e-02,  -4.925000e-04, 1.042000e-06,
     1.506000e+00, //����� �������� ����������� n1
     5.411000e-01,  -7.238000e+00, 1.203000e-01,  -6.450000e-04, 1.208000e-06,
     -1.200000e-01, 5.000000e-03,  1.500000e-02,  -1.197500e-02, 9.983000e-05,
     0.000000e+00,

     -1.556050e+01, 8.248000e-01,  7.691320e+01,  -1.721000e-01, 5.756000e-03,
     -3.635000e-06, -8.607000e-01, 7.861000e-03,  -5.711000e-06, 1.279100e+00,
     -1.576000e-02, 6.499000e-05,  -5.145000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.411000e-01,  -2.152000e-01, 4.167000e-03,  1.587000e-06,  -1.651000e-09,
     -1.200000e-01, 5.000000e-03,  1.500000e-02,  -1.698000e-02, 1.448000e-04,
     -9.535000e-08,

     -3.322830e+01, 1.784000e-01,  5.550636e+02,  1.020400e+00,  2.499000e-03,
     -1.519000e-06, 7.833000e-01,  2.861000e-03,  -1.944000e-06, -4.400000e+00,
     3.024000e-02,  -3.283000e-05, 1.012000e-08,  1.506000e+00,  5.411000e-01,
     -3.800000e+00, 1.972000e-02,  -1.833000e-05, 4.938000e-09,  -1.200000e-01,
     5.000000e-03,  1.500000e-02,  1.083000e-02,  6.694000e-05,  -4.277000e-08,

     -1.819080e+01, 7.000000e-01,  1.146386e+02,  -5.101900e+00, 6.258000e-02,
     -1.672200e-04, -7.804000e-01, 7.173000e-03,  -5.578000e-06, -4.384000e+00,
     8.063000e-02,  -4.925000e-04, 1.042000e-06,
     1.506000e+00, //����� �������� ����������� n1
     5.515000e-01,  -6.683000e+00, 1.101400e-01,  -5.837500e-04, 1.083000e-06,
     -1.200000e-01, 2.500000e-02,  7.500000e-03,  -9.900000e-03, 8.212000e-05,
     3.125000e-09,

     -1.564080e+01, 7.754000e-01,  6.791620e+01,  -1.721000e-01, 5.756000e-03,
     -3.635000e-06, -7.540000e-01, 6.850000e-03,  -4.600000e-06, 1.279100e+00,
     -1.576000e-02, 6.499000e-05,  -5.145000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.515000e-01,  -2.162000e-01, 4.086000e-03,  1.270000e-06,  -1.587000e-09,
     -1.200000e-01, 2.500000e-02,  7.500000e-03,  -1.249000e-02, 1.111000e-04,
     -7.706000e-08,

     -3.277310e+01, 1.899000e-01,  5.842450e+02,  1.020400e+00,  2.499000e-03,
     -1.519000e-06, 7.250000e-01,  2.675000e-03,  -1.750000e-06, -4.400000e+00,
     3.024000e-02,  -3.283000e-05, 1.012000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.515000e-01,  -3.700000e+00, 1.783000e-02,  -1.506000e-05, 3.580000e-09,
     -1.200000e-01, 2.500000e-02,  7.500000e-03,  8.317000e-03,  4.837000e-05,
     -3.039000e-08,

     -1.852090e+01, 6.419000e-01,  1.159569e+02,  -5.101900e+00, 6.258000e-02,
     -1.672200e-04, -8.220000e-01, 8.330000e-03,  -1.233000e-05, 9.776000e-01,
     -2.570000e-02, 2.027000e-04,  -4.708000e-07,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -5.352000e+00, 8.615000e-02,  -4.437500e-04, 8.125000e-07,
     -1.000000e-01, 2.083000e-02,  6.251000e-03,  -7.680000e-03, 6.362000e-05,
     3.125000e-09,

     -1.522290e+01, 7.569000e-01,  5.581650e+01,  -1.721000e-01, 5.756000e-03,
     -3.635000e-06, -5.700000e-01, 5.250000e-03,  -3.000000e-06, 1.290300e+00,
     -1.547000e-02, 5.964000e-05,  -4.503000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -1.486000e-01, 3.263000e-03,  3.143000e-06,  -3.429000e-09,
     -1.000000e-01, 2.083000e-02,  6.250000e-03,  -7.879000e-03, 7.258000e-05,
     -3.658000e-08,

     -3.167150e+01, 2.265000e-01,  5.715408e+02,  1.020400e+00,  2.499000e-03,
     -1.519000e-06, 6.100000e-01,  2.343000e-03,  -1.433000e-06, -8.980000e+00,
     4.087000e-02,  -3.950000e-05, 1.123000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -3.700000e+00, 1.750000e-02,  -1.500000e-05, 3.704000e-09,
     -1.000000e-01, 2.083000e-02,  6.251000e-03,  4.667000e-03,  4.606000e-05,
     -2.722000e-08,

     -1.865220e+01, 6.124000e-01,  1.164154e+02,  -5.101900e+00, 6.258000e-02,
     -1.672200e-04, -7.376000e-01, 7.597000e-03,  -1.209000e-05, 9.776000e-01,
     -2.570000e-02, 2.027000e-04,  -4.708000e-07,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -4.799000e+00, 7.779200e-02,  -4.075000e-04, 7.708300e-07,
     -1.000000e-01, 2.750000e-02,  3.800000e-03,  -5.600000e-03, 4.667000e-05,
     0.000000e+00,

     -1.697520e+01, 6.736000e-01,  8.544400e+01,  -1.721000e-01, 5.756000e-03,
     -3.635000e-06, -4.760000e-01, 4.400000e-03,  -2.400000e-06, 1.290300e+00,
     -1.547000e-02, 5.964000e-05,  -4.503000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -1.495000e-01, 3.182000e-03,  2.825000e-06,  -3.365000e-09,
     -1.000000e-01, 2.750000e-02,  3.750000e-03,  -4.882000e-03, 4.692000e-05,
     -1.742000e-08,

     -2.975920e+01, 2.948000e-01,  5.283389e+02,  1.020400e+00,  2.499000e-03,
     -1.519000e-06, 9.333000e-02,  3.038000e-03,  -1.711000e-06, -8.980000e+00,
     4.087000e-02,  -3.950000e-05, 1.123000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -4.400000e+00, 1.981000e-02,  -1.806000e-05, 4.938000e-09,
     -1.000000e-01, 2.750000e-02,  3.750000e-03,  -1.333000e-02, 7.167000e-05,
     -3.518000e-08,

     -1.865860e+01, 6.038000e-01,  1.163531e+02,  -5.101900e+00, 6.258000e-02,
     -1.672200e-04, -3.150000e-01, 2.325000e-03,  2.500000e-06,  -5.632000e-01,
     5.743000e-03,  -9.250000e-06, 4.167000e-09,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -4.903000e+00, 8.270800e-02,  -4.587500e-04, 9.167000e-07,
     -1.200000e-01, 4.116000e-02,  1.433000e-03,  -4.963000e-03, 4.136000e-05,
     0.000000e+00,

     -1.730450e+01, 6.382000e-01,  8.195960e+01,  -1.721000e-01, 5.756000e-03,
     -3.635000e-06, -2.920000e-01, 2.800000e-03,  -8.000000e-07, 2.057000e-01,
     -2.912000e-03, 1.739000e-05,  -8.565000e-09,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -8.190000e-02, 2.358000e-03,  4.698000e-06,  -5.206000e-09,
     -1.300000e-01, 4.389000e-02,  1.821000e-03,  -5.195000e-03, 4.664000e-05,
     -2.164000e-08,

     -2.884630e+01, 3.140000e-01,  5.097230e+02,  1.020400e+00,  2.499000e-03,
     -1.519000e-06, -3.333000e-01, 3.522000e-03,  -1.889000e-06, -1.578000e+01,
     5.757000e-02,  -5.322000e-05, 1.512000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -3.600000e+00, 1.694000e-02,  -1.556000e-05, 4.321000e-09,
     -1.200000e-01, 4.116000e-02,  1.433000e-03,  -3.500000e-03, 4.317000e-05,
     -2.056000e-08,

     -1.864950e+01, 5.974000e-01,  1.162144e+02,  -5.101900e+00, 6.258000e-02,
     -1.672200e-04, -5.161000e-01, 5.341000e-03,  -8.672000e-06, -5.632000e-01,
     5.743000e-03,  -9.250000e-06, 4.167000e-09,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -5.115000e+00, 8.507500e-02,  -4.587500e-04, 8.750000e-07,
     -1.100000e-01, 3.810000e-02,  1.178000e-03,  -4.110000e-03, 3.463000e-05,
     -3.125000e-09,

     -1.826600e+01, 5.797000e-01,  1.009417e+02,  -1.721000e-01, 5.756000e-03,
     -3.635000e-06, -3.113000e-01, 2.839000e-03,  -1.089000e-06, 2.057000e-01,
     -2.911000e-03, 1.739000e-05,  -8.565000e-09,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -8.286000e-02, 2.278000e-03,  4.381000e-06,  -5.143000e-09,
     -1.100000e-01, 3.810000e-02,  1.178000e-03,  -5.017000e-03, 4.282000e-05,
     -2.132000e-08,

     -2.629940e+01, 3.817000e-01,  4.341220e+02,  1.020400e+00,  2.499000e-03,
     -1.519000e-06, -4.333000e-01, 3.522000e-03,  -1.889000e-06, -1.578000e+01,
     5.757000e-02,  -5.322000e-05, 1.512000e-08,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -3.600000e+00, 1.653000e-02,  -1.528000e-05, 4.321000e-09,
     -1.100000e-01, 3.810000e-02,  1.178000e-03,  -1.500000e-03, 3.250000e-05,
     -1.389000e-08,

     -1.870740e+01, 5.772000e-01,  1.163395e+02,  -5.101900e+00, 6.258000e-02,
     -1.672200e-04, -2.531000e-01, 1.929000e-03,  1.495000e-06,  4.842000e-01,
     -1.604000e-02, 1.405000e-04,  -3.375000e-07,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -3.137000e+00, 4.774200e-02,  -2.275000e-04, 3.958300e-07,
     -9.000000e-02, 3.118000e-02,  9.662000e-04,  -3.030000e-03, 2.532000e-05,
     -5.556000e-10,

     -1.927820e+01, 5.118000e-01,  1.165792e+02,  -1.721000e-01, 5.756000e-03,
     -3.635000e-06, -3.307000e-01, 2.878000e-03,  -1.378000e-06, 1.499000e-03,
     -2.399000e-04, 7.006000e-06,  -5.999000e-10,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  -2.048000e-01, 3.596000e-03,  -1.587000e-06, 3.175000e-10,
     -9.000000e-02, 3.117000e-02,  9.662000e-04,  -5.455000e-03, 4.273000e-05,
     -2.273000e-08,

     -2.366270e+01, 4.231000e-01,  3.364318e+02,  1.020400e+00,  2.499000e-03,
     -1.519000e-06, -1.750000e-01, 2.642000e-03,  -1.417000e-06, -9.750000e+00,
     3.383000e-02,  -2.694000e-05, 6.481000e-09,
     1.506000e+00, //����� �������� ����������� n1
     5.585000e-01,  1.000000e-01,  2.639000e-03,  -2.778000e-07, -6.173000e-10,
     -9.000000e-02, 3.118000e-02,  9.662000e-04,  -8.333000e-04, 2.817000e-05,
     -1.129000e-08}
/*% ***************                       */
,
       AD[38] = {-0.028, -0.045, -0.047, -0.035, -0.011, 0.022,  0.057,  0.090,
                 0.114,  0.125,  0.118,  0.096,  0.060,  0.013,  -0.037, -0.086,
                 -0.128, -0.162, -0.185, -0.199, -0.202, -0.193, -0.173, -0.140,
                 -0.096, -0.042, 0.015,  0.070,  0.115,  0.144,  0.155,  0.145,
                 0.120,  0.084,  0.044,  0.006,  -0.023, -0.040};

int MF0[7] = {75, 100, 125, 150, 175, 200, 250};

/*########################## ������ ATM   ###########################*/
/* ======================================================================
����������: ������ ��������� ��������� �� ������ ���� 25645.115-84 (version 90)
           ��� ����� "���������� �������".
     ����� ���������: ATM(H1,F135,F,F0,KP,T0).
������� ����������:
H1 - �������������� ������ (��.),
F135 - ������� �� 135 ����� �������� ������� ��������� ����������
      (10e-22*��/(�*�*��)),
F - �������������� �������� ������� (10e-22*��/(�*�*��)),
F0 - ������������� �������� ������� (10e-22*��/(�*�*��)),
KP - ������ ������������ �������������,
T0 - ����� ������� �� 0-�� ������ 1958 �. ( �����),
AT1[1:75] - ������ �������� ��������� ( � ���*����),
AD[1:38] - ������ ������������ ������������ ������� .
�������� ����������:
 RO - ��������� ��������� (��/(�*�*��)).
����������� �������:
������� ����������:
�=365.88115, F135=75, F=75.5, F0=75, KP=3.3, T0=9938.95330996.
------------------------------------------------------------------- */
/* ======================================================================
   Destination: Calculation of atmosphere density value without bulge.
   Call: RO=ATM(H1,F81,F,F0,KP,T0).
   Input data: H1 is altitude in km,
               F81 is is the weighting-averaged value of solar activity index
                  for preceding 81 days,
               F is daily averaged values of solar activity index,
               F0 is model's fixed level of solar activity index,
               KP is daily averaged planetary index of geomagnetic activity,
               T0 is time in days count of from 0h UTC Dec 31 1957.
   Output data: RO is atmospheric density value in kg/(m*m*km).
  ------------------------------------------------------------------- */

;
/**
 * @brief Calculation of atmosphere density value without bulge.
 * @param  H1  altitude in km.
 * @param  F81 weighting-averaged value of solar activity index for preceding 81
 * days.
 * @param  F    daily averaged values of solar activity index.
 * @param  F0  model's fixed level of solar activity index
 * @param  KP  daily averaged planetary index of geomagnetic activity.
 * @param  T0  time in days count of from 0h UTC Dec 31 1957.
 * @return    double atmospheric density value in kg(m*m*km)
 */
double JVC_DLL ATM(double H1, double F81, double F, int F0, double KP,
                   double T0)
//		; double ATM(double H1,double F81,double F,int F0,double
// KP,double  T0)
{
  double DJ1;
  double E, E1, E2, E3, E4, E7, E8, E9;
  int I, DJ0;
  ;
  if (H1 < 120.0) {
    if (H1 < 20) {
      RO = 1228.0;
      E2 = 0.090764;
      E3 = -0.0020452;
      E4 = 0;
    } else if (H1 < 60) {
      RO = 90.13;
      E2 = 0.16739;
      E3 = 0.00062669;
      E4 = 20;
    } else if (H1 < 100) {
      RO = 0.31043;
      E2 = 0.12378;
      E3 = -8.700E-4;
      E4 = 60;
    } else {
      RO = 0.53675E-3;
      E2 = 0.17527;
      E3 = 0.1287E-2;
      E4 = 100;
    };
    RO *= exp(-E2 * (H1 - E4) + E3 * sqr(H1 - E4));
  } else {
    E1 = (E = H1 * H1) * H1;
    RO = exp(AT1[0] - AT1[1] * sqrt(H1 - AT1[2])) *
         (1 + (AT1[22] + AT1[23] * H1 + AT1[24] * E) * (F135 - F0));
    DJ1 = T0 + 0.125;
    I = DJ1;
    I /= 365;
    DJ0 = I * 365 + (I + 1) / 4;
    DJ1 -= DJ0;
    if (DJ1 > 365)
      DJ1 -= 365;
    if (DJ1 < 0)
      DJ1 += 365;
    I = DJ1 / 10;
    E9 = AD[I] + (DJ1 - 10 * I) * 0.1 * (AD[I + 1] - AD[I]);
    E7 = AT1[3] + AT1[4] * H1 + AT1[5] * E;
    E8 = AT1[6] + AT1[7] * H1 + AT1[8] * E;
    RO *= (1 + E7 * E9) * (1 + E8 * (F - F81) / F);
    E9 = AT1[15] + AT1[16] * H1 + AT1[17] * E + AT1[18] * E1;
    RO *= (1 + E9 * (AT1[19] + AT1[20] * KP + AT1[21] * KP * KP)) * 9.80665E3;
  }
  return (RO);
}

/*
 double
AT[450]={-18.299,0.7009,115.3429,-2.6122,0.02935,-0.6318E-4,-0.15661,-0.8904E-3,0.18296E-4,
     1.2249,
       -0.04205,0.3833E-3,-0.9827E-6,3,0.541052,-3.9877,0.03831,0.1059E-4,-0.4408E-6,-0.17,0.04746,0.006107,-
       0.01188,0.9983E-4,0,-16.1768,0.8016,86.3329,-0.7379,0.008524,-0.5328E-5,-0.9068,0.007446,-0.4864E-5,3.4183
       ,-0.03872,0.1432E-3,-0.1274E-6,3,0.541052,3.6776,-0.03663,0.1334E-3,-0.1247E-6,-0.17,0.04746,0.006107,
      -0.8571E-2,0.7952E-4,0.7936E-8,-31.9086,0.2336,491.2201,1.3952,0.2956E-2,-0.1968E-5,2.95,-0.2122E-2,
       0.3704E-6,38.0659,-0.1029,0.9405E-4,-0.2764E-7,3,0.541052,24.5361,-0.06672,0.6173E-4,-1.818E-8,-0.17,
       0.04746,0.6107E-2,0.142,-2.289E-4,0.10368E-6,-18.1908,0.7,114.6386,-2.6122,0.02935,-0.6318E-4,0.28435,-
       0.6976E-2,0.3839E-4,-1.0495,0.6964E-2,0.3783E-4,-0.1915E-6,3,0.551524,-1.722,-0.21025E-2,0.2415E-3,-
       0.8699E-6,-0.16,0.04933,0.004,-0.008798,0.7332E-4,0,-16.0652,0.7675,77.1052,-0.7379,0.008524,-0.5328E-5,-
       0.7321,0.006244,-0.3681E-5,2.012,-0.02331,0.9072E-4,-0.759E-7,3,0.551524,2.6952,-0.2571E-1,0.9481E-4,-
       0.8389E-7,-0.16,0.04933,0.004,-0.8114E-2,0.6924E-4,0.1587E-8,-31.0887,0.2417,490.7284,1.3952,0.2956E-2,-
       0.1968E-5,0.39,0.3473E-2,-0.2178E-5,18.1358,-0.036,0.2471E-4,-0.5258E-8,3,0.551524,21.462,-0.05357,
       0.4649E-4,-0.1287E-7,-0.16,0.04933,0.004,0.118,-0.1933E-3,0.8889E-7,-18.5209,0.6419,115.9569,-2.6122,
       0.02935,-0.6318E-4,0.3204,-0.6899E-2,0.3525E-4,0.1608,-0.01437,0.159E-3,-0.42E-6,3,0.558505,-1.4436,-
       0.3861E-2,0.2271E-3,-0.7887E-6,-0.14,0.04507,0.002786,-0.007196,0.5997E-4,0,-15.9559,0.7362,70.5386,-
       0.7379,0.008524,-0.5328E-5,-0.8203,0.6628E-2,-0.4695E-5,1.544,-0.01795,0.7066E-4,-0.56E-7,3,0.558505,
       1.9076,-0.01738,0.669E-4,-0.5614E-7,-0.14,0.04507,0.2786E-2,-0.6E-2,0.5332E-4,0,-29.9882,0.2654,479.9537,
       1.3952,0.2956E-2,-0.1968E-5,0.1238,0.3439E-2,-0.1992E-5,-3.642,0.323E-1,-0.412E-4,0.1478E-7,3,0.558505,
       17.75,-0.04072,0.3333E-4,-0.8704E-8,-0.14,0.04507,0.2786E-2,0.4E-2,0.5224E-4,-0.2592E-7,-18.6522,0.6124,
       116.4154,-2.6122,0.02935,-0.6318E-4,0.58624,-0.10503E-1,0.4682E-4,1.0526,-0.03015,0.2517E-3,-0.6128E-6,3.2
       ,0.558505,0.1718,-0.03416,0.4124E-3,-0.1164E-5,-0.13,0.04389,0.001821,-0.0056,0.4667E-4,0,-16.746,0.6805,
       80.4406,-0.7379,0.008524,-0.5328E-5,-0.599,0.5078E-2,-0.316E-5,1.466,-0.01629,0.5941E-4,-0.432E-7,3.2,
       0.558505,1.5609,-0.01323,0.5113E-4,-0.4051E-7,-0.13,0.04389,0.1821E-2,-0.6968E-2,0.5684E-4,-0.1428E-7,-
       28.8362,0.2911,461.9691,1.3952,0.2956E-2,-0.1968E-5,-0.45,0.4289E-2,-0.2259E-5,-20.17,0.082,-0.8718E-4,
       0.2824E-7,3.2,0.558505,-11.456,0.05307,-0.5963E-4,0.202E-7,-0.13,0.04389,0.1821E-2,-0.014,0.8668E-4,-
       0.4444E-7,-18.6495,0.5974,116.2144,-2.6122,0.02935,-0.6318E-4,0.4325,-0.78E-2,0.35E-4,-1.707,+0.03118,-
       0.1926E-3,0.4277E-6,3.2,0.558505,1.1137,-0.05049,0.5047E-3,-0.1344E-5,-0.11,0.0381,0.001178,-0.0048,0.4E-4
       ,0,-18.266,0.5797,100.9417,-0.7379,0.008524,-0.5328E-5,-0.4309,0.353E-2,-0.1326E-5,-0.3247,0.2102E-2,
       0.2471E-5,0.4509E-8,3.2,0.558505,1.5767,-0.01334,0.4906E-4,-0.3848E-7,-0.11,0.0381,0.1178E-2,-0.4629E-2,
       0.4048E-4,-0.7936E-8,-26.2658,0.3492,399.0603,1.3952,0.2956E-2,-0.1968E-5,-0.553,0.4163E-2,-0.204E-5,-
       37.009,0.128,-0.1259E-3,0.3842E-7,3.2,0.558505,-15.5637,0.06105,-0.6206E-4,0.1941E-7,-0.11,0.0381,
       0.1178E-2,-0.022,0.9176E-4,-0.452E-7,-18.7074,0.5772,116.3395,-2.6122,0.02935,-0.6318E-4,0.415,-0.7108E-2,
       0.3042E-4,0.2253,-0.6403E-2,0.4728E-4,-0.7963E-7,3.2,0.558505,-3.6114,0.04682,-0.1505E-3,0.9242E-7,-0.09,
       0.03118,0.9662E-3,-0.321E-2,0.2675E-4,0,-19.2125,0.5095,115.2277,-0.7379,0.008524,-0.5328E-5,-0.3518,
       0.2727E-2,-0.5682E-6,-0.2987,0.2074E-2,0.1273E-5,0.4149E-8,3.2,0.558505,0.78,-0.5567E-2
       ,0.25E-4,-0.1852E-7,-0.09,0.03117,0.9662E-3,-0.3714E-2,0.3124E-4,-0.9524E-8,-23.0783,0.414,284.6955,1.3952
       ,0.2956E-2,-0.1968E-5,-1.1467,0.5124E-2,-0.2356E-5,-32.8571,0.107,-0.9706E-4,0.2734E-7,3.2,0.558505,-18.45
       ,0.06637,-0.641E-4,0.192E-7,-0.09,0.03117,0.9662E-3,-0.032,0.10244E-3,-0.4964E-7}

          ,
AD[38]={-0.06,-0.075,-0.077,-6.5E-2,-0.04,0,0.04,0.077,0.11,0.122,0.122,0.1,0.065,2.5E-2,
          -0.01,-0.05,-8.2E-2,-0.11,-0.14,-0.157,-0.166,-0.17,-0.16,-0.14,-0.11,-0.055,0,6.2E-2,1.22E-1,0.161,0.17,0.156,
       0.12,0.073,0.027,-0.02,-0.045,-0.07};
int MF0[6]={75,100,125,150,200,250} ;
*/
/*########################## ������ ATM   ###########################*/
/* ======================================================================
 ����������: ������ ��������� ��������� �� ������ ���� 25645.115-84 (version 84)
             ��� ����� "���������� �������".
       ����� ���������: ATM(H1,F135,F,F0,KP,T0).
 ������� ����������:
  H1 - �������������� ������ (��.),
  F135 - ������� �� 135 ����� �������� ������� ��������� ����������
        (10e-22*��/(�*�*��)),
  F - �������������� �������� ������� (10e-22*��/(�*�*��)),
  F0 - ������������� �������� ������� (10e-22*��/(�*�*��)),
  KP - ������ ������������ �������������,
  T0 - ����� ������� �� 0-�� ������ 1958 �. ( �����),
  AT1[1:75] - ������ �������� ��������� ( � ���*����),
  AD[1:38] - ������ ������������ ������������ ������� .
 �������� ����������:
   RO - ��������� ��������� (��/(�*�*��)).
 ����������� �������:
 ������� ����������:
  �=365.88115, F135=75, F=75.5, F0=75, KP=3.3, T0=9938.95330996.
 �������� ����������:
  RO=1.82115027151e-9.
 ------------------------------------------------------------------- */
/*
        ; double ATM(double H1,double F135,double F,int F0,double KP,double T0)
                {
                  double  DJ1 ;
                  double
                         E
                      ,  E1
                      ,  E2
                      ,  E3
                      ,  E4
                      ,  E7
                      ,  E8
              ,  EPS=0
              ,  E9;
                    int  I
                      ,  DJ0;

                  if  (H1<185.0&&H1>85.0)
                    {  if  (H1<90)
                           {
                              E2=0.0
                            ; E3=0.055
                            ; E4=85.0;
                           }
                        else if  (H1<95)
                           {
                               E2=0.055
                            ;  E3=0.071
                            ;  E4=90;
                           }
                        else if  (H1<100)
                           {
                               E2=0.071
                            ;  E3=0.039
                            ;  E4=95;
                           }
                        else if  (H1<105)
                           {
                               E2=0.039
                            ;  E3=0.117
                            ;  E4=100;
                           }
                        else if  (H1<110)
                           {
                               E2=0.117
                            ;  E3=0.219
                            ;  E4=105;
                           }
                        else if  (H1<115)
                           {
                               E2=0.219
                            ;  E3=0.309
                            ;  E4=110;
                           }
                        else if  (H1<120)
                           {
                               E2=0.309
                            ;  E3=0.275
                            ;  E4=115;
                           }
                        else if  (H1<125)
                           {
                               E2=0.275
                            ;  E3=0.197
                            ;  E4=120;
                           }
                        else if  (H1<130)
                           {
                            ;  E2=0.197
                            ;  E3=0.169
                            ;  E4=125;
                           }
                        else if  (H1<135)
                           {
                               E2=0.169
                            ;  E3=0.152
                            ;  E4=130;
                           }
                        else if  (H1<140)
                           {
                               E2=0.152
                            ;  E3=0.143
                            ;  E4=135;
                           }
                        else if  (H1<145)
                           {
                               E2=0.143
                            ;  E3=0.136
                            ;  E4=140;
                           }
                        else if  (H1<150)
                           {
                               E2=0.136
                            ;  E3=0.131
                            ;  E4=145;
                           }
                        else if  (H1<155)
                           {
                               E2=0.131
                            ;  E3=0.118
                            ;  E4=150;
                           }
                        else if  (H1<160)
                           {
                               E2=0.118
                            ;  E3=0.115
                            ;  E4=155;
                           }
                        else if  (H1<165)
                           {
                               E2=0.115
                            ;  E3=0.105
                            ;  E4=160;
                           }
                        else if  (H1<170)
                           {
                               E2=0.105
                            ;  E3=0.097
                            ;  E4=165;
                           }
                        else if  (H1<175)
                           {
                               E2=0.097
                            ;  E3=0.043
                            ;  E4=170
                            ;
                           }
                        else if  (H1<180)
                           {
                               E2=0.043
                            ;  E3=0.015
                            ;  E4=175
                            ;
                           }
                        else
                           {
                               E2=0.015
                            ;  E3=0.0
                            ;  E4=180
                            ;
                           }
                      EPS=E2+0.2*(E3-E2)*(H1-E4);
                    }

      if  (H1<120.0)
                    {  if  (H1<20)
                           {   RO=1228.0
                            ;  E2=0.090764
                            ;  E3=-0.0020452
                            ;  E4=0;
                           }
                        else if  (H1<60)
                                {   RO=90.13
                                 ;  E2=0.16739
                                 ;  E3=0.00062669
                                 ;  E4=20;
                                }
                        else if  (H1<100)
                                {  RO=0.31043
                                 ;  E2=0.12378
                            ;  E3=-8.700E-4
                                 ;  E4=60;
                                }
                        else
                                {  RO=0.53675E-3
                                 ;  E2=0.17527
                                 ;  E3=0.1287E-2
                                 ;  E4=100;
                           }
          ;  RO*=(1.0-EPS)*exp(-E2*(H1-E4)+E3*sqr(H1-E4));
                    }   else
                         {
                         E1=(E=H1*H1)*H1
                      ;
   RO=exp(AT1[0]-AT1[1]*sqrt(H1-AT1[2]))*(1+(AT1[22]+AT1[23]*H1+AT1[24]*E)*(F135-F0))
                      ;  DJ1=T0+0.125
                                ;  I=DJ1
                      ;  I/=365
                      ;  DJ0=I*365+(I+1)/4
                      ;  DJ1-=DJ0
                                ;  if (DJ1>365) DJ1-=365
                      ;  if (DJ1<0) DJ1+=365
                                ;  I=DJ1/10
                                ;  E9=AD[I]+(DJ1-10*I)*0.1*(AD[I+1]-AD[I])
                                ;  E7=AT1[3]+AT1[4]*H1+AT1[5]*E
                      ;  E8=AT1[6]+AT1[7]*H1+AT1[8]*E
                      ;  RO*=(1+E7*E9)*(1+E8*(F-F135)/F135)
                                ;  E9=AT1[15]+AT1[16]*H1+AT1[17]*E+AT1[18]*E1
        ;  RO*=(1.0-EPS)*(1+E9*(AT1[19]+AT1[20]*KP+AT1[21]*KP*KP))*9.80665E3
                      ;
                    }
                    return(RO);
                 }
*/
/*########################## ����� TAGOCT ###########################*/
/*########################## ������ BESSEL0   ###########################*/
/* ======================================================================
   ����������: ���������� �������� ���������������� ������� ������� �������
               ���� �������� ������� �� ������� ���������.
   ����� ���������: I0(X).
   I0(X) - �������� ����������-��������.
   ����������� �������:
         ������� ����������:
    X=.232361931547.
   �������� ����������:
    I0=1.01354363652.
   ------------------------------------------------------------------- */
/* ======================================================================
   Destination: Calculation of Bessel function I0.
   Call: Y=I0(X).
   Input data: X.
   Output data: Y is value of Bessel function I0.
  ------------------------------------------------------------------- */

;
/**
 * @brief Calculation of Bessel function I0.
 * @param  X
 * @return   double value of Bessel function I0.
 */
double JVC_DLL I0(double X)
//        ;  double I0(double X)
{
  double T;

  if (X <= 10) {
    T = 0.01 * X * X;
    T = (((((((((((0.668831515 * T + 0.1439945401) * T + 9.6664625536) * T +
                 26.2421075769) *
                    T +
                95.9019523737) *
                   T +
               239.2391892805) *
                  T +
              471.3119673418) *
                 T +
             678.08454241) *
                T +
            678.1808514731) *
               T +
           434.026680086) *
              T +
          156.2500502179) *
             T +
         24.9999990983) *
            T +
        1.0000000027;
  } else {
    T = 10 / X;
    T = exp(X) *
        ((((((0.107871E-5 * T - 0.72753E-6) * T + 0.580397E-5) * T +
            0.286977E-4) *
               T +
           0.2806006E-3) *
              T +
          0.00498677224) *
             T +
         0.39894228047) /
        sqrt(X);
  }
  return (T);
}
/*########################## ����� BESSEL0 ###########################*/
/*########################## ������ BESSEL1   ###########################*/
/* ======================================================================
   ����������: ���������� �������� ���������������� ������� ������� �������
               ���� ������� ������� �� ������� ���������.
   ����� ���������: I1(X).
   I1(X) - �������� ����������-��������.
   ����������� �������:
   ������� ����������:
    X=.232361931547.
   �������� ����������:
    I1=0.11696683807.
   ------------------------------------------------------------------- */
/* ======================================================================
   Destination: Calculation of Bessel function I1.
   Call: Y=I1(X).
   Input data: X.
   Output data: Y is value of Bessel function I1.
  ------------------------------------------------------------------- */

;
/**
 * @brief Calculation of Bessel function I1.
 * @param  X
 * @return   double value of Bessel function I1.
 */
double JVC_DLL I1(double X)
//        ;  double I1(double X)
{
  double T;
  if (X <= 10) {
    T = 0.01 * X * X;
    T = X *
        ((((((((((((0.02411174292 * T + 0.01604459097) * T + 0.41187059355) *
                      T +
                  1.35566931981) *
                     T +
                 5.28380368332) *
                    T +
                14.98220505757) *
                   T +
               33.65156773642) *
                  T +
              56.51119244161) *
                 T +
             67.81726179804) *
                T +
            54.25343507056) *
               T +
           26.04166836566) *
              T +
          6.2499999695) *
             T +
         0.50000000009);
  } else {
    T = 10 / X;
    T = exp(X) *
        ((((((-0.123545E-5 * T + 0.74632E-6) * T - 0.726227E-5) * T -
            0.4031439E-4) *
               T -
           0.00046761773) *
              T -
          0.01496032838) *
             T +
         0.39894228032) /
        sqrt(X) * sign(X);
  }
  return (T);
}
/*########################## ����� BESSEL1 ###########################*/

/*########################## ������ K��� ###########################*/
/* ======================================================================
   ����������: ���������� �������������, ������������ ���������� ��
               ������������� ���������.
   ����� ���������: KATM(H,B1,B2).

   ������� ����������:
    � - ������� ������� (��), E -��������������,CI - ������� ����������,
                H - ������ ���������� ��������� (��), B1 - ��������, �����������
   ��- ��������� ������ ���������� ��������� �� ������ (1/��), B2 - ��������,
                ����������� ����������� �������� "���������� �������" ���������
   �� ������ (1/(��*��)), ��������������� �������������� ���������� GME,
                �������� �������� ����� CBZ, ����������� ��� ������  ���������
   �����- ���� C20, �������������� ������ ����� RZ. �������� ����������: ������
   �������������, ������������ ���������� �� ������������� ����- ����� KA[1:6,
   0:3]. ����������� �������: ������� ����������: �=6786., E=.01, CI=cos(PI/4),
   H=36, GME=398600., CBZ=6.30038738, C20=-.001082626848, RZ=6378.137. ��������
   ����������:  KA { 3.5770269e-1 -0.18407400e-2 3.78030160e-1 1.128539666e-2
    2.1078500e-1 -1.36138900e-2 2.05982779e-1 1.520200000e-2
    6.6528000e-3  1.80230000e-3
                 -1.27576000e-2 1.15962900e-1  1.15184560e-2
                  1.52232809e-3
                                5.79814500e-4  5.75922800e-5
                            }.
   ------------------------------------------------------------------- */
/* ======================================================================
         Destination: Calculation of slow-changing functions for
     taken into account pertubations due to atmospheric drag.
     (see formulas (1.2.10), (1.2.13) Topic 5.5, Report 1).
         Call: KATM(H,B1,B2).
     Input data: H is scale height,
     B1, B2 (see foormula (1.2.5) Topic 5.5, Report 1).
         It is used variables, declared as global:
         �(semimajor axis, km), E (eccentricity),
     SI2(sin(I/2)),CI(cos(I)), where I is inclination.
         Otput data: static array  KA[0:5,0:4].
   ------------------------------------------------------------------- */

;
/**
 * Calculation of slow-changing functions for taken into account pertubations
 * due to atmospheric drag.
 * @param  H  scale height.
 * @param  B1
 * @param  B2
 * @return    void
 */
void JVC_DLL KATM(double H, double B1, double B2)
//         ;void  KATM(double H,double B1,double B2)
{
  double E2, E3, E4, C, S2, J0, J1, J2, J3, J4, J5, J6, J7, WN, A1, X, Z, Z2,
      GR, PAB;
  ;
  E2 = E * E;
  E3 = E2 * E;
  E4 = E3 * E;
  Z = A * E / H;
  Z2 = Z * Z;
  S2 = 1 - CI * CI;
  WN = CBZ * A * sqrt(A / GME) / 86400.0;
  X = sqr(A * E) * B2;
  GR = 1 - WN * sqrt((1 - E) / (1 + E)) * (1 - E) * CI
      //���� ��� 	  ;  GR=1-WN*sqrt((1-E)/(1+E))*(1-E)/(1+E)*CI
      ;
  A1 = -0.375 * C20 * (3 * S2 - 2) * E * RZ * RZ / (A * (1 - E2) * H);
  Z = Z + A1;
  Z2 = Z * Z;
  A1 = B1 * A * E;
  if (Z < 20) {
    C = GR * GR * exp(-Z);
    J0 = C * I0(Z);
    J1 = C * I1(Z);
    J2 = J0 - 2 * J1 / Z;
    J3 = J1 - 4 * J2 / Z;
    J4 = J2 - 6 * J3 / Z;
    J5 = J3 - 8 * J4 / Z;
    J6 = J4 - 10 * J5 / Z;
    KA[0][0] =
        0.25 *
        ((4 + 3 * E2 + 21.0 / 16.0 * E4 + WN * WN * S2 + 12 * WN * E2 * CI) *
             J0 +
         (8 * E + 3 * E3 + 16 * E * WN * CI) * J1 +
         (3 * E2 + 7.0 / 4.0 * E4 + 8 * WN * E2 * CI) * J2 + E3 * J3 +
         4 * A1 * (J0 - J1 + E * (-J0 + 2 * J1 - J2)) +
         2 * X * (3 * J0 - 4 * J1 + J2 - E * (4 * J0 - 7 * J1 + 4 * J2 - J3)));
    KA[0][1] = 0.125 * (4 * E * J0 + 8 * J1 + 12 * E * J2 +
                        X * (-8 * J0 + 14 * J1 - 8 * J2 + 2 * J3 +
                             E * (9 * J0 - 20 * (J1 - J2) - 12 * J3 + 3 * J4)));
    KA[0][2] = 0.25 * (4 * J2 + 8 * E * J3 +
                       X * (J0 - 4 * J1 + 6 * J2 - 4 * J3 + J4 +
                            E * (2 * J1 - 8 * J2 + 12 * J3 - 8 * J4 + 2 * J5)))
        // ����		    ;  KA[0][3]=-0.5*E*J2+J3+2.5*E*J4
        // ����		    ;  KA[0][4]=J4+2*E*J5
        ;
    KA[0][3] = -0.5 * E * J2 + J3 + 2.5 * E * J4 +
               0.125 * X * (2 * J1 - 8 * J2 + 12 * J3 - 8 * J4 + 2 * J5);
    KA[0][4] = J4 + 2 * E * J5 + 0.5 * X * (J2 - 4 * J3 + 6 * J4 - 4 * J5 + J6);
    KA[1][0] =
        0.0625 * (1 - E2) *
        ((8 * E + 3 * E3 + 36 * E * WN * CI + X * (14 * E - 16)) * J0 +
         (16 + 6 * E2 + X * (28 - 24 * E)) * J1 +
         (8 * E + 4 * E3 + 28 * E * WN * CI + 16 * X * (E - 1)) * J2 +
         4 * A1 *
             (-2 * J0 + 4 * J1 - 2 * J2 + E * (2 * J0 - 3 * J1 + 2 * J2 - J3)) +
         (2 * E2 + 4 * X * (1 - 2 * E)) * J3 + 2 * X * E * J4);
    KA[1][1] = 0.0625 * (1 - E2) *
               ((8 + X * (14 - 8 * E)) * J0 + (8 * E - X * (24 - 16 * E)) * J1 +
                (8 + 16 * X * (1 - E)) * J2 + (8 * E - X * (8 - 14 * E)) * J3 +
                X * (2 - 8 * E) * J4 + 2 * X * E * J5);
    KA[1][2] =
        0.0625 * (1 - E2) *
        (-4 * E * J0 + 8 * (J1 + J3) + 12 * E * J4 -
         X * (8 * J0 - 16 * (J1 - J2) - 14 * J3 + 8 * J4 - 2 * J5 +
              E * (4 * J0 - 13 * J2 + 20 * J3 - 20 * J4 + 12 * J5 - 3 * J6)));
    KA[1][3] = 0.125 * (1 - E2) *
               (-4 * E * J1 + 4 * J2 + 4 * E * J3 + 4 * J4 + 8 * E * J5)
        // ����			    ;
        // KA[1][4]=0.0625*(1-E2)*(4*(J3-J5)-8*E*(J2-J6))
        ;
    KA[1][4] = 0.125 * (1 - E2) * (4 * (J3 - J5) - 8 * E * (J2 - J6));
    KA[2][0] =
        1 / (32 * GR) * CI2 * SI * WN * ((4 + E2) * J0 - 8 * E * J1 + E2 * J2);
    KA[2][1] = 1 / (64 * GR) * CI2 * SI * WN *
               (11 * E2 * J0 - 16 * E * J1 + (8 - 6 * E2) * J2);
    KA[2][2] = 0;
    KA[2][3] = 0;
    KA[2][4] = 0;
    PAB = 1 / (32 * GR) * WN * SI2 *
          (11 * E2 * J0 - 16 * E * J1 + (8 - 6 * E2) * J2);
    KA[3][0] = 2 * E * GR * SI2 * PAB;
    KA[3][1] = 0.125 * (-4 * J0 - 6 * E * J1 + 4 * J2 + 6 * E * J3 -
                        X * (5 * J0 - 4 * (J1 + J2 - J3) - J4));
    KA[3][2] = 0.125 * (-4 * J1 - 8 * E * J2 + 4 * J3 + 8 * E * J4 +
                        X * (4 * (J0 - J4) - 6 * J1 + 5 * J3 + J5));
    KA[3][3] =
        0.125 * (2 * E * J1 - 4 * J2 - 12 * E * J3 + 4 * J4 + 10 * E * J5);
    KA[3][4] =
        0.125 * (4 * E * J2 - 4 * J3 - 16 * E * J4 + 4 * J5 + 12 * E * J6);
    KA[4][0] = PAB;
    KA[4][1] = 0;
    KA[4][2] = 0;
    KA[4][3] = 0;
    KA[4][4] = 0;
    KA[5][0] = 2 * GR * SI2 * PAB;
    KA[5][1] = 0.375 * E * (J0 - J2);
    KA[5][2] = 0.375 * E * (J1 - J3);
    KA[5][3] = 0.375 * E * (J2 - J4);
    KA[5][4] = 0.375 * E * (J3 - J5)

        ;
  } else {
    C = GR * GR / sqrt(PI2 * Z);
    J1 = (1 - 8 * E + 3 * E2) / (4 * (1 - E2));
    J2 = (3 - 16 * E + 50 * E2 + 16 * E3 - 5 * E4) / (32 * sqr(1 - E2));
    J3 = (5 - 24 * E + 45 * E2 - 80 * E3 - 249 * E4 - 24 * E4 * E +
          7 * E3 * E3) /
         (128 * sqr(1 - E2) * (1 - E2))
        //���� ���	    ;  J0=C*sqrt((1+E)/(1-E))*(1+E)/(1-E)
        ;
    J0 = C * sqrt((1 + E) / (1 - E)) * (1 + E);
    J4 = -(1 + E) / (1 - E);
    J5 = -E * J4 / (1 - E);
    J6 = 4 * J4;
    J7 = 2 * (1 + 4 * E + 3 * E * E) / sqr(1 - E);
    KA[0][0] = J0 * (1 + 0.5 * J1 / Z + 0.75 * J2 / Z2 + 1.375 * J3 / (Z * Z2) +
                     0.75 * X * (1 + 2.5 * J1 / Z) / Z2 +
                     0.5 * A1 * (1 + 1.5 * J1 / Z) / Z);
    KA[0][1] = J0 * (1 + 0.5 * (J1 + J4) / Z + 0.75 * (J5 + J4 * J1 + J2) / Z2);
    KA[0][2] = J0 * (1 + 0.5 * (J1 + J6) / Z + 0.75 * (J7 + J6 * J1 + J2) / Z2);
    KA[0][3] = 0;
    KA[0][4] = 0;
    ;
    J1 = -(3 + 4 * E - 3 * E2) / (4 * (1 - E2));
    J2 = (-5 + 24 * E + 26 * E2 + 8 * E3 - 5 * E4) / (32 * sqr(1 - E2));
    J3 = -(7 - 20 * E + 27 * E2 + 200 * E3 + 101 * E4 + 12 * E4 * E -
           7 * E3 * E3) /
         (128 * sqr(1 - E2) * (1 - E2));
    J0 = J0 * (1 - E);
    KA[1][0] = J0 * (1 + 0.5 * J1 / Z + 0.75 * J2 / Z2 + 1.375 * J3 / (Z * Z2) +
                     0.75 * X * (1 + 2.5 * J1 / Z) / Z2 +
                     0.5 * A1 * (1 + 1.5 * J1 / Z) / Z);
    KA[1][1] = J0 * (1 + 0.5 * (J1 + J4) / Z + 0.75 * (J5 + J4 * J1 + J2) / Z2);
    KA[1][2] = J0 * (1 + 0.5 * (J1 + J6) / Z + 0.75 * (J7 + J6 * J1 + J2) / Z2);
    KA[1][3] = KA[1][4] = 0;
    J1 = (1 + 8 * E + 11 * E2) / (4 * (1 - E2));
    J0 = 0.125 * C * CI2 * sqr(1 - E) * SI * WN / GR;
    KA[2][0] = J0 * (1 + 0.5 * J1 / Z + 0.75 * X * (1 + 2.5 * J1 / Z) / Z2);
    J1 = -(15 + 24 * E + 5 * E2) / (4 * (1 - E2))
        //����			    ;  J2=(1+0.5*J1/Z+0.75*X*(1+2.5*J1/Z)/Z2)
        ;
    J2 = (1 + 0.5 * J1 / Z);
    KA[2][1] = J0 * J2;
    KA[2][2] = KA[2][3] = KA[2][4] = 0
        // J2 ��� ����� �������� �� ��������� ????
        ;
    J2 = 0.125 * C * WN * sqr(1 - E) * J2 * SI2;
    KA[3][0] = 2 * E * SI2 * J2;
    ;
    KA[3][1] = KA[3][2] = KA[3][3] = KA[3][4] = 0;
    KA[4][0] = J2 / GR;
    KA[4][1] = KA[4][2] = KA[4][3] = KA[4][4] = 0;
    KA[5][0] = 2 * SI2 * J2;
    KA[5][1] = KA[5][2] = KA[5][3] = KA[5][4] = 0;
  };
}
/*########################## ����� KATM ###########################*/
/*########################## ������ F107kp ###########################*/
;
/**
 * @brief do nothing.
 * @param  t    [description]
 * @param  f107 [description]
 * @param  f81  [description]
 * @param  kp   [description]
 * @return      always -1.
 */
int JVC_DLL F107KP(double t, double *f107, double *f81, double *kp)
//		;  int F107KP(double t,double *f107,double *f81,double *kp)
{
  double ti, ti1, fi1, ki1, fai1, seek;
  int db, i_bd_fkp = -1;
  /*
            if(i_bd_fkp==-1) return(i_bd_fkp);

                          if(i_bd_fkp<-1)
                           {
                            if((i_bd_fkp=d4use("f107kp"))<0)
                                   {i_bd_fkp=-1;
                                    return(i_bd_fkp);}
                                   if((i4index("jd1957","JD1957",0,0))<0)
                                   {i_bd_fkp=-1;
                                    return(i_bd_fkp);}
                           }
                           db=d4select(i_bd_fkp);
                          if (d4reccount()>0)
                            {
  //��� F10.7 ������ � �� ���������� �� 20.00 UT, �.�. t+5/6    � ������������
  ���������� 1.7 ����� seek=t+5.0/6.0-1.7; if(d4seek((char*)&seek)!=3)
                                   {
                                    ti = f4value(f4ref("JD1957"));
                                    *f107=f4value(f4ref("f107"));
                                    *f81 =f4value(f4ref("f81"));
                                    if (d4recno()>1) d4skip(-1);
                                          ti1= f4value(f4ref("JD1957"));
                                          if (ti!=ti1)
                                                  {
                                                    fi1 =f4value(f4ref("f107"));
                                                    fai1=f4value(f4ref("f81"));
                                                    *f107+=(fi1-*f107)/(ti1-ti)*(t-ti);
                                                    *f81
  +=(fai1-*f81)/(ti1-ti)*(t-ti);
                                                  }

                                   }
  //��� ������������� KP ������ � �� ���������� �� 12.00 UT, �.�. t+0.5  �
  ������������ ���������� 0.6 ����� seek=t-0.1; if(d4seek((char*)&seek)!=3)
                                   {
                                    ti = f4value(f4ref("JD1957"));
                                    *kp  =f4value(f4ref("kp"));
                                    if (d4recno()>1) d4skip(-1);
                                          ti1= f4value(f4ref("JD1957"));
                                          if (ti!=ti1)
                                                  {
                                                    ki1 =f4value(f4ref("kp"));
                                                    *kp
  +=(ki1-*kp)/(ti1-ti)*(t-ti);
                                                  }

                                   }

                            }else i_bd_fkp=-1;
                                  d4select(db);
  */
  return (i_bd_fkp);
}
/*########################## ����� F107KP ###########################*/

#endif
