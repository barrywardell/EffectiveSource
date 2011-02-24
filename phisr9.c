#include <math.h>

static int I0_0;
static int I0_1;
static int I0_2;
static int I0_3;
static int I0_4;
static int I0_5;
static int I0_6;
static int I0_7;
static int I0_8;
static int I0_9;
static int I0_10;
static double R0_17;
static double R0_18;
static double R0_20;
static double R0_22;
static double R0_24;
static double R0_26;
static double R0_28;
static double R0_31;
static double R0_32;
static double R0_33;
static double R0_36;
static double R0_37;
static double R0_38;
static double R0_40;
static double R0_41;
static double R0_42;
static double R0_44;
static double R0_47;
static double R0_48;
static double R0_49;
static double R0_50;
static double R0_51;
static double R0_52;
static double R0_53;
static double R0_54;
static double R0_55;
static double R0_56;
static double R0_57;
static double R0_58;
static double R0_59;
static double R0_60;
static double R0_61;
static double R0_62;
static double R0_64;
static int initialize = 1;

#include "phisr9.h"

int Initialize_phisr9()
{
  if(initialize)
  {
    R0_37 = (double) 0.04743083003952569;
    R0_41 = (double) 102.57470363819434;
    R0_64 = (double) 1.5;
    R0_51 = (double) 15.540224983755685;
    I0_9 = (int) 81;
    I0_1 = (int) -16;
    R0_58 = (double) 0.021573322379181633;
    R0_36 = (double) 30.73517786561265;
    R0_40 = (double) -30.486842105263158;
    R0_33 = (double) -0.021066676098596592;
    R0_53 = (double) -0.019684731834585762;
    R0_61 = (double) 1.9655620532813516;
    R0_57 = (double) -0.021739130434782608;
    I0_5 = (int) 5;
    R0_60 = (double) 2.909605009203791;
    R0_18 = (double) 0.06453000909551339;
    R0_52 = (double) -8.149310819916733;
    R0_50 = (double) 17.083154652734223;
    I0_7 = (int) 1944;
    R0_28 = (double) -210.49398960363874;
    R0_26 = (double) 372.96539961013644;
    R0_31 = (double) -0.47243356403005826;
    I0_10 = (int) 9;
    R0_38 = (double) 0.630323062817342;
    I0_0 = (int) 2;
    R0_56 = (double) 1.2806324110671936;
    R0_22 = (double) -172.3421052631579;
    R0_48 = (double) -0.00018693342388201837;
    R0_32 = (double) -0.08047388072894354;
    I0_8 = (int) 216;
    I0_4 = (int) -3;
    R0_24 = (double) 493.1895635676843;
    R0_62 = (double) -1.1640309347232585;
    R0_44 = (double) -33.027636352495925;
    I0_3 = (int) 4;
    R0_49 = (double) -5.1875;
    R0_55 = (double) 0.00018133307696307566;
    R0_20 = (double) -0.006844259178682439;
    R0_59 = (double) -0.6666666666666666;
    R0_54 = (double) 0.0014320902789711864;
    I0_6 = (int) 3;
    R0_42 = (double) 47.173489278752434;
    R0_17 = (double) 0.041666666666666664;
    I0_2 = (int) 15;
    R0_47 = (double) 0.001815517914137198;
    initialize = 0;
  }
  return 0;
}

void Uninitialize_phisr9()
{
  if( !initialize)
  {
    initialize = 1;
  }
}

int phisr9(double dr, double dth, double dph, double *Res)
{
  double R0_0;
  double R0_1;
  double R0_2;
  double R0_3;
  double R0_4;
  double R0_5;
  double R0_6;
  double R0_7;
  double R0_8;
  double R0_9;
  double R0_10;
  double R0_11;
  double R0_12;
  double R0_13;
  double R0_14;
  double R0_15;
  double R0_16;
  double R0_19;
  double R0_21;
  double R0_23;
  double R0_25;
  double R0_27;
  double R0_29;
  double R0_30;
  double R0_34;
  double R0_35;
  double R0_39;
  double R0_43;
  double R0_45;
  double R0_46;
  double R0_63;
  int err = 0;
  R0_0 = dr;
  R0_1 = dth;
  R0_2 = dph;
  R0_3 = cos(R0_2);
  R0_4 = (double) I0_0;
  R0_4 = R0_4 * R0_2;
  R0_5 = cos(R0_4);
  R0_6 = (double) I0_1;
  R0_6 = R0_6 * R0_3;
  R0_7 = (double) I0_2;
  R0_7 = R0_7 + R0_6 + R0_5;
  R0_8 = R0_1 * R0_1;
  if( I0_3 == 0)
  {
    if( R0_1 == 0)
    {
      err = 1;
      goto error_label;
    }
    else
    {
      R0_9 = 1;
    }
  }
  else
  {
    int S0 = I0_3;
    double S1 = R0_1;
    int S2 = 0;
    if( S0 < 0)
    {
      S2 = 1;
      S0 = -S0;
    }
    R0_9 = 1;
    while( S0)
    {
      if( S0 & 1)
      {
        R0_9 = S1 * R0_9;
      }
      S1 = S1 * S1;
      S0 = S0 >> 1;
    }
    if( S2)
    {
      R0_9 = 1 / R0_9;
    }
  }
  R0_10 = (double) I0_3;
  R0_10 = R0_10 * R0_3;
  R0_11 = -R0_5;
  R0_12 = (double) I0_4;
  R0_12 = R0_12 + R0_10 + R0_11;
  if( I0_3 == 0)
  {
    if( R0_0 == 0)
    {
      err = 1;
      goto error_label;
    }
    else
    {
      R0_13 = 1;
    }
  }
  else
  {
    int S0 = I0_3;
    double S1 = R0_0;
    int S2 = 0;
    if( S0 < 0)
    {
      S2 = 1;
      S0 = -S0;
    }
    R0_13 = 1;
    while( S0)
    {
      if( S0 & 1)
      {
        R0_13 = S1 * R0_13;
      }
      S1 = S1 * S1;
      S0 = S0 >> 1;
    }
    if( S2)
    {
      R0_13 = 1 / R0_13;
    }
  }
  if( I0_5 == 0)
  {
    if( R0_0 == 0)
    {
      err = 1;
      goto error_label;
    }
    else
    {
      R0_14 = 1;
    }
  }
  else
  {
    int S0 = I0_5;
    double S1 = R0_0;
    int S2 = 0;
    if( S0 < 0)
    {
      S2 = 1;
      S0 = -S0;
    }
    R0_14 = 1;
    while( S0)
    {
      if( S0 & 1)
      {
        R0_14 = S1 * R0_14;
      }
      S1 = S1 * S1;
      S0 = S0 >> 1;
    }
    if( S2)
    {
      R0_14 = 1 / R0_14;
    }
  }
  if( I0_6 == 0)
  {
    if( R0_0 == 0)
    {
      err = 1;
      goto error_label;
    }
    else
    {
      R0_15 = 1;
    }
  }
  else
  {
    int S0 = I0_6;
    double S1 = R0_0;
    int S2 = 0;
    if( S0 < 0)
    {
      S2 = 1;
      S0 = -S0;
    }
    R0_15 = 1;
    while( S0)
    {
      if( S0 & 1)
      {
        R0_15 = S1 * R0_15;
      }
      S1 = S1 * S1;
      S0 = S0 >> 1;
    }
    if( S2)
    {
      R0_15 = 1 / R0_15;
    }
  }
  R0_16 = R0_0 * R0_0;
  R0_19 = R0_18 * R0_13;
  R0_21 = R0_20 * R0_14;
  R0_23 = R0_22 * R0_9;
  R0_25 = R0_24 * R0_12;
  R0_27 = R0_26 * R0_7;
  R0_29 = R0_28 * R0_7;
  R0_30 = (double) I0_7;
  R0_30 = R0_30 + R0_29;
  R0_29 = R0_8 * R0_30;
  R0_30 = R0_32 * R0_8;
  R0_34 = R0_33 * R0_7;
  R0_35 = R0_31 + R0_30 + R0_34;
  R0_30 = R0_15 * R0_35;
  R0_35 = R0_37 * R0_8;
  R0_34 = R0_38 * R0_7;
  R0_39 = R0_36 + R0_35 + R0_34;
  R0_35 = R0_16 * R0_39;
  R0_34 = R0_40 * R0_9;
  R0_39 = R0_41 * R0_12;
  R0_43 = R0_42 * R0_7;
  R0_45 = R0_44 * R0_7;
  R0_46 = (double) I0_8;
  R0_46 = R0_46 + R0_45;
  R0_45 = R0_8 * R0_46;
  R0_34 = R0_34 + R0_39 + R0_43 + R0_45;
  R0_39 = R0_0 * R0_34;
  R0_19 = R0_19 + R0_21 + R0_23 + R0_25 + R0_27 + R0_29 + R0_30 + R0_35 + R0_39;
  R0_23 = R0_47 * R0_13;
  R0_27 = R0_48 * R0_14;
  R0_30 = R0_49 * R0_9;
  R0_39 = R0_50 * R0_12;
  R0_43 = R0_51 * R0_7;
  R0_46 = R0_52 * R0_7;
  R0_21 = (double) I0_9;
  R0_21 = R0_21 + R0_46;
  R0_46 = R0_8 * R0_21;
  R0_29 = R0_54 * R0_8;
  R0_34 = R0_55 * R0_7;
  R0_45 = R0_53 + R0_29 + R0_34;
  R0_29 = R0_15 * R0_45;
  R0_21 = R0_57 * R0_8;
  R0_35 = R0_58 * R0_7;
  R0_45 = R0_56 + R0_21 + R0_35;
  R0_21 = R0_16 * R0_45;
  R0_35 = R0_59 * R0_9;
  R0_25 = R0_60 * R0_12;
  R0_34 = R0_61 * R0_7;
  R0_45 = R0_62 * R0_7;
  R0_63 = (double) I0_10;
  R0_63 = R0_63 + R0_45;
  R0_45 = R0_8 * R0_63;
  R0_35 = R0_35 + R0_25 + R0_34 + R0_45;
  R0_25 = R0_0 * R0_35;
  R0_23 = R0_23 + R0_27 + R0_30 + R0_39 + R0_43 + R0_46 + R0_29 + R0_21 + R0_25;
  R0_30 = pow(R0_23, R0_64);
  R0_23 = 1 / R0_30;
  R0_30 = R0_17 * R0_19 * R0_23;
  *Res = R0_30;
  error_label:
  return err;
}

