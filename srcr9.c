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
static int I0_11;
static int I0_12;
static int I0_13;
static int I0_14;
static int I0_15;
static double R0_18;
static double R0_21;
static double R0_24;
static double R0_26;
static double R0_28;
static double R0_30;
static double R0_32;
static double R0_34;
static double R0_38;
static double R0_40;
static double R0_42;
static double R0_45;
static double R0_47;
static double R0_49;
static double R0_52;
static double R0_54;
static double R0_56;
static double R0_58;
static double R0_65;
static double R0_67;
static double R0_69;
static double R0_71;
static double R0_73;
static double R0_76;
static double R0_78;
static double R0_80;
static double R0_83;
static double R0_85;
static double R0_87;
static double R0_90;
static double R0_92;
static double R0_94;
static double R0_100;
static double R0_102;
static double R0_104;
static double R0_107;
static double R0_113;
static double R0_117;
static double R0_120;
static double R0_121;
static double R0_126;
static double R0_129;
static double R0_130;
static double R0_135;
static double R0_137;
static double R0_139;
static double R0_142;
static double R0_174;
static double R0_182;
static double R0_185;
static double R0_187;
static double R0_192;
static double R0_193;
static double R0_194;
static double R0_195;
static double R0_198;
static double R0_199;
static double R0_201;
static double R0_202;
static double R0_205;
static double R0_206;
static double R0_207;
static int initialize = 1;

#include "srcr9.h"

int Initialize_srcr9()
{
  if(initialize)
  {
    R0_45 = (double) 0.04743083003952569;
    R0_54 = (double) 102.57470363819434;
    R0_185 = (double) -0.0625;
    R0_198 = (double) 0.021786214969646375;
    R0_126 = (double) 1.5;
    R0_202 = (double) -0.1368851835736488;
    I0_7 = (int) 81;
    R0_73 = (double) 15.540224983755685;
    I0_1 = (int) -16;
    R0_85 = (double) 0.021573322379181633;
    I0_14 = (int) -8;
    R0_135 = (double) 0.09486166007905138;
    R0_49 = (double) 30.73517786561265;
    R0_117 = (double) 3.5;
    R0_52 = (double) -30.486842105263158;
    R0_113 = (double) 2.5;
    R0_139 = (double) -689.3684210526316;
    R0_40 = (double) -0.021066676098596592;
    R0_80 = (double) -0.019684731834585762;
    R0_94 = (double) 1.9655620532813516;
    R0_100 = (double) -0.043478260869565216;
    R0_83 = (double) -0.021739130434782608;
    R0_120 = (double) 0.007262071656548792;
    I0_12 = (int) -4;
    R0_192 = (double) -0.125;
    I0_15 = (int) 6;
    I0_6 = (int) 5;
    R0_92 = (double) 2.909605009203791;
    R0_24 = (double) 0.06453000909551339;
    R0_121 = (double) -0.0009346671194100919;
    R0_201 = (double) 0.7743601091461606;
    R0_18 = (double) -8.149310819916733;
    R0_71 = (double) 17.083154652734223;
    R0_207 = (double) 0.001322314049586777;
    I0_9 = (int) 1944;
    R0_194 = (double) -2068.1052631578946;
    I0_13 = (int) -2;
    R0_182 = (double) 0.15625;
    R0_34 = (double) -210.49398960363874;
    R0_32 = (double) 372.96539961013644;
    R0_42 = (double) -0.47243356403005826;
    I0_8 = (int) 9;
    I0_11 = (int) 16;
    R0_174 = (double) 0.25;
    R0_199 = (double) -0.0037386684776403675;
    R0_47 = (double) 0.630323062817342;
    I0_0 = (int) 2;
    R0_87 = (double) 1.2806324110671936;
    R0_107 = (double) -2.6666666666666665;
    R0_142 = (double) -121.94736842105263;
    R0_28 = (double) -172.3421052631579;
    R0_67 = (double) -0.00018693342388201837;
    R0_38 = (double) -0.08047388072894354;
    I0_10 = (int) 216;
    I0_4 = (int) -3;
    R0_137 = (double) -0.16094776145788708;
    R0_30 = (double) 493.1895635676843;
    R0_104 = (double) -20.75;
    R0_21 = (double) -1.1640309347232585;
    R0_129 = (double) 0.25812003638205355;
    R0_195 = (double) -365.8421052631579;
    R0_205 = (double) -0.001322314049586777;
    R0_206 = (double) 0.9636363636363636;
    R0_58 = (double) -33.027636352495925;
    I0_3 = (int) 4;
    R0_69 = (double) -5.1875;
    R0_78 = (double) 0.00018133307696307566;
    R0_130 = (double) -0.0342212958934122;
    R0_26 = (double) -0.006844259178682439;
    R0_90 = (double) -0.6666666666666666;
    R0_187 = (double) -62.25;
    R0_76 = (double) 0.0014320902789711864;
    I0_5 = (int) 3;
    R0_56 = (double) 47.173489278752434;
    R0_193 = (double) 0.041666666666666664;
    I0_2 = (int) 15;
    R0_65 = (double) 0.001815517914137198;
    R0_102 = (double) 0.002864180557942373;
    initialize = 0;
  }
  return 0;
}

void Uninitialize_srcr9()
{
  if( !initialize)
  {
    initialize = 1;
  }
}

int srcr9(double dr, double dth, double dph, double *Res)
{
  if(sqrt(dr*dr + dth*dth + 2.*(1-cos(dph))) < 0.02)
  {
    *Res = 0.0;
    return 0;
  }

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
  double R0_17;
  double R0_19;
  double R0_20;
  double R0_22;
  double R0_23;
  double R0_25;
  double R0_27;
  double R0_29;
  double R0_31;
  double R0_33;
  double R0_35;
  double R0_36;
  double R0_37;
  double R0_39;
  double R0_41;
  double R0_43;
  double R0_44;
  double R0_46;
  double R0_48;
  double R0_50;
  double R0_51;
  double R0_53;
  double R0_55;
  double R0_57;
  double R0_59;
  double R0_60;
  double R0_61;
  double R0_62;
  double R0_63;
  double R0_64;
  double R0_66;
  double R0_68;
  double R0_70;
  double R0_72;
  double R0_74;
  double R0_75;
  double R0_77;
  double R0_79;
  double R0_81;
  double R0_82;
  double R0_84;
  double R0_86;
  double R0_88;
  double R0_89;
  double R0_91;
  double R0_93;
  double R0_95;
  double R0_96;
  double R0_97;
  double R0_98;
  double R0_99;
  double R0_101;
  double R0_103;
  double R0_105;
  double R0_106;
  double R0_108;
  double R0_109;
  double R0_110;
  double R0_111;
  double R0_112;
  double R0_114;
  double R0_115;
  double R0_116;
  double R0_118;
  double R0_119;
  double R0_122;
  double R0_123;
  double R0_124;
  double R0_125;
  double R0_127;
  double R0_128;
  double R0_131;
  double R0_132;
  double R0_133;
  double R0_134;
  double R0_136;
  double R0_138;
  double R0_140;
  double R0_141;
  double R0_143;
  double R0_144;
  double R0_145;
  double R0_146;
  double R0_147;
  double R0_148;
  double R0_149;
  double R0_150;
  double R0_151;
  double R0_152;
  double R0_153;
  double R0_154;
  double R0_155;
  double R0_156;
  double R0_157;
  double R0_158;
  double R0_159;
  double R0_160;
  double R0_161;
  double R0_162;
  double R0_163;
  double R0_164;
  double R0_165;
  double R0_166;
  double R0_167;
  double R0_168;
  double R0_169;
  double R0_170;
  double R0_171;
  double R0_172;
  double R0_173;
  double R0_175;
  double R0_176;
  double R0_177;
  double R0_178;
  double R0_179;
  double R0_180;
  double R0_181;
  double R0_183;
  double R0_184;
  double R0_186;
  double R0_188;
  double R0_189;
  double R0_190;
  double R0_191;
  double R0_196;
  double R0_197;
  double R0_200;
  double R0_203;
  double R0_204;
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
  R0_13 = R0_0 * R0_0;
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
  if( I0_5 == 0)
  {
    if( R0_1 == 0)
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
    int S0 = I0_5;
    double S1 = R0_1;
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
  if( I0_3 == 0)
  {
    if( R0_0 == 0)
    {
      err = 1;
      goto error_label;
    }
    else
    {
      R0_16 = 1;
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
    R0_16 = 1;
    while( S0)
    {
      if( S0 & 1)
      {
        R0_16 = S1 * R0_16;
      }
      S1 = S1 * S1;
      S0 = S0 >> 1;
    }
    if( S2)
    {
      R0_16 = 1 / R0_16;
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
      R0_17 = 1;
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
    R0_17 = 1;
    while( S0)
    {
      if( S0 & 1)
      {
        R0_17 = S1 * R0_17;
      }
      S1 = S1 * S1;
      S0 = S0 >> 1;
    }
    if( S2)
    {
      R0_17 = 1 / R0_17;
    }
  }
  R0_19 = R0_18 * R0_7;
  R0_20 = (double) I0_7;
  R0_20 = R0_20 + R0_19;
  R0_22 = R0_21 * R0_7;
  R0_23 = (double) I0_8;
  R0_23 = R0_23 + R0_22;
  R0_25 = R0_24 * R0_16;
  R0_27 = R0_26 * R0_17;
  R0_29 = R0_28 * R0_9;
  R0_31 = R0_30 * R0_12;
  R0_33 = R0_32 * R0_7;
  R0_35 = R0_34 * R0_7;
  R0_36 = (double) I0_9;
  R0_36 = R0_36 + R0_35;
  R0_37 = R0_8 * R0_36;
  R0_39 = R0_38 * R0_8;
  R0_41 = R0_40 * R0_7;
  R0_43 = R0_42 + R0_39 + R0_41;
  R0_44 = R0_14 * R0_43;
  R0_46 = R0_45 * R0_8;
  R0_48 = R0_47 * R0_7;
  R0_50 = R0_49 + R0_46 + R0_48;
  R0_51 = R0_13 * R0_50;
  R0_53 = R0_52 * R0_9;
  R0_55 = R0_54 * R0_12;
  R0_57 = R0_56 * R0_7;
  R0_59 = R0_58 * R0_7;
  R0_60 = (double) I0_10;
  R0_60 = R0_60 + R0_59;
  R0_61 = R0_8 * R0_60;
  R0_62 = R0_53 + R0_55 + R0_57 + R0_61;
  R0_63 = R0_0 * R0_62;
  R0_64 = R0_25 + R0_27 + R0_29 + R0_31 + R0_33 + R0_37 + R0_44 + R0_51 + R0_63;
  R0_66 = R0_65 * R0_16;
  R0_68 = R0_67 * R0_17;
  R0_70 = R0_69 * R0_9;
  R0_72 = R0_71 * R0_12;
  R0_74 = R0_73 * R0_7;
  R0_75 = R0_8 * R0_20;
  R0_77 = R0_76 * R0_8;
  R0_79 = R0_78 * R0_7;
  R0_81 = R0_80 + R0_77 + R0_79;
  R0_82 = R0_14 * R0_81;
  R0_84 = R0_83 * R0_8;
  R0_86 = R0_85 * R0_7;
  R0_88 = R0_87 + R0_84 + R0_86;
  R0_89 = R0_13 * R0_88;
  R0_91 = R0_90 * R0_9;
  R0_93 = R0_92 * R0_12;
  R0_95 = R0_94 * R0_7;
  R0_96 = R0_8 * R0_23;
  R0_97 = R0_91 + R0_93 + R0_95 + R0_96;
  R0_98 = R0_0 * R0_97;
  R0_99 = R0_66 + R0_68 + R0_70 + R0_72 + R0_74 + R0_75 + R0_82 + R0_89 + R0_98;
  R0_101 = R0_100 * R0_13 * R0_1;
  R0_103 = R0_102 * R0_14 * R0_1;
  R0_105 = R0_104 * R0_15;
  R0_106 = (double) I0_0;
  R0_106 = R0_106 * R0_1 * R0_20;
  R0_108 = R0_107 * R0_15;
  R0_109 = (double) I0_0;
  R0_109 = R0_109 * R0_1 * R0_23;
  R0_110 = R0_108 + R0_109;
  R0_111 = R0_0 * R0_110;
  R0_112 = R0_101 + R0_103 + R0_105 + R0_106 + R0_111;
  R0_114 = pow(R0_99, R0_113);
  R0_115 = 1 / R0_114;
  R0_114 = (double) I0_8;
  R0_114 = R0_114 + R0_0;
  R0_116 = R0_114 * R0_114;
  R0_118 = pow(R0_99, R0_117);
  R0_119 = 1 / R0_118;
  R0_118 = R0_120 * R0_14;
  R0_122 = R0_121 * R0_16;
  R0_123 = (double) I0_5;
  R0_123 = R0_123 * R0_13 * R0_81;
  R0_124 = (double) I0_0;
  R0_124 = R0_124 * R0_0 * R0_88;
  R0_125 = R0_118 + R0_122 + R0_91 + R0_93 + R0_95 + R0_96 + R0_123 + R0_124;
  R0_127 = pow(R0_99, R0_126);
  R0_128 = 1 / R0_127;
  R0_127 = R0_129 * R0_14;
  R0_131 = R0_130 * R0_16;
  R0_132 = (double) I0_5;
  R0_132 = R0_132 * R0_13 * R0_43;
  R0_133 = (double) I0_0;
  R0_133 = R0_133 * R0_0 * R0_50;
  R0_134 = R0_127 + R0_131 + R0_53 + R0_55 + R0_57 + R0_61 + R0_132 + R0_133;
  R0_136 = R0_135 * R0_13 * R0_1;
  R0_138 = R0_137 * R0_14 * R0_1;
  R0_140 = R0_139 * R0_15;
  R0_141 = (double) I0_0;
  R0_141 = R0_141 * R0_1 * R0_36;
  R0_143 = R0_142 * R0_15;
  R0_144 = (double) I0_0;
  R0_144 = R0_144 * R0_1 * R0_60;
  R0_145 = R0_143 + R0_144;
  R0_146 = R0_0 * R0_145;
  R0_147 = R0_136 + R0_138 + R0_140 + R0_141 + R0_146;
  R0_148 = (double) I0_11;
  R0_148 = R0_148 * R0_3;
  R0_149 = (double) I0_12;
  R0_149 = R0_149 * R0_5;
  R0_150 = R0_148 + R0_149;
  R0_151 = (double) I0_12;
  R0_151 = R0_151 * R0_3;
  R0_152 = (double) I0_3;
  R0_152 = R0_152 * R0_5;
  R0_153 = R0_151 + R0_152;
  R0_154 = sin(R0_2);
  R0_155 = (double) I0_11;
  R0_155 = R0_155 * R0_154;
  R0_156 = sin(R0_4);
  R0_157 = (double) I0_13;
  R0_157 = R0_157 * R0_156;
  R0_158 = R0_155 + R0_157;
  R0_159 = (double) I0_12;
  R0_159 = R0_159 * R0_154;
  R0_160 = (double) I0_0;
  R0_160 = R0_160 * R0_156;
  R0_161 = R0_159 + R0_160;
  R0_162 = R0_73 * R0_158;
  R0_163 = R0_85 * R0_13 * R0_158;
  R0_164 = R0_78 * R0_14 * R0_158;
  R0_165 = R0_18 * R0_8 * R0_158;
  R0_166 = R0_71 * R0_161;
  R0_167 = R0_94 * R0_158;
  R0_168 = R0_21 * R0_8 * R0_158;
  R0_169 = R0_92 * R0_161;
  R0_170 = R0_167 + R0_168 + R0_169;
  R0_171 = R0_0 * R0_170;
  R0_172 = R0_162 + R0_163 + R0_164 + R0_165 + R0_166 + R0_171;
  R0_173 = (double) I0_13;
  R0_173 = R0_173 * R0_114;
  R0_175 = R0_174 + R0_173 + R0_116;
  R0_176 = sin(R0_1);
  R0_177 = R0_176 * R0_176;
  R0_178 = R0_174 * R0_177;
  R0_179 = R0_116 + R0_178;
  R0_180 = 1 / R0_179;
  R0_181 = tan(R0_1);
  R0_183 = R0_112 * R0_112;
  R0_184 = R0_182 * R0_64 * R0_183 * R0_119;
  R0_183 = R0_100 * R0_13;
  R0_186 = R0_102 * R0_14;
  R0_188 = R0_187 * R0_8;
  R0_189 = (double) I0_0;
  R0_189 = R0_189 * R0_20;
  R0_190 = (double) I0_14;
  R0_190 = R0_190 * R0_8;
  R0_191 = (double) I0_0;
  R0_191 = R0_191 * R0_23;
  R0_190 = R0_190 + R0_191;
  R0_191 = R0_0 * R0_190;
  R0_183 = R0_183 + R0_186 + R0_188 + R0_189 + R0_191;
  R0_186 = R0_185 * R0_64 * R0_183 * R0_115;
  R0_188 = R0_192 * R0_147 * R0_112 * R0_115;
  R0_191 = R0_135 * R0_13;
  R0_190 = R0_137 * R0_14;
  R0_189 = R0_194 * R0_8;
  R0_183 = (double) I0_0;
  R0_183 = R0_183 * R0_36;
  R0_196 = R0_195 * R0_8;
  R0_197 = (double) I0_0;
  R0_197 = R0_197 * R0_60;
  R0_196 = R0_196 + R0_197;
  R0_197 = R0_0 * R0_196;
  R0_191 = R0_191 + R0_190 + R0_189 + R0_183 + R0_197;
  R0_190 = R0_193 * R0_191 * R0_128;
  R0_191 = R0_125 * R0_125;
  R0_189 = R0_182 * R0_191 * R0_64 * R0_119;
  R0_191 = R0_192 * R0_125 * R0_134 * R0_115;
  R0_197 = R0_198 * R0_13;
  R0_183 = R0_199 * R0_14;
  R0_196 = (double) I0_15;
  R0_196 = R0_196 * R0_0 * R0_81;
  R0_200 = (double) I0_0;
  R0_200 = R0_200 * R0_88;
  R0_197 = R0_197 + R0_183 + R0_196 + R0_200;
  R0_183 = R0_185 * R0_197 * R0_64 * R0_115;
  R0_196 = R0_201 * R0_13;
  R0_197 = R0_202 * R0_14;
  R0_200 = (double) I0_15;
  R0_200 = R0_200 * R0_0 * R0_43;
  R0_203 = (double) I0_0;
  R0_203 = R0_203 * R0_50;
  R0_196 = R0_196 + R0_197 + R0_200 + R0_203;
  R0_197 = R0_193 * R0_196 * R0_128;
  R0_189 = R0_189 + R0_191 + R0_183 + R0_197;
  R0_191 = R0_175 * R0_189;
  R0_189 = -R0_0;
  R0_183 = (double) I0_14;
  R0_183 = R0_183 + R0_189;
  R0_189 = R0_185 * R0_125 * R0_64 * R0_115;
  R0_197 = R0_193 * R0_134 * R0_128;
  R0_189 = R0_189 + R0_197;
  R0_197 = (double) I0_13;
  R0_197 = R0_197 * R0_183 * R0_189;
  R0_183 = R0_185 * R0_64 * R0_112 * R0_115;
  R0_189 = R0_193 * R0_147 * R0_128;
  R0_183 = R0_183 + R0_189;
  R0_183 = R0_183 * R0_181;
  R0_189 = -R0_183;
  R0_184 = R0_184 + R0_186 + R0_188 + R0_190 + R0_191 + R0_197 + R0_189;
  R0_186 = R0_180 * R0_184;
  R0_184 = R0_73 * R0_150;
  R0_188 = R0_85 * R0_13 * R0_150;
  R0_190 = R0_78 * R0_14 * R0_150;
  R0_191 = R0_18 * R0_8 * R0_150;
  R0_197 = R0_71 * R0_153;
  R0_189 = R0_94 * R0_150;
  R0_183 = R0_21 * R0_8 * R0_150;
  R0_196 = R0_92 * R0_153;
  R0_189 = R0_189 + R0_183 + R0_196;
  R0_183 = R0_0 * R0_189;
  R0_184 = R0_184 + R0_188 + R0_190 + R0_191 + R0_197 + R0_183;
  R0_188 = R0_185 * R0_184 * R0_64 * R0_115;
  R0_184 = R0_32 * R0_150;
  R0_190 = R0_47 * R0_13 * R0_150;
  R0_191 = R0_40 * R0_14 * R0_150;
  R0_197 = R0_34 * R0_8 * R0_150;
  R0_183 = R0_30 * R0_153;
  R0_189 = R0_56 * R0_150;
  R0_196 = R0_58 * R0_8 * R0_150;
  R0_200 = R0_54 * R0_153;
  R0_189 = R0_189 + R0_196 + R0_200;
  R0_196 = R0_0 * R0_189;
  R0_184 = R0_184 + R0_190 + R0_191 + R0_197 + R0_183 + R0_196;
  R0_190 = R0_193 * R0_184 * R0_128;
  R0_184 = R0_172 * R0_172;
  R0_191 = R0_182 * R0_64 * R0_119 * R0_184;
  R0_184 = R0_32 * R0_158;
  R0_197 = R0_47 * R0_13 * R0_158;
  R0_183 = R0_40 * R0_14 * R0_158;
  R0_196 = R0_34 * R0_8 * R0_158;
  R0_189 = R0_30 * R0_161;
  R0_200 = R0_56 * R0_158;
  R0_203 = R0_58 * R0_8 * R0_158;
  R0_204 = R0_54 * R0_161;
  R0_200 = R0_200 + R0_203 + R0_204;
  R0_203 = R0_0 * R0_200;
  R0_184 = R0_184 + R0_197 + R0_183 + R0_196 + R0_189 + R0_203;
  R0_197 = R0_192 * R0_115 * R0_172 * R0_184;
  R0_188 = R0_188 + R0_190 + R0_191 + R0_197;
  R0_191 = 1 / R0_175;
  R0_183 = R0_174 + R0_116;
  R0_196 = R0_207 * R0_183;
  R0_183 = R0_206 + R0_196;
  R0_196 = (double) I0_13;
  R0_196 = R0_196 * R0_183;
  R0_183 = (double) I0_8;
  R0_183 = R0_183 + R0_0 + R0_196;
  R0_196 = R0_114 * R0_191 * R0_183 * R0_180;
  R0_191 = R0_181 * R0_181;
  R0_183 = R0_180 * R0_191;
  R0_191 = R0_205 + R0_196 + R0_183;
  R0_188 = R0_188 * R0_191;
  R0_186 = R0_186 + R0_188;
  *Res = R0_186;
  error_label:
  return err;
}

