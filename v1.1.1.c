/* Wed 22 Sep 2010 22:42:52 CEST
   4th order singular field and corresponding effective source in Schwarzschild
   for a particle in a circular orbit at radius 10M 
   
   Variables: dr, dth, dph
   Function: pow, cos, sin, tan, sec
   */
#include <math.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
  double phis, phinum, s2, src, srcnum, dr, dth, dph, srcapprox;

  dth = 0.0;
  
  for(dr=-1.; dr<1.; dr+=0.1)
  {
    for(dph=-0.5; dph<0.5; dph+=0.05)
    {
      s2 = 7*pow(dr,2) + 80*(8*pow(dph,2) + 7*pow(dth,2));
     
      phis = (5186160*pow(dr,8) - 295323*pow(dr,9) - 3292800*pow(dr,6)*(-31360 + 3688*pow(dph,2) + 1463*pow(dth,2)) - 
           20580*pow(dr,7)*(-31360 + 34872*pow(dph,2) + 23443*pow(dth,2)) + 
           51380224e5*pow(8*pow(dph,2) + 7*pow(dth,2),2)*(pow(dph,2)*(240 + 7*pow(dph,2)) + 6*(35 + 9*pow(dph,2))*pow(dth,2) + 
              2*pow(dth,4)) + 6272e4*pow(dr,4)*(64*pow(dph,2)*(7056 + 743*pow(dph,2)) + 672*(588 + 167*pow(dph,2))*pow(dth,2) + 
              37387*pow(dth,4)) - 156800*pow(dr,5)*(480*pow(dph,2)*(5488 - 1201*pow(dph,2)) - 252*(-7840 + 7397*pow(dph,2))*pow(dth,2) + 
              178409*pow(dth,4)) + 1605632e4*pow(dr,2)*(160*pow(dph,4)*(1008 + 157*pow(dph,2)) + 
              28*pow(dph,2)*(10080 + 2951*pow(dph,2))*pow(dth,2) + 7*(17640 + 9893*pow(dph,2))*pow(dth,4) + 14161*pow(dth,6)) - 
           200704e3*pow(dr,3)*(48*pow(dph,4)*(9520 + 3907*pow(dph,2)) + 6*pow(dph,2)*(125440 + 68567*pow(dph,2))*pow(dth,2) + 
              7*(44100 + 49673*pow(dph,2))*pow(dth,4) + 126469*pow(dth,6)) - 
           131072e4*dr*(96*pow(dph,6)*(35280 + 2549*pow(dph,2)) + 588*pow(dph,4)*(14560 + 2341*pow(dph,2))*pow(dth,2) + 
              49*pow(dph,2)*(147e3 + 38021*pow(dph,2))*pow(dth,4) + 686*(2940 + 1339*pow(dph,2))*pow(dth,6) + 156065*pow(dth,8)))/
           (2.1504e7*sqrt(35)*pow(s2,3.5));

      if(s2>1e-16)
        src = -(1293025674264576e9*pow(dph,8) + 4729480547401728e7*pow(dph,12) + 888955151056896e8*pow(dph,8)*dr + 
          1116682907025408e7*pow(dph,12)*dr + 35356170780672e9*pow(dph,6)*pow(dr,2) - 392706039742464e7*pow(dph,8)*pow(dr,2) + 
          21008127080005632e6*pow(dph,10)*pow(dr,2) + 65687229825024e7*pow(dph,12)*pow(dr,2) + 
          7734162358272e9*pow(dph,6)*pow(dr,3) + 23005187014656e7*pow(dph,8)*pow(dr,3) + 
          5435219238715392e6*pow(dph,10)*pow(dr,3) + 23202487074816e7*pow(dph,4)*pow(dr,4) - 
          286459957346304e6*pow(dph,6)*pow(dr,4) - 17378924241616896e5*pow(dph,8)*pow(dr,4) + 
          351150419214336e6*pow(dph,10)*pow(dr,4) + 117462590816256e6*pow(dph,4)*pow(dr,5) + 8231728447488e6*pow(dph,6)*pow(dr,5) - 
          35554501117083648e4*pow(dph,8)*pow(dr,5) - 845924007936e6*pow(dph,2)*pow(dr,6) + 9509092196352e6*pow(dph,4)*pow(dr,6) + 
          4029122022998016e4*pow(dph,6)*pow(dr,6) - 1809088335839232e4*pow(dph,8)*pow(dr,6) + 2590642274304e5*pow(dph,2)*pow(dr,7) - 
          17923756720128e5*pow(dph,4)*pow(dr,7) + 4660418854060032e3*pow(dph,6)*pow(dr,7) - 92522938368e5*pow(dr,8) + 
          11678424662016e4*pow(dph,2)*pow(dr,8) - 299017947820032e3*pow(dph,4)*pow(dr,8) + 139846422331392e3*pow(dph,6)*pow(dr,8) - 
          237090029568e4*pow(dr,9) + 1181887601664e4*pow(dph,2)*pow(dr,9) - 14411308520448e3*pow(dph,4)*pow(dr,9) - 
          178850144256e3*pow(dr,10) + 505753908307200*pow(dph,2)*pow(dr,10) - 251048182425600*pow(dph,4)*pow(dr,10) - 
          4603235616e3*pow(dr,11) + 21395496165120*pow(dph,2)*pow(dr,11) - 196084243740*pow(dr,12) + 474266264640*pow(dph,2)*pow(dr,12) - 
          4891442052*pow(dr,13) + 64387617*pow(dr,14) + 282849366245376e90*pow(dph,6)*pow(dth,2) + 
          21001565443719168e8*pow(dph,8)*pow(dth,2) + 3488084675002368e8*pow(dph,10)*pow(dth,2) + 
          2368863442305024e8*pow(dph,6)*dr*pow(dth,2) + 13742817310408704e7*pow(dph,8)*dr*pow(dth,2) + 
          823575548264448e8*pow(dph,10)*dr*pow(dth,2) + 371239793197056e8*pow(dph,4)*pow(dr,2)*pow(dth,2) + 
          460182660317184e8*pow(dph,6)*pow(dr,2)*pow(dth,2) + 38761135758901248e6*pow(dph,8)*pow(dr,2)*pow(dth,2) + 
          48445620486144e8*pow(dph,10)*pow(dr,2)*pow(dth,2) + 1148523110203392e7*pow(dph,4)*pow(dr,3)*pow(dth,2) + 
          12790410460004352e6*pow(dph,6)*pow(dr,3)*pow(dth,2) + 15586746678902784e6*pow(dph,8)*pow(dr,3)*pow(dth,2) - 
          20302176190464e7*pow(dph,2)*pow(dr,4)*pow(dth,2) + 291481243877376e6*pow(dph,4)*pow(dr,4)*pow(dth,2) - 
          42237344247447552e5*pow(dph,6)*pow(dr,4)*pow(dth,2) + 109846670384037888e4*pow(dph,8)*pow(dr,4)*pow(dth,2) + 
          39335466369024e6*pow(dph,2)*pow(dr,5)*pow(dth,2) + 1330925474414592e5*pow(dph,4)*pow(dr,5)*pow(dth,2) - 
          106836963938009088e4*pow(dph,6)*pow(dr,5)*pow(dth,2) - 2960734027776e6*pow(dr,6)*pow(dth,2) + 
          12355588718592e6*pow(dph,2)*pow(dr,6)*pow(dth,2) + 9770391372005376e4*pow(dph,4)*pow(dr,6)*pow(dth,2) - 
          57354959166898176e3*pow(dph,6)*pow(dr,6)*pow(dth,2) - 6337821278208e5*pow(dr,7)*pow(dth,2) - 
          5637290459136e5*pow(dph,2)*pow(dr,7)*pow(dth,2) + 11941199189090304e3*pow(dph,4)*pow(dr,7)*pow(dth,2) - 
          2136908132352e4*pow(dr,8)*pow(dth,2) - 287074234116096e3*pow(dph,2)*pow(dr,8)*pow(dth,2) + 
          363129071920742400*pow(dph,4)*pow(dr,8)*pow(dth,2) + 1861262575872e3*pow(dr,9)*pow(dth,2) - 
          21314424268032e3*pow(dph,2)*pow(dr,9)*pow(dth,2) + 83517626757600*pow(dr,10)*pow(dth,2) - 
          320366113804800*pow(dph,2)*pow(dr,10)*pow(dth,2) + 5439280334880*pow(dr,11)*pow(dth,2) + 178660964664*pow(dr,12)*pow(dth,2) + 
          1484959172788224e9*pow(dph,4)*pow(dth,4) + 40687376246243328e8*pow(dph,6)*pow(dth,4) + 
          86185874999476224e7*pow(dph,8)*pow(dth,4) + 1670579069386752e8*pow(dph,4)*dr*pow(dth,4) + 
          3132848735256576e8*pow(dph,6)*dr*pow(dth,4) + 20349442708209664e7*pow(dph,8)*dr*pow(dth,4) - 
          162417409523712e8*pow(dph,2)*pow(dr,2)*pow(dth,4) + 6347813755551744e7*pow(dph,4)*pow(dr,2)*pow(dth,4) + 
          4485880756568064e7*pow(dph,6)*pow(dr,2)*pow(dth,4) + 1197026041659392e7*pow(dph,8)*pow(dr,2)*pow(dth,4) + 
          131964145238016e7*pow(dph,2)*pow(dr,3)*pow(dth,4) + 12848454904578048e6*pow(dph,4)*pow(dr,3)*pow(dth,4) + 
          214397245053730816e5*pow(dph,6)*pow(dr,3)*pow(dth,4) - 35528808333312e7*pow(dr,4)*pow(dth,4) + 
          46420079935488e7*pow(dph,2)*pow(dr,4)*pow(dth,4) - 52947491234512896e5*pow(dph,4)*pow(dr,4)*pow(dth,4) + 
          15973123188850688e5*pow(dph,6)*pow(dr,4)*pow(dth,4) - 6106513932288e7*pow(dr,5)*pow(dth,4) + 
          302741663121408e5*pow(dph,2)*pow(dr,5)*pow(dth,4) - 129556442190249984e4*pow(dph,4)*pow(dr,5)*pow(dth,4) + 
          55348543488e7*pow(dr,6)*pow(dth,4) + 560096074911744e5*pow(dph,2)*pow(dr,6)*pow(dth,4) - 
          70184990552162304e3*pow(dph,4)*pow(dr,6)*pow(dth,4) + 25169027315712e4*pow(dr,7)*pow(dth,4) + 
          8649760935574528e3*pow(dph,2)*pow(dr,7)*pow(dth,4) - 3590577025056e4*pow(dr,8)*pow(dth,4) + 
          253941108762521600*pow(dph,2)*pow(dr,8)*pow(dth,4) - 3328802863836800*pow(dr,9)*pow(dth,4) - 6951275702560*pow(dr,10)*pow(dth,4) - 
          433113092063232e9*pow(dph,2)*pow(dth,6) + 24432218889781248e8*pow(dph,4)*pow(dth,6) + 
          101135083527733248e7*pow(dph,6)*pow(dth,6) - 13534784126976e9*pow(dph,2)*dr*pow(dth,6) + 
          17981651262898176e7*pow(dph,4)*dr*pow(dth,6) + 23879116944048128e7*pow(dph,6)*dr*pow(dth,6) - 
          189486977777664e8*pow(dr,2)*pow(dth,6) + 2482907809579008e7*pow(dph,2)*pow(dr,2)*pow(dth,6) + 
          49392296518483968e6*pow(dph,4)*pow(dr,2)*pow(dth,6) + 1404653937885184e7*pow(dph,6)*pow(dr,2)*pow(dth,6) - 
          245740924305408e7*pow(dr,3)*pow(dth,6) + 1506167696130048e6*pow(dph,2)*pow(dr,3)*pow(dth,6) + 
          196239609869893632e5*pow(dph,4)*pow(dr,3)*pow(dth,6) + 150111858720768e6*pow(dr,4)*pow(dth,6) - 
          3581413892603904e6*pow(dph,2)*pow(dr,4)*pow(dth,6) + 143546207664340992e4*pow(dph,4)*pow(dr,4)*pow(dth,6) + 
          1092927209472e7*pow(dr,5)*pow(dth,6) - 74587057627815936e4*pow(dph,2)*pow(dr,5)*pow(dth,6) + 
          1108569446186496e4*pow(dr,6)*pow(dth,6) - 39156082447663104e3*pow(dph,2)*pow(dr,6)*pow(dth,6) + 
          1926630653553536e3*pow(dr,7)*pow(dth,6) + 495890400544e5*pow(dr,8)*pow(dth,6) - 378973955555328e9*pow(dth,8) + 
          5580198141493248e8*pow(dph,2)*pow(dth,8) + 61652633546391552e7*pow(dph,4)*pow(dth,8) - 
          331602211110912e8*dr*pow(dth,8) + 1735836064284672e7*pow(dph,2)*dr*pow(dth,8) + 
          14556871809564672e7*pow(dph,4)*dr*pow(dth,8) + 625560803868672e7*pow(dr,2)*pow(dth,8) + 
          33626227051855872e6*pow(dph,2)*pow(dr,2)*pow(dth,8) + 856286577033216e7*pow(dph,4)*pow(dr,2)*pow(dth,8) + 
          490477313851392e6*pow(dr,3)*pow(dth,8) + 111375127168155648e5*pow(dph,2)*pow(dr,3)*pow(dth,8) - 
          8062028778725376e5*pow(dr,4)*pow(dth,8) + 76433864040382464e4*pow(dph,2)*pow(dr,4)*pow(dth,8) - 
          16587064542418944e4*pow(dr,5)*pow(dth,8) - 8272992437047296e3*pow(dr,6)*pow(dth,8) + 1123387082539008e8*pow(dth,10) + 
          18870872769036288e7*pow(dph,2)*pow(dth,10) + 683506598412288e7*dr*pow(dth,10) + 
          4455622737133568e7*pow(dph,2)*dr*pow(dth,10) + 971305806987264e7*pow(dr,2)*pow(dth,10) + 
          262095455125504e7*pow(dph,2)*pow(dr,2)*pow(dth,10) + 27347767106207744e5*pow(dr,3)*pow(dth,10) + 
          17605229465829376e4*pow(dr,4)*pow(dth,10) + 230937254166528e8*pow(dth,12) + 54526851678208e8*dr*pow(dth,12) + 
          3207461863424e8*pow(dr,2)*pow(dth,12) + 588e3*(8 + dr)*
           (5751424*pow(dr,10) + 259651*pow(dr,11) - 62720*pow(dr,8)*(-31360 + 241272*pow(dph,2) + 43211*pow(dth,2)) - 
             392*pow(dr,9)*(-658560 + 1584040*pow(dph,2) + 360759*pow(dth,2)) - 
             50176e3*pow(dr,6)*(-32*pow(dph,2)*(112 + 205*pow(dph,2)) + 4*(-3136 + 2279*pow(dph,2))*pow(dth,2) + 12733*pow(dth,4)) + 
             4480*pow(dr,7)*(960*pow(dph,2)*(-18032 + 79029*pow(dph,2)) + 1120*(11172 + 70951*pow(dph,2))*pow(dth,2) + 587069*pow(dth,4)) - 
             16777216e5*pow(8*pow(dph,2) + 7*pow(dth,2),2)*
              (2560*pow(dph,4) + 2*pow(dph,2)*(560 + 239*pow(dph,2))*pow(dth,2) - (980 + 2301*pow(dph,2))*pow(dth,4) + 413*pow(dth,6)) + 
             409600*pow(dr,5)*(-960*pow(dph,4)*(47824 + 51809*pow(dph,2)) - 7056*pow(dph,2)*(4760 + 7723*pow(dph,2))*pow(dth,2) + 
                7203*(1200 + 4411*pow(dph,2))*pow(dth,4) + 1689275*pow(dth,6)) - 
             8192e3*pow(dr,4)*(512*pow(dph,4)*(11760 - 35519*pow(dph,2)) - 37632*pow(dph,2)*(140 + 829*pow(dph,2))*pow(dth,2) - 
                27440*(336 + 437*pow(dph,2))*pow(dth,4) + 5336051*pow(dth,6)) - 
             262144e4*pow(dr,2)*(1024*pow(dph,6)*(2800 + 409*pow(dph,2)) + 13440*pow(dph,4)*(224 + 3*pow(dph,2))*pow(dth,2) - 
                1568*pow(dph,2)*(840 + 1877*pow(dph,2))*pow(dth,4) - 196*(7840 + 8627*pow(dph,2))*pow(dth,6) + 555317*pow(dth,8)) - 
             32768e3*pow(dr,3)*(-7680*pow(dph,6)*(-2800 + 351*pow(dph,2)) + 768*pow(dph,4)*(57820 + 20627*pow(dph,2))*pow(dth,2) - 
                56*pow(dph,2)*(-388080 + 127951*pow(dph,2))*pow(dth,4) - 784*(735 + 29809*pow(dph,2))*pow(dth,6) + 966231*pow(dth,8)) + 
             2097152e4*dr*(737280*pow(dph,8) + 192*pow(dph,6)*(6160 + 711*pow(dph,2))*pow(dth,2) - 
                40*pow(dph,4)*(-4704 + 12217*pow(dph,2))*pow(dth,4) - 56*pow(dph,2)*(7350 + 8557*pow(dph,2))*pow(dth,6) + 
                49*(-2940 + 2089*pow(dph,2))*pow(dth,8) + 49049*pow(dth,10)))*pow(cos(dth),-2) - 
          24500*(8 + dr)*dth*(62694912*pow(dr,10) + 3997665*pow(dr,11) - 2744*pow(dr,9)*(24*(-74480 + 62647*pow(dph,2)) + 31199*pow(dth,2)) - 
             439040*pow(dr,8)*(-94080 + 214680*pow(dph,2) + 117341*pow(dth,2)) - 
             33554432e5*pow(8*pow(dph,2) + 7*pow(dth,2),3)*
              (3*pow(dph,2)*(-560 + 239*pow(dph,2)) - 2*(735 + 157*pow(dph,2))*pow(dth,2) + 14*pow(dth,4)) + 
             7340032e4*pow(dr,2)*pow(8*pow(dph,2) + 7*pow(dth,2),2)*
              (672*pow(dph,2)*(40 - 13*pow(dph,2)) + (23520 + 5423*pow(dph,2))*pow(dth,2) + 1519*pow(dth,4)) - 
             10035200*pow(dr,6)*(768*pow(dph,2)*(-1960 + 2659*pow(dph,2)) + 56*(-23520 + 45583*pow(dph,2))*pow(dth,2) + 85e33*pow(dth,4)) + 
             31360*pow(dr,7)*(1920*pow(dph,2)*(18032 + 6659*pow(dph,2)) + 392*(82320 + 133009*pow(dph,2))*pow(dth,2) + 22177351*pow(dth,4)) + 
             4194304e4*dr*pow(8*pow(dph,2) + 7*pow(dth,2),2)*
              (12*pow(dph,4)*(-6160 + 2133*pow(dph,2)) + 5*pow(dph,2)*(-22344 + 1357*pow(dph,2))*pow(dth,2) + 
                14*(-2940 + 481*pow(dph,2))*pow(dth,4) + 3185*pow(dth,6)) - 
             172032e3*pow(dr,4)*(1024*pow(dph,4)*(-11760 + 7603*pow(dph,2)) + 1792*pow(dph,2)*(-11760 + 6799*pow(dph,2))*pow(dth,2) + 
                392*(-23520 + 16223*pow(dph,2))*pow(dth,4) + 1451919*pow(dth,6)) + 
             2150400*pow(dr,5)*(1024*pow(dph,4)*(23520 + 28313*pow(dph,2)) + 10976*pow(dph,2)*(4560 + 11099*pow(dph,2))*pow(dth,2) + 
                588*(43120 + 198559*pow(dph,2))*pow(dth,4) + 26163697*pow(dth,6)) - 
             131072e3*pow(dr,3)*(1536*pow(dph,6)*(13720 + 1651*pow(dph,2)) - 112*pow(dph,4)*(-388080 + 64189*pow(dph,2))*pow(dth,2) - 
                294*pow(dph,2)*(-94080 + 85319*pow(dph,2))*pow(dth,4) - 1029*(-4900 + 13687*pow(dph,2))*pow(dth,6) + 578641*pow(dth,8)))*tan(dth)
          )/(5.376e7*sqrt(35)*(8 + dr)*pow(10 + dr,2)*pow(7*pow(dr,2) + 80*(8*pow(dph,2) + 7*pow(dth,2)),5.5));
      else
        src = 0.0;

      printf("%g\t%g\t%g\t%g\n", dr, dph, phis, src);
    }
  }
}

static double phisnum060000, phisnum070000, phisnum080000, phisnum090000, phisnum040200, phisnum050200, phisnum060200, phisnum070200, phisnum020400, phisnum030400, phisnum040400, phisnum050400, phisnum000600, phisnum010600, phisnum020600, phisnum030600, phisnum000800, phisnum010800, phisnum040002, phisnum050002, phisnum060002, phisnum070002, phisnum020202, phisnum030202, phisnum040202, phisnum050202, phisnum000402, phisnum010402, phisnum020402, phisnum030402, phisnum000602, phisnum010602, phisnum020004, phisnum030004, phisnum040004, phisnum050004, phisnum000204, phisnum010204, phisnum020204, phisnum030204, phisnum000404, phisnum010404, phisnum000006, phisnum010006, phisnum020006, phisnum030006, phisnum000206, phisnum010206, phisnum000008, phisnum010008;
static double sqrt35 = sqrt(35.);

static void phisnum_init(const double r1)
{
  phisnum060000=103262208000*sqrt35;
  phisnum070000=645388800*sqrt35;
  phisnum080000=5186160*sqrt35;
  phisnum090000=-295323*sqrt35;
  phisnum040200=24782929920000*sqrt35;
  phisnum050200=-309786624000*sqrt35;
  phisnum060200=-4817366400*sqrt35;
  phisnum070200=-482456940*sqrt35;
  phisnum020400=1982634393600000*sqrt35;
  phisnum030400=-61957324800000*sqrt35;
  phisnum040400=2344912640000*sqrt35;
  phisnum050400=-27974531200*sqrt35;
  phisnum000600=52870250496000000*sqrt35;
  phisnum010600=-2643512524800000*sqrt35;
  phisnum020600=227373547520000*sqrt35;
  phisnum030600=-25382834176000*sqrt35;
  phisnum000800=503526195200000*sqrt35;
  phisnum010800=-204557516800000*sqrt35;
  phisnum040002=28323348480000*sqrt35;
  phisnum050002=-413048832000*sqrt35;
  phisnum060002=-12143846400*sqrt35;
  phisnum070002=-717665760*sqrt35;
  phisnum020202=4531735756800000*sqrt35;
  phisnum030202=-151057858560000*sqrt35;
  phisnum040202=7038689280000*sqrt35;
  phisnum050202=292282099200*sqrt35;
  phisnum000402=181269430272000000*sqrt35;
  phisnum010402=-9441116160000000*sqrt35;
  phisnum020402=1111916216320000*sqrt35;
  phisnum030402=-69786988544000*sqrt35;
  phisnum000602=14746124288000000*sqrt35;
  phisnum010602=-1203967098880000*sqrt35;
  phisnum020004=2589563289600000*sqrt35;
  phisnum030004=-91713699840000*sqrt35;
  phisnum040004=2982461440000*sqrt35;
  phisnum050004=90392064000*sqrt35;
  phisnum000204=207165063168000000*sqrt35;
  phisnum010204=-11221440921600000*sqrt35;
  phisnum020204=1326701608960000*sqrt35;
  phisnum030204=-82570027008000*sqrt35;
  phisnum000404=33494768025600000*sqrt35;
  phisnum010404=-2441909370880000*sqrt35;
  phisnum000006=78920024064000000*sqrt35;
  phisnum010006=-4439251353600000*sqrt35;
  phisnum020006=403334758400000*sqrt35;
  phisnum030006=-37639225344000*sqrt35;
  phisnum000206=21785214976000000*sqrt35;
  phisnum010206=-1804216565760000*sqrt35;
  phisnum000008=2301834035200000*sqrt35;
  phisnum010008=-320738426880000*sqrt35;
}

inline double phisnum(const double dr, const double dtheta, const double dphi)
{
  double phisnum;

  phisnum = dphi6*phisnum000006 + dphi8*phisnum000008 + dphi4*dtheta2*phisnum000204 + dphi6*dtheta2*phisnum000206 + dphi2*dtheta4*phisnum000402 + dphi4*dtheta4*phisnum000404 + dtheta6*phisnum000600 + dphi2*dtheta6*phisnum000602 + dtheta8*phisnum000800 + dphi6*dr1*phisnum010006 + dphi8*dr1*phisnum010008 + dphi4*dr1*dtheta2*phisnum010204 + dphi6*dr1*dtheta2*phisnum010206 + dphi2*dr1*dtheta4*phisnum010402 + dphi4*dr1*dtheta4*phisnum010404 + dr1*dtheta6*phisnum010600 + dphi2*dr1*dtheta6*phisnum010602 + dr1*dtheta8*phisnum010800 + dphi4*dr2*phisnum020004 + dphi6*dr2*phisnum020006 + dphi2*dr2*dtheta2*phisnum020202 + dphi4*dr2*dtheta2*phisnum020204 + dr2*dtheta4*phisnum020400 + dphi2*dr2*dtheta4*phisnum020402 + dr2*dtheta6*phisnum020600 + dphi4*dr3*phisnum030004 + dphi6*dr3*phisnum030006 + dphi2*dr3*dtheta2*phisnum030202 + dphi4*dr3*dtheta2*phisnum030204 + dr3*dtheta4*phisnum030400 + dphi2*dr3*dtheta4*phisnum030402 + dr3*dtheta6*phisnum030600 + dphi2*dr4*phisnum040002 + dphi4*dr4*phisnum040004 + dr4*dtheta2*phisnum040200 + dphi2*dr4*dtheta2*phisnum040202 + dr4*dtheta4*phisnum040400 + dphi2*dr5*phisnum050002 + dphi4*dr5*phisnum050004 + dr5*dtheta2*phisnum050200 + dphi2*dr5*dtheta2*phisnum050202 + dr5*dtheta4*phisnum050400 + dr6*phisnum060000 + dphi2*dr6*phisnum060002 + dr6*dtheta2*phisnum060200 + dr7*phisnum070000 + dphi2*dr7*phisnum070002 + dr7*dtheta2*phisnum070200 + dr8*phisnum080000 + dr9*phisnum090000;

  return phisnum;
}

static double phisden000000;

static void phisden_init(const double r1)
{
  phisden000000=Power(Power(7*Power(\[CapitalDelta]r,2) + 80*(7*Power(\[CapitalDelta]\[Theta],2) + 8*Power(\[CapitalDelta]\[Phi],2)),3.5),0.2857142857142857);
}

inline double phisden(const double dr, const double dtheta, const double dphi)
{
  double phisden;

  phisden = phisden000000;

  return phisden;
}


