#include "gnss.h"

namespace {
#define DIM 11
#define EDIM 11
#define MEDIM 11
typedef void (*Hfun)(double *, double *, double *);
const static double MAHA_THRESH_6 = 3.8414588206941227;
const static double MAHA_THRESH_20 = 3.8414588206941227;
const static double MAHA_THRESH_7 = 3.8414588206941227;
const static double MAHA_THRESH_21 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4085845534313808973) {
   out_4085845534313808973[0] = delta_x[0] + nom_x[0];
   out_4085845534313808973[1] = delta_x[1] + nom_x[1];
   out_4085845534313808973[2] = delta_x[2] + nom_x[2];
   out_4085845534313808973[3] = delta_x[3] + nom_x[3];
   out_4085845534313808973[4] = delta_x[4] + nom_x[4];
   out_4085845534313808973[5] = delta_x[5] + nom_x[5];
   out_4085845534313808973[6] = delta_x[6] + nom_x[6];
   out_4085845534313808973[7] = delta_x[7] + nom_x[7];
   out_4085845534313808973[8] = delta_x[8] + nom_x[8];
   out_4085845534313808973[9] = delta_x[9] + nom_x[9];
   out_4085845534313808973[10] = delta_x[10] + nom_x[10];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_9087850996606955167) {
   out_9087850996606955167[0] = -nom_x[0] + true_x[0];
   out_9087850996606955167[1] = -nom_x[1] + true_x[1];
   out_9087850996606955167[2] = -nom_x[2] + true_x[2];
   out_9087850996606955167[3] = -nom_x[3] + true_x[3];
   out_9087850996606955167[4] = -nom_x[4] + true_x[4];
   out_9087850996606955167[5] = -nom_x[5] + true_x[5];
   out_9087850996606955167[6] = -nom_x[6] + true_x[6];
   out_9087850996606955167[7] = -nom_x[7] + true_x[7];
   out_9087850996606955167[8] = -nom_x[8] + true_x[8];
   out_9087850996606955167[9] = -nom_x[9] + true_x[9];
   out_9087850996606955167[10] = -nom_x[10] + true_x[10];
}
void H_mod_fun(double *state, double *out_7803289466902203486) {
   out_7803289466902203486[0] = 1.0;
   out_7803289466902203486[1] = 0;
   out_7803289466902203486[2] = 0;
   out_7803289466902203486[3] = 0;
   out_7803289466902203486[4] = 0;
   out_7803289466902203486[5] = 0;
   out_7803289466902203486[6] = 0;
   out_7803289466902203486[7] = 0;
   out_7803289466902203486[8] = 0;
   out_7803289466902203486[9] = 0;
   out_7803289466902203486[10] = 0;
   out_7803289466902203486[11] = 0;
   out_7803289466902203486[12] = 1.0;
   out_7803289466902203486[13] = 0;
   out_7803289466902203486[14] = 0;
   out_7803289466902203486[15] = 0;
   out_7803289466902203486[16] = 0;
   out_7803289466902203486[17] = 0;
   out_7803289466902203486[18] = 0;
   out_7803289466902203486[19] = 0;
   out_7803289466902203486[20] = 0;
   out_7803289466902203486[21] = 0;
   out_7803289466902203486[22] = 0;
   out_7803289466902203486[23] = 0;
   out_7803289466902203486[24] = 1.0;
   out_7803289466902203486[25] = 0;
   out_7803289466902203486[26] = 0;
   out_7803289466902203486[27] = 0;
   out_7803289466902203486[28] = 0;
   out_7803289466902203486[29] = 0;
   out_7803289466902203486[30] = 0;
   out_7803289466902203486[31] = 0;
   out_7803289466902203486[32] = 0;
   out_7803289466902203486[33] = 0;
   out_7803289466902203486[34] = 0;
   out_7803289466902203486[35] = 0;
   out_7803289466902203486[36] = 1.0;
   out_7803289466902203486[37] = 0;
   out_7803289466902203486[38] = 0;
   out_7803289466902203486[39] = 0;
   out_7803289466902203486[40] = 0;
   out_7803289466902203486[41] = 0;
   out_7803289466902203486[42] = 0;
   out_7803289466902203486[43] = 0;
   out_7803289466902203486[44] = 0;
   out_7803289466902203486[45] = 0;
   out_7803289466902203486[46] = 0;
   out_7803289466902203486[47] = 0;
   out_7803289466902203486[48] = 1.0;
   out_7803289466902203486[49] = 0;
   out_7803289466902203486[50] = 0;
   out_7803289466902203486[51] = 0;
   out_7803289466902203486[52] = 0;
   out_7803289466902203486[53] = 0;
   out_7803289466902203486[54] = 0;
   out_7803289466902203486[55] = 0;
   out_7803289466902203486[56] = 0;
   out_7803289466902203486[57] = 0;
   out_7803289466902203486[58] = 0;
   out_7803289466902203486[59] = 0;
   out_7803289466902203486[60] = 1.0;
   out_7803289466902203486[61] = 0;
   out_7803289466902203486[62] = 0;
   out_7803289466902203486[63] = 0;
   out_7803289466902203486[64] = 0;
   out_7803289466902203486[65] = 0;
   out_7803289466902203486[66] = 0;
   out_7803289466902203486[67] = 0;
   out_7803289466902203486[68] = 0;
   out_7803289466902203486[69] = 0;
   out_7803289466902203486[70] = 0;
   out_7803289466902203486[71] = 0;
   out_7803289466902203486[72] = 1.0;
   out_7803289466902203486[73] = 0;
   out_7803289466902203486[74] = 0;
   out_7803289466902203486[75] = 0;
   out_7803289466902203486[76] = 0;
   out_7803289466902203486[77] = 0;
   out_7803289466902203486[78] = 0;
   out_7803289466902203486[79] = 0;
   out_7803289466902203486[80] = 0;
   out_7803289466902203486[81] = 0;
   out_7803289466902203486[82] = 0;
   out_7803289466902203486[83] = 0;
   out_7803289466902203486[84] = 1.0;
   out_7803289466902203486[85] = 0;
   out_7803289466902203486[86] = 0;
   out_7803289466902203486[87] = 0;
   out_7803289466902203486[88] = 0;
   out_7803289466902203486[89] = 0;
   out_7803289466902203486[90] = 0;
   out_7803289466902203486[91] = 0;
   out_7803289466902203486[92] = 0;
   out_7803289466902203486[93] = 0;
   out_7803289466902203486[94] = 0;
   out_7803289466902203486[95] = 0;
   out_7803289466902203486[96] = 1.0;
   out_7803289466902203486[97] = 0;
   out_7803289466902203486[98] = 0;
   out_7803289466902203486[99] = 0;
   out_7803289466902203486[100] = 0;
   out_7803289466902203486[101] = 0;
   out_7803289466902203486[102] = 0;
   out_7803289466902203486[103] = 0;
   out_7803289466902203486[104] = 0;
   out_7803289466902203486[105] = 0;
   out_7803289466902203486[106] = 0;
   out_7803289466902203486[107] = 0;
   out_7803289466902203486[108] = 1.0;
   out_7803289466902203486[109] = 0;
   out_7803289466902203486[110] = 0;
   out_7803289466902203486[111] = 0;
   out_7803289466902203486[112] = 0;
   out_7803289466902203486[113] = 0;
   out_7803289466902203486[114] = 0;
   out_7803289466902203486[115] = 0;
   out_7803289466902203486[116] = 0;
   out_7803289466902203486[117] = 0;
   out_7803289466902203486[118] = 0;
   out_7803289466902203486[119] = 0;
   out_7803289466902203486[120] = 1.0;
}
void f_fun(double *state, double dt, double *out_4982008935116428445) {
   out_4982008935116428445[0] = dt*state[3] + state[0];
   out_4982008935116428445[1] = dt*state[4] + state[1];
   out_4982008935116428445[2] = dt*state[5] + state[2];
   out_4982008935116428445[3] = state[3];
   out_4982008935116428445[4] = state[4];
   out_4982008935116428445[5] = state[5];
   out_4982008935116428445[6] = dt*state[7] + state[6];
   out_4982008935116428445[7] = dt*state[8] + state[7];
   out_4982008935116428445[8] = state[8];
   out_4982008935116428445[9] = state[9];
   out_4982008935116428445[10] = state[10];
}
void F_fun(double *state, double dt, double *out_7365073110426550313) {
   out_7365073110426550313[0] = 1;
   out_7365073110426550313[1] = 0;
   out_7365073110426550313[2] = 0;
   out_7365073110426550313[3] = dt;
   out_7365073110426550313[4] = 0;
   out_7365073110426550313[5] = 0;
   out_7365073110426550313[6] = 0;
   out_7365073110426550313[7] = 0;
   out_7365073110426550313[8] = 0;
   out_7365073110426550313[9] = 0;
   out_7365073110426550313[10] = 0;
   out_7365073110426550313[11] = 0;
   out_7365073110426550313[12] = 1;
   out_7365073110426550313[13] = 0;
   out_7365073110426550313[14] = 0;
   out_7365073110426550313[15] = dt;
   out_7365073110426550313[16] = 0;
   out_7365073110426550313[17] = 0;
   out_7365073110426550313[18] = 0;
   out_7365073110426550313[19] = 0;
   out_7365073110426550313[20] = 0;
   out_7365073110426550313[21] = 0;
   out_7365073110426550313[22] = 0;
   out_7365073110426550313[23] = 0;
   out_7365073110426550313[24] = 1;
   out_7365073110426550313[25] = 0;
   out_7365073110426550313[26] = 0;
   out_7365073110426550313[27] = dt;
   out_7365073110426550313[28] = 0;
   out_7365073110426550313[29] = 0;
   out_7365073110426550313[30] = 0;
   out_7365073110426550313[31] = 0;
   out_7365073110426550313[32] = 0;
   out_7365073110426550313[33] = 0;
   out_7365073110426550313[34] = 0;
   out_7365073110426550313[35] = 0;
   out_7365073110426550313[36] = 1;
   out_7365073110426550313[37] = 0;
   out_7365073110426550313[38] = 0;
   out_7365073110426550313[39] = 0;
   out_7365073110426550313[40] = 0;
   out_7365073110426550313[41] = 0;
   out_7365073110426550313[42] = 0;
   out_7365073110426550313[43] = 0;
   out_7365073110426550313[44] = 0;
   out_7365073110426550313[45] = 0;
   out_7365073110426550313[46] = 0;
   out_7365073110426550313[47] = 0;
   out_7365073110426550313[48] = 1;
   out_7365073110426550313[49] = 0;
   out_7365073110426550313[50] = 0;
   out_7365073110426550313[51] = 0;
   out_7365073110426550313[52] = 0;
   out_7365073110426550313[53] = 0;
   out_7365073110426550313[54] = 0;
   out_7365073110426550313[55] = 0;
   out_7365073110426550313[56] = 0;
   out_7365073110426550313[57] = 0;
   out_7365073110426550313[58] = 0;
   out_7365073110426550313[59] = 0;
   out_7365073110426550313[60] = 1;
   out_7365073110426550313[61] = 0;
   out_7365073110426550313[62] = 0;
   out_7365073110426550313[63] = 0;
   out_7365073110426550313[64] = 0;
   out_7365073110426550313[65] = 0;
   out_7365073110426550313[66] = 0;
   out_7365073110426550313[67] = 0;
   out_7365073110426550313[68] = 0;
   out_7365073110426550313[69] = 0;
   out_7365073110426550313[70] = 0;
   out_7365073110426550313[71] = 0;
   out_7365073110426550313[72] = 1;
   out_7365073110426550313[73] = dt;
   out_7365073110426550313[74] = 0;
   out_7365073110426550313[75] = 0;
   out_7365073110426550313[76] = 0;
   out_7365073110426550313[77] = 0;
   out_7365073110426550313[78] = 0;
   out_7365073110426550313[79] = 0;
   out_7365073110426550313[80] = 0;
   out_7365073110426550313[81] = 0;
   out_7365073110426550313[82] = 0;
   out_7365073110426550313[83] = 0;
   out_7365073110426550313[84] = 1;
   out_7365073110426550313[85] = dt;
   out_7365073110426550313[86] = 0;
   out_7365073110426550313[87] = 0;
   out_7365073110426550313[88] = 0;
   out_7365073110426550313[89] = 0;
   out_7365073110426550313[90] = 0;
   out_7365073110426550313[91] = 0;
   out_7365073110426550313[92] = 0;
   out_7365073110426550313[93] = 0;
   out_7365073110426550313[94] = 0;
   out_7365073110426550313[95] = 0;
   out_7365073110426550313[96] = 1;
   out_7365073110426550313[97] = 0;
   out_7365073110426550313[98] = 0;
   out_7365073110426550313[99] = 0;
   out_7365073110426550313[100] = 0;
   out_7365073110426550313[101] = 0;
   out_7365073110426550313[102] = 0;
   out_7365073110426550313[103] = 0;
   out_7365073110426550313[104] = 0;
   out_7365073110426550313[105] = 0;
   out_7365073110426550313[106] = 0;
   out_7365073110426550313[107] = 0;
   out_7365073110426550313[108] = 1;
   out_7365073110426550313[109] = 0;
   out_7365073110426550313[110] = 0;
   out_7365073110426550313[111] = 0;
   out_7365073110426550313[112] = 0;
   out_7365073110426550313[113] = 0;
   out_7365073110426550313[114] = 0;
   out_7365073110426550313[115] = 0;
   out_7365073110426550313[116] = 0;
   out_7365073110426550313[117] = 0;
   out_7365073110426550313[118] = 0;
   out_7365073110426550313[119] = 0;
   out_7365073110426550313[120] = 1;
}
void h_6(double *state, double *sat_pos, double *out_3494771982687573302) {
   out_3494771982687573302[0] = sqrt(pow(-sat_pos[0] + state[0], 2) + pow(-sat_pos[1] + state[1], 2) + pow(-sat_pos[2] + state[2], 2)) + state[6];
}
void H_6(double *state, double *sat_pos, double *out_129949877751880805) {
   out_129949877751880805[0] = (-sat_pos[0] + state[0])/sqrt(pow(-sat_pos[0] + state[0], 2) + pow(-sat_pos[1] + state[1], 2) + pow(-sat_pos[2] + state[2], 2));
   out_129949877751880805[1] = (-sat_pos[1] + state[1])/sqrt(pow(-sat_pos[0] + state[0], 2) + pow(-sat_pos[1] + state[1], 2) + pow(-sat_pos[2] + state[2], 2));
   out_129949877751880805[2] = (-sat_pos[2] + state[2])/sqrt(pow(-sat_pos[0] + state[0], 2) + pow(-sat_pos[1] + state[1], 2) + pow(-sat_pos[2] + state[2], 2));
   out_129949877751880805[3] = 0;
   out_129949877751880805[4] = 0;
   out_129949877751880805[5] = 0;
   out_129949877751880805[6] = 1;
   out_129949877751880805[7] = 0;
   out_129949877751880805[8] = 0;
   out_129949877751880805[9] = 0;
   out_129949877751880805[10] = 0;
}
void h_20(double *state, double *sat_pos, double *out_3726883287393159815) {
   out_3726883287393159815[0] = sqrt(pow(-sat_pos[0] + state[0], 2) + pow(-sat_pos[1] + state[1], 2) + pow(-sat_pos[2] + state[2], 2)) + sat_pos[3]*state[10] + state[6] + state[9];
}
void H_20(double *state, double *sat_pos, double *out_3913751108416185265) {
   out_3913751108416185265[0] = (-sat_pos[0] + state[0])/sqrt(pow(-sat_pos[0] + state[0], 2) + pow(-sat_pos[1] + state[1], 2) + pow(-sat_pos[2] + state[2], 2));
   out_3913751108416185265[1] = (-sat_pos[1] + state[1])/sqrt(pow(-sat_pos[0] + state[0], 2) + pow(-sat_pos[1] + state[1], 2) + pow(-sat_pos[2] + state[2], 2));
   out_3913751108416185265[2] = (-sat_pos[2] + state[2])/sqrt(pow(-sat_pos[0] + state[0], 2) + pow(-sat_pos[1] + state[1], 2) + pow(-sat_pos[2] + state[2], 2));
   out_3913751108416185265[3] = 0;
   out_3913751108416185265[4] = 0;
   out_3913751108416185265[5] = 0;
   out_3913751108416185265[6] = 1;
   out_3913751108416185265[7] = 0;
   out_3913751108416185265[8] = 0;
   out_3913751108416185265[9] = 1;
   out_3913751108416185265[10] = sat_pos[3];
}
void h_7(double *state, double *sat_pos_vel, double *out_8983076725070806951) {
   out_8983076725070806951[0] = (sat_pos_vel[0] - state[0])*(sat_pos_vel[3] - state[3])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2)) + (sat_pos_vel[1] - state[1])*(sat_pos_vel[4] - state[4])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2)) + (sat_pos_vel[2] - state[2])*(sat_pos_vel[5] - state[5])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2)) + state[7];
}
void H_7(double *state, double *sat_pos_vel, double *out_6658591309417436314) {
   out_6658591309417436314[0] = pow(sat_pos_vel[0] - state[0], 2)*(sat_pos_vel[3] - state[3])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + (sat_pos_vel[0] - state[0])*(sat_pos_vel[1] - state[1])*(sat_pos_vel[4] - state[4])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + (sat_pos_vel[0] - state[0])*(sat_pos_vel[2] - state[2])*(sat_pos_vel[5] - state[5])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) - (sat_pos_vel[3] - state[3])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[1] = (sat_pos_vel[0] - state[0])*(sat_pos_vel[1] - state[1])*(sat_pos_vel[3] - state[3])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + pow(sat_pos_vel[1] - state[1], 2)*(sat_pos_vel[4] - state[4])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + (sat_pos_vel[1] - state[1])*(sat_pos_vel[2] - state[2])*(sat_pos_vel[5] - state[5])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) - (sat_pos_vel[4] - state[4])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[2] = (sat_pos_vel[0] - state[0])*(sat_pos_vel[2] - state[2])*(sat_pos_vel[3] - state[3])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + (sat_pos_vel[1] - state[1])*(sat_pos_vel[2] - state[2])*(sat_pos_vel[4] - state[4])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + pow(sat_pos_vel[2] - state[2], 2)*(sat_pos_vel[5] - state[5])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) - (sat_pos_vel[5] - state[5])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[3] = -(sat_pos_vel[0] - state[0])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[4] = -(sat_pos_vel[1] - state[1])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[5] = -(sat_pos_vel[2] - state[2])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[6] = 0;
   out_6658591309417436314[7] = 1;
   out_6658591309417436314[8] = 0;
   out_6658591309417436314[9] = 0;
   out_6658591309417436314[10] = 0;
}
void h_21(double *state, double *sat_pos_vel, double *out_8983076725070806951) {
   out_8983076725070806951[0] = (sat_pos_vel[0] - state[0])*(sat_pos_vel[3] - state[3])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2)) + (sat_pos_vel[1] - state[1])*(sat_pos_vel[4] - state[4])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2)) + (sat_pos_vel[2] - state[2])*(sat_pos_vel[5] - state[5])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2)) + state[7];
}
void H_21(double *state, double *sat_pos_vel, double *out_6658591309417436314) {
   out_6658591309417436314[0] = pow(sat_pos_vel[0] - state[0], 2)*(sat_pos_vel[3] - state[3])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + (sat_pos_vel[0] - state[0])*(sat_pos_vel[1] - state[1])*(sat_pos_vel[4] - state[4])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + (sat_pos_vel[0] - state[0])*(sat_pos_vel[2] - state[2])*(sat_pos_vel[5] - state[5])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) - (sat_pos_vel[3] - state[3])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[1] = (sat_pos_vel[0] - state[0])*(sat_pos_vel[1] - state[1])*(sat_pos_vel[3] - state[3])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + pow(sat_pos_vel[1] - state[1], 2)*(sat_pos_vel[4] - state[4])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + (sat_pos_vel[1] - state[1])*(sat_pos_vel[2] - state[2])*(sat_pos_vel[5] - state[5])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) - (sat_pos_vel[4] - state[4])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[2] = (sat_pos_vel[0] - state[0])*(sat_pos_vel[2] - state[2])*(sat_pos_vel[3] - state[3])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + (sat_pos_vel[1] - state[1])*(sat_pos_vel[2] - state[2])*(sat_pos_vel[4] - state[4])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) + pow(sat_pos_vel[2] - state[2], 2)*(sat_pos_vel[5] - state[5])/pow(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2), 3.0/2.0) - (sat_pos_vel[5] - state[5])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[3] = -(sat_pos_vel[0] - state[0])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[4] = -(sat_pos_vel[1] - state[1])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[5] = -(sat_pos_vel[2] - state[2])/sqrt(pow(sat_pos_vel[0] - state[0], 2) + pow(sat_pos_vel[1] - state[1], 2) + pow(sat_pos_vel[2] - state[2], 2));
   out_6658591309417436314[6] = 0;
   out_6658591309417436314[7] = 1;
   out_6658591309417436314[8] = 0;
   out_6658591309417436314[9] = 0;
   out_6658591309417436314[10] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void gnss_update_6(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_6, H_6, NULL, in_z, in_R, in_ea, MAHA_THRESH_6);
}
void gnss_update_20(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_20, H_20, NULL, in_z, in_R, in_ea, MAHA_THRESH_20);
}
void gnss_update_7(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_7, H_7, NULL, in_z, in_R, in_ea, MAHA_THRESH_7);
}
void gnss_update_21(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_21, H_21, NULL, in_z, in_R, in_ea, MAHA_THRESH_21);
}
void gnss_err_fun(double *nom_x, double *delta_x, double *out_4085845534313808973) {
  err_fun(nom_x, delta_x, out_4085845534313808973);
}
void gnss_inv_err_fun(double *nom_x, double *true_x, double *out_9087850996606955167) {
  inv_err_fun(nom_x, true_x, out_9087850996606955167);
}
void gnss_H_mod_fun(double *state, double *out_7803289466902203486) {
  H_mod_fun(state, out_7803289466902203486);
}
void gnss_f_fun(double *state, double dt, double *out_4982008935116428445) {
  f_fun(state,  dt, out_4982008935116428445);
}
void gnss_F_fun(double *state, double dt, double *out_7365073110426550313) {
  F_fun(state,  dt, out_7365073110426550313);
}
void gnss_h_6(double *state, double *sat_pos, double *out_3494771982687573302) {
  h_6(state, sat_pos, out_3494771982687573302);
}
void gnss_H_6(double *state, double *sat_pos, double *out_129949877751880805) {
  H_6(state, sat_pos, out_129949877751880805);
}
void gnss_h_20(double *state, double *sat_pos, double *out_3726883287393159815) {
  h_20(state, sat_pos, out_3726883287393159815);
}
void gnss_H_20(double *state, double *sat_pos, double *out_3913751108416185265) {
  H_20(state, sat_pos, out_3913751108416185265);
}
void gnss_h_7(double *state, double *sat_pos_vel, double *out_8983076725070806951) {
  h_7(state, sat_pos_vel, out_8983076725070806951);
}
void gnss_H_7(double *state, double *sat_pos_vel, double *out_6658591309417436314) {
  H_7(state, sat_pos_vel, out_6658591309417436314);
}
void gnss_h_21(double *state, double *sat_pos_vel, double *out_8983076725070806951) {
  h_21(state, sat_pos_vel, out_8983076725070806951);
}
void gnss_H_21(double *state, double *sat_pos_vel, double *out_6658591309417436314) {
  H_21(state, sat_pos_vel, out_6658591309417436314);
}
void gnss_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
}

const EKF gnss = {
  .name = "gnss",
  .kinds = { 6, 20, 7, 21 },
  .feature_kinds = {  },
  .f_fun = gnss_f_fun,
  .F_fun = gnss_F_fun,
  .err_fun = gnss_err_fun,
  .inv_err_fun = gnss_inv_err_fun,
  .H_mod_fun = gnss_H_mod_fun,
  .predict = gnss_predict,
  .hs = {
    { 6, gnss_h_6 },
    { 20, gnss_h_20 },
    { 7, gnss_h_7 },
    { 21, gnss_h_21 },
  },
  .Hs = {
    { 6, gnss_H_6 },
    { 20, gnss_H_20 },
    { 7, gnss_H_7 },
    { 21, gnss_H_21 },
  },
  .updates = {
    { 6, gnss_update_6 },
    { 20, gnss_update_20 },
    { 7, gnss_update_7 },
    { 21, gnss_update_21 },
  },
  .Hes = {
  },
  .sets = {
  },
  .extra_routines = {
  },
};

ekf_init(gnss);
