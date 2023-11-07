#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6690086502958356170) {
   out_6690086502958356170[0] = delta_x[0] + nom_x[0];
   out_6690086502958356170[1] = delta_x[1] + nom_x[1];
   out_6690086502958356170[2] = delta_x[2] + nom_x[2];
   out_6690086502958356170[3] = delta_x[3] + nom_x[3];
   out_6690086502958356170[4] = delta_x[4] + nom_x[4];
   out_6690086502958356170[5] = delta_x[5] + nom_x[5];
   out_6690086502958356170[6] = delta_x[6] + nom_x[6];
   out_6690086502958356170[7] = delta_x[7] + nom_x[7];
   out_6690086502958356170[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1138017061106975654) {
   out_1138017061106975654[0] = -nom_x[0] + true_x[0];
   out_1138017061106975654[1] = -nom_x[1] + true_x[1];
   out_1138017061106975654[2] = -nom_x[2] + true_x[2];
   out_1138017061106975654[3] = -nom_x[3] + true_x[3];
   out_1138017061106975654[4] = -nom_x[4] + true_x[4];
   out_1138017061106975654[5] = -nom_x[5] + true_x[5];
   out_1138017061106975654[6] = -nom_x[6] + true_x[6];
   out_1138017061106975654[7] = -nom_x[7] + true_x[7];
   out_1138017061106975654[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_8731402594447162530) {
   out_8731402594447162530[0] = 1.0;
   out_8731402594447162530[1] = 0;
   out_8731402594447162530[2] = 0;
   out_8731402594447162530[3] = 0;
   out_8731402594447162530[4] = 0;
   out_8731402594447162530[5] = 0;
   out_8731402594447162530[6] = 0;
   out_8731402594447162530[7] = 0;
   out_8731402594447162530[8] = 0;
   out_8731402594447162530[9] = 0;
   out_8731402594447162530[10] = 1.0;
   out_8731402594447162530[11] = 0;
   out_8731402594447162530[12] = 0;
   out_8731402594447162530[13] = 0;
   out_8731402594447162530[14] = 0;
   out_8731402594447162530[15] = 0;
   out_8731402594447162530[16] = 0;
   out_8731402594447162530[17] = 0;
   out_8731402594447162530[18] = 0;
   out_8731402594447162530[19] = 0;
   out_8731402594447162530[20] = 1.0;
   out_8731402594447162530[21] = 0;
   out_8731402594447162530[22] = 0;
   out_8731402594447162530[23] = 0;
   out_8731402594447162530[24] = 0;
   out_8731402594447162530[25] = 0;
   out_8731402594447162530[26] = 0;
   out_8731402594447162530[27] = 0;
   out_8731402594447162530[28] = 0;
   out_8731402594447162530[29] = 0;
   out_8731402594447162530[30] = 1.0;
   out_8731402594447162530[31] = 0;
   out_8731402594447162530[32] = 0;
   out_8731402594447162530[33] = 0;
   out_8731402594447162530[34] = 0;
   out_8731402594447162530[35] = 0;
   out_8731402594447162530[36] = 0;
   out_8731402594447162530[37] = 0;
   out_8731402594447162530[38] = 0;
   out_8731402594447162530[39] = 0;
   out_8731402594447162530[40] = 1.0;
   out_8731402594447162530[41] = 0;
   out_8731402594447162530[42] = 0;
   out_8731402594447162530[43] = 0;
   out_8731402594447162530[44] = 0;
   out_8731402594447162530[45] = 0;
   out_8731402594447162530[46] = 0;
   out_8731402594447162530[47] = 0;
   out_8731402594447162530[48] = 0;
   out_8731402594447162530[49] = 0;
   out_8731402594447162530[50] = 1.0;
   out_8731402594447162530[51] = 0;
   out_8731402594447162530[52] = 0;
   out_8731402594447162530[53] = 0;
   out_8731402594447162530[54] = 0;
   out_8731402594447162530[55] = 0;
   out_8731402594447162530[56] = 0;
   out_8731402594447162530[57] = 0;
   out_8731402594447162530[58] = 0;
   out_8731402594447162530[59] = 0;
   out_8731402594447162530[60] = 1.0;
   out_8731402594447162530[61] = 0;
   out_8731402594447162530[62] = 0;
   out_8731402594447162530[63] = 0;
   out_8731402594447162530[64] = 0;
   out_8731402594447162530[65] = 0;
   out_8731402594447162530[66] = 0;
   out_8731402594447162530[67] = 0;
   out_8731402594447162530[68] = 0;
   out_8731402594447162530[69] = 0;
   out_8731402594447162530[70] = 1.0;
   out_8731402594447162530[71] = 0;
   out_8731402594447162530[72] = 0;
   out_8731402594447162530[73] = 0;
   out_8731402594447162530[74] = 0;
   out_8731402594447162530[75] = 0;
   out_8731402594447162530[76] = 0;
   out_8731402594447162530[77] = 0;
   out_8731402594447162530[78] = 0;
   out_8731402594447162530[79] = 0;
   out_8731402594447162530[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_5880714312440560991) {
   out_5880714312440560991[0] = state[0];
   out_5880714312440560991[1] = state[1];
   out_5880714312440560991[2] = state[2];
   out_5880714312440560991[3] = state[3];
   out_5880714312440560991[4] = state[4];
   out_5880714312440560991[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5880714312440560991[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5880714312440560991[7] = state[7];
   out_5880714312440560991[8] = state[8];
}
void F_fun(double *state, double dt, double *out_6668270723209344525) {
   out_6668270723209344525[0] = 1;
   out_6668270723209344525[1] = 0;
   out_6668270723209344525[2] = 0;
   out_6668270723209344525[3] = 0;
   out_6668270723209344525[4] = 0;
   out_6668270723209344525[5] = 0;
   out_6668270723209344525[6] = 0;
   out_6668270723209344525[7] = 0;
   out_6668270723209344525[8] = 0;
   out_6668270723209344525[9] = 0;
   out_6668270723209344525[10] = 1;
   out_6668270723209344525[11] = 0;
   out_6668270723209344525[12] = 0;
   out_6668270723209344525[13] = 0;
   out_6668270723209344525[14] = 0;
   out_6668270723209344525[15] = 0;
   out_6668270723209344525[16] = 0;
   out_6668270723209344525[17] = 0;
   out_6668270723209344525[18] = 0;
   out_6668270723209344525[19] = 0;
   out_6668270723209344525[20] = 1;
   out_6668270723209344525[21] = 0;
   out_6668270723209344525[22] = 0;
   out_6668270723209344525[23] = 0;
   out_6668270723209344525[24] = 0;
   out_6668270723209344525[25] = 0;
   out_6668270723209344525[26] = 0;
   out_6668270723209344525[27] = 0;
   out_6668270723209344525[28] = 0;
   out_6668270723209344525[29] = 0;
   out_6668270723209344525[30] = 1;
   out_6668270723209344525[31] = 0;
   out_6668270723209344525[32] = 0;
   out_6668270723209344525[33] = 0;
   out_6668270723209344525[34] = 0;
   out_6668270723209344525[35] = 0;
   out_6668270723209344525[36] = 0;
   out_6668270723209344525[37] = 0;
   out_6668270723209344525[38] = 0;
   out_6668270723209344525[39] = 0;
   out_6668270723209344525[40] = 1;
   out_6668270723209344525[41] = 0;
   out_6668270723209344525[42] = 0;
   out_6668270723209344525[43] = 0;
   out_6668270723209344525[44] = 0;
   out_6668270723209344525[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6668270723209344525[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6668270723209344525[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6668270723209344525[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6668270723209344525[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6668270723209344525[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6668270723209344525[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6668270723209344525[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6668270723209344525[53] = -9.8000000000000007*dt;
   out_6668270723209344525[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6668270723209344525[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6668270723209344525[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6668270723209344525[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6668270723209344525[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6668270723209344525[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6668270723209344525[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6668270723209344525[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6668270723209344525[62] = 0;
   out_6668270723209344525[63] = 0;
   out_6668270723209344525[64] = 0;
   out_6668270723209344525[65] = 0;
   out_6668270723209344525[66] = 0;
   out_6668270723209344525[67] = 0;
   out_6668270723209344525[68] = 0;
   out_6668270723209344525[69] = 0;
   out_6668270723209344525[70] = 1;
   out_6668270723209344525[71] = 0;
   out_6668270723209344525[72] = 0;
   out_6668270723209344525[73] = 0;
   out_6668270723209344525[74] = 0;
   out_6668270723209344525[75] = 0;
   out_6668270723209344525[76] = 0;
   out_6668270723209344525[77] = 0;
   out_6668270723209344525[78] = 0;
   out_6668270723209344525[79] = 0;
   out_6668270723209344525[80] = 1;
}
void h_25(double *state, double *unused, double *out_6309852841180432979) {
   out_6309852841180432979[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6040165476284428233) {
   out_6040165476284428233[0] = 0;
   out_6040165476284428233[1] = 0;
   out_6040165476284428233[2] = 0;
   out_6040165476284428233[3] = 0;
   out_6040165476284428233[4] = 0;
   out_6040165476284428233[5] = 0;
   out_6040165476284428233[6] = 1;
   out_6040165476284428233[7] = 0;
   out_6040165476284428233[8] = 0;
}
void h_24(double *state, double *unused, double *out_7716672714826793939) {
   out_7716672714826793939[0] = state[4];
   out_7716672714826793939[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2482239877482105097) {
   out_2482239877482105097[0] = 0;
   out_2482239877482105097[1] = 0;
   out_2482239877482105097[2] = 0;
   out_2482239877482105097[3] = 0;
   out_2482239877482105097[4] = 1;
   out_2482239877482105097[5] = 0;
   out_2482239877482105097[6] = 0;
   out_2482239877482105097[7] = 0;
   out_2482239877482105097[8] = 0;
   out_2482239877482105097[9] = 0;
   out_2482239877482105097[10] = 0;
   out_2482239877482105097[11] = 0;
   out_2482239877482105097[12] = 0;
   out_2482239877482105097[13] = 0;
   out_2482239877482105097[14] = 1;
   out_2482239877482105097[15] = 0;
   out_2482239877482105097[16] = 0;
   out_2482239877482105097[17] = 0;
}
void h_30(double *state, double *unused, double *out_2554258583712588207) {
   out_2554258583712588207[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7878882267297515185) {
   out_7878882267297515185[0] = 0;
   out_7878882267297515185[1] = 0;
   out_7878882267297515185[2] = 0;
   out_7878882267297515185[3] = 0;
   out_7878882267297515185[4] = 1;
   out_7878882267297515185[5] = 0;
   out_7878882267297515185[6] = 0;
   out_7878882267297515185[7] = 0;
   out_7878882267297515185[8] = 0;
}
void h_26(double *state, double *unused, double *out_7555931851373433852) {
   out_7555931851373433852[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8665075278551067159) {
   out_8665075278551067159[0] = 0;
   out_8665075278551067159[1] = 0;
   out_8665075278551067159[2] = 0;
   out_8665075278551067159[3] = 0;
   out_8665075278551067159[4] = 0;
   out_8665075278551067159[5] = 0;
   out_8665075278551067159[6] = 0;
   out_8665075278551067159[7] = 1;
   out_8665075278551067159[8] = 0;
}
void h_27(double *state, double *unused, double *out_3273066789736100294) {
   out_3273066789736100294[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5704118955497090274) {
   out_5704118955497090274[0] = 0;
   out_5704118955497090274[1] = 0;
   out_5704118955497090274[2] = 0;
   out_5704118955497090274[3] = 1;
   out_5704118955497090274[4] = 0;
   out_5704118955497090274[5] = 0;
   out_5704118955497090274[6] = 0;
   out_5704118955497090274[7] = 0;
   out_5704118955497090274[8] = 0;
}
void h_29(double *state, double *unused, double *out_3548260852020606183) {
   out_3548260852020606183[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8389113611611907369) {
   out_8389113611611907369[0] = 0;
   out_8389113611611907369[1] = 1;
   out_8389113611611907369[2] = 0;
   out_8389113611611907369[3] = 0;
   out_8389113611611907369[4] = 0;
   out_8389113611611907369[5] = 0;
   out_8389113611611907369[6] = 0;
   out_8389113611611907369[7] = 0;
   out_8389113611611907369[8] = 0;
}
void h_28(double *state, double *unused, double *out_668661468927402643) {
   out_668661468927402643[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3306714594542376795) {
   out_3306714594542376795[0] = 1;
   out_3306714594542376795[1] = 0;
   out_3306714594542376795[2] = 0;
   out_3306714594542376795[3] = 0;
   out_3306714594542376795[4] = 0;
   out_3306714594542376795[5] = 0;
   out_3306714594542376795[6] = 0;
   out_3306714594542376795[7] = 0;
   out_3306714594542376795[8] = 0;
}
void h_31(double *state, double *unused, double *out_6034658778895927090) {
   out_6034658778895927090[0] = state[8];
}
void H_31(double *state, double *unused, double *out_8038867176317715683) {
   out_8038867176317715683[0] = 0;
   out_8038867176317715683[1] = 0;
   out_8038867176317715683[2] = 0;
   out_8038867176317715683[3] = 0;
   out_8038867176317715683[4] = 0;
   out_8038867176317715683[5] = 0;
   out_8038867176317715683[6] = 0;
   out_8038867176317715683[7] = 0;
   out_8038867176317715683[8] = 1;
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

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_6690086502958356170) {
  err_fun(nom_x, delta_x, out_6690086502958356170);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1138017061106975654) {
  inv_err_fun(nom_x, true_x, out_1138017061106975654);
}
void car_H_mod_fun(double *state, double *out_8731402594447162530) {
  H_mod_fun(state, out_8731402594447162530);
}
void car_f_fun(double *state, double dt, double *out_5880714312440560991) {
  f_fun(state,  dt, out_5880714312440560991);
}
void car_F_fun(double *state, double dt, double *out_6668270723209344525) {
  F_fun(state,  dt, out_6668270723209344525);
}
void car_h_25(double *state, double *unused, double *out_6309852841180432979) {
  h_25(state, unused, out_6309852841180432979);
}
void car_H_25(double *state, double *unused, double *out_6040165476284428233) {
  H_25(state, unused, out_6040165476284428233);
}
void car_h_24(double *state, double *unused, double *out_7716672714826793939) {
  h_24(state, unused, out_7716672714826793939);
}
void car_H_24(double *state, double *unused, double *out_2482239877482105097) {
  H_24(state, unused, out_2482239877482105097);
}
void car_h_30(double *state, double *unused, double *out_2554258583712588207) {
  h_30(state, unused, out_2554258583712588207);
}
void car_H_30(double *state, double *unused, double *out_7878882267297515185) {
  H_30(state, unused, out_7878882267297515185);
}
void car_h_26(double *state, double *unused, double *out_7555931851373433852) {
  h_26(state, unused, out_7555931851373433852);
}
void car_H_26(double *state, double *unused, double *out_8665075278551067159) {
  H_26(state, unused, out_8665075278551067159);
}
void car_h_27(double *state, double *unused, double *out_3273066789736100294) {
  h_27(state, unused, out_3273066789736100294);
}
void car_H_27(double *state, double *unused, double *out_5704118955497090274) {
  H_27(state, unused, out_5704118955497090274);
}
void car_h_29(double *state, double *unused, double *out_3548260852020606183) {
  h_29(state, unused, out_3548260852020606183);
}
void car_H_29(double *state, double *unused, double *out_8389113611611907369) {
  H_29(state, unused, out_8389113611611907369);
}
void car_h_28(double *state, double *unused, double *out_668661468927402643) {
  h_28(state, unused, out_668661468927402643);
}
void car_H_28(double *state, double *unused, double *out_3306714594542376795) {
  H_28(state, unused, out_3306714594542376795);
}
void car_h_31(double *state, double *unused, double *out_6034658778895927090) {
  h_31(state, unused, out_6034658778895927090);
}
void car_H_31(double *state, double *unused, double *out_8038867176317715683) {
  H_31(state, unused, out_8038867176317715683);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);