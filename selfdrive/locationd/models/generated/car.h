#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3699073491920453798);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5512259151403388105);
void car_H_mod_fun(double *state, double *out_6004474647322099520);
void car_f_fun(double *state, double dt, double *out_5290008103764663250);
void car_F_fun(double *state, double dt, double *out_2929614667293202327);
void car_h_25(double *state, double *unused, double *out_2763657023161734409);
void car_H_25(double *state, double *unused, double *out_5353771594371325908);
void car_h_24(double *state, double *unused, double *out_732815993837596622);
void car_H_24(double *state, double *unused, double *out_7526421193376825474);
void car_h_30(double *state, double *unused, double *out_3038851085446240298);
void car_H_30(double *state, double *unused, double *out_2835438635864077281);
void car_h_26(double *state, double *unused, double *out_3982360938347902733);
void car_H_26(double *state, double *unused, double *out_9095274913245382132);
void car_h_27(double *state, double *unused, double *out_7765330290822580054);
void car_H_27(double *state, double *unused, double *out_5010201947664502192);
void car_h_29(double *state, double *unused, double *out_2065034001907613);
void car_H_29(double *state, double *unused, double *out_2325207291549685097);
void car_h_28(double *state, double *unused, double *out_762567587029514275);
void car_H_28(double *state, double *unused, double *out_7407606308619215671);
void car_h_31(double *state, double *unused, double *out_6861960321441106273);
void car_H_31(double *state, double *unused, double *out_8725261058230818008);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}