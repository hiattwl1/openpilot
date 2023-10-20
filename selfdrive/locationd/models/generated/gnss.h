#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void gnss_update_6(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void gnss_update_20(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void gnss_update_7(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void gnss_update_21(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void gnss_err_fun(double *nom_x, double *delta_x, double *out_4085845534313808973);
void gnss_inv_err_fun(double *nom_x, double *true_x, double *out_9087850996606955167);
void gnss_H_mod_fun(double *state, double *out_7803289466902203486);
void gnss_f_fun(double *state, double dt, double *out_4982008935116428445);
void gnss_F_fun(double *state, double dt, double *out_7365073110426550313);
void gnss_h_6(double *state, double *sat_pos, double *out_3494771982687573302);
void gnss_H_6(double *state, double *sat_pos, double *out_129949877751880805);
void gnss_h_20(double *state, double *sat_pos, double *out_3726883287393159815);
void gnss_H_20(double *state, double *sat_pos, double *out_3913751108416185265);
void gnss_h_7(double *state, double *sat_pos_vel, double *out_8983076725070806951);
void gnss_H_7(double *state, double *sat_pos_vel, double *out_6658591309417436314);
void gnss_h_21(double *state, double *sat_pos_vel, double *out_8983076725070806951);
void gnss_H_21(double *state, double *sat_pos_vel, double *out_6658591309417436314);
void gnss_predict(double *in_x, double *in_P, double *in_Q, double dt);
}