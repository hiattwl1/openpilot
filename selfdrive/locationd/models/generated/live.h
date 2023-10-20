#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_6199050048544017636);
void live_err_fun(double *nom_x, double *delta_x, double *out_7089086926548306499);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_3683298852657733819);
void live_H_mod_fun(double *state, double *out_4137737081612590829);
void live_f_fun(double *state, double dt, double *out_5505404667903329923);
void live_F_fun(double *state, double dt, double *out_3042775597903365000);
void live_h_4(double *state, double *unused, double *out_5960744769219280775);
void live_H_4(double *state, double *unused, double *out_4751454028888311057);
void live_h_9(double *state, double *unused, double *out_9134496086248111623);
void live_H_9(double *state, double *unused, double *out_4992643675517901702);
void live_h_10(double *state, double *unused, double *out_6171174241093153049);
void live_H_10(double *state, double *unused, double *out_5524634505072389290);
void live_h_12(double *state, double *unused, double *out_7141019455285610596);
void live_H_12(double *state, double *unused, double *out_8675833636789278764);
void live_h_35(double *state, double *unused, double *out_7558707577539954684);
void live_H_35(double *state, double *unused, double *out_5930270604464265055);
void live_h_32(double *state, double *unused, double *out_4655941239594316770);
void live_H_32(double *state, double *unused, double *out_3225246091562515118);
void live_h_13(double *state, double *unused, double *out_3156527716099568111);
void live_H_13(double *state, double *unused, double *out_8448013961875848770);
void live_h_14(double *state, double *unused, double *out_9134496086248111623);
void live_H_14(double *state, double *unused, double *out_4992643675517901702);
void live_h_33(double *state, double *unused, double *out_6640895078276883960);
void live_H_33(double *state, double *unused, double *out_2779713599825407451);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}