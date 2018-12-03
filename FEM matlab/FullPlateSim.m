
clc; close all; clear all
i = 9
    nc = i*2;
    no = nc;

[res,E_unitCell,nu_unitCell,sigma_max_Q,n_elems_Q,eta] = FEM_master(0.2,'LSTiso',nc,no,true);

sigma_max_Q


sigma_max_Q_1010 = 6.3436e+03

  