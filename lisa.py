'''
Lisa is a code to model the performance of solar cells using an extended Hovel 
model.

Copyright 2016-2018 Piotr Kowalczewski

Reference technical paper:
Piotr Kowalczewski, Lisa Redorici, Angelo Bozzola, Lucio Claudio Andreani,
"Silicon solar cells reaching the efficiency limits: from simple to complex 
modelling," Journal of Optics 18, 054001 (2016)

If you use this code for your research, please cite the paper above.

License: GNU General Public License v3.0 (see: LICENSE file)

This file is part of Lisa.

Lisa is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

NOTATION:
 - index "B" refers to the base, index "E" refers to the emitter

'''

import numpy as np
import scipy.optimize as opt
import sys
#from parameters import *
import parameters as params

if __name__ == "__main__":
    """ Usage: python p.py [thickness in um] [SRV in cm/s] """

    # ======================== PHYSICAL CONSTANTS ========================
    # Planck's constant, (eV*s)
    h = 4.135667516 * 1E-15      
    # speed of light (m/s)
    c = 299792458
    # electron charge (C)
    q = 1.602 * 1E-19          
    # Boltzmann constant (eV/K)
    kB = 8.617332478*1E-5
    # Boltzmann constant (J/K)
    kB_J = 1.380648813*1E-23

    if ( len(sys.argv)>1 ):
        th = float(sys.argv[1])*1E-4
        S_B = float(sys.argv[2])

    # base thickness
    th_base = th - params.th_emitter
    
    # set dopings
    # base
    N_B = max(params.N_a_B,params.N_d_B)
    N_E = max(params.N_a_E,params.N_d_E)

    # Band-gap narrowing is defined according to Schenk.
    # See: Schenk, J. Appl. Phys. 84, 3684-95,1998.
    # Returns delta_Eg in eV.
    # We use our own parametrization prepared for params.ni_0 = 9.65E9.
    # Can be checked with 
    # http://www.pvlighthouse.com.au/calculators/Band%20gap%20calculator/Band%20gap%20calculator.aspx
    def bgn(N_dop):
        A1 = -1.63096E-5;
        A2 = 0.10876;
        x0 = 18.82367;
        dx = 0.88411;
        
        delta_Eg = A2 + (A1-A2)/(1.0 + np.exp( (np.log10(N_dop)-x0)/dx));
        
        return delta_Eg

    ni_B = np.sqrt(params.ni_0**2 *np.exp(bgn(N_B)/(kB*params.T)) )
    ni_E = np.sqrt(params.ni_0**2 *np.exp(bgn(N_E)/(kB*params.T)) )

    # Built-in potential of the junction.
    # Eq. (6.2) from Nelson, The Physics of Solar Cells 
    # (London: Imperial College Press), 2003
    Vbi=kB_J*params.T/q * np.log(N_B*N_E/(params.ni_0**2))

    # Width of SCR in n-type and p-type
    # Eq. (6.10) and (6.11) from Nelson
    W_E_scr=(1/N_E)*np.sqrt(2*Vbi*12*8.85E-12*1E-2/((1/N_E+1/N_B)*q))
    W_B_scr=(1/N_B)*np.sqrt(2*Vbi*12*8.85E-12*1E-2/((1/N_E+1/N_B)*q)) 
    # total width of the SCR
    w_scr=W_E_scr+W_B_scr

    # width of the quasi-neutral n-type region (cm)
    w_n = params.th_emitter-W_E_scr

    if (w_n<0):
        print('\nWarning: thickness of the quasi-neutral n-type region \
              is negative! Change the doping.\n')

    # ======================== OUTPUT ========================
    if (params.params_flag):
        print('\n============ Parameters of the structure ============\n')
        print('Built-in potential of the junction (V): \t\t {:.3f}'
              .format(Vbi))
        print('Cell thickness (um): \t\t\t\t\t {:.3f}'.format(th*1E4))
        print('Emitter thickness (um): \t\t\t\t {:.3f}'
              .format(params.th_emitter*1E4))
        print('SCR thickness (um): \t\t\t\t\t {:.3f}'.format(w_scr*1E4))
        print('Thickness of the quasi-neutral n-type region (um): \t {:.3f}'
              .format(w_n*1E4))

    en_vec = np.linspace(params.en_start,params.en_stop,params.en_points)

    # ======================== EXTERNAL DATA FILES ========================
    # Optical functions of silicon. Default optical functions taken from 
    # Green, Sol. Energy Mat. Sol. Cells 92 1305-10, 2008.
    # The data file has the following format:
    # Energy (eV) Wavelength (nm) eps1 eps2 n k Absorption coefficient (cm-1)
    en_vec_TEMP = []
    n_vec_TEMP = []
    abs_coeff_vec_TEMP = []                        
    for line in open("./Data/"+params.si_f):
        line = line.split()
        en_vec_TEMP.append(float(line[0]))
        n_vec_TEMP.append(float(line[4]))
        abs_coeff_vec_TEMP.append(float(line[6]))
    
    en_vec_TEMP.reverse()
    n_vec_TEMP.reverse()
    abs_coeff_vec_TEMP.reverse()

    n_vec = np.interp(en_vec, en_vec_TEMP, n_vec_TEMP)
    abs_coeff_vec = np.interp(en_vec, en_vec_TEMP, abs_coeff_vec_TEMP)

    # Solar spectrum. Units:
    # Energy (eV)  AM1.5G (Wm-2eV-1)
    en_vec_TEMP = []
    am15g_vec_TEMP = []
    for line in open("./Data/"+params.solar_f):
        line = line.split()
        en_vec_TEMP.append(float(line[0]))
        am15g_vec_TEMP.append(float(line[1]))
        
    en_vec_TEMP.reverse()
    am15g_vec_TEMP.reverse()
    am15g_vec = np.interp(en_vec, en_vec_TEMP, am15g_vec_TEMP)
    
    # Integrated solar spectrum/energy -- for spect_fact
    # (such a simple notation works).
    solar_spect_int = np.trapz(am15g_vec/en_vec,en_vec)

    # For a given thickness, prepare alpha_LT -- effective absorption 
    # coefficient (including light trapping). Approximation from Green 
    # (instead of integrating), see: Green, Prog. Photovolt.,
    # Res. Appl. 10 235-41 (2002).
    def alpha_LT_func(i):
        a = 0.935
        b = 0.67        

        x = a * (abs_coeff_vec[i]*th)**b
        eff_LPE = (2.0+x)/(1.0+x)
        return abs_coeff_vec[i] * eff_LPE

    # Diffusion length calculated for Auger recombination
    # Model taken from: Richter A et al., Phys. Rev. B 86, 165202 (2012)
    # Notation:
    # D -- diffusion coefficient (cm2/s)
    # v -- voltage (V)
    # Nd -- donor doping (cm-3)
    # Na -- acceptor doping (cm-3)
    # Approximation: n and p are not calculated from d-d equations; 
    # they are taken from approximate formulas. We have checked it and it is 
    # a pretty good approximation (the self-consisten approach is not 
    # necessary in this case).
    def L_Auger(D,v,Nd,Na):
        
        # Approximation for voltage close to 0 (check if this can be done 
        # better). If voltage is lower than 0.01, there is "divide by zero 
        # encountered in double_scalars" exception thrown.
        # dn term from illumination should be added at low voltage
        if (v<0.01):
            v=0.01
        
        N0_eeh = 3.3E+17
        N0_ehh = 7.0E+17
        
        # p-type material
        if (Na>Nd):
            ni = np.sqrt(params.ni_0**2 *np.exp(bgn(Na)/(kB*params.T)) )
            p0 = 0.5*(Na - Nd + np.sqrt((Na - Nd)*(Na - Nd) + 4 * ni**2))
            n0 = ni**2 / p0
        
        # n-type material
        else:
            ni = np.sqrt(params.ni_0**2 *np.exp(bgn(Nd)/(kB*params.T)) )
            n0 = 0.5*(Nd - Na + np.sqrt((Nd - Na)*(Nd - Na) + 4 * ni**2));
            p0 = ni**2 / n0
            
        g_eeh = 1.0 + 13.0 * (1.0 - np.tanh(pow((n0 / N0_eeh),0.66)))
        g_ehh = 1.0 + 7.5 * (1.0 - np.tanh(pow((p0 / N0_ehh),0.63)))
        
        
        # We calcualte dn from Eq. (2) from Richter et al., 
        # IEEE J. Photovolt. 3 1184-91 (2013).
        eq_a = 1.0
        eq_b = n0+p0
        eq_c = n0*p0-ni**2 *np.exp(q*v/(kB_J*params.T))
        eq_delta = eq_b**2 -4*eq_a*eq_c
        dn = (-eq_b+np.sqrt(eq_delta))/(2*eq_a)
        
        # Auger recombination rate (1/cm3s), see Eq. (2) from Richter2013 
        # and Eq. (18) from Richter2012 (full citations above).
        R_Aug = ni**2 * (np.exp(q*v/(kB_J*params.T)) - 1 ) * ( 2.5e-31 * g_eeh * n0 + 8.5E-32 * g_ehh * p0 + 3E-29 * dn**0.92)
        
        tau = dn/R_Aug
        
        # Diffusion length (cm).
        L = np.sqrt(D*tau)    

        return L

    # Concentration of holes as a f of voltage.
    def p(v):
        ni = np.sqrt(params.ni_0**2 *np.exp(bgn(N_B)/(kB*params.T)) )
        p0 = 0.5*(N_B + np.sqrt((N_B)*(N_B) + 4 * ni**2))
        n0 = ni**2 / p0
        
        # We need to calcualte dn from Eq. (2) Richter 2013 
        # (full citation above).
        eq_a = 1.0
        eq_b = n0+p0
        eq_c = n0*p0-ni**2 *np.exp(q*v/(kB_J*params.T))
        eq_delta = eq_b**2 -4*eq_a*eq_c
        dn = (-eq_b+np.sqrt(eq_delta))/(2*eq_a)
        
        p = dn+p0
        
        return p

    def n(v):
        ni = np.sqrt(params.ni_0**2 *np.exp(bgn(N_E)/(kB*params.T)) )
        n0 = 0.5*(N_E + np.sqrt((N_E)*(N_E) + 4 * ni**2))
        p0 = ni**2 / n0
        
        # We need to calcualte dn from Eq. (2) Richter 2013 
        # (full citation above).
        eq_a = 1.0
        eq_b = n0+p0
        eq_c = n0*p0-ni**2 *np.exp(q*v/(kB_J*params.T))
        eq_delta = eq_b**2 -4*eq_a*eq_c
        dn = (-eq_b+np.sqrt(eq_delta))/(2*eq_a)
        
        n = dn+n0
        
        return n

    # Returns total current as a function of voltage.
    def J(v):
        
        # Changes of SCR with voltage. This change is important below cell 
        # thickness of around 10 um.
        dV = Vbi-v
        if (dV<0):
            dV=0
        
        W_E_scr=(1/N_E)*np.sqrt(2*dV*12*8.85E-12*1E-2/((1/N_E+1/N_B)*q)) 
        W_B_scr=(1/N_B)*np.sqrt(2*dV*12*8.85E-12*1E-2/((1/N_E+1/N_B)*q))

        # Base.
        L_B = L_Auger(params.D_B,v,0,N_B)
        t_B = L_B**2 / params.D_B
        
        # ==========================================================
        # TEST
        # Inlucde SRH recombination
        # L_B_Aug = L_B
        # t_B_Aug = t_B
        
        # diffusion length in cm
        # L_B_SRH = 200E-4
        # t_B_SRH = L_B_SRH**2 / params.D_B
        
        # t_B = t_B_SRH * t_B_Aug / (t_B_SRH + t_B_Aug)
        # L_B = np.sqrt(t_B*params.D_B)
        
        # ==========================================================
        
        # Emitter.
        L_E = L_Auger(params.D_E,v,N_E,0)
        t_E = L_E**2 / params.D_E

        # Currents as a function of energy.
        J_B_vec = []
        J_E_vec = []
        J_scr_vec = []

        for i in range(len(en_vec)):
            
            alpha_LT = alpha_LT_func(i)

            # Reflection at the back interface.
            Rb = 1
            
            z_B = (th_base - W_B_scr) / L_B
            
            # Base -- current at x = W_B_scr.
            # Spectral factor to make A_B integrable. Integrating this quantity 
            # with respect to energy gives 1.
            spect_fact = am15g_vec[i]/en_vec[i] / solar_spect_int
            
            gamma_B = L_B **2 / (params.D_B * (1-alpha_LT**2*L_B**2) )
            
            F1 = Rb * np.exp(-2*alpha_LT*th)
            
            # 1/q: W --> eV
            # 1E-4: 1/m2 --> 1/cm2 
            a_minus = alpha_LT * am15g_vec[i]/en_vec[i] / (1 - F1 * (1-1/n_vec[i]**2) ) * 1/q * 1E-4 
            a_plus = a_minus * Rb

            A_B = ni_B**2/p(v) * (np.exp(q*v/(kB_J*params.T)) - 1) * spect_fact - gamma_B * ( a_minus * np.exp(-alpha_LT*(W_B_scr + params.th_emitter) )  + a_plus * np.exp(-alpha_LT*(2*th - (W_B_scr + params.th_emitter) ))   ) 
            
            B_B = ( - A_B * ( params.D_B/L_B * np.sinh(z_B) + S_B * np.cosh(z_B) ) - gamma_B * ( a_minus * np.exp(-alpha_LT*th ) * (S_B - alpha_LT * params.D_B) + a_plus *  np.exp(-alpha_LT*(2*th- (th_base + params.th_emitter) )) * (S_B + alpha_LT * params.D_B) ) ) / (params.D_B/L_B*np.cosh(z_B) + S_B * np.sinh(z_B))
            
            J_B = q*params.D_B * ( B_B/L_B + gamma_B * ( -alpha_LT * a_minus * np.exp(-alpha_LT*(W_B_scr + params.th_emitter) )  + alpha_LT * a_plus * np.exp(-alpha_LT*(2*th- (W_B_scr + params.th_emitter) ))   )  ) 
        
            # Units '1E3': A --> mA.
            J_B = J_B * 1E3
        
            J_B_vec.append(J_B)
        
            # Emitter -- current at x = W_E_scr.
            gamma_E = L_E**2 / (params.D_E * ( 1-alpha_LT**2*L_E**2 ) )
            
            z_E = (params.th_emitter - W_E_scr) / L_E
            
            A_E = ni_E**2/n(v) * (np.exp(q*v/(kB_J*params.T)) - 1) * spect_fact - gamma_E*( a_minus * np.exp(-alpha_LT*(params.th_emitter - W_E_scr) )  + a_plus * np.exp(-alpha_LT*(2*th- (params.th_emitter-W_E_scr) ))   ) 
            
            B_E = (-A_E * ( params.D_E/L_E * np.sinh(z_E)  + params.S_E *np.cosh(z_E) )  -gamma_E * (a_minus*(params.S_E+alpha_LT*params.D_E) + a_plus*np.exp(-2*alpha_LT*th)*(params.S_E-alpha_LT*params.D_E)   ) )  / (params.D_E/L_E*np.cosh(z_E) + params.S_E*np.sinh(z_E) )
            
            J_E = q*params.D_E* ( -B_E/L_E + gamma_E * ( -alpha_LT * a_minus * np.exp(-alpha_LT*(-W_E_scr + params.th_emitter) )  + alpha_LT * a_plus * np.exp(-alpha_LT*(2*th - (-W_E_scr + params.th_emitter) ))   ) )
            
            # units '1E3': A --> mA
            J_E = - J_E * 1E3
        
            J_E_vec.append(J_E)
            
            # SCR.
            J_gen = q/alpha_LT * ( a_minus * ( np.exp(-alpha_LT*(-W_E_scr + params.th_emitter) ) - np.exp(-alpha_LT*(W_B_scr + params.th_emitter) )  ) + a_plus * ( np.exp(-alpha_LT*(2*th - (W_B_scr + params.th_emitter) )) - np.exp(-alpha_LT*(2*th - (-W_E_scr + params.th_emitter) )) )  )
            
            J_rec = q*params.ni_0*(W_E_scr + W_B_scr) / (t_E+t_B) * 2*np.sinh(q*v/(2*kB_J*params.T)) * (np.pi/2)    /   (  q*(Vbi-v)/ (kB_J*params.T)  )
            
            J_scr = J_gen - J_rec
            
            # units '1E3': A --> mA
            J_scr = J_scr * 1E3
            
            J_scr_vec.append(J_scr)
        
        # Integrate with respect to energy to get the total current.
        # For V = 0 it should be close to Jsc.
        
        J_B_tot = np.trapz(J_B_vec,en_vec)
        J_E_tot = np.trapz(J_E_vec,en_vec)
        J_scr_tot = np.trapz(J_scr_vec,en_vec)
        
        return J_B_tot + J_scr_tot + J_E_tot
        
    if (params.jv_flag):

        v_vec = np.concatenate((np.linspace(0.0,0.59,10),
                                np.linspace(0.6,0.85,250)))
        
        j_vec = []

        for v in v_vec:
            J_V = J(v)
            j_vec.append(J_V)
            
        # Calculating Voc.
        j_buff = j_vec[::-1]
        v_buff = v_vec[::-1]
        Voc = np.interp(0,j_buff,v_buff)

        # Truncate jv vectors above Voc.
        i_max = np.where(v_vec>Voc)[0][0]
        v_vec = v_vec[:i_max+1]
        j_vec = j_vec[:i_max+1]

        power = [v * j for v, j in zip(v_vec, j_vec)]
        p_max = max(power)
        i_mpp = np.where(power==p_max)[0][0]
        v_m = v_vec[i_mpp]

        Jsc = j_vec[0]

        # Save JV to a file.
        JV_f = open("./Results/JV.dat","w")
        for i in range(len(v_vec)):
            print('{:.3f} {:.3f}'.format(v_vec[i], j_vec[i]), file = JV_f)

        print("\nJV characteristic saved to JV.dat file.")
        JV_f.close()
        
        eff = (max(j_vec*v_vec))
    
    else:
        
        # Optimization: looking for maximum (minimum).
        # https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.optimize.minimize_scalar.html#scipy.optimize.minimize_scalar
        optimization = opt.minimize_scalar(lambda v: -J(v)*v)
        eff = -optimization['fun']
        # The spectrum is normalized, so eff=p_max.
        p_max = -optimization['fun']
        Vm = optimization['x']
        Jsc = J(0)
        # Find Voc which is root.
        Voc = opt.brenth(lambda v: J(v), Vm, 1)
    
    if (params.results_flag):
        print('\n============ Results ============\n')
        print("Thickness (um): {0:.3f}".format(th*1E4)) 
        print("Bottom SRV (cm/s): {0:.3f}".format(S_B))
        print("Open-circuit voltage (V): {0:.3f}".format(Voc))
        print("Short-circuit current density (mA/cm2): {0:.3f}".format(Jsc))
        print("Fill Factor: {0:.3f}".format((p_max/(Jsc*Voc))))
        print("Efficiency (%): {0:.3f}".format(eff))
    else:
        print("{:.3f} \t {:.3f} \t {:.3f} \t {:.3f} \t {:.3f} \t {:.3f} \t"
              .format( (th*1E4), S_B, Voc, Jsc, p_max/(Jsc*Voc), eff) )
        
