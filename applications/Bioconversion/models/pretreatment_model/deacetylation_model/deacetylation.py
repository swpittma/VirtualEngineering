# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 14:33:24 2025

@author: spittman
"""

# Set up Excel Relation and Data Frame
# knowns should be in excel file at this point
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
import math
import matplotlib.pyplot as plt

###########################################
# PHYSICAL AND CHEMICAL CONSTANTS
###########################################
###
# Equation Constants
F = 96485 # Faraday Constant -Amps*Sec/mol- -C/mol-
R = 0.008314 # Ideal gas constant -8.314 J/mol*K- -kJ/mol-
R_gas = 0.08206 # -L*atm/mol*K-
rho_h2o = 1000 # density of water -g/L-
mrho_h2o = 55.56 # density of water -mol/m3-
Ea = 38 # Arrhenius activation energy for acetyl groups -J/mol-
Arc = 1.1e4 # Arrhenius pre-exponent factor  for acetyl groups -L/(mol-s)-
Ea_L = 55 # Arrhenius activation energy for lignin -J/mol-
Arc_L = 1.5e6 # Arrhenius pre-exponent factor for lignin -L/(mol-s)-
Ea_X = 45 # Arrhenius activation energy for xylan -J/mol-
Arc_X = 7e3 # Arrhenius pre-exponent factor for xylan -L/(mol-s)-
###
# Stoichiometric Equilibrium Constants
pKa_ace = 4.76
pK_Naace = 4.76


###
# Typical Composition of Cob Types
# Corn Husk
Husk_C = 0.38 # fraction of cellulose -%/100-
Husk_H = 0.35 # fraction of hemicellulose -%/100- 
Husk_L = 0.15 # fraction of lignin -%/100-
Husk_Ep = 0.126 # effective porosity
Husk_T = 5.5 # tortuosity
Husk_v = 0.65 # void fraction
# Corn Cob
Cob_C = 0.33 # fraction of cellulose -%/100-
Cob_H = 0.36 # fraction of hemicellulose -%/100- 
Cob_L = 0.17 # fraction of lignin -%/100-
Cob_Ep = 0.126 # effective porosity
Cob_T = 10.3 # tortuosity
Cob_v = 0.65 # void fraction
# Corn Stalk
Stalk_C = 0.40 # fraction of cellulose -%/100-
Stalk_H = 0.28 # fraction of hemicellulose -%/100- 
Stalk_L = 0.21 # fraction of lignin -%/100-
Stalk_Ep = 0.671 # effective porosity
Stalk_T = 1.36 # tortuosity
Stalk_v = 0.65 # void fraction
# Acetyl Groups
Acy_P = 0.023 #0.023 #0.0223 # weight percentage of acetate in biomass
# Acy_P = 0.027 # switchgrass
# Lignin Groups
Lig_P = 0.153 # weight percentage of lignin in biomass
# Lig_P = 0.192 # switchgrass
# Xylan Groups
Xy_P = 0.213 # corn stover
# Xy_P = 0.218 # switchgrass

loops = (3,6)
corn = np.zeros(loops)
corn[0,0] = Husk_Ep
corn[0,1] = Husk_T
corn[0,2] = Husk_v
corn[0,3] = Husk_C
corn[0,4] = Husk_H
corn[0,5] = Husk_L
corn[1,0] = Cob_Ep
corn[1,1] = Cob_T
corn[1,2] = Cob_v
corn[1,3] = Cob_C
corn[1,4] = Cob_H
corn[1,5] = Cob_L
corn[2,0] = Stalk_Ep
corn[2,1] = Stalk_T
corn[2,2] = Stalk_v
corn[2,3] = Stalk_C
corn[2,4] = Stalk_H
corn[2,5] = Stalk_L

###
# Molar Mass - g/mol - 
M_NaOH = 40  
M_OH = 17.01
M_Na = 22.99
M_ace = 59.04 # acetate - C2H3O2
M_hace = 60.05 # acetic acid - C2H4O2
M_Naace = 82.03 # sodium acetate
M_Cel = 162.14 # cellulose - C6H10O5
M_Hem = 174.15 # hemicellulose - modeled as acetylated xylan
M_acy = 43.045 # acetyl group - C2H3O
M_X = 132.11 # deacetylated xylan
M_Lig = 188 # lignin - common approximation
M_H = 1

# Create matrix -g/mol-
loops = (3,12)
const = np.zeros(loops)
const[0,0] = M_ace
const[0,1] = M_hace
const[0,2] = M_Naace
const[0,3] = M_Cel
const[0,4] = M_Hem
const[0,5] = M_Lig
const[0,6] = M_acy
const[0,7] = M_X
const[0,8] = M_NaOH
const[0,9] = M_Na
const[0,10] = M_OH
const[0,11] = M_H


#########################################
# INITIAL CONDITIONS
#########################################
###
# VARIABLE SET INITIAL CONDITIONS 
t_final = 7200 # Total reaction time -Sec- -120 min- #7200 #5400
Tf = 365.15 # Operating Fluid Temperature -K- -363.15|90 C- -323.15|50 deg C-
V_DA = 55 # Volume of Deacetylation Reactor -L- #30 #20
pHi_DA = 7 # Inital pH
pOHi_DA = 14-pHi_DA # Initial pOH
Part = 0.75 # particle size -inch-
Decomp = 2 # storage time -unitless- -(months)-

###
# Initial Influent Content -g/L-
C_Husk = 10 # corn husk
C_Cob = 10 # corn cob
C_Stalk = 10 # corn stalk
C_NaOH = 9.09 # sodium hydoxide # 4.55 6.36 9.09 # 9.47 5.35 # 2.61 2.01 1.41
# Biomass
W_BM = 5000 # biomass -g- #600
# Acetyl Groups
C_acy = W_BM*Acy_P/V_DA # acetyl groups -g/L-2
Molar_acy = C_acy/M_acy
C_ace_max = C_acy*(M_ace/M_acy)
# Lignin Groups
C_lig = W_BM*Lig_P/V_DA # -g/L-
C_lig_max = C_lig
Molar_lig = C_lig/M_Lig
# Xylan Groups
C_xy = W_BM*Xy_P/V_DA # acetyl groups -g/L-2
C_xyo_max = C_xy*(M_X/M_Hem)
Molar_xy = C_xy/M_X

#######################################
# Set UP Matrices
#######################################

### =
# Calculation Matrix -row 0 g/L- -row 1 mol/L-
loops=(1,12)
DAC = np.zeros(loops)
M_DAC = np.zeros(loops)
dM_DAC = np.zeros(loops)
DAC[0,0] = 0                                               # Acetate
DAC[0,1] = 0                                               # Acetic Acid
DAC[0,2] = 0                                               # Sodium Acetate
DAC[0,3] = (C_Husk*Husk_C)+(C_Cob*Cob_C)+(C_Stalk*Stalk_C) # Cellulose
DAC[0,4] = C_xy # Hemicellulose
DAC[0,5] = C_lig # Lignin
DAC[0,6] = C_acy                                           # Acetyl Groups
DAC[0,7] = 0                                               # Deacetylated Xylan
DAC[0,8] = C_NaOH                                          # Sodium Hydroxide
DAC[0,9] = 0                                               # Sodium
DAC[0,10] = 10**(-pOHi_DA)*M_OH                            # Hydroxide
DAC[0,11] = 10**(-pHi_DA)*M_H                              # Hydrogen
# Convert
M_DAC[0,:] = DAC[0,:]/const[0,:] # -mol/L-
# NaOH Dissociates Completely
M_DAC[0,9] = M_DAC[0,9] + (1*M_DAC[0,8]) # Sodium
M_DAC[0,10] = M_DAC[0,10] + (1*M_DAC[0,8]) # Hydroxide
M_DAC[0,8] =M_DAC[0,8] - M_DAC[0,8] # Sodium Hydroxide
# Convert
DAC[0,:] = M_DAC[0,:]*const[0,:] # -g/L-

# Free Acids
M_DAC[0,10] = M_DAC[0,10] - 0.02

# Arrehnius Rate Kinetics
# Lignin
kT_arh_L = Arc_L*np.exp(-Ea_L/(R*Tf)) # -L/(mol-s)-
# Acetate
kT_arh = Arc*np.exp(-Ea/(R*Tf)) # -L/(mol-s)-
# Xylan
kT_arh_X = Arc_X*np.exp(-Ea_X/(R*Tf)) # -L/(mol-s)-

###############################################
### DEFINE FUNCTION
def DACy(t, y):
    # Arrehnius Rate Kinetics
    # Lignin
    kT_arh_L = Arc_L*np.exp(-Ea_L/(R*Tf)) # -L/(mol-s)-
    # Acetate
    kT_arh = Arc*np.exp(-Ea/(R*Tf)) # -L/(mol-s)-
    # Xylan
    kT_arh_X = Arc_X*np.exp(-Ea_X/(R*Tf)) # -L/(mol-s)-
    # Define Rate Equation
    # Lignin
    dligdt = kT_arh_L*y[3]*y[1] # -mol/L-s-
    # Change value
    dClig_dt = -dligdt
    dCslig_dt = +dligdt
    # Acetate
    dacedt = kT_arh*y[0]*y[1] # -mol/L-s-
    # Change value
    dCacy_dt = -dacedt
    dCace_dt = +dacedt
    # Xylan
    dxydt = kT_arh_X*y[5]*y[1] # -mol/L-s-
    dCxyo_dt = +dxydt
    dCxy_dt = -dxydt
    dCOH_dt  =-(dacedt+(1.2*dligdt)+(dxydt))
    return [dCacy_dt, dCOH_dt, dCace_dt, dClig_dt, dCslig_dt, dCxy_dt, dCxyo_dt]

###############################################
### INTEGRATE
def deacetylate(ve_params,verbose=True,show_plots=True):
    print("doing deacetylation..")
    Tf=ve_params.pt_in['deacetylation temperature']
    Acy_P = ve_params.pt_in['acetylfrac']
    t_span = (0,t_final)
    t_span = (0,t_final)
    Cacetyl = M_DAC[0,6]
    COH = M_DAC[0,10]
    Cace = M_DAC[0,0]
    Cxy = M_DAC[0,4]
    Cxyo = 0
    Cxyoh = 0
    Clig = M_DAC[0,5]
    Cslig = 0
    y0 = [Cacetyl, COH, Cace, Clig, Cslig, Cxy, Cxyo]
    sol_DA = solve_ivp(DACy,t_span, y0, method='BDF', rtol=1e-6, atol=1e-9)

    ############################################
    # EXPERIMENTAL TIME POINTS
    ############################################
    times = np.array([0,5,10,15,20,25,30,40,50,60,70,80,90,100,110,120], dtype=int)

    # Acetate Experimental Time Points
    acetate_interp = np.interp(times, sol_DA.t/60, sol_DA.y[2]*M_ace)
    acetyl_interp = np.interp(times, sol_DA.t/60, sol_DA.y[0]*M_acy)
    DAC_exp_ace = np.vstack((times, acetate_interp, acetyl_interp)).T

    Macetate_interp = np.interp(times, sol_DA.t/60, sol_DA.y[2])
    Macetyl_interp = np.interp(times, sol_DA.t/60, sol_DA.y[0])

    # Lignin Experimental Time Points
    lignin_interp = np.interp(times, sol_DA.t/60, sol_DA.y[3]*M_Lig)
    Slignin_interp = np.interp(times, sol_DA.t/60, sol_DA.y[4]*M_Lig)
    DAC_exp_lig = np.vstack((times, Slignin_interp, lignin_interp)).T

    # Xylan Experimental Time Points
    Xylan_interp = np.interp(times, sol_DA.t/60, sol_DA.y[5]*M_Hem)
    Xylose_interp = np.interp(times, sol_DA.t/60, sol_DA.y[6]*M_X)
    OH_interp = np.interp(times, sol_DA.t/60, sol_DA.y[1]*M_NaOH)
    DAC_exp_xyo = np.vstack((times, Xylan_interp, Xylose_interp)).T

    # pH Experimental Time Points
    # Calculate pKw
    Hnot0 = 55.83
    Kw2 = 10**(-14)*((Hnot0/R)*((1/289.15)-(1/Tf)))
    pKw = -np.log10(Kw2)
    OH_interp = np.interp(times, sol_DA.t/60, sol_DA.y[1]*M_OH)
    pH = pKw+np.log10(OH_interp)
    DAC_exp_pH = np.vstack((times, pH)).T


    ############################################
    # CREATE GRAPHS
    ############################################
    # Acetate
    '''plt.plot(sol_DA.t/60, sol_DA.y[2]*M_ace,color=(1, 0.5, 0.8))
    plt.title('Deacetylation Time Profile')
    plt.xlabel('Time (min)')
    plt.ylabel('Species Concentration (g/L)')
    plt.show() 

    # Lignin
    plt.plot(sol_DA.t/60, sol_DA.y[4]*M_Lig,color=(0.1, 0.3, 0.7))
    plt.title('Lignin Solubilization Time Profile')
    plt.xlabel('Time (min)')
    plt.ylabel('Species Concentration (g/L)')
    plt.show() 

    # Xylose
    plt.plot(times, DAC_exp_xyo[:,2],color=(0, 0.7, 0.3))
    plt.title('Xylose Time Profile')
    plt.xlabel('Time (min)')
    plt.ylabel('Species Concentration (g/L)')
    plt.show() 

    # pH
    plt.plot(times, DAC_exp_pH[:,1],color=(1, 0.5, 0.5))
    plt.title('pH Time Profile')
    plt.xlabel('Time (min)')
    plt.ylabel('Species Concentration (g/L)')
    plt.show()'''
