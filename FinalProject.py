# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 12:20:31 2025

@author: Hunter

TODO Graph Ultimate strength sigma_f against loading rate d_sigma/d_t

How to get ultimate strength:
    - Get MAX LOAD before failure and convert to stress using sigma_max equation in paper
        - Calculate D_eff
        - Calculate f
        - Calculate sigma_max for max load

How to get Loading Rate:
    -Graph force against time
    -Make a linear fit for an appropriate portion of the curve
    -Calculate the derivative of the fitted line (=m of linear curve)
    
I think there is a Weibull distribution here somewhere?

Plot ultimate strength against loading rate, then fit a curve (try to?)    

"""
import glob
import numpy as np
import matplotlib.pyplot as plt



nu = 0.25

def retrieve_sample_data(directory):
    """
    Given a parent folder, this function returns the characteristics of the sample test
    given in LTCC_033_sandgestrahlt, and the load data from the respective dat file,
    for each sample in a single loading rate. 
    
    Load is in Newtons(N) and lengths are in millimeters (mm).
    
    I ran into issues with the original xlsx file. I could not use it as xlsx 
    because I got an encoding error. I then had to change all of the commas to 
    decimals so that np.loadtxt would stop interpreting them as strings.
    
    """
    DAT_files = sorted(glob.glob(f"{directory}LTCC_*.dat")) #Provide file name for for loop.
    props = np.loadtxt(f"{directory}LTCC_033_sandgestrahlt.csv",delimiter=";",skiprows=2,usecols=(1,2,3))
    sample_data = {}
    for i,file_name in enumerate(DAT_files):
        
        sample_data[f"Sample{i}"] = {}
        Loading = np.loadtxt(file_name,skiprows=2)
        sample_data[f"Sample{i}"]["Loading"] = Loading
        D_ball =  props[i,0]
        sample_data[f"Sample{i}"]["D_ball"] = D_ball
        h = props[i,1]
        sample_data[f"Sample{i}"]["Thickness"] = h
        L = props[i,2]
        sample_data[f"Sample{i}"]["Side Length"] = L
        
    return sample_data

def get_f(h,nu,r,r_supp):
    """
    Function to get prefactor f from sample thickness, poisson's ratio, effective radius,
    and support radius. This is dimensionless.
    """
    m_1 = 0.697
    m_2 = -0.118
    m_3 = -0.728
    term1 = m_1*(1+nu)
    #print(f"{term1} is term1")
    term2 = m_2*np.log(h/r_supp)
    #print(f"{term2} is term2")
    term3 = m_3*((r*h**2)/r_supp**3)**(1/4)
    #print(f"{term3} is term3")
    inner_f = term1 + term2 + term3
    factor = np.exp(inner_f)
    return factor

def get_r_eff(L,h,r_supp):
    """
    Fucntion to get effective radius for a square sample using side length,
    thickness, and support radius. This is in units of mm.
    """
    D_eff = L*(1.053 - 0.017*(h*L/(r_supp)**2))
    r_eff = D_eff/2
    
    return r_eff

def graph_measurements(DAT_file_np):
    """Currentlty unused function to be repurposed for line fitting.
    """
    Zeit = DAT_file_np[0]
    Kraft = DAT_file_np[2]
    
    plt.figure(figsize=(7, 5)) 
    plt.plot(Zeit,Kraft)
    
    plt.show()
    return DAT_file_np

def calculate_sigma_max(sample_data_dict):
    """
    Takes in the sample data dictionary from retrieve_sample_data, and then uses
    the max load, side length, thickness, and ball diameter to calculate the
    ultimate stress before failure. It then returns all stresses for this 
    loading rate as a list ready to graph.
    
    Because load parameter is in N and lengths are in mm, stress is in units of MPa.
    
    """
    #If stresses seem off, check first that B3B is the ball diameter and not the radius. Currently assuming it is the diameter
    sigma_max_list=[]
    for i in range(len(sample_data_dict)):
        #Extract values from dictionary
        max_load = np.max(sample_data_dict[f"Sample{i}"]["Loading"])
        L = sample_data_dict[f"Sample{i}"]["Side Length"]
        h = sample_data_dict[f"Sample{i}"]["Thickness"]/1000
        D_ball = sample_data_dict[f"Sample{i}"]["D_ball"]
        r_supp = D_ball/2
        
        #Call previously defined functions
        r_eff = get_r_eff(L,h,r_supp)
        f_factor = get_f(h,nu,r_eff,r_supp)
        
        #Calculate sigma_max and store for this sample
        sigma_max = f_factor * (max_load/(h)**2)
        sigma_max_list.append(sigma_max)
        
    return sigma_max_list

#Function calls
sample_data_dict = retrieve_sample_data("CuratedRawData_Sand_033_Box1_RT_schnell(1)/")
sigma_max_list = calculate_sigma_max(sample_data_dict)

