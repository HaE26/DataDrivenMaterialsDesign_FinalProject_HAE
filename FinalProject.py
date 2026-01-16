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


Progress:
    Almost all steps done. Just need to calculate the propertoes from load rate
    versus max stress, and then something about that Weibull distribution. 
    The code should then get cleaned up.
    
NOTES for Presentaiton:
    - Make a block diagram of processing steps
    - Share what what worked, what was hard, how you worked around it
    - Optimize code to be faster and/or easier to read/fewer lines.
"""
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 13:34:52 2025

@author: Hunter
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize



nu = 0.25

def graphing(
    x, y, title, xaxislabel, yaxislabel,
    legendlabel=None,
    font_size=16,
    color='blue',
    line_style='-',
    plot_type='plot',
    marker_type=None):
    """
    I don't know if this actually makes things easier, but I put most of the graphing into a function. I don't know.
    """
    if plot_type == 'plot':
        plt.plot(
            x, y,
            color=color,
            linestyle=line_style,
            marker=marker_type,
            label=legendlabel)
        
    elif plot_type == 'scatter':
        plt.scatter(
            x, y,
            color=color,
            marker=marker_type,
            label=legendlabel
        )
    else:
        raise ValueError(f"Unsupported plot_type for this funny function: {plot_type}")

    plt.title(title, fontsize=font_size)
    plt.xlabel(xaxislabel, fontsize=font_size)
    plt.ylabel(yaxislabel, fontsize=font_size)

def retrieve_sample_data(directory):
    """
    Given a parent folder, this function returns the characteristics of the sample test
    given in LTCC_033_sandgestrahlt, and the load data from the respective dat file,
    for each sample in a single loading rate. 
    
    Load is in Newtons(N) and lengths are in millimeters (mm).
    
    I ran into issues with the original xlsx file. I could not use it as xlsx 
    because I got an encoding error. I then had to change all of the commas to 
    decimals so that np.loadtxt would stop interpreting them as strings.
    This also means that the units of h have to be divided by 100
    
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

def get_linear_func(x,m,b):
    return m*x + b

def calculate_sigma(sample_data_dict):
    """
    Takes in the sample data dictionary from retrieve_sample_data, and then uses
    the max load, side length, thickness, and ball diameter to calculate the stresses and
    ultimate stress before failure. It then edits the loading data to add a 4th column containing stress,
    and adds a key value pair in the dictionary for max stress before failure. and then returns
    the entire dictionary with the added details.
    
    Because load parameter is in N and lengths are in mm, stress is in units of MPa.
    
    THE UNITS OF H ARE BEING CORRECTED HERE!!
    
    """
    #If stresses seem off, check first that B3B is the ball diameter and not the radius. Currently assuming it is the diameter
    for i in range(len(sample_data_dict)):
        #Extract values from dictionary
        max_load = np.max(sample_data_dict[f"Sample{i}"]["Loading"])
        load_data = sample_data_dict[f"Sample{i}"]["Loading"]
        L = sample_data_dict[f"Sample{i}"]["Side Length"]
        h = sample_data_dict[f"Sample{i}"]["Thickness"]/1000
        D_ball = sample_data_dict[f"Sample{i}"]["D_ball"]
        r_supp = D_ball/2
        
        #Call previously defined functions
        r_eff = get_r_eff(L,h,r_supp)
        f_factor = get_f(h,nu,r_eff,r_supp)
        
        #Calculate sigma_max and store for this sample
        sigma = f_factor * (load_data[:, 2]/(h)**2)
        load_data = np.column_stack((load_data, sigma))
        sample_data_dict[f"Sample{i}"]["Loading"] = load_data
        
        sigma_max = np.max(sigma)
        sample_data_dict[f"Sample{i}"]["Max Stress"] = sigma_max
        
    return sample_data_dict


def graph_load_measurements(sample_data_dict):
    """
    Uses the calculate_sigma function to graph stress over time for each sample in
    sample_data_dict. Also performs a linear fit of the loading data for stress over time
    to find loading rate in MPa/s, and then returns the slope in the dictionary
    under key "Loading Rate".
    
    CHANGES TO MAKE:
        - Improve appearance of graphs with fonts and titles
        - Add legend contianing max stress
    """
    sample_data_dict = calculate_sigma(sample_data_dict)
    
    for i in range(len(sample_data_dict)):
        #Extracting data to be plotted
        sample = sample_data_dict[f"Sample{i}"]
        load_data = sample["Loading"]
        Zeit = load_data[:, 0]
        Spannung = load_data[:, 3]
        
        #Splitting the data for a linear fit. We are taking the top 75% of data while slicing off the last value.
        Zeit_75 = Zeit[int(len(Zeit)*0.25):int(len(Zeit)-1)]
        Spannung_75 = Spannung[int(len(Spannung)*0.25):int(len(Spannung)-1)]
        
        #Calculating a linear fit.
        param, cov = optimize.curve_fit(get_linear_func,
                                        Zeit_75,
                                        Spannung_75, 
                                        p0=[Zeit_75[-1], Spannung_75[0]], 
                                        maxfev=50000)
        linear_fit = get_linear_func(Zeit_75, *param)
        """
        #Graphing Code
        plt.figure(figsize=(7, 5)) 
        graphing(Zeit, Spannung, title='Load over Time', xaxislabel='Time t (s)',
                 yaxislabel='Stress sigma (MPa)', legendlabel='Stress')
        plt.plot(Zeit_75,linear_fit,
                color='black',
                label = 'Linear Function fit')
        
        plt.legend(prop={'size': 15})
        """
        
        #Storing the slope of the fitted line as the loading rate.
        m, b = param
        sample_data_dict[f"Sample{i}"]["Loading Rate"] = m
        
    return sample_data_dict
        
def graph_sigma_max_vs_loading_rate(*sample__data_dicts):
    """
    This function is kind of a mess. I was just trying to get it to work.
    
    This function takes each dictionary for each loading rate, and then 
    graphs the max stress against loading rate for each sample. It then makes a linear fit.
    
    NEXT STEP:
        Need to calculate N and that other coefficient. If forgot what
        it's called because I don't have the paper with me rn.
    """
    
    
    plt.figure(figsize=(7, 5))
    loading_rate_list=[]
    max_stress_list =[]
    for sample_data_dict in sample__data_dicts:
        for i in range(len(sample_data_dict)):
            sample = sample_data_dict[f"Sample{i}"]
            
            loading_rate = sample["Loading Rate"]
            loading_rate_list.append(float(loading_rate))
            
            max_stress = sample["Max Stress"]
            max_stress_list.append(float(max_stress))
            
            
    loading_rate_array = np.asarray(loading_rate_list, dtype=float)
    max_stress_array   = np.asarray(max_stress_list,   dtype=float)

    log_loading_rate_array = np.log10(loading_rate_array)
    log_max_stress_array = np.log10(max_stress_array)
    plt.scatter(log_loading_rate_array, log_max_stress_array)
    param, cov = optimize.curve_fit(get_linear_func,
                                    log_loading_rate_array,
                                    log_max_stress_array, 
                                    p0=[loading_rate_array[-1], max_stress_array[0]], 
                                    maxfev=50000)

    linear_fit = get_linear_func(log_loading_rate_array, *param)
    m, b = param
    n = 1/m - 1 
    log10_B = b * (n + 1) - np.log10(n + 1)
    
    plt.plot(log_loading_rate_array, linear_fit)
    plt.xlabel("Log (dσ/dt) MPa/s",
               fontsize=17)
    plt.ylabel("Log (σ_max) (MPa)",
               fontsize=17)
    plt.title("Ultimate Stress vs Loading Rate",
               fontsize=17)
    plt.text(1.325, 2.38, f'The n quantity is {np.round(n,3)} \nThe lg(B*sigma_c^(n-2) quantity is {np.round(log10_B,3)}', 
             fontsize=17, color='black', 
             bbox=dict(facecolor='white', alpha=0.5)) # Alpha is the transparency
    plt.show()


    #print(f'the n value is {n}.\n The Bsigma^(n-2) quantity is: {y}')
    print(n)
    print(log10_B)
    
#Function calls
sample_data_dict_schnell = retrieve_sample_data("CuratedRawData_Sand_033_Box1_RT_schnell(1)/")
sample_data_dict_schnell = graph_load_measurements(sample_data_dict_schnell)

sample_data_dict_langsam = retrieve_sample_data("CuratedRawData_Sand_033_Box2_RT_langsam(1)/")
sample_data_dict_langsam = graph_load_measurements(sample_data_dict_langsam)

graph_sigma_max_vs_loading_rate(sample_data_dict_schnell,sample_data_dict_langsam)
