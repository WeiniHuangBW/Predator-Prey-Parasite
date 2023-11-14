##### THIS SCRIPT DOES STOCHASTIC SIMULATIONS USING THE GILLESPIE ALGORITHM #####
##### THIS SCRIPT RECORDS THE COEXISTENCE OF PREDATOR, PREY AND PARASITES IN A PREDATOR-PREY-PARASITE SYSTEM WHERE THE PARASITE IS TRANSMITTED TROPHICALLY FROM INFECTED PREY TO PREDATORS #####
##### THIS SCRIPT RECORDS THE FREQUENCIES OF INFECTED INDIVIDUALS IN A PREDATOR-PREY-PARASITE SYSTEM WHERE THE PARASITE IS TRANSMITTED TROPHICALLY FROM INFECTED PREY TO PREDATORS #####

## Create a .csv file called "output.csv" (13 columns)
## Assign the following names to the columns: Qx, Qy, rx rp, coexistence, coexistence_predator_and_prey, coexistence_prey, extinction, cyan, yellow, black, magenta, white
## Import relevant packages and modules ##
import pandas as pd, math, statistics, random, numpy.random as nura, numpy as np, array as arr, matplotlib.pyplot as plt, matplotlib.patches as mpatches, sys, getopt, time
from csv import writer

# initial population sizes
uninfected_prey = 800 # uninfected prey
parasite = 1000 # free-living parasites
uninfected_predator = 100 # uninfected predator
infected_prey = 0 # at the beginning of the simulation, all predator are uninfected
infected_predator = 0 # at the beginning of the simulation, all predator are uninfected
sum_prey = uninfected_prey + infected_prey # total number of prey individuals
sum_predator = uninfected_predator + infected_predator # total number of predator individuals
sum_parasite = parasite + infected_prey + infected_predator # total number of parasite individuals (this number is to keep track of the extinction of the parasite)

# prey parameters 
gx = 2.0 # growth rate of the prey
dx = 0.1 # intrinsic death rate of the prey
pop_limit = 2000 # population limit (competition between prey for carrying capacity)

# parasite parameters
n_z = 6 # number of parasite offspring per reproduction
dz = 0.09 # intrinsic death rate of the free-living parasite
S = 0.0005 # scaling factor for prey-parasite population sizes

#predator parameters
fy = 0.01 # predation rate
ky = 0.2 # reproduction rate of the predator
re = 1.0 # reproduction cost of the predator due to parasite exposure
dy = 1.0 # intrinsic death rate of the predator

### PARAMETERS THAT VARY BETWEEN SIMULATIONS (CAN BE CHANGED TO TRY OTHER PARAMETERS OF INTEREST)###
# infection probability prey
Qx = float(sys.argv[1]) # first position in command line set in shell script to vary the parameter values in for loop; this parameter can be changed to fixed value simply by adding a number instear of assigning a command line
# infection probability predator
Qy = float(sys.argv[2]) # second position in command line set in shell script to vary the parameter values in for loop; this parameter can be changed to fixed value simply by adding a number instear of assigning a command line
# reproduction cost of the prey due to parasite infection
rx = float(sys.argv[2]) # third position in command line set in shell script to vary the parameter values in for loop; this parameter can be changed to fixed value simply by adding a number instear of assigning a command line
# reproduction cost of the predator due to parasite infection
rp = float(sys.argv[2]) # third position in command line set in shell script to vary the parameter values in for loop; this parameter can be changed to fixed value simply by adding a number instear of assigning a command line

# Empty arrays to store the abundances of subpopulations
store_uninfected_prey = []
store_infected_prey = []
store_uninfected_predator = []
store_infected_predator = []

# Variables to record extinctions (simulation puporses)
extinction_prey = False
extinction_predator = False
extinction_parasite = False

# Assign extinctions (output purposes; this is done at the end of the script after the algorithm finishes)
## if there is no coexistence, then we assign "extinction result"
coexistence = "0"
coexistence_predator_and_prey = "0"
coexistence_prey = "0"
extinction = "1"

### Assign color for outcome of the frequencies of infected prey and infected predators (output purposes; this is done at the end of the script after the algorithm finishes)
## if there is no coexistence, then we assign a "white color"
cyan = "0"
yellow = "0"
black = "0"
magenta = "0"
white = "1"

# Time variables
recording_time = 130 # start recording the abundance of host subpopulations in time point 130 (for a total of 20 time points as the final time is 150, i.e., population sizes are stable and less computer memory is required)
max_time = 150 # determine final time (the algorithm will run until reaching max time)
Time = 0 # initial Gillespie time (we add dt_next_event and run until reaching max_time)
dt_next_event = 0 # random time step after event occurs (following the Gillespie algorithm). This quantity is summed to the total time (continuos time simulation)
n = 0 # number of time points across simulations in which we record abundances of subpopulations (in units of one, i.e., n+1)

### GILLESPIE ALGORITHM ###
while Time < max_time and (uninfected_prey + infected_prey > 0): # ALGORITHM STARTS: repeat simulation until reaching max time

## STEP 1 ##
## CALCULATE VALUES FOR EVERY POSSIBLE EVENT ACCORDING TO CURRENT CONDITIONS ##      
###### events uninfected prey ######
    prey_growth = uninfected_prey * gx # uninfected prey reproduuces
    prey_death = uninfected_prey * dx # uninfected prey dies
    prey_competition = uninfected_prey * (uninfected_prey + infected_prey) * (1 /pop_limit) # uninfected prey dies due to competition

###### events infected prey ######
    infected_prey_growth = infected_prey * gx * rx # infected prey reproduces
    infected_prey_death = infected_prey * dx # infected prey dies
    infected_prey_competition = infected_prey * (infected_prey + uninfected_prey) * (1 /pop_limit) # infected prey dies due to competition

###### events free-living parasite ######
    infection_prey = parasite * uninfected_prey * Qx * S # free-living parasite infects prey
    non_infection_prey = parasite * uninfected_prey * (1-Qx) * S # free-living parasite fails infecting prey
    parasite_death = parasite * dz # free-living parasite dies

###### events uninfected predator ######
    predator_growth = uninfected_predator * uninfected_prey * fy * ky # predator reproduces after feeding
    predator_non_growth = uninfected_predator * uninfected_prey * fy * (1-ky) # predator does not reproduce after feeding

    predator_exposure_growth = uninfected_predator * infected_prey * fy * (1-Qy) * re * ky # uninfected predator exposed to parasite reproduces
    predator_exposure_non_growth = uninfected_predator * infected_prey * fy * (1-Qy) * (1 - (re * ky)) # uninfected predator exposed to parasite does not reproduce
    predator_infection_growth = uninfected_predator * infected_prey * fy * Qy * rp * ky # uninfected predator infected by parasite reproduces
    predator_infection_non_growth = uninfected_predator * infected_prey * fy * Qy * (1 - (rp * ky)) # uninfected predator infected by parasite does not reproduce
    predator_death = uninfected_predator * dy # predator dies
    
    # events infected predator
    infected_predator_growth = infected_predator * uninfected_prey * fy * rp * ky # infected predator reproduces after feeding
    infected_predator_non_growth = infected_predator * uninfected_prey * fy * (1 - (rp * ky)) # infected predator does not reproduce after feeding
    
    infected_predator_exposure_growth = infected_predator * infected_prey * fy * (1-Qy) * re * rp * ky # infected predator exposed to the parasite reproduces
    infected_predator_exposure_non_growth = infected_predator * infected_prey * fy * (1-Qy) * (1 - (re * rp * ky))  # infected predator exposed to the parasite does not reproduce
    infected_predator_infection_growth = infected_predator * infected_prey * fy * Qy * rp * rp * ky # infected predator infected by parasite reproduces
    infected_predator_infection_non_growth = infected_predator * infected_prey * fy * Qy * (1 - (rp * rp * ky)) # infected predator infected by parasite does not reproduce
    infected_predator_death = infected_predator * dy # infected predator dies
    
    # Store the sum all events
    sum_events = (prey_growth + prey_death + prey_competition + 
    infected_prey_growth + infected_prey_death + infected_prey_competition + 
    infection_prey + non_infection_prey + parasite_death + predator_growth + predator_non_growth + 
    predator_exposure_growth + predator_exposure_non_growth + 
    predator_infection_growth + predator_infection_non_growth + predator_death + 
    infected_predator_growth + infected_predator_non_growth + 
    infected_predator_exposure_growth + infected_predator_exposure_non_growth  + 
    infected_predator_infection_growth + infected_predator_infection_non_growth + infected_predator_death)

    ## STEP 2 ##
    ## CALCULATE NEXT TIME STEP AND NEXT EVENT ##
    ## GENERATE RANDOM NUMBERS ##

    # Next time step
    dt_next_event = np.random.exponential(scale=1/sum_events)

    # Next event
    URN = random.uniform(0,1) # unit-interval uniform random number generator for next event
    P = 0 # variable for doing cummulative sum in picking the next event

    ## STEP 3 ##
    ## THE NEXT EVENT HAPPENS, UPDATE POPULATION SIZES AND ADD TIME STEP TO TOTAL TIME ##
    
    ############### Uninfected prey ####################
    occurrence = False
    while not occurrence:
        if URN > P and URN <= P + prey_growth/sum_events:
            bx = 1 # uninfected prey increases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += prey_growth/sum_events #! use += to modify in place
        if occurrence:
            break
            
        if URN > P and URN <= P + prey_death/sum_events:
            bx = 0 # uninfected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += prey_death/sum_events
        if occurrence:
            break
    
        if URN > P and URN <= P + prey_competition/sum_events:
            bx = 0 # uninfected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += prey_competition/sum_events
        if occurrence:
            break

    ############### Infected prey #################
        if URN > P and URN <= P + infected_prey_growth/sum_events:
            bx = 6 # infected prey reproduces
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += infected_prey_growth/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + infected_prey_death/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += infected_prey_death/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + infected_prey_competition/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += infected_prey_competition/sum_events
        if occurrence:
            break
                                        
    ################ Free-living parasite ####################
        if URN > P and URN <= P + infection_prey/sum_events:
            bx = 3 # prey is carrying a parasite (now it is infected)
            bz = 0 # parasite gets in the prey
            by = 2 # nothing happens to predator
            occurrence = True
        P += infection_prey/sum_events
        if occurrence:
            break
        
        if URN > P and URN <= P + non_infection_prey/sum_events:
            bx = 2 # nothing happens to prey (it is not infected)
            bz = 2 #  nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += non_infection_prey/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + parasite_death/sum_events:
            bx = 2 # nothing happens to prey
            bz = 0 # free-living parasite decreases by one
            by = 2 # nothing happens to predator
            occurrence = True
        P += parasite_death/sum_events
        if occurrence:
            break

    ################ Uninfected predator #####################
        if URN > P and URN <= P + predator_growth/sum_events:
            bx = 0 # uninfected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 1 # predator increases by one
            occurrence = True
        P += predator_growth/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + predator_non_growth/sum_events:
            bx = 0 # uninfected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += predator_non_growth/sum_events
        if occurrence:
            break
        
        if URN > P and URN <= P + predator_death/sum_events:
            bx = 2 # nothing happens to prey
            bz = 2 # nothing happens to free-living parasite
            by = 0 # predator decreases by one
            occurrence = True
        P += predator_death/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + predator_exposure_growth/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 1 # predator increases by one
            occurrence = True
        P += predator_exposure_growth/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + predator_exposure_non_growth/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += predator_exposure_non_growth/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + predator_infection_growth/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 1 # parasite reproduces in the predator
            by = 3 # predator is infected and increases by one
            occurrence = True
        P += predator_infection_growth/sum_events
        if occurrence:
            break
                        
        if URN > P and URN <= P + predator_infection_non_growth/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 1 # parasite reproduces in the predator
            by = 4 # predator is infected and does not reproduce
            occurrence = True
        P += predator_infection_non_growth/sum_events
        if occurrence:
            break

    ################ Infected predator #####################
        if URN > P and URN <= P + infected_predator_growth/sum_events:
            bx = 0 # uninfected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 6 # predator increases by one
            occurrence = True
        P += infected_predator_growth/sum_events
        if occurrence:
            break
                            
        if URN > P and URN <= P + infected_predator_non_growth/sum_events:
            bx = 0 # uninfected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += infected_predator_non_growth/sum_events
        if occurrence:
            break
    
        if URN > P and URN <= P + infected_predator_exposure_growth/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 6 # predator increases by one
            occurrence = True
        P += infected_predator_exposure_growth/sum_events
        if occurrence:
            break
                            
        if URN > P and URN <= P + infected_predator_exposure_non_growth/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 2 # nothing happens to free-living parasite
            by = 2 # nothing happens to predator
            occurrence = True
        P += infected_predator_exposure_non_growth/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + infected_predator_infection_growth/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 1 # parasite reproduces in predator
            by = 6 # uninfected predator increases by one
            occurrence = True
        P += infected_predator_infection_growth/sum_events
        if occurrence:
            break
                                
        if URN > P and URN <= P + infected_predator_infection_non_growth/sum_events:
            bx = 5 # infected prey decreases by one
            bz = 1 # parasite reproduces in the predator
            by = 2 # nothing happens to predator
            occurrence = True
        P += infected_predator_infection_non_growth/sum_events
        if occurrence:
            break

        if URN > P and URN <= P + infected_predator_death/sum_events:
            bx = 2 # nothing happens to prey
            bz = 2 # nothing happens to free-living parasite
            by = 5 # infected predator decreases by one
            occurrence = True
        P += infected_predator_death/sum_events
        if occurrence:
            break

##################### PREY EVENTS #########################
    if(bx == 1): # uninfected prey reproduces
        uninfected_prey += 1

    if(bx == 0): # uninfected prey dies
        uninfected_prey -= 1
    
    if(bx == 3):# free-living parasite infects prey
        uninfected_prey -= 1
        infected_prey += 1

    if(bx == 5): # infected prey dies
        infected_prey -= 1

    if(bx == 6): # infected prey reproduces
        uninfected_prey += 1

################### PARASITE EVENTS #######################
    if(bz == 0): # free-living parasite dies
        parasite -= 1
    
    if(bz == 1): # parasite reproduces in the predator
        parasite += n_z
    
################### PREDATOR EVENTS ########################
    if(by == 1): # uninfected predator reproduces
        uninfected_predator += 1

    if(by == 0): # uninfected predator dies
        uninfected_predator -= 1

    if(by == 3): # uninfected predator gets infected and reproduces
        infected_predator += 1

    if(by == 4): # uninfected predator gets infected and does not reproduce
        uninfected_predator -= 1
        infected_predator += 1

    if(by == 5): # infected predator dies
        infected_predator -= 1
        
    if(by == 6): # infected predator reproduces
        uninfected_predator += 1

    # Advance a step in time
    Time += dt_next_event # continuous time simulation
    
    # Sum uninfected and infected hosts
    sum_prey = uninfected_prey + infected_prey # for extinction purposes and calculate total abundance of prey
    sum_predator = uninfected_predator + infected_predator # for extinction purposes and calculate total abundance of predator
    sum_parasite = parasite + infected_prey + infected_predator # for extinction purposes

    # Record extinctions
    if(sum_prey <= 0):
        extinction_prey = True
    if(sum_predator <= 0):
        extinction_predator = True
    if(sum_parasite <= 0):
        extinction_parasite = True

    # Record abundance subpopulations
    if Time > n and not extinction_parasite and not extinction_predator and not extinction_prey:
        if n > recording_time:
            store_uninfected_prey.append(uninfected_prey)
            store_infected_prey.append(infected_prey)
            store_uninfected_predator.append(uninfected_predator)
            store_infected_predator.append(infected_predator)
        n += 1

# Record average of abundances of subpopulations
if store_uninfected_prey != []:
    av_uninfected_prey = sum(store_uninfected_prey[:]) / len(store_uninfected_prey)
else:
    av_uninfected_prey = 0
if store_infected_prey != []:
    av_infected_prey = sum(store_infected_prey[:]) / len(store_infected_prey)
else:
    av_infected_prey = 0
if store_uninfected_predator != []:
    av_uninfected_predator = sum(store_uninfected_predator[:]) / len(store_uninfected_predator)
else:
    av_uninfected_predator = 0
if store_infected_predator != []:
    av_infected_predator = sum(store_infected_predator[:]) / len(store_infected_predator)
else:
    av_infected_predator = 0

# Record coexistence/extinctions
if(sum_prey <= 0):
    extinction_prey = True
if(sum_predator <= 0):
    extinction_predator = True
if(sum_parasite <= 0):
    extinction_parasite = True

if not extinction_prey and not extinction_parasite and not extinction_predator:
    coexistence = "1"
    coexistence_predator_and_prey = "0"
    coexistence_prey = "0"
    extinction = "0"
if not extinction_prey and extinction_parasite and not extinction_predator:
    coexistence = "0"
    coexistence_predator_and_prey = "1"
    coexistence_prey = "0"
    extinction = "0"
if not extinction_prey and extinction_parasite and extinction_predator:
    coexistence = "0"
    coexistence_predator_and_prey = "0"
    coexistence_prey = "1"
    extinction = "0"
if extinction_prey and extinction_parasite and extinction_predator:
    coexistence = "0"
    coexistence_predator_and_prey = "0"
    coexistence_prey = "0"
    extinction = "1"

### Assign color for outcome of the frequencies of infected prey and infected predators
## if there is no coexistence, then we assign a white color
if av_uninfected_prey < av_infected_prey and av_uninfected_predator < av_infected_predator and coexistence == "1":
    ## higher frequency of infected prey and infected predators
    cyan = "1"
    yellow = "0"
    black = "0"
    magenta = "0"
    white = "0"

if av_uninfected_prey >= av_infected_prey and av_uninfected_predator < av_infected_predator and coexistence == "1":
    ## higher frequency of uninfected prey and infected predators
    cyan = "0"
    yellow = "1"
    black = "0"
    magenta = "0"
    white = "0"

if av_uninfected_prey < av_infected_prey and av_uninfected_predator >= av_infected_predator and coexistence == "1":
    ## higher frequency of infected prey and uninfected predators
    cyan = "0"
    yellow = "0"
    black = "1"
    magenta = "0"
    white = "0"

if av_uninfected_prey >= av_infected_prey and av_uninfected_predator >= av_infected_predator and coexistence == "1":
    ## higher frequency of uninfected prey and uninfected predators
    cyan = "0"
    yellow = "0"
    black = "0"
    magenta = "1"
    white = "0"

# Save output (a row will be appended in output file every time a simulation fnishes until "for loop" finishes running in the shell script for every parameter value)
List = [str(Qx),str(Qy),str(rx),str(rp),str(coexistence),str(coexistence_predator_and_prey),str(coexistence_prey),str(extinction),str(cyan),str(yellow),str(black),str(magenta),str(white)]

# Open our existing CSV file in append mode (a row will be appended in output file every time a simulation fnishes until "for loop" finishes running in the shell script for every parameter value)
# Create a file object for this file
with open("output.csv", 'a') as f_object:
    writer_object = writer(f_object)
    writer_object.writerow(List)
    f_object.close()

### When all the rows have been appended run the following code in a separate python script called "average_timepoints.py"###
#import sys
#import os
#import pandas as pd

#m_df = pd.read_csv("output.csv")
#df = m_df.groupby(['Qx','Qy','rp','rx']).sum()
#df.to_csv("Average_output.csv", index = True)

### These lines will sum all the values from the simulation outputs. You will end up with a .csv file called "Average_output.csv" that contains the total number of occorence of coexistence events and frequencies of subpopulations
### The frequencies of subpopulations are assigned to a color which represents the percentage of color for each panel in our plots