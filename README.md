# Healthcare_NaturalDisaster_Pandemic
This repository set is for the case study discussed in the manuscript titled "Orchestrating Performance of Healthcare Networks Subjected to the Compound Events of Natural Disasters and Pandemic"

Orchestrating Healthcare Networks Subjected to the Compound Events of Natural Disasters and COVID-19 Pandemic
Emad M. Hassan, Hussam Mahmoud*
Department of Civil and Environmental Engineering, Colorado State University, Fort Collins, CO, USA

Contents
Overview
Code Contents
System Requirements
       Hardware Requirements
       Software Requirements
Installation Guide
Demo
Run using your data
Results

Overview
We investigate the impact of the combined effects of natural disasters, focusing on wildfires, and disease outbreak on the operation of healthcare systems. We use Butte County, California as a testbed and utilize the Camp Fire data combined with different spread courses of COVID-19, using a modified SEIR disease transmission model, to evaluate the demand on the healthcare system. The capacity of the healthcare system is quantified using a fault tree describing the availability of staffed beds over time. The quality of the hospitalization service is quantified using the accessibility and effectiveness of the offered service in terms of patients’ waiting time and treatment time. The total functionality is described as a combination of quantity and quality of the hospitalization service. Patients are distributed using the patient-driven model. The interaction between the hospitals is simulated using the hospital interaction model and includes patient, staff, and resources transfer. Here we present a sample scenario, which assumes the wildfire occurred 40 days after the epidemic outbreak. 
Code Contents
The main code file
Healthcare_Wildfire_Epedimic_Model.mat: the main code file developed using Matlab software and including input, analysis, and embedded function sections.

The input files
model_input.mat: to present the main model parameters including a) healthcare system input (number of beds, number of staff, rating, ambulance availability, treatment ability, agreement, etc.), b) community-related input (population data, insurance, private transportation, travel time, and housing functionality as well as the functionality of power, water, wastewater, transportation, telecommunications, fuel), and c) the total number of regular patients at each census tract. 

demand.mat: to introduce the number of disease-related patients at each census tract and each age group. This demand classified based on the type of hospitalization service needed for each patient, which includes ER, inpatients, ICU, and M. ventilators.

PDM_input.mat: to present is the input for the patient-driven model, including P1: Case criticality, P2: Insurance availability, P3: Media effect “word of mouth”, P4: Positive previous experience availability, P5: Brand availability, P6: Transportation functionality, P7: Detour functionality, P8: Less travel time availability, P9: Less waiting time, P10: Hospital ability to treat specific injuries, P11: Ambulance availability, and P12: Private transportation availability.

Patient_dist0.mat: to distribute the patient demand on each healthcare facility on the first day.

HIM_input.mat: to present the input for the hospitals' interaction model including I1: Case criticality, I2: Insurance availability, I3_4: Transportation, I5: Travel time, I6: Agreement, I7: Telecommunication, I8: Waiting time, I9: Treatment ability, I10: Ambulance, and I11: Air transportation.

functionality_input_em.mat: to introduce the parameters of the success tress used to calculate the total number of ER beds including R1_em: Physicians availability, R2_em: Nurses availability, R3_em: Supporting staff availability, R4_em: Alternative staffing availability, R5_em: Corridors functionality, R6_em: Elevator functionality, R7_em: Stairs functionality, R8_em: Municipal water functionality, R9_em: Backup water system functionality, R10_em: Municipal power functionality, R11_em: Backup power system functionality, R12_em: Transportation network functionality, R13_em: Transportation detours availability, R14_em: Ambulance service functionality, R15_em: Telecommunication service functionality, R16_em: Backup telecom service functionality, R17_em: Municipal wastewater functionality, R18_em: Backup wastewater functionality, R19_em: Drinking water system functionality, R20_em: Backup drinking water functionality, R21_em: Structural components functionality, R22_em: Non-structural components functionality, R23_em: Contents functionality, R24_em: Backup space functionality, R25_em: Oxygen availability, R26_em: Surgical supply availability, R27_em: Rx availability, R28_em: Fuel supply availability, R29_em: Food supply availability, and R30_em: Other supply availability.

functionality_input_in.mat: to introduce the parameters of the success tress used to calculate the total number of inpatient beds including R1_in: Physicians availability, R2_in: Nurses availability, R3_in: Supporting staff availability, R4_in: Alternative staffing availability, R5_in: Corridors functionality, R6_in: Elevator functionality, R7_in: Stairs functionality, R8_in: Municipal water functionality, R9_in: Backup water system functionality, R10_in: Municipal power functionality, R11_in: Backup power system functionality, R12_in: Transportation network functionality, R13_in: Transportation detours availability, R14_in: Ambulance service functionality, R15_in: Telecommunication service functionality, R16_in: Backup telecom service functionality, R17_in: Municipal wastewater functionality, R18_in: Backup wastewater functionality, R19_in: Drinking water system functionality, R20_in: Backup drinking water functionality, R21_in: Structural components functionality, R22_in: Non-structural components functionality, R23_in: Contents functionality, R24_in: Backup space functionality, R25_in: Oxygen availability, R26_in: Surgical supply availability, R27_in: Rx availability, R28_in: Fuel supply availability, R29_in: Food supply availability, and R30_in: Other supply availability.

functionality_input_ICU.mat: to introduce the parameters of the success tress used to calculate the total number of ICU beds including R1_ICU: Physicians availability, R2_ICU: Nurses availability, R3_ICU: Supporting staff availability, R4_ICU: Alternative staffing availability, R5_ICU: Corridors functionality, R6_ICU: Elevator functionality, R7_ICU: Stairs functionality, R8_ICU: Municipal water functionality, R9_ICU: Backup water system functionality, R10_ICU: Municipal power functionality, R11_ICU: Backup power system functionality, R12_ICU: Transportation network functionality, R13_ICU: Transportation detours availability, R14_ICU: Ambulance service functionality, R15_ICU: Telecommunication service functionality, R16_ICU: Backup telecom service functionality, R17_ICU: Municipal wastewater functionality, R18_ICU: Backup wastewater functionality, R19_ICU: Drinking water system functionality, R20_ICU: Backup drinking water functionality, R21_ICU: Structural components functionality, R22_ICU: Non-structural components functionality, R23_ICU: Contents functionality, R24_ICU: Backup space functionality, R25_ICU: Oxygen availability, R26_ICU: Surgical supply availability, R27_ICU: Rx availability, R28_ICU: Fuel supply availability, R29_ICU: Food supply availability, and R30_ICU: Other supply availability.

functionality_input_V.mat: to introduce the parameters of the success tress used to calculate the total number of M. ventilator beds including R1_V: Physicians availability, R2_V: Nurses availability, R3_V: Supporting staff availability, R4_V: Alternative staffing availability, R5_V: Corridors functionality, R6_V: Elevator functionality, R7_V: Stairs functionality, R8_V: Municipal water functionality, R9_V: Backup water system functionality, R10_V: Municipal power functionality, R11_V: Backup power system functionality, R12_V: Transportation network functionality, R13_V: Transportation detours availability, R14_V: Ambulance service functionality, R15_V: Telecommunication service functionality, R16_V: Backup telecom service functionality, R17_V: Municipal wastewater functionality, R18_V: Backup wastewater functionality, R19_V: Drinking water system functionality, R20_V: Backup drinking water functionality, R21_V: Structural components functionality, R22_V: Non-structural components functionality, R23_V: Contents functionality, R24_V: Backup space functionality, R25_V: Oxygen availability, R26_V: Surgical supply availability, R27_V: Rx availability, R28_V: Fuel supply availability, R29_V: Food supply availability, and R30_V: Other supply availability.

PDM_mod.mat: to estimate the change in patient-driven model parameters including P9: Waiting time and P10_em, P10_in, P10_ICU, P10_V: Different hospital bed availability.

HIM_mod.mat: to estimate the change in hospital interaction model parameters including I2: Insurance availability and I8: Waiting time.

Embedded Functions
Function (functionality_tree): to calculate the probability of occurrence of main and intermediate events in the success tree used to estimate hospitals’ capacity.
Function (Patient_dist): to distribute the total patient demand on the available hospital beds.
Function (Overflow): to calculate the patient overflow at each hospital and for each hospitalization service. 
Function (HIM_redist): to redistribute the patients to consider the patient transfer between the investigated hospitals.

System Requirements
Hardware Requirements
The code package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:
RAM: 8+ GB
CPU:  4+ cores, 3.3+ GHz/core
The runtimes below are generated using a computer with the recommended specs (32 GB RAM, 8 cores@3.6 GHz).
10 sec.

Software Requirements
The code developed and run using Matlab R2019b operating systems. 
Installation Guide
No installation required.
Demo
The code is divided into three main sections:
1.	Input section: in this section, the main model parameters are presented including: 
•	Healthcare providers' annual reports.
•	2018 Census reports.
•	Campfire data.
•	COVID-19 disease transmission data based on the SEIR model.
•	Other soci-technical data from OnTheMap website, Google maps, OSHPD, California Department of Health, insurance, etc.

2.	Analysis section: this section uses the input parameters to 
•	Calculation of quantity, quality, and total functionality for different bed types including ER, inpatients, ICU, and M. ventilators. 
•	Calculation of the available number of each type of staffed beds.
•	Distribution of patients using Patient-driven model.
•	Transferring of patients from saturated hospitals using the hospital interaction model.
•	Calculation of the daily number of treated and untreated number of patients as well as the healthcare functionality indicators. 

3.	Embedded functions:  the impeded functions include:
•	Success tree analysis using functionality_tree 
•	Probability tree analysis using Patient_dist
•	Calculation of patient overflow using Overflow
•	Redistribution of patients using HIM_redist

Run using your data
The following input files can be replaced to run the analysis using your data:
model_input.mat, demand.mat, PDM_input.mat, Patient_dist1.mat, HIM_input.mat, functionality_input_em.mat, functionality_input_in.mat, functionality_input_ICU.mat, functionality_input_V.mat, PDM_mod.mat, and HIM_mod.mat.   
Results
The provided sample analysis investigated the patients’ distribution and healthcare functionality indicators over time during the epidemic and wildfire in Butte County, CA. The epidemic outbreak took place on day 10 and the wildfire occurred on day 10. The epidemic data is based on CoVID-19 disease transmission rates in the US and the wildfire data is collected from the 2018 Campfire in Paradise city. 
