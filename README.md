# WhittJansen_public
Matlab code for the box model reported in Whitt and Jansen (2020) "Slower nutrient stream suppresses SubarcticAtlantic biological productivity in global warming"
Submitted to PNAS. January 15, 2020.

This code should not be treated as a black box. Contact us with questions: dwhitt@ucar.edu, mfj@uchicago.edu

Running scripts starting "plot" will generate several of the figure files in the paper.

We suggest starting with plot_figS1S2S3.m, which will run and visualize a single scenario,
beginning with a 200-year "early-21st-century equilibrium" scenario with climatechangeflag=0.
 
To run the transient global warming scenarios in Fig. 2, switch the climatechangeflag as follows:
climatechangeflag = 1 (MLD and PSI changing together); fig S3
climatechangeflag = 2 (only MLD changing, fixed PSI)
climatechangeflag = 3 (only PSI changing, fixed MLD)

