# WhittJansen_public
Matlab code for the box model reported in Whitt and Jansen (2020) "Slower nutrient stream suppresses SubarcticAtlantic biological productivity in global warming"
Submitted to PNAS. January 15, 2020.

This code should not be treated as a black box. Contact us with questions.

Running scripts in the top directory with names starting "plot" will generate several of the figure files in the paper.

Scripts with names starting "run" are helper scripts to run various scenarios and save output in matfiles for quick visualization.

The subdirectory 1Dmodel contains the main model source code.

The subdirectory utility contains various utility scripts.

The subdirectory parameter_selection includes scripts to generate 248,832 sensitivity simulations.
This takes several hours on a 36-core node on the Casper system at NCAR using the matlab parallel toolbox.
The output data is saved already, so relevant analysis scripts can be run quickly.

A good starting point is: plot_figS1S2S3.m, which will run and visualize a single scenario,
beginning by default with a 200-year "early-21st-century equilibrium" scenario with climatechangeflag=0.
 
To run the transient global warming scenarios in Fig. 2, switch the climatechangeflag as follows:
climatechangeflag = 1 (MLD and PSI changing together); fig S3
climatechangeflag = 2 (only MLD changing, fixed PSI)
climatechangeflag = 3 (only PSI changing, fixed MLD)

The results of these simulations can also be visualized from existing output in mat files using plot_fig2.m.
And, the more extensive results of figure 3 can be visualized using plot_fig3.m using existing output in mat files.

