# TUNLtasks

## 1) Prcoess all recording sessions with CNMF-E
## 2)	Register cells across recordings. This is particularly important since I want to see how individual cells LEARN the contingencies of the TUNL task. This means registering the cell across ALL sessions, from the first to final recording session. I realize this won’t be possible for all cells, and many cells will likely pop in and out.
	Figure 1 - Cells tracked. Show which cells successfully tracked across days. I am imagining a checkerboard plot with cell ID on the Y-axis, recording session on the X-axis. Successful registration would be filled by a black square, unsuccessful with a white square. 
	
## 3) Reliability of behavioral correlate
	Step 1 - align data to each of the anchoring points of the task. (trial initiation, touch sample, touch choice, reward). This should only be done on the actual TUNL task. Later we will look at what cells are doing during all the training prior to the TUNL task. 
	Step 2 - similar to split half reliability, randomly divide the trials in half and compute the pearson correlation of the mean activity of each half. Repeat this process 100 times. Return the average correlation value of these 100 iterations. Compare these values to shuffled data (circular shift) and find the 99th percentile. I am curious what percentage of these cells pass this threshold for each anchoring point. This should be done on the entire trial (NOT just the delay period!) and should be done on all anchoring points (trial initiation, touch sample, touch choice, reward) independently.  All analysis will be done on cells that pass this criteria unless otherwise specified. 
		Figure 1 - Sequence of population data. Y-axis is the cell ID, X-axis is the duration of the trial. Plot mean fluorescence of the cells that pass criteria, ordered by their peak fluorescence. Plots for each anchoring point can be made. 
		Figure 2 - Population representation of task contingencies. Normalize the fluorescence of each neuron between 0-1, average these values across all cells to create a single line plot of the activity across the trial. This plot can be made for each anchoring point, with a label to show where in the task the anchoring point is. 

## 4) Memory related firing
	Step 1. Of cells that pass the reliability threshold create separate plots based on the location of the sample location (position 1 or 5 in the current data, but this will include positions 1-5 in the future). 
		Figure 1. Sample selective activity of individual cells. Plot rasters for each neuron based on the location of the sample. 
		Figure 2. Sample selective activity in the population. Plot mean activity of all cells ordered by peak fluorescence based on the location of the sample. 

Step 2. What is happening during error trials? For each cell compare the activity during correct trials to error trials. This analysis should be done independently for each anchoring point. 
	Figure 1. Visual comparison of correct vs error trials. In one figure, show rasters of each correct trial and each type of error trials. 
	Figure 2. Population comparison of correct vs error trials. In one figure, show population sequence (mean fluorescence, order by peak) of each correct trial and each type of error trials. 
	Metric 1. Compute correlation between mean fluorescence activity of correct versus incorrect trials. 
	Metric 2. Compute correlation between mean fluorescence activity of correct versus subsampled correct trials. Based on the number of error trials, you will randomly select and remove that number of correct trials. Computer the correlation between the mean fluorescence of these trials. Repeat this random selection and correlation process 100 times. Metric 2 is the average of these 100 correlation values. 
	Figure 3. Comparison of correlation values between correct and incorrect trials. CDF plot of this Metric 1 and Metric 2 of each cell. 
	Figure 4. Locating the source of the error signal. Plot Metric 1 and Metric 2 as a function of time during the trials. 

	
## 5) Influence of spatial location on activity patterns during the TUNL task. 
Step 1) Create spatial firing maps of all cells (even those that do not pass the reliability test described above). Compute spatial information and split half reliability of the spatial map for each cells. 
Figure 1 - Comparison of the TUNL task reliability score to the spatial reliability score. This will be two scatterplots. TUNL reliability versus spatial information, and TUNL reliability versus split-half spatial map correlation. Here will see if there are dedicated spatially tuned neurons or TUNL task related neurons.

Step 2) Create a generalized linear model to determine if spatial activity explains TUNL task reliability. This will be done only on neurons that pass the TUNL task reliability score. See this Eichenbaum paper for how to do this analysis. https://www.ncbi.nlm.nih.gov/pubmed/23707613


## 6) Evolution of TUNL task coding during learning 
Step 1) For each cell registered across days…
1)	Plot TUNL reliability across days (line plot for each cell) this will give us indication of whether the code was stable across days 
2)	For each cell, plot mean fluorescence for each day. In this case Y-axis is the day, X-axis is the time within the trial.

Here, we will want to re-calculate the TUNL task reliability across all tria
Figure 1 - Single cell reliability as a function of task performance. Scatterplot of TUNL task reliability scores of all cells across all sessions as a function of the behavioral performance. 

Step 2) Determine how time cells (delay activity cells) change when the duration of the delay increases?
Figure 1 - Time cell activity across multiple delay periods. Find cells that pass the reliability threshold and have a peak of mean activity during the delay period. For each cell, show the rasters of all trials across sessions. Since the max delay is 8 seconds, the X axis should be 8 seconds. Trials that lasted only 2 second will be empty of course and displayed in WHITE from seconds 2-8. 


## 7) Population level analysis
Step 1 - Can we predict correct versus incorrect trials based on population activity? Could use a support vector machine (SVM) analysis - train the data one session and apply to the next session. Of course this means only using cells that are registered across the two sessions. At what point within a trial can we estimate the likelihood of success in selecting the correct choice???? This would be amazing and give insight into which aspects of population activity are required for task performance. Does activity during the sample, prior to the delay predict behavioral performance? Or does it even predict the fidelity of the delay activity???
Step 2 - PCA exploratory analysis. Examine PCA state space as a function of task contingencies and performance. I am imagining that we can look at progression through state space on a trial by trial basis. How consistent is movement through state space? How does this relate to behavioral performance? How does this influence the fidelity of time cells? 





