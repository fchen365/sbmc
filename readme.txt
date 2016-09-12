Filename	--	Content

initLable.r	--	Generate initial labeling
irlbaMod.R	--	modified version of the irlba function where A is a list instead of a matrix
VEM_SBM  	--	Variational EM regarding or disregarding covariates

Simulation related work
job.R   	--	Work required by Hoffman2
gen_demo	--	Generating a random positive-definite matrix with user-specified positive eigenvalues
MR_cal.R	--	Function to Compute misclassified rate
plot_MR.R	--	Plot the misclassified rate vs. number of blocks

Implement on real dataset (Facebook, Amazon)
getRes.R	--	Function to return blockmodel results
EgoReader.R	--	Function to read input network dataset
CirGraph.R	--	Function to generate circle graph. V1: Split the overlaping circles; V2: Directly remove the overlaping node
casc.R  	--	Covariate Assisted Spectral Clustering
facebook2.R	--	Main flow of analysis on Facebook dataset
amazon.R	--	Main flow of analysis on Amazon dataset