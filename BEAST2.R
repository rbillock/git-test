#Download and install Rtools on your machine before running this if trying to do BEAST2 in R. 
#If just making basic xml file, don't bother with devtools and Rtools

#Need to do this with R closed and then reopen
#install.packages("devtools", repos = 'http://cran.us.r-project.org')
library(devtools)
find_rtools()

devtools::install_github("richelbilderbeek/babette")
library(babette)
library(ape)
library(geiger)
library(devtools)
library(ggplot2)
library(knitr)
library(phangorn)
library(rmarkdown)
library(seqinr)
library(stringr)
library(testit)
library(tracerer)

memory.limit()
memory.limit(size=20000)

#assign working library
setwd("J:/Transmission Networks/Dynamics Paper/BEAST_2")
set.seed(12345)
install_beast2()
create_beast2_input_file((input_filenames = "cluster265.fasta"), output_filename = "beast265.xml")
create_beast2_input_file((input_filenames = "cluster325.fasta"), output_filename = "beast325.xml")
create_beast2_input_file((input_filenames = "cluster355.fasta"), output_filename = "beast355.xml")
create_beast2_input_file((input_filenames = "cluster618.fasta"), output_filename = "beast618.xml")
create_beast2_input_file((input_filenames = "cluster652.fasta"), output_filename = "beast652.xml")
create_beast2_input_file((input_filenames = "cluster715.fasta"), output_filename = "beast715.xml")
create_beast2_input_file((input_filenames = "cluster872.fasta"), output_filename = "beast872.xml")
create_beast2_input_file((input_filenames = "cluster891.fasta"), output_filename = "beast891.xml")
create_beast2_input_file((input_filenames = "cluster1076.2.fasta"), output_filename = "beast1076.2.xml")
create_beast2_input_file((input_filenames = "cluster1926.fasta"), output_filename = "beast1926.xml")

#BEAST2 in R - can't achieve necessary parameters
mcmc <- create_mcmc(chain_length = 10000, store_every = 1000)
out <- bbt_run(fasta_filenames = "cluster265.fasta", site_models = create_gtr_site_model(gamma_site_model = create_gamma_site_model(prop_invariant=0.5)), 
       clock_models = create_rln_clock_model(),tree_priors = create_bd_tree_prior(),mcmc=mcmc)
traces <- remove_burn_ins(traces = out$estimates, burn_in_fraction = 0.2)
sum_stats <- calc_summary_stats(traces = traces, sample_interval = 1000)

