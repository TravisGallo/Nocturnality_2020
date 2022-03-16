# R script, JAGS models, and data used to assess how urbanization alters diel activity of mammals

# Mammals adjust diel activity across gradients of urbanization

# **Authors**

Travis Gallo^1^,^2^, Mason Fidino^2^, Brian Gerber^3^, Adam A. Ahlers^4^, Julia L. Angstmann^5^, Max Amaya^6^, Amy L. Concilio^7^, David Drake^8^, Danielle Gray^9^, Elizabeth W. Lehrer^2^, Maureen H. Murray^2&, Travis J. Ryan^5^, Colleen Cassady St. Clair^10^, Carmen M. Salsbury^5^, Heather A. Sander^11^, Theodore Stankowich^6^, Jacque Williamson^12^, J. Amy Belaire^13^, Kelly Simon^14^, Seth B. Magle^2^

1.  Environmental Science and Policy, College of Science, George Mason University, Fairfax, VA 22030 USA

2.  Urban Wildlife Institute, Conservation and Science Department, Lincoln Park Zoo, Chicago, IL 60614 USA

3.  Department of Natural Resource Science, The University of Rhode Island, Kingston, RI 02881, USA

4.  Department of Horticulture and Natural Resources, Kansas State University, Manhattan, KS 66502 USA

5.  Department of Biological Sciences and Center for Urban Ecology and Sustainability, Butler University, Indianapolis, IN 46208 USA

6.  Department of Biological Sciences, California State University Long Beach, Long Beach, CA 90840 USA

7.  Department of Environmental Science and Policy, St. Edward’s University, Austin, TX 78704 USA

8.  Austin Parks and Recreation, City of Austin, TX 78704 USA

9.  Department of Forest and Wildlife Ecology, University of Wisconsin-Madison, Madison, WI 53706, USA

10. Department of Biological Sciences, University of Alberta, Edmonton, Canada

11. Department of Geographical and Sustainability Sciences, University of Iowa, Iowa City, IA 52242 USA

12. Department of Education & Conservation, Brandywine Zoo, Wilmington, Delaware, 19802 USA

13. The Nature Conservancy in Texas, San Antonio, Texas 78215 USA

14. Texas Parks and Wildlife, Austin, Texas 78774 USA

# **File Descriptions**

## **Scripts:**

**2021_Gallo_etal_DielActivity_UtilityScripts.R** - a script that contains all the functions used in this analysis. Source at very beginning to load needed packages and functions.

Technically, you only need three scripts to run this analysis. These scripts are located in the `Analysis` subfolder:

**Scripts/Analysis/2021_Gallo_etal_DielActivity_FiguresScript.R** – the only file that you need to open. Sources the two files below to load data sets, format data for JAGS model, run JAGS model, and create figures from manuscript.

**Scripts/Analysis/2021_Gallo_etal_DielActivity_fitJAGS.R** - script that loads data sets, formats data to be fit in JAGS model, and fits the models using JAGS.

**Scripts/Analysis/2021_Gallo_etal_DielActivity_SoftmaxModel** - multinomial model with categorical LASSO regularization used to estimate temporal resource selection. Read in `2021_Gallo_etal_DielActivity_fitJAGS.R`.

We have also include scripts for processing photo data and calculating landcape covariates. These scripts are located in the `DataProcessing` subfolder:

**Scripts/DataProcessing/2021_Gallo_etal_DielActivity_CovariateSetup.R** - script used to calculate landscape covariates at species-specific scales. Note that this script requires several spatial datasets that are publicly avaliable, but we do not provide in this repository. Sources are indicated in the script and/or the manuscript. This script produces `2021_Gallo_etal_DielActivity_data_list_8sp.rds` which is read in `021_Gallo_etal_DielActivity_fitJAGS.R`.

**Scripts/DataProcessing/2021_Gallo_etal_DielActivity_photo_data_setup_10cities.R** - script used to clean and organize photo observations from 10 U.S. cities. This script produces `2021_Gallo_etal_DielActivity_cleaned_photo_dataset_9_species.rds`, which is read into `2021_Gallo_etal_DielActivity_CovariateSetup.R`.

## **Data Files:**

**2021_Gallo_etal_DielActivity_data_list_8sp.rds** - This is a list that contains data for each species that we ultimately used. Each list element contains observation data and site-level covariate data for a single species. This data is read into `2021_Gallo_etal_DielActivity_fitJAGS.R`.

**2021_Gallo_etal_DielActivity_sitecovsBuffered.rds** - These are site level covariates calculated at three scales - 500-m, 1000-m, and 1500-m buffers.

**2021_Gallo_etal_DielActivity_cleaned_photo_dataset_9\_species.rds** - A cleaned dataset of photo observations for only the species and cities that we analyzed.

**CityFiles/** - this folder contains the raw observation data for each city. Read in `2021_Gallo_etal_DielActivity_photo_data_setup_10cities.R`.

**group_mean_centered_variables.rds** - covariate values that are group mean centered by city and scaled by global standard deviation. These are separated from our original data, because the procedure was performed post-hoc to satisfy peer-review.

**Note:** All of these files must be within your working directory for the analysis to work. Several analysis were done in parallel. Therefore, you will need to adjust the settings accordingly.
