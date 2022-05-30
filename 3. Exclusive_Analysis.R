

# Clean up the R envrionment 
rm(list=ls())

# Do you want the Inclusive or Exclusive analyses
data_use<-"Exclusive"

# Do you want to impute the income data for missing values as the average of an individuals income (note income post retirement is not predecited to change much)
impute_income<-FALSE

# By what method do you want to GAMs - default is "GCV.Cp", but "REML" is also very comoonly used
fit_with<-"REML"

# Set the directories for the data and the folder to write output too
home_dir<-"/Users/asenior/Dropbox (Sydney Uni)/Alan Cohen/Analyses"
home_dir<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Alan Cohen/Analyses"
write_dir<-paste0(home_dir, "/", paste0(data_use, "_ImpIncome_", impute_income))

# Set to the home dir to read data and funtions
setwd(home_dir)

# Load the packages
library(openssl)
library(plyr)
library(mixexp)
library(sp)
library(mgcv)
library(geometry)
library(MASS)
library(gplots)
library(Vennerable)
library(stringr)
library(ggplot2)
# Function for generatig colors for surfaces
source("0. Header_Functions.R")

# Read in the encrypted data
encrypted_data<-readRDS("Encrypted_Data.rds")

# Unencrypt
passphrase<-charToRaw("NuAgeOtag0Syd")
key<-sha256(passphrase)
data<-unserialize(aes_cbc_decrypt(encrypted_data, key=key))

# Remove the encrypted data
rm(encrypted_data)

# That is everything we need read in. From here we can work in the write directory
setwd(write_dir)

# Calculate BMI
data$BMI<-data$POIMES / data$TAIMES^2

# What outcomes do we want to model - here the physiological dysregulation scores
traits<-c("blood_2", "immune_2", "liver_2", "lipid_2", "vit_2", "global_2", "pheno_age", "BA10")
trait_label<-c("Oxygen Transport", "Leukopoiesis", "Liver/Kidney", "Lipids", "Micronutrients", "Global", "Phenotypic Age", "Biological Age")

# We will recode smoking status to current/not current
data$FUMER<-as.numeric((is.na(match(data$FUMER, c(1, 2))) == F))

##################################################
########### DROPPING MISSING DATA ################
##################################################

# Ditch any missing macronutrient intake data
data<-drop.missing(data=data, check=c("PROTEIN_kJ.d", "CARBO_kJ.d", "LIPID_kJ.d"))

# Pull out the names of the data on micronutrients 
names(data)[c(448:479)]
# Lets give those a nicer name
names(data)[c(448:479)]<-c("Calcium", "Iron", "Magnesium", "Phosphorous", "Potassium", "Sodium", "Zinc", "Selenium", "Vitamin_A", "A_Carotenoid", "A_Tocopherol", "B_Cryptoxan.", "Lycopene", "Lutein_Zeaxanthin", "Vitamin_D", "Vitamin_C", "Thiamine", "Riboflavin", "Niacin_Equiv.", "Pantothenic_Acid", "Vitamin_B6", "Vitamin_B12", "Vitamin_K", "Folic_Acid", "Folic_Acid_Equiv.", "Cholesterol", "Trans._F._Acids", "Sat._F._Acids", "Linoleic_2", "Linoleic_3", "Mono._F._Acids", "Poly._F._Acids")

# Which ones do we want to analyse - note we are leaving out a few based on Valerie's advice as they are high resolution decompostions of other more inclusive variables
micros<-names(data)[c(448:456, 458, 460:470, 472:475, 478:479)]

# Drop the missing micronutrient data
data<-drop.missing(data=data, check=micros)

# There is something a bit weird going on with the eduction data - there are missing values for some visits - I am going to fill these in with the average reported level for a subject (this should not change over the time of the followups)
edu<-ddply(data, .(sujetno), summarise, sum(is.na(SCOLAR)), mean(SCOLAR, na.rm=T), sd(SCOLAR, na.rm=T))
# Note SD is NA for all individuals suggesting education does not vary within individuals so this is fine
data$SCOLAR<-edu$..2[match(data$sujetno, edu$sujetno)]

# Impute the missing income data if possible for some individuals
inc<-ddply(data, .(sujetno), summarise, sum(is.na(REVFAM)), mean(REVFAM, na.rm=T), sd(REVFAM, na.rm=T))
# OK so there is some within-individual variance in income so we are making assumptions when we impute these missing values at their average. We need to know the sensitivity of our analyses to this so will do analysis with and without the imputation
if(impute_income == T){
	tag<-which(is.na(data$REVFAM) == T)
	data$REVFAM[tag]<-inc$..2[match(data$sujetno[tag], inc$sujetno)]
}


# These are the confounders we are interested in
confounders<-c("age_visit", "REVFAM", "SCOLAR", "ALCH_mass_d", "PASE", "com_b_nb", "POIMES", "TAIMES", "sex", "FUMER")

# Drop missing confounder data
data<-drop.missing(data=data, check=confounders)
# lost more - note we were using a lot here until the education thing above

##################################################
########### SUBSET ANY DATA IF WE WANT ###########
##################################################

# Do you want to subset to get the more exlusive dataset
if(data_use == "Exclusive"){

	# Also lets drop anyone who is DIAB = 1 or 2
	data<-data[-which(data$DIAB == 1 | data$DIAB == 2),]
	
	# Drop anyone who has BMI out side the range of 22-29.9
	data<-data[-which(data$BMI < 22 | data$BMI > 29.9),]
	
	# Drop individuals on prescribed diets
	data<-data[-which(data$PRESCRIB == 1),]
		
	# Work out who has a variable weight over observations
	weights<-ddply(data, .(sujetno), summarise, mu_mass = mean(POIMES, na.rm=T), sd_mass = sd(POIMES, na.rm=T), n_obs = sum((is.na(POIMES) == F)))
	weights$CV_mass<-weights$sd_mass / weights$mu_mass
	
	# Lets drop anyone who has a CV in mass > 0.04
	use<-weights$sujet[which(weights$CV_mass < 0.04)]
	
	# Now find those individuals that are OK
	data<-data[-which(is.na(match(data$sujet, use)) == T),]

}
dim(data)

# Generate a summary of the dataset in question
means<-as.data.frame(round(cbind(apply(data[,c(confounders, "PROTEIN_kJ.d", "CARBO_kJ.d", "LIPID_kJ.d", "A_Tocopherol", "Vitamin_C", "Trans._F._Acids", traits)], 2, mean, na.rm=T), apply(data[,c(confounders, "PROTEIN_kJ.d", "CARBO_kJ.d", "LIPID_kJ.d", "A_Tocopherol", "Vitamin_C", "Trans._F._Acids", traits)], 2, sd, na.rm=T)), 3))
means<-cbind(rownames(means), means)
names(means)<-c("variable", "mean", "SD")
means<-rbind(data.frame(variable="n_obs", mean=dim(data)[1], SD=NA), data.frame(variable="n_ind", mean=length(unique(data$sujetno)), SD=NA), means)
write.table(means, file="summary_stats.csv", sep=",", row.names=F, col.names=names(means))

##################################################
############### Z-TRANSFORMATIONS ################
##################################################

# Scale the micronutrient variables
for(i in 1:length(micros)){
	data$new<-c(scale(data[,micros[i]]))
	names(data)[dim(data)[2]]<-paste0("z.", micros[i])
}

# I'm going to include z transformed income, age, alcohol, as well as factor(sex) - all additive.
data$z.age<-c(scale(data$age_visit))
data$z.income<-c(scale(data$REVFAM))
data$z.education<-c(scale(data$SCOLAR))
data$z.alcohol<-c(scale(data$ALCH_mass_d))
data$z.PASE<-c(scale(data$PASE))
data$z.comorb<-c(scale(data$com_b_nb)) # I have real concerns about including this - n.comorb could be a result of the diet, and inturn dysregulation - we could be "throwing the baby out with the bathwater"
data$z.mass<-c(scale(data$POIMES))
data$z.height<-c(scale(data$TAIMES))
data$sex<-as.factor(data$sex)
data$smoke<-as.factor(data$FUMER)

# Let's Z-transform them to help interpretation 
for(i in 1:length(traits)){
	data$new<-c(scale(data[, traits[i]]))
	names(data)[dim(data)[2]]<-paste0("z.", traits[i])
}

# Rename - going forwards we will use these outcomes
traits<-paste0("z.", traits)

#######################################################
########### NORMALISING FOR REQUIREMENTS ##############
#######################################################

# Fit a multi-response model for the intake of each macronutrient as a function of age, mass height sex and PASE
MR_intake_model<-gam(list(PROTEIN_kJ.d ~ s(z.age, z.mass, z.height, by = sex) + sex + s(z.PASE), CARBO_kJ.d ~ s(z.age, z.mass, z.height, by = sex) + sex + s(z.PASE), LIPID_kJ.d ~ s(z.age, z.mass, z.height, by = sex) + sex + s(z.PASE)), data=data, family=mvn(d=3), method=fit_with)
summary(MR_intake_model)
save(MR_intake_model, file="MR_intake_model.Rdata")
load("MR_intake_model.Rdata")

# Write the output from the GAM to a csv file
csv.file<-"Intake_model.csv"
write.table(Sys.time(), file=csv.file, sep=",", row.names=F, col.names=F)
write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)

# Write the linear terms
p.table<-round(summary(MR_intake_model)$p.table, 4)
p.table<-as.data.frame(cbind(row.names(p.table), p.table))
names(p.table)[1]<-"Coef."
suppressWarnings(write.table(p.table, file=csv.file, sep=",", row.names=F, col.names=names(p.table), append=T))
write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)

# Write the smooth terms
s.table<-round(summary(MR_intake_model)$s.table, 4)
s.table<-as.data.frame(cbind(row.names(s.table), s.table))
names(s.table)[1]<-"Coef."
suppressWarnings(write.table(s.table, file=csv.file, sep=",", row.names=F, col.names=names(s.table), append=T))	
write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)

# Write the n and deviance explained
dev.expl<-paste0("n = ", summary(MR_intake_model)$n, ": % Dev. Explained = ", round(summary(MR_intake_model)$dev.expl * 100, 2), ": AIC = ", round(AIC(MR_intake_model)))
write.table(dev.expl, file=csv.file, sep=",", row.names=F, col.names=F, append=T)

# Calculate the residuals
Res<-data[,c("PROTEIN_kJ.d", "CARBO_kJ.d", "LIPID_kJ.d")] - predict(MR_intake_model) 

# Calculate each individuals intake relative to their predicted value and add to the dataframe
rel_intake<-as.data.frame(Res / predict(MR_intake_model) * 100)
names(rel_intake)<-c("rel_PROTEIN", "rel_CARBO", "rel_LIPID")
data<-cbind(data, rel_intake)

#######################################################
######## CREATE THE DATA LIST FOR THE TRAITS ##########
#######################################################

# Drop any missing dysregulation scores, from this new dataset add the individual complete datasets in to a list
data_list<-list()
for(i in 1:length(traits)){
	print(traits[i])
	data_list[[i]]<-drop.missing(data=data, check=traits[i])
}

#######################################################
###################### MODEL 1 ########################
#######################################################

# Formulas to fit to each trait - here the same one for all traits
model1<-" ~ s(PROTEIN_kJ.d, CARBO_kJ.d, LIPID_kJ.d) + s(sujetno, bs=\"re\")"
formula_list<-as.list(rep(model1, length(traits)))

# Fit the models and get the results - subject ID is the random effect and will be ignored, but needs to be specified
predict_val<-data.frame(sujetno=0)

# The order for plotting nutrients
XYZ_list<-list()
for(i in 1:length(traits)){
	XYZ_list[[i]]<-c("PROTEIN_kJ.d", "CARBO_kJ.d", "LIPID_kJ.d")
}
# labels for plotting nutrients
labels_list<-list()
for(i in 1:length(traits)){
	labels_list[[i]]<-c("Protein (kJ/day)", "Carbohydrate (kJ/day)", "Lipid (kJ/day)")
}

# Run the my.gam function to get the output
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_1.csv", pdf.file="Model_1.pdf", XYZ_list=XYZ_list, scale_surface=0.8, labels_list=labels_list, method.fit=fit_with)

# Check the model residuals
test.resid(models, pdf.file="Model_1_Resid.pdf", traits=traits)

#######################################################
###################### MODEL 2 ########################
#######################################################

# Formulas to fit to each trait - here the same one for all traits
model2<-" ~ s(rel_PROTEIN, rel_CARBO, rel_LIPID) + s(sujetno, bs=\"re\")"
formula_list<-as.list(rep(model2, length(traits)))

# Fit the models and get the results - subject ID is the random effect and will be ignored, but needs to be specified
predict_val<-data.frame(sujetno=0)

# The order for plotting nutrients
XYZ_list<-list()
for(i in 1:length(traits)){
	XYZ_list[[i]]<-c("rel_PROTEIN", "rel_CARBO", "rel_LIPID")
}

# labels for plotting nutrients
labels_list<-list()
for(i in 1:length(traits)){
	labels_list[[i]]<-c("Protein (% Relative)", "Carbohydrate (% Relative)", "Lipid (% Relative)")
}

# Lets also slice at +/- 20% 
slice_at<-list()
for(i in 1:length(traits)){
	slice_at[[i]]<-c(-20, 0, 20)
}

# Run the my.gam function to get the output
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_2.csv", pdf.file="Model_2.pdf", XYZ_list=XYZ_list, scale_surface=0.8, labels_list=labels_list, slice_at=slice_at, method.fit=fit_with)

# Check the model residuals
test.resid(models, pdf.file="Model_2_Resid.pdf", traits=traits)

#######################################################
###################### MODEL 3 ########################
#######################################################

# Formulas to fit to each trait - here the same one for all traits
model3<-" ~ s(rel_PROTEIN, rel_CARBO, rel_LIPID) + s(z.alcohol) + s(z.income) + s(z.education) + s(z.age) + s(z.PASE) + sex + smoke + s(sujetno, bs=\"re\")"
formula_list<-as.list(rep(model3, length(traits)))

# Fit the models and get the results - note numeric predictors are set as 0, which is the mean, and sex male (1), smoking no (0) (i.e. at te intercept) - subject ID is the random effect and will be ignored, but needs to be specified
predict_val<-data.frame(z.alcohol=0, z.income=0, z.education=0, z.age=0, z.PASE=0, sex=as.factor(1), smoke=as.factor(0), sujetno=0)

# Run the my.gam function to get the output
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_3.csv", pdf.file="Model_3.pdf", XYZ_list=XYZ_list, scale_surface=0.8, labels_list=labels_list, slice_at=slice_at, method.fit=fit_with)

# Check the model residuals
test.resid(models, pdf.file="Model_3_Resid.pdf", traits=traits)

#######################################################
###################### MODEL 4 ########################
#######################################################

# Formulas to fit to each trait - here the same one for all traits - same as model 3 but including comobidities
model4<-paste0(model3, " + s(z.comorb, k=3)")
formula_list<-as.list(rep(model4, length(traits)))

# Fit the models and get the results - note numeric predictors are set as 0, which is the mean, and sex male (1), smoking no (0) (i.e. at te intercept) - subject ID is the random effect and will be ignored, but needs to be specified
predict_val<-data.frame(z.alcohol=0, z.income=0, z.education=0, z.age=0, z.PASE=0, sex=as.factor(1), smoke=as.factor(0), sujetno=0, z.comorb=0)

# Run the my.gam function to get the output - note I aam requesting the SE here 
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_4.csv", pdf.file="Model_4.pdf", XYZ_list=XYZ_list, scale_surface=0.8, labels_list=labels_list, slice_at=slice_at, method.fit=fit_with, include_se = T)

# Check the model residuals
test.resid(models, pdf.file="Model_4_Resid.pdf", traits=traits)

#######################################################
################ MICRONUTRIENT CLUSTERING #############
#######################################################

# Starting with an evalution of the covariance among the micronutrients

# Correlation matrix
cor_mat<-cor(data[,micros])

# Cluster on the correlations them selves - makes defining the cutoff much more simple than going throgh Euclidean distance first
hc<-hclust(as.dist(1-cor_mat))
# Cut the tree at anything above cut_at
cut_at<-0.65
tag<-cutree(hc, h=1-cut_at)

# Which have > 1 within the cluster
multi_group<-sort(unique(tag[which(duplicated(tag) == T)]))

# Create a pallette for each multi-nutrient cluster
my.pallette<-rainbow(length(multi_group))
row.col<-rep("white", length(tag))
for(i in 1:length(multi_group)){
	row.col[which(is.na(match(tag, multi_group[i])) == F)]<-my.pallette[i]
}

# Visualise as a heatmap
pdf("heatmap_micro.pdf", height=10, width=10)
pal<-colorRampPalette(c("Blue", "lightgreen", "Red"))
map<-heatmap.2(cor_mat, trace="none", col=pal(length(seq(-1, 1, 0.05))-1), Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc), keysize=1.5, denscol="black", key.title="", key.ylab="", key.xlab="Correlation (Pearson)", dendrogram="row", cexRow=1.5, cexCol=1.5, key.par=list(mar=c(6,1,2,7), cex.lab=2.25, cex.axis=2), RowSideColors = row.col, lmat = rbind(c(0,0,4),c(3,1,2),c(0,0,5)), lwid = c(2,0.4,7), lhei = c(0.1,5.5,2), margins=c(12,12), key.xtickfun=function(){
	x<-sort(c(seq(-1, 1, 0.5), cut_at))
	return(list(at=(x+1)/2, labels=x))
	})
legend("bottomleft", legend=paste0("Cluster ", seq(1:length(multi_group))), pch=15, col=my.pallette, cex=1.3)
dev.off()

# Get the actual cluster names assigned in the heatmap for use later
cluster.IDs<-match(row.col, my.pallette)
names(cluster.IDs)<-names(tag)

########################################################
################## MICRONUTRIENT PCA ###################
########################################################

# There are some strongly correlated micronutrients - lets do PCA within cluster
write.table(Sys.time(), file="micro_PCA.csv", sep=",", row.names=F, col.names=F)

# object to hold the model - i.e. PC1 where relevant
micro_model<-list()

# For each cluster in the tag do a pca
for(i in 1:max(tag)){
	
	# pull out the cluster and keep the name
	cluster.i<-names(which(tag == i))
	names.i<-cluster.i
	
	# If there are more than 1 nutrients in the cluster do the pca
	if(length(cluster.i) > 1){
		
		# The name of the cluster in the heatmap above
		cluster.name.i<-unique(cluster.IDs[match(names.i, names(cluster.IDs))])
		
		# Add a correlation matrix to the .csv file
		write.table(" ", file="micro_PCA.csv", sep=",", append=T, row.names=F, col.names=F)
		write.table(paste0("Cluster ", cluster.name.i), file="micro_PCA.csv", sep=",", append=T, row.names=F, col.names=F)
		
		# Do the PCA
		pca<-prcomp(data[,cluster.i])
		
		write.table("", file="micro_PCA.csv", sep=",", append=T, row.names=F, col.names=F)
		suppressWarnings(write.table(round(summary(pca)$importance, 2), file="micro_PCA.csv", sep=",", row.names=row.names(summary(pca)$importance), col.names=F, append=T))
		
		# Add the PC1 to the dataframe name the column
		names.i<-paste0("Cluster_", cluster.name.i, "_PC1")
		data$pc<-c(scale(pca$x[,1]))
		write.table(" ", file="micro_PCA.csv", sep=",", append=T, row.names=F, col.names=F)
		suppressWarnings(write.table(round(cor(data[,c(cluster.i, "pc")]), 2), file="micro_PCA.csv", sep=",", row.names=F, col.names=c(cluster.i, "PC1"), append=T))
		names(data)[dim(data)[2]]<-paste0("z.", names.i)
		
			
	}
	
	# Save the name for modelling
	micro_model[[i]]<-names.i
}

# Get the new list of micro nutrients that includes the PC1 for clusters
micro_model<-unlist(micro_model)

#######################################################
### RECREATE THE DATA LIST FOR THE TRAITS WITH PC1s ###
#######################################################

# Drop any missing dysregulation scores, from this new dataset add the individual complete datasets in to a list
data_list<-list()
for(i in 1:length(traits)){
	print(traits[i])
	data_list[[i]]<-drop.missing(data=data, check=traits[i])
}

############################################################
########### 3-WAY MICRONUTRIENT COMBINATIONS ###############
############################################################

# Check we are doing the inclusive dataset
if(data_use == "Inclusive" & impute_income == FALSE){
		
	# Combinations of smoothers
	micro_cmbn<-combn(paste0("z.", micro_model), 3)
	smoothers<-paste0("s(", apply(micro_cmbn, 2, paste, collapse=", "), ")")
		
	# Some tables to hold results
	results_dev<-as.data.frame(array(0, c(length(smoothers), length(traits)+1)))
	names(results_dev)<-c("smoother", traits)
	results_dev$smoother<-smoothers
	results_p<-results_dev
	results_F<-results_dev
	
	save_at<-100
	check_save<-round(seq(0, length(smoothers)+save_at, save_at))
	count<-2
	
	# Loop for the gam to do the micronutrient smoothers
	for(j in 1:length(smoothers)){
		
		# Drop any missing dysregulation scores, from this new dataset add the individual complete datasets in to a list
		data_list<-list()
		for(i in 1:length(traits)){
			data_list[[i]]<-drop.missing(data=data, check=traits[i], verbose=F)
		}
		
		# Print the progress - stops you going insane
		print(paste0("smoother ", j, " / ", length(smoothers)))
		#Sys.sleep(2)
		
		formula_fit<-paste0(" ~ ", smoothers[j], " + s(sujetno, bs=\"re\")")
		formula_list<-as.list(rep(formula_fit, length(traits)))
		
		# Run the models - note I am not plotting or getting any output - just checking the difference in deviance
		models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, method.fit=fit_with)
		
		# Pull out the deviance for micro models
		deviance_micro<-list()
		p_micro<-list()
		f_micro<-list()
		for(i in 1:length(traits)){
			deviance_micro[[i]]<-summary(models[[i]])$dev.expl * 100
			p_micro[[i]]<-summary(models[[i]])$s.table[1,4]
			f_micro[[i]]<-summary(models[[i]])$s.table[1,3]
		}
		
		# Add results to the tables
		results_dev[j,-1]<-unlist(deviance_micro)
		results_p[j,-1]<-unlist(p_micro)
		results_F[j,-1]<-unlist(f_micro)
		
		if(j == check_save[count]){
			write.table(results_dev, file="Deviance.csv", sep=",", row.names=F, col.names=names(results_dev))
			write.table(results_p, file="pval.csv", sep=",", row.names=F, col.names=names(results_p))
			write.table(results_F, file="Fval.csv", sep=",", row.names=F, col.names=names(results_F))	
			count<-count+1
		}
		
	}
	
	write.table(results_dev, file="Deviance.csv", sep=",", row.names=F, col.names=names(results_dev))
	write.table(results_p, file="pval.csv", sep=",", row.names=F, col.names=names(results_p))
	write.table(results_F, file="Fval.csv", sep=",", row.names=F, col.names=names(results_F))
	
	# Read in p values and deviance explained
	results_p<-read.csv("pval.csv")
	results_dev<-read.csv("Deviance.csv")
	
	# Frequency histograms of p-values
	pdf("Histogram_p.pdf", height=10, width=20)
	par(mfrow=c(2,4), cex.axis=1.5, cex.lab=1.5, mar=c(6,6,5.5,1))
	
	panel<-c("A)", "B)", "C)", "D)", "E)", "F)", "G)", "H)")
	
	for(i in 1:length(traits)){
		hi<-hist(results_p[,i+1], main="", xlab="", ylab="", pch=16, cex=1.25, breaks=20, col="lightgrey", border="white", cex.axis=2, xlim=c(0,1))
		mtext(trait_label[i], cex=2.5, font=2, line=2.2)
		mtext(panel[i], cex=2.5, at=-0.1, line=2.2)
		
		# Expected number of false values 
		exp_T1<-dim(results_p)[1]*0.05
		abline(h=exp_T1, col=2, lwd=1.5)
		abline(v=0.2, col=4, lwd=1.5)
		mtext("p", side=1, line=3, cex=2)
		mtext("Count", side=2, line=3.2, cex=2)
		
		# Find the counts within the upper quadrant
		discrep_i<-hi$counts[which(hi$breaks < 0.2)]-exp_T1
		n_in_i<-sum(discrep_i * (discrep_i > 0))
		ith_percent<-round(n_in_i / dim(results_p)[1] * 100, 0)
		mtext(paste0("", ith_percent, "%"), at=0.1, line=-1, cex=2)
	}
	
	dev.off()
	
	
	# Get p values adjusted for fdr and ask if they are significant
	adjusted<-apply(results_p[,-1], 2, p.adjust, method="fdr")
	for(i in 1:length(traits)){
		adjusted[,i]<-(adjusted[,i] < (0.05))
	}
	
	# Now only do this next 'discovery step' for the inclusive dataset - for the exclusive dataset we will select the same 'target' combination identified in the inclusive analysis and test it alone

	# Lets draw a venn diagram for the overlap amongst systems
	
	# Overall areas for the scores
	areas<-apply(adjusted, 2, sum)
	
	# Find those with > 0 hits
	venn.plot<-names(areas)[which(areas > 0)]
		
	list_smoothers<-list()
	for(k in 1:length(venn.plot)){
		list_smoothers[[k]]<-results_dev[which(adjusted[, venn.plot[k]] == 1),1]
		#names(list_smoothers)[k]<-venn.plot[k]
	}
	
	# Work out which 
	venn_micro<-Venn(list_smoothers, SetNames=trait_label[match(venn.plot, traits)])
	print(venn_micro)
	
	# Draw the venn diagram
	pdf("Venn.signif.pdf")
	plot(venn_micro, type="battle")
	dev.off()
	
	# Save the venn table - might be a better way of displaying
	write.table(venn_micro@IndicatorWeight, file="Venn_table.csv", sep=",", row.names=F, col.names=colnames(venn_micro@IndicatorWeight))
	
	# Which is that model that appears for all 4?
	cross_system_significants<-apply(adjusted, 1, sum)
	across_4<-which(cross_system_significants == max(cross_system_significants))
	
	# Find the smoother with highest deviance explained across all traits
	mean_dev<-apply(results_dev[across_4,-1], 1, mean)
	
	# Print the significant smoothers
	print(cbind(results_dev[across_4,1], mean_dev))
	
	micro_smooth<-results_dev[across_4[which(mean_dev == max(mean_dev))],1]
	
}else{
	
	# Otherwise use the micro smnoother selected by the other analysis
	micro_smooth<-"s(z.A_Tocopherol, z.Vitamin_C, z.Trans._F._Acids)"	
}

# What do the predictions for these micronutrients look like - need to refit the model to find out

# Formulas to fit to each trait - here the same one for all traits
model_micro_k<-paste0(" ~ ", micro_smooth, " + s(sujetno, bs=\"re\")")
formula_list<-as.list(rep(model_micro_k, length(traits)))

# Fit the models and get the results - subject ID is the random effect and will be ignored, but needs to be specified
predict_val<-data.frame(sujetno=0)

# The order for plotting nutrients
XYZ_list_micro<-list()
for(i in 1:length(traits)){
	XYZ_list_micro[[i]]<-c("z.A_Tocopherol", "z.Vitamin_C", "z.Trans._F._Acids")
}

# labels for plotting nutrients
labels_list_micro<-list()
for(i in 1:length(traits)){
	labels_list_micro[[i]]<-c("A Tocopherol Intake (Z)", "Vitamin C Intake (Z)", "Trans. F. Acids (Z)")
}

# Lets also slice at +/- 1/2SD 
slice_at_micro<-list()
for(i in 1:length(traits)){
	slice_at_micro[[i]]<-c(-0.5, 0, 0.5)
}

# Run the my.gam function to get the output
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_micro.csv", pdf.file="Model_micro.pdf", XYZ_list=XYZ_list_micro, scale_surface=0.8, labels_list=labels_list_micro, slice_at=slice_at_micro, method.fit=fit_with)

# Check the model residuals
test.resid(models, pdf.file="Model_micro_Resid.pdf", traits=traits)

#######################################################
###################### MODEL 5 ########################
#######################################################

# Lets fit models for that smoother with correction for confounders

# Add the micronutrient smoother in to model 3 and create the list
model5<-paste0(" ~ ", micro_smooth, " + s(z.alcohol) + s(z.income) + s(z.education) + s(z.age) + s(z.PASE) + sex + smoke + s(sujetno, bs=\"re\")")
formula_list<-as.list(rep(model5, length(traits)))

# For the predicted values add in 0's for the micronutrients 
predict_val<-data.frame(z.alcohol=0, z.income=0, z.education=0, z.age=0, z.PASE=0, sex=as.factor(1), smoke=as.factor(0), sujetno=0)
m<-as.data.frame(array(0, c(1, length(micro_model))))
names(m)<-paste0("z.", micro_model)
predict_val<-cbind(predict_val, m)

# Run the models and obtain the output with the additional 3D micronutrient term
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_5.csv", pdf.file="Model_5.pdf", XYZ_list=XYZ_list_micro, scale_surface=0.8, labels_list=labels_list_micro, slice_at=slice_at_micro, method.fit=fit_with)

# Test the residuals
test.resid(models, pdf.file="Model_5_Resid.pdf", traits=traits)

#######################################################
###################### MODEL 6 ########################
#######################################################

# Lets fit models for that smoother with correction for confounders and n comobidities

# Add the micronutrient smoother in to model 3 and create the list
model6<-paste0(model5, " + s(z.comorb, k=3)")
formula_list<-as.list(rep(model6, length(traits)))

# For the predicted values add in 0's for the micronutrients 
predict_val<-data.frame(z.alcohol=0, z.income=0, z.education=0, z.age=0, z.PASE=0, sex=as.factor(1), smoke=as.factor(0), z.comorb=0, sujetno=0)
m<-as.data.frame(array(0, c(1, length(micro_model))))
names(m)<-paste0("z.", micro_model)
predict_val<-cbind(predict_val, m)

# Run the models and obtain the output with the additional 3D micronutrient term
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_6.csv", pdf.file="Model_6.pdf", XYZ_list=XYZ_list_micro, scale_surface=0.8, labels_list=labels_list_micro, slice_at=slice_at_micro, method.fit=fit_with)

# Test the residuals
test.resid(models, pdf.file="Model_6_Resid.pdf", traits=traits)

#######################################################
###### MIXTURE MODELS AND RMTS FOR MICRO INTAKE #######
#######################################################

# Steve suggests that we check how diet composotion in terms of macronutrients associated with micronutrient intake

# Create a datalist for the micronutrients
micro_data_list<-list()
for(i in 1:length(micro_model)){
	micro_data_list[[i]]<-drop.missing(data=data, check=paste0("z.", micro_model[i]))	
}

# Run the mixture model function 
mixture.models(data_list=micro_data_list, traits=paste0("z.", micro_model), titles=micro_model, p_P="PROTEIN_pct", p_C="CARBO_pct", p_F="LIPID_pct", pdf.file="micro_RMT.pdf")

#######################################################
###################### MODEL 7 ########################
#######################################################

# Add the micronutrient smoother in to model 3 so now we have macronutrients and micronutrients in the same model and create the list
model7<-paste0(model3, " + ", micro_smooth)
formula_list<-as.list(rep(model7, length(traits)))

# For the predicted values add in 0's for the micronutrients 
predict_val<-data.frame(z.alcohol=0, z.income=0, z.education=0, z.age=0, z.PASE=0, sex=as.factor(1), smoke=as.factor(0), sujetno=0)
m<-as.data.frame(array(0, c(1, length(micro_model))))
names(m)<-paste0("z.", micro_model)
predict_val<-cbind(predict_val, m)

# Run the models and obtain the output with the additional 3D micronutrient term
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_7.csv", pdf.file="Model_7_Macro.pdf", XYZ_list=XYZ_list, scale_surface=0.8, labels_list=labels_list, slice_at=slice_at, method.fit=fit_with)

# Test the residuals
test.resid(models, pdf.file="Model_7_Resid.pdf", traits=traits)

# Lets plot the effects of the micronutrients

# We need to redefine the predictors list to include the macro nutrients - I will set as 0, which is average for someone you age sex weight and PASE level
predict_val<-data.frame(rel_PROTEIN=0, rel_CARBO=0, rel_LIPID=0, z.alcohol=0, z.income=0, z.education=0, z.age=0, z.PASE=0, sex=as.factor(1), smoke=as.factor(0), sujetno=0)

# Re run the my.gam function to get the surfaces - note the function actually re-runs the models - some redundancy here but what the hey
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", pdf.file="Model_7_Micro.pdf", XYZ_list=XYZ_list_micro, labels_list=labels_list_micro, slice_at=slice_at_micro, scale_surface=0.8, method.fit=fit_with)

#######################################################
###################### MODEL 8 ########################
#######################################################

# Add the micronutrient smoother and coomorbidities in to model 3 so now we have macronutrients and micronutrients in the same model and create the list
model8<-paste0(model7, " + s(z.comorb, k=3)")
formula_list<-as.list(rep(model8, length(traits)))

# For the predicted values add in 0's for the micronutrients 
predict_val<-data.frame(z.alcohol=0, z.income=0, z.education=0, z.age=0, z.PASE=0, sex=as.factor(1), smoke=as.factor(0), sujetno=0, z.comorb=0)
m<-as.data.frame(array(0, c(1, length(micro_model))))
names(m)<-paste0("z.", micro_model)
predict_val<-cbind(predict_val, m)

# Run the models and obtain the output with the additional 3D micronutrient term
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", csv.file="Model_8.csv", pdf.file="Model_8_Macro.pdf", XYZ_list=XYZ_list, scale_surface=0.8, labels_list=labels_list, slice_at=slice_at, method.fit=fit_with)

# Test the residuals
test.resid(models, pdf.file="Model_8_Resid.pdf", traits=traits)

# We also need to redefine the predictors list to include the macro nutrients - I will set as 0, which is average for someone you age sex weight and PASE level
predict_val<-data.frame(rel_PROTEIN=0, rel_CARBO=0, rel_LIPID=0, z.alcohol=0, z.income=0, z.education=0, z.age=0, z.PASE=0, sex=as.factor(1), smoke=as.factor(0), sujetno=0, z.comorb=0)

# Re run the my.gam function to get the surfaces - note the function actually re-runs the models - some redundancy here but what the hey
models<-my.gam(data_list=data_list, traits=traits, formula_list=formula_list, predict_val=predict_val, exclude="s(sujetno)", pdf.file="Model_8_Micro.pdf", XYZ_list=XYZ_list_micro, slice_at=slice_at_micro, scale_surface=0.8, method.fit=fit_with)
