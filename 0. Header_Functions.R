
# FUnction to work out whether points in nutrient spaace fall within the hull given by the raw data

inhull <- function(testpts, calpts, hull=convhulln(calpts), tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)) { 

	# https://tolstoy.newcastle.edu.au/R/e8/help/09/12/8784.html
	calpts <- as.matrix(calpts) 
	testpts <- as.matrix(testpts) 
	p <- dim(calpts)[2] 
	cx <- dim(testpts)[1] # rows in testpts
	nt <- dim(hull)[1] # number of simplexes in hull 
	nrmls <- matrix(NA, nt, p)
	
	degenflag <- matrix(TRUE, nt, 1) 
	for (i in 1:nt){ 
		nullsp<-t(Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))
		if (dim(nullsp)[1] == 1){
			nrmls[i,]<-nullsp
			degenflag[i]<-FALSE
		}
	}

	if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")
	nrmls <- nrmls[!degenflag,] 
	nt <- dim(nrmls)[1] 
	
	center = apply(calpts, 2, mean) 
	a<-calpts[hull[!degenflag,1],] 
	nrmls<-nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)
	
	dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
	nrmls <- nrmls*matrix(dp, nt, p)
	
	aN <- diag(a %*% t(nrmls)) 
	val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min) 
	
	val[abs(val) < tol] <- 0 
	as.integer(sign(val)) 
}

# A function created to find the outer perimeter over which the surface should be fitted for proportional data
findConvex.prop<-function(x,y,rgnames,res=101){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	x.new<-seq(0,1,len=res)
	y.new<-seq(0,1,len=res)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}

# A function to run the compositional analyses
# Arguments
# csv.file = csv file to write results to
# pdf.file = pdf file to write figures to
# traits = a character vector of the traits we are interested in
# titles = a character vector of titles for the traits we are interested in
# data = the dataset for plotting
# pX = propotion of energy coming from nutrient X

mixture.models<-function(data_list, traits, titles=traits, p_P="p_P", p_C="p_C", p_F="p_F", csv.file="MM_results.csv", pdf.file="RMT.pdf"){
			
	# Set the resolution of the surface
	surface.resolution<-501
	
	# How many values to round surface
	round.surf<-3
	
	# This specifies the color scheme for surface - it is actually a function that returns a function
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	
	# How many different colours should we use on the plot
	no.cols<-256
	
	# Get the colors to use from the pallette specified above
	map<-rgb.palette(no.cols)
	
	# How many levels should there be on the surface
	nlev<-3
		
	# Labels for each
	labels<-c("Protein (%)", "Carbohydrate (%)", "Fat (%)")
	
	# results file for the jth cutoff
	res.file<-csv.file
		
	## plot the RMT surface
	iso.lines<-seq(1, 0, -0.2)
	
	# List to hold the models
	AIC_models<-list()
	
	# Open the pdf file for plotting
	pdf(pdf.file, height=5, width=5)

	# Set the layout
	par(mfrow=c(1,1), mar=c(5,5,5,1))	
	
	for(k in 1:length(traits)){
		
		# Find the kth dataset
		data<-data_list[[k]]
		
		# Make sure the proportions are closed off to 1
		data$p_P<-data[,p_P] / (data[, p_P] + data[, p_C] + data[, p_F])
		data$p_C<-data[,p_C] / (data[, p_P] + data[, p_C] + data[, p_F])
		data$p_F<-data[,p_F] / (data[, p_P] + data[, p_C] + data[, p_F])
		
		## estimate convex hull and predict
		mdff2<-findConvex.prop(data$p_P, data$p_C, c("p_P","p_C"), surface.resolution)
		mdff2$p_F<-with(mdff2, 1 - p_P - p_C)

		# Create table for results
		if(k == 1){write.table("", file=res.file, sep=",", row.names=F, col.names=F)}
		
		# Find the right outcome
		data$this.outcome<-data[,traits[k]]
		
		# Variable to hold the four models of he scheffe's polynomials
		mmods<-list()
		
		# Fit the intercept model
		mmods[[1]]<-lm(this.outcome ~ 1, data=data)
		
		# Fit Scheffes polynomials
		for(i in 1:4){
			model<-i
			mmods[[i+1]]<-MixModel(frame=data, response="this.outcome", mixcomps=c("p_P","p_C","p_F"), model=model)
		}
		
		# Find minimal model based on AIC
		AICs<-unlist(lapply(mmods, AIC))
		deltas<-AICs - min(AICs)
		options<-which(deltas <= 2)
		min.model<-min(options)
		model.AIC<-mmods[[min.model]]
		AIC_models[[k]]<-model.AIC
		
		# Write the results to the table
		write.table("", file=res.file, sep=",", row.names=F, col.names=F, append=T)
		write.table(traits[k], file=res.file, sep=",", row.names=F, col.names=F, append=T)
		write.table(cbind(seq(1, 5, 1), round(AICs, 2)), file=res.file, sep=",", row.names=F, col.names=c("Model", "AIC"), append=T)
		write.table("", file=res.file, sep=",", row.names=F, col.names=F, append=T)
		write.table(paste("Model ", min.model, " favoured by AIC."), file=res.file, sep=",", row.names=F, col.names=F, append=T)	
		write.table("", file=res.file, sep=",", row.names=F, col.names=F, append=T)
		res.k<-as.data.frame(round(summary(model.AIC)$coef[,c(1:3)], 4))
		res.k$df<-(dim(data)[1]) - (dim(res.k)[1])
		res.k$p<-(round(summary(model.AIC)$coef[,4], 4))
		res.k<-cbind(row.names(res.k), res.k)
		colnames(res.k)[1]<-"Coef."
		write.table(res.k, file=res.file, sep=",", row.names=F, col.names=colnames(res.k), append=T)
		
		# Get the predicted surface
		mdff2$fit<-predict(model.AIC, newdata=mdff2)
		NAs<-which(apply((is.na(mdff2[,c(1:3)]) == F), 1, prod) == 0)
		mdff2$fit[NAs]<-NA
		surf<-matrix(mdff2$fit, nrow=sqrt(dim(mdff2)[1]))
		surf<-round(surf, round.surf)
		
		# Find minimal and maximal values so as to scale sensibly
		mn<-min(surf, na.rm=TRUE)
		mx<-max(surf, na.rm=TRUE)
		null<-0
		if(mn == mx){
			null<-1
			mn<-mn - mn*0.025
			mx<-mx + mx*0.025
		}
		locs<-(range(surf, na.rm=TRUE) - mn) / (mx-mn) * no.cols
		
		# Actually plots the surface using all of the above info above
		plot(-10, -10, bty="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xlab="", ylab="", xaxt="n", yaxt="n")
		
		# Adds some axes
		axis(1, at = seq(0, 1, 0.2), labels=seq(0, 1, 0.2)*100)
		axis(2, at = seq(0, 1, 0.2), labels=seq(0, 1, 0.2)*100)
			
		# Add the Isolines
		for(i in 1:length(iso.lines)){
			abline(a = iso.lines[i], b=-1)
		}
		
		# Add the surface
		image(seq(0, 1, length.out = surface.resolution), seq(0, 1, length.out = surface.resolution), surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE, main="", add=T)
		# Adds a contour over the top (can add a title using main)
		if(sd(surf, na.rm=T) > 0){
			contour(seq(0, 1, length.out = surface.resolution), seq(0, 1, length.out = surface.resolution), surf, add=TRUE, levels=pretty(range(mn,mx), nlev), labcex=0.8)
		}
		
		# Add the axes labels
		mtext(labels[1], side=1, line=2.25)
		mtext(labels[2], side=2, line=2.25)
		text(0.55, 0.55, labels[3], srt=-45, cex=1.1)
		mtext(titles[k], cex=1.5, line=2)
			
		if(null == 1){
			text(mean(data$p_P), mean(data$p_C), round(mean(mdff2$fit, na.rm=T), 2), cex=0.85, srt=-75)
		}
	
	}

	# Close the plotting file
	dev.off()

	# Spit out the AIC favoured models 
	return(AIC_models)
	
}


# A function to fit gams for intake for a suite of traits, produce model output ans generate surfaces

my.gam<-function(data_list, traits, formula_list, XYZ_list=NA, predict_val=NA, exclude=NULL, csv.file=NA, pdf.file=NA, slice_at=NA, fit.resolution=101, no.cols=256, nlev=8, include_se=F, markers=NA, scale_surface=NA, cex.axis=2, labels_list=XYZ_list, cex.lab=2, direction=1, method.fit="GCV.Cp"){
	
	# If we want to plot the surfaces
	if(is.na(pdf.file) == F){
		# Set the layout
		# Open the pdf file for plotting
		if(include_se == T){
			# DO you want surfaces for SE
			if(direction == 1){
				# Set the layout for direction 1 - reading left to right: lower, middle, upper slice
				pdf(pdf.file, height=6 * 5, width=3 * 5)
				par(mfrow=c(6, 3), mar=c(6,6,5,1))
			}else{
				# Set the layout for direction 2 - reading top to bottom: lower, middle, upper slice
				pdf(pdf.file, height=3 * 5, width=6 * 5)
				par(mar=c(6,6,5,1))
				layout(as.matrix(array(seq(1, 6*3, 1), c(3, 6))))
			}	
		}else{
			pdf(pdf.file, height=3 * 5, width=3 * 5)
			if(direction == 1){
				par(mfrow=c(3, 3), mar=c(6,6,5,1))
			}else{
				par(mar=c(6,6,5,1))
				layout(as.matrix(array(seq(1, 3*3, 1), c(3, 3))))
			}		
		}
	}
	
	# Order the markers
	markers_list<-list()
	if(is.list(markers) == T){
		markers_list[[1]]<-list(markers[[1]], markers[[2]])
		markers_list[[2]]<-list(markers[[1]], markers[[3]])
		markers_list[[3]]<-list(markers[[2]], markers[[3]])
	}
	
	# This specifies the color scheme for surface
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	map<-rgb.palette(no.cols)
	
	# create the csv file to write to, if you want one
	if(is.na(csv.file) == F){
		write.table(Sys.time(), file=csv.file, sep=",", row.names=F, col.names=F)
		write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
	}
		
	# List to hold the fitted models
	models.list<-list()
	
	# Formatting the progress bar
	 pb <- txtProgressBar(min = 0, max = length(traits), style = 3)
	 progress<-0

	# Loop for the proteins
	for(k in 1:length(traits)){
		
		# Pull out the kth dataset
		data<-data_list[[k]]
		
		# write the trait to the table
		if(is.na(csv.file) == F){
			write.table(traits[k], file=csv.file, sep=",", row.names=F, col.names=F, append=T)
		}
				
		# Get the kth trait as the outcome
		data$outcome<-data[,traits[k]]
		
		# Format the formula 
		formula_fit<-as.formula(paste0("outcome ", formula_list[[k]]))
		
		# Fit the model	
		GAM<-gam(formula_fit, data=data, method=method.fit)
		
		# Save the GAM
		models.list[[k]]<-GAM
		names(models.list)[k]<-traits[k]
		
		# Writing output
		if(is.na(csv.file) == F){
			# Write the linear terms
			p.table<-round(summary(GAM)$p.table, 4)
			p.table<-as.data.frame(cbind(row.names(p.table), p.table))
			names(p.table)[1]<-"Coef."
			suppressWarnings(write.table(p.table, file=csv.file, sep=",", row.names=F, col.names=names(p.table), append=T))
			write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
			
			# Write the smooth terms
			s.table<-round(summary(GAM)$s.table, 4)
			s.table<-as.data.frame(cbind(row.names(s.table), s.table))
			names(s.table)[1]<-"Coef."
			suppressWarnings(write.table(s.table, file=csv.file, sep=",", row.names=F, col.names=names(s.table), append=T))	
			write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
			
			# Write the n and deviance explained
			dev.expl<-paste0("n = ", summary(GAM)$n, ": % Dev. Explained = ", round(summary(GAM)$dev.expl * 100, 2), ": AIC = ", round(AIC(GAM)))
			write.table(dev.expl, file=csv.file, sep=",", row.names=F, col.names=F, append=T)		
			write.table(" ", file=csv.file, sep=",", row.names=F, col.names=F, append=T)
		}
		
		# If we are plotting the surface make the surfaces
		if(is.na(pdf.file) == F){
			
			# List for the order of plots
			list.order<-list()
			list.order[[1]]<-XYZ_list[[k]][c(1,2,3)]
			list.order[[2]]<-XYZ_list[[k]][c(1,3,2)]
			list.order[[3]]<-XYZ_list[[k]][c(2,3,1)]
			
			# List for the labels
			labels.order<-list()
			labels.order[[1]]<-labels_list[[k]][c(1,2,3)]
			labels.order[[2]]<-labels_list[[k]][c(1,3,2)]
			labels.order[[3]]<-labels_list[[k]][c(2,3,1)]
			
			# Lists to hold the predictions for the nth combinations
			predictors.list<-list()
			predictions.list<-list()
			xyz.list<-list() 
			cv.list<-list()
			
			# Go through the nutrient orders and get the predicted values, note we will do the predictions for all three combinations, then find the minimla and maximal values, then go back though the nutrients orders and plot out
			for(n in 1:3){
			
				# Order to plot the nutrients
				nutrient.order<-list.order[[n]]
				
				# Values to predict over
				x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
				y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
				
				# If we do not specify values to slice at, use the 25, 50, and 75 %ile
				if(is.list(slice_at) == F){
					z.vals<-round(quantile(data[,nutrient.order[3]])[c(2:4)])
				}else{
					z.vals<-slice_at[[k]]	
				}
				
				# Fitted list to hold some results for later
				x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=fit.resolution)
				y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=fit.resolution)
				z.new<-z.vals
				predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
				names(predictors)<-nutrient.order
				in.poly<-as.numeric(inhull(predictors[,c(1:3)], data[,names(predictors)]) != -1)
				
				# Add the predictors for the additional 'confounders'
				predictors<-cbind(predictors, predict_val)
				
				# Do the predictions
				predictions<-predict(GAM, newdata=predictors, type="response", exclude=exclude, se.fit=T)
				
				# Edit out based on the marker list
				predictions$fit[which(in.poly == 0)]<-NA
				predictions$se.fit[which(in.poly == 0)]<-NA
				
				# Save the nth set of predictions
				predictions.list[[n]]<-predictions$fit
				predictors.list[[n]]<-predictors
				xyz.list[[n]]<-list(x.new, y.new, z.new)
				cv.list[[n]]<-predictions$se.fit
			}
			
			# Find the min and max values across all predictions
			mn<-min(unlist(predictions.list), na.rm=T)
			mx<-max(unlist(predictions.list), na.rm=T)
			
			# If no color scale is specified for the outcomes scale by the predicted values
			if(sum(is.na(scale_surface)) < length(scale_surface)){
				# Find the absolute max values across all predictions, and the scale
				upp_abs<-max(abs(c(scale_surface, mn, mx)))	
				mn<-(-upp_abs)
				mx<-upp_abs
			}
			
			# now do the coefficient of the error
			mn.cv<-min(unlist(cv.list), na.rm=T)
			mx.cv<-max(unlist(cv.list), na.rm=T)

			# Now go back though the predictions and plot
			for(n in 1:3){
				
				# Order to plot the nutrients
				nutrient.order<-list.order[[n]]
				labs<-labels.order[[n]]
				
				# Pull out the nth set of predictors and predictions
				predictors<-predictors.list[[n]]
				predictions<-predictions.list[[n]]
				x.new<-xyz.list[[n]][[1]]
				y.new<-xyz.list[[n]][[2]]
				z.new<-xyz.list[[n]][[3]]
					
				# Do the 3 quantiles for the predictions
				for(i in 1:length(z.new)){
								
					# Subset for the ith quantile
					ith_Quantile<-predictions[which(predictors[, nutrient.order[3]] == z.new[i])]
											
					surf<-matrix(ith_Quantile, nrow=fit.resolution)
								
					locs<-round((range(surf, na.rm=TRUE) - mn) / (mx-mn) * no.cols)
					image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE)
					mtext(paste0(labs[3], " = ", z.new[i]), line=1, cex=cex.lab)
					mtext(labs[1], side=1, line=4, cex=cex.lab)
					mtext(labs[2], side=2, line=4, cex=cex.lab)
					if(i == 3 & n == 1){
						mtext(traits[k], line=2, font=1, cex=1, at = max(x.new)*0.9)
					}
					axis(1, cex.axis=cex.axis)
					axis(2, cex.axis=cex.axis)
					contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn, mx), nlev), labcex=1, lwd=3)			
					# Add any markers
					if(is.list(markers) == T){
						abline(v=markers_list[[n]][[1]], col="grey")
						abline(h=markers_list[[n]][[2]], col="grey")
					}
								
				}
				
				# Now the 3 quantiles for the errors if we want those
				if(include_se == T){
					
					# Pull out the nth set of errors
					predictions<-cv.list[[n]]
					
					# GO through each quantile
					for(i in 1:length(z.new)){
									
						# Subset for the ith quantile
						ith_Quantile<-predictions[which(predictors[, nutrient.order[3]] == z.new[i])]
												
						surf<-matrix(ith_Quantile, nrow=fit.resolution)
									
						locs<-round((range(surf, na.rm=TRUE) - mn.cv) / (mx.cv-mn.cv) * no.cols)
						image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE)
						mtext(paste0("se"), line=1, cex=cex.lab)
						mtext(labs[1], side=1, line=4, cex=cex.lab)
						mtext(labs[2], side=2, line=4, cex=cex.lab)
						axis(1, cex.axis=cex.axis)
						axis(2, cex.axis=cex.axis)
						contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn.cv, mx.cv), nlev), labcex=1, lwd=3)
									
					}
				}
			}
		}		
		
		# Update the progress
		progress<-progress + 1
		setTxtProgressBar(pb, progress)
	}
	
	if(is.na(pdf.file) == F){
		# Close the file 
		dev.off()
	}
	
	# Return the models listed
	return(models.list)

}

# Function to drop all missing data for a set of variables from a predictor
drop.missing<-function(data, check, verbose=T){
	
	if(verbose == T){
		print("n.rows before:")
		print(dim(data)[1])
	}
		
	for(i in 1:length(check)){
		missing<-which(is.na(data[,check[i]]) == T)
		if(length(missing) > 0){
			data<-data[-missing,]
		}
	}
	
	if(verbose == T){
		print("n.rows after:")
		print(dim(data)[1])
	}
	
	return(data)
}

# Function to test residuals for a list of models via GAM
test.resid<-function(models, pdf.file, traits){
	
	# Open the file for plotting
	pdf(pdf.file)
	par(mfrow=c(1,1))
	
	# For each model in the list		
	for(i in 1:length(models)){
		
		# Get the predictions and pearson residuals	
		x<-predict(models[[i]])
		y<-resid(models[[i]], type="pearson")
		
		# Plot them
		plot(x, y, xlab="Predicted Values", ylab="Pearson Residuals", main=traits[i])
		abline(h=0, col="gray")
		
		# Fit a GAM to check
		if(sd(x) > 0){
			test<-gam(y ~ s(x))
			new.x<-seq(min(x), max(x), len=1000)
			new.y<-predict(test, newdata=data.frame(x=new.x), se=T)
			lines(new.x, new.y$fit, col=2, lwd=2)
			lines(new.x, new.y$fit + new.y$se.fit*1.96, col=2, lwd=2, lty=2)
			lines(new.x, new.y$fit - new.y$se.fit*1.96, col=2, lwd=2, lty=2)
			mtext(paste0("p = ", as.character(round(summary(test)$s.table), 3)[4]))	
		}	
	}
	
	# Close the file
	dev.off()

}

# Function to plot models in 3D
surface_3d<-function(data_list, traits, models, XYZ_list, predict_val, slice_at=NA, fit.resolution=101, no.cols=256, folder_name="3d_interactive", scale_surface=NA, exclude=NULL){
	
	# Load plotly
	require(plotly)
	
	# Formatting the progress bar
	pb <- txtProgressBar(min = 0, max = length(traits), style = 3)
	progress<-0
	
	# Get the directory
	directory<-getwd()
	
	# Create the new directory
	if(file.exists(folder_name)){
		unlink(folder_name, recursive=T)
	}
	dir.create(folder_name)
	new_dir<-paste0(directory, "/", folder_name)
	setwd(new_dir)
	
	# This specifies the color scheme for surface - the usual colours for GF surfaces
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	
	# Get the colors to use from the pallette specified above
	map<-rgb.palette(no.cols)

	# Loop for the traits
	for(k in 1:length(traits)){
		
		# Pull out the kth dataset
		data<-data_list[[k]]
				
		# Find the model
		GAM<-models[[k]]

		# Order to plot the nutrients
		nutrient.order<-XYZ_list[[k]]
		
		# Values to predict over
		x.limits<-c(floor(min(data[,nutrient.order[1]])), ceiling(max(data[,nutrient.order[1]])))
		y.limits<-c(floor(min(data[,nutrient.order[2]])), ceiling(max(data[,nutrient.order[2]])))
		
		# If we do not specify values to slice at, use the 25, 50, and 75 %ile
		if(is.list(slice_at) == F){
			z.vals<-round(quantile(data[,nutrient.order[3]])[c(2:4)])
		}else{
			z.vals<-slice_at	[[k]]
		}
		
		# Fitted list to hold some results for later
		x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=fit.resolution)
		y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=fit.resolution)
		z.new<-z.vals
		predictors<-as.data.frame(expand.grid(x.new, y.new, z.new))
		names(predictors)<-nutrient.order
		in.poly<-as.numeric(inhull(predictors[,c(1:3)], data[,names(predictors)]) != -1)
		
		# Add the predictors for the additional 'confounders'
		predictors<-as.data.frame(cbind(predictors, predict_val))
		
		# Do the predictions
		predictions<-predict(GAM, newdata=predictors, type="response", exclude=exclude, se.fit=T)
		
		# Edit out based on the marker list
		fitted_values<-as.data.frame(cbind(predictions$fit[-which(in.poly == 0)], predictors[-which(in.poly == 0),c(1:3)]))
		names(fitted_values)<-c("fit", "x", "y", "z")
		
		# Find the min and max values across all predictions
		mn<-min(unlist(fitted_values$fit), na.rm=T)
		mx<-max(unlist(fitted_values$fit), na.rm=T)
		
		# If a color scale is specified for the outcomes double check it is within the range, if not expand the range
		if(is.na(scale_surface) == F){
			# Find the absolute max values across all predictions, and the scale
			upp_abs<-max(abs(c(scale_surface, mn, mx)))	
			mn<-(-upp_abs)
			mx<-upp_abs
		}
		
		# Folder for interactive plots if we want those
		int_dir<-paste0(traits[k])
		dir.create(int_dir, showWarnings = F)
		setwd(int_dir)
			
		# Scale the colours appropriately	
		locsplt<-round((range(fitted_values$fit, na.rm=TRUE) - mn) / (mx-mn) * no.cols)
		col.pal<-map[locsplt[1]:locsplt[2]]
		
		# Make the plotly object
		plotly_obj<-plot_ly(fitted_values, x=~x, y=~y, z=~z, color=~fit, colors=map, type="scatter3d", mode="markers") %>% 
		layout(title = traits[k], 
				scene=list(
				xaxis=list(title=nutrient.order[1]),
				yaxis=list(title=nutrient.order[2]),
				zaxis=list(title=nutrient.order[3])))
		plotly_obj<-colorbar(plotly_obj, limits=c(mn, mx))		
		
		# Save it
		htmlwidgets::saveWidget(plotly_obj, "index.html", selfcontained=F)
		
		# Reset the wd
		setwd(new_dir)
										
		# Update the progress
		progress<-progress + 1
		setTxtProgressBar(pb, progress)

	# Close loop for traits
	}
	
	# Make sure to reset the wd
	setwd(directory)

}




