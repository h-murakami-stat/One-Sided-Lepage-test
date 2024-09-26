###	Test the hypothesis that is the equality of location and scale parameters in the first and second samples

library(mvtnorm)
library(MASS)
library(medicaldata)
library(survival)

###	"LEPAGE" denotes the Lepage statistic for the one-sided alternative
###	"WMD" denotes the combination of Wilcoxon and Mood statistics for the one-sided alternative
###	"ADP" denotes the adaptive statistic based on "LEPAGE" and "WMD" for the one-sided alternative
###	"MAX" denotes the maximum statistic based on "LEPAGE" and "WMD" for the one-sided alternative
### "YES" denotes Interim Analysis 
### "NO" denotes NO Interim Analysis 
### "WITH" denotes the statistic with weighted constant in two-stage
### "WITHOUT" denotes the statistic with unweighted constant in two-stage
### "BK" denotes the adaptive design of Bauer and Kohne (1994)
### "LW" denotes the adaptive design of Lehmacher and Wassmer (1999)



##################################
##### Linear Rank Statistics #####
##################################

LRS = function(x, y){
	m1 = length(x)  ## Size of the first sample
	n1 = length(y)  ## Size of the second sample
	N = m1+n1  ## Total sample sizes

	x1 = sort(x)
	y1 = sort(y)
	z1 = rank(c(x1, y1))
	zy1 = z1[(m1 + 1):N]  ## Rank of the second sample

	u0 <- c(x1, y1)
	u1 <- sort(unique(u0))
	t <- tabulate(match(u0, u1), length(u1))  

## Score functions
	u2 = unique(sort(rank(u0)))/N  ## Wilcoxon score
	u3 = (N + 1)/(2*N) - abs(u2 - (N + 1)/(2*N))  ## Ansari-Bradley score
	u4 = (u2 - (N + 1)/(2*N))^2  ## Mood score

## Wilcoxon rank sum statistic
	T1 = sum(zy1/N)
	mu1 = n1/N*sum(u2*t)
	sd1 = sqrt(m1*n1/(N^2*(N-1))*(N*sum(t*u2^2) - (sum(t*u2))^2))

## Ansari-Bradley statistic
	T2 = n1*(N + 1)/(2*N) - sum(abs(zy1/N - (N + 1)/(2*N)))
	mu2 = n1/N*sum(u3*t)
	sd2 = sqrt(m1*n1/(N^2*(N-1))*(N*sum(t*u3^2) - (sum(t*u3))^2))

## Mood statistic
	T3 = sum( (zy1/N - (N + 1)/(2*N))^2 )
	mu3 = n1/N*sum(u4*t)
	sd3 = sqrt(m1*n1/(N^2*(N-1))*(N*sum(t*u4^2) - (sum(t*u4))^2))

## Covariances
## Wilcoxon-Ansari-Bradley
	covWAB = m1*n1/(N^2*(N-1))*(N*sum(t*u2*u3)-(sum(t*u2))*(sum(t*u3)))  
## Wilcoxon-Mood
	covWM  = m1*n1/(N^2*(N-1))*(N*sum(t*u2*u4)-(sum(t*u2))*(sum(t*u4)))
## Ansari-Bradley-Mood
	covABMD = m1*n1/(N^2*(N-1))*(N*sum(t*u3*u4)-(sum(t*u3))*(sum(t*u4)))

	return(list(
		T.W = T1, T.W.EX = mu1, T.W.SD = sd1, T.W.STND = (T1 - mu1)/sd1, 
		T.AB = T2, T.AB.EX = mu2, T.AB.SD = sd2, T.AB.STND = (T2 - mu2)/sd2, 
		T.MD = T3, T.MD.EX = mu3, T.MD.SD = sd3, T.MD.STND = (T3 - mu3)/sd3, 
		T.COV.WAB = covWAB, T.COV.WMD = covWM, T.COV.ABMD = covABMD
		))
	}



####################
##### Selector #####
####################

lower.mean = function(z, p){
	len.ind = p * length(z)
	Lp0 = floor(len.ind)

	if( Lp0 < 1 ){
		Lp = min(z)
	} else{
		z = sort(z, partial = 1:(Lp0 + 1))
		Lp1 = len.ind - Lp0
		Lp = (sum(z[1:Lp0]) + Lp1*z[Lp0 + 1])/len.ind
	}
	return(Lp)
}

upper.mean <- function(z, p){ -lower.mean(-z, p) }


Q1 = function(z){
	U05 = upper.mean(z, 0.05)
	M50 = mean(z, trim = 0.25)
	L05 = lower.mean(z, 0.05)
	return((U05 - M50)/(M50 - L05))
}

Q2 = function(z){
	L05 = lower.mean(z, 0.05)
	L50 = lower.mean(z, 0.5)
	U05 = upper.mean(z, 0.05)
	U50 = upper.mean(z, 0.5)
	return((U05 - L05)/(U50 - L50))
}



##################################################
##### First Stage One-Sided Lepage-type Test #####	
##################################################

ONE.first.LEPAGE.Stat = function(x,y){
	z0 = c(x, y)
	m1 = length(x)
	n1 = length(y)
	N = m1 + n1

## Wilcoxon rank sum statistic
	T1 = LRS(x, y)$T.W
	mu1 = LRS(x, y)$T.W.EX
	sd1 = LRS(x, y)$T.W.SD
	T1S = LRS(x, y)$T.W.STND

## Ansari-Bradley statistic
	T2 = LRS(x, y)$T.AB
	mu2 = LRS(x, y)$T.AB.EX
	sd2 = LRS(x, y)$T.AB.SD
  	T2S = LRS(x, y)$T.AB.STND
	COVWAB = LRS(x, y)$T.COV.WAB
	
## Mood statistic
	T3 = LRS(x, y)$T.MD
	mu3 = LRS(x, y)$T.MD.EX
	sd3 = LRS(x, y)$T.MD.SD
  	T3S = LRS(x, y)$T.MD.STND
  	COVWMD = LRS(x, y)$T.COV.WMD

  	COVABMD = LRS(x, y)$T.COV.ABMD

## First Stage Statistics
	first.WAB = ( (T1 - mu1) - (T2 - mu2) )/sqrt(sd1^2 + sd2^2 - 2*COVWAB)
	first.WMD = ( (T1 - mu1) + (T3 - mu3) )/sqrt(sd1^2 + sd3^2 + 2*COVWMD)

## Weights of Statistics
	WGT1 = 1 - pnorm(T1S, lower.tail = FALSE)
	WGT2 = 1 - pnorm(T2S, lower.tail = TRUE)
	WGT3 = 1 - pnorm(T3S, lower.tail = FALSE)

	WGT.WAB1 = WGT1/(WGT1 + WGT2)  
	WGT.WAB2 = 1 - WGT.WAB1  

	WGT.WMD1 = WGT1/(WGT1 + WGT3)
	WGT.WMD2 = 1 - WGT.WMD1

## P-value of Statistic in the First Stage
	pval.first.WAB = pnorm(first.WAB, lower.tail = FALSE)
	pval.first.WMD = pnorm(first.WMD, lower.tail = FALSE)

	if(Q2(z0) > 7){
		first.ADP = first.WAB
		pval.first.ADP = pval.first.WAB
		WGT.ADP1 = WGT.WAB1
		WGT.ADP2 = WGT.WAB2
		STATISTIC = "OLS1.WAB"
	} else {
		first.ADP = first.WMD
		pval.first.ADP = pval.first.WMD
		WGT.ADP1 = WGT.WMD1
		WGT.ADP2 = WGT.WMD2
		STATISTIC = "OLS2.WMD"
	}

	first.MAX = max(first.WAB, first.WMD)
	meanvec = c(0, 0)
	a12 = a21 = (sd1^2 + COVWMD - COVWAB - COVABMD)/sqrt( (sd1^2 + sd2^2 - 2*COVWAB)*(sd1^2 + sd3^2 + 2*COVWMD) )
	asymcv = rep(first.MAX, 2)
	corrmtrx = matrix (c(1.0, a12, a21, 1.0), 2, 2)
	pval.first.MAX = 1- pmvnorm(meanvec, corrmtrx, lower = -Inf, upper = asymcv)[1]

return(list(
	T.first.WAB = first.WAB, 
	T.first.WMD = first.WMD, 
	T.first.ADP = first.ADP, 
	T.first.MAX = first.MAX, 
	pvalue.first.WAB = pval.first.WAB, 
	pvalue.first.WMD = pval.first.WMD, 
	pvalue.first.ADP = pval.first.ADP, 
	pvalue.first.MAX = pval.first.MAX, 
	Weight.WAB1 = WGT.WAB1, 
	Weight.WAB2 = WGT.WAB2, 
	Weight.WMD1 = WGT.WMD1, 
	Weight.WMD2 = WGT.WMD2,
	STAT = STATISTIC
))
}


###################################################
##### Second Stage One-Sided Lepage-type Test #####	
###################################################
ONE.second.LEPAGE.Stat = function(x, y, wab1, wab2, wmd1, wmd2, TEST){
	m2 = length(x)
	n2 = length(y)
	N = m2 + n2

## Wilcoxon rank sum statistic
	T1 = LRS(x, y)$T.W
	mu1 = LRS(x, y)$T.W.EX
	sd1 = LRS(x, y)$T.W.SD
	T1S = LRS(x, y)$T.W.STND

## Ansari-Bradley statistic
	T2 = LRS(x, y)$T.AB
	mu2 = LRS(x, y)$T.AB.EX
	sd2 = LRS(x, y)$T.AB.SD
  	T2S = LRS(x, y)$T.AB.STND
	COVWAB = LRS(x, y)$T.COV.WAB
	
## Mood statistic
	T3 = LRS(x, y)$T.MD
	mu3 = LRS(x, y)$T.MD.EX
	sd3 = LRS(x, y)$T.MD.SD
  	T3S = LRS(x, y)$T.MD.STND
  	COVWMD = LRS(x, y)$T.COV.WMD

  	COVABMD = LRS(x, y)$T.COV.ABMD

	second.WAB = ( wab1*(T1 - mu1) - wab2*(T2 - mu2) )/sqrt(wab1^2*sd1^2 + wab2^2*sd2^2 - 2*wab1*wab2*COVWAB)
	second.WMD = ( wmd1*(T1 - mu1) + wmd2*(T3 - mu3) )/sqrt(wmd1^2*sd1^2 + wmd2^2*sd3^2 + 2*wmd1*wmd2*COVWMD)

	pval.second.WAB = pnorm(second.WAB, lower.tail = FALSE)
	pval.second.WMD = pnorm(second.WMD, lower.tail = FALSE)

	if(TEST == "OLS1.WAB"){
	second.ADP = second.WAB	
	pval.second.ADP = pnorm(second.ADP, lower.tail = FALSE)
	} else if(TEST == "OLS2.WMD"){
	second.ADP = second.WMD
	pval.second.ADP = pnorm(second.ADP, lower.tail = FALSE)
	} 		
	
	second.MAX = max(second.WAB, second.WMD)
	meanvec = c(0, 0)
	a12 = a21 = (wab1*wmd1*sd1^2 + wab1*wmd2*COVWMD - wab2*wmd1*COVWAB - wab2*wmd2*COVABMD)/sqrt( (wab1^2*sd1^2 + wab2^2*sd2^2 - 2*wab1*wab2*COVWAB)*(wmd1^2*sd1^2 + wmd2^2*sd3^2 + 2*wmd1*wmd2*COVWMD)	)
	
	corrmtrx = matrix (c(1.0, a12, a21, 1.0), 2, 2)
	asymcv = rep(second.MAX, 2)
	pval.second.MAX = 1- pmvnorm(meanvec, corrmtrx, lower = -Inf, upper = asymcv)[1]


	return(list(
	T.second.WAB = second.WAB, 
	T.second.WMD = second.WMD, 
	T.second.ADP = second.ADP,
	T.second.MAX = second.MAX,
	pvalue.second.WAB = pval.second.WAB, 
	pvalue.second.WMD = pval.second.WMD, 
	pvalue.second.ADP = pval.second.ADP,
	pvalue.second.MAX = pval.second.MAX
	))
}

##################################
##### Bauer and Kohne (1994) #####
##################################

p.val.form.BK = function(alpha){
	if(alpha == 0.05){
		alpha1 = 0.0233
 		calpha0 = 0.0087
	} else if(alpha == 0.025){
		alpha1 = 0.0102
		calpha0 = 0.0038
	}
return(list(cr.prob.alpha1 = alpha1, cr.prob.final.dicision = calpha0))
}


########################################
##### Lehmacher and Wassmer (1999) #####
#####   Weighted Inverse Normal    #####
########################################

p.val.form.LW = function(alpha){
	if(alpha == 0.05){
		alpha1 = 0.00443356
		calpha0 = 0.05
	} else if(alpha == 0.025){
		alpha1 = 0.00107783
		calpha0 = 0.025
	}
return(list(cr.prob.alpha1 = alpha1, cr.prob.final.dicision = calpha0))
}    
    

###	Test the hypothesis that is the equality of location and scale parameters in the first and second samples

One.Sided.Lepage.type.statistics = function(x1, y1, x2, y2, xx, yy, significance.level, method = c("LEPAGE", "WMD", "ADP", "MAX"), Interim = c("YES", "NO"), weighted = c("WITH", "WITHOUT"), design = c("BK", "LW")){
###	"LEPAGE" denotes the Lepage statistic for the one-sided alternative
###	"WMD" denotes the combination of Wilcoxon and Mood statistics for the one-sided alternative
###	"ADP" denotes the adaptive statistic based on "LEPAGE" and "WMD" for the one-sided alternative
###	"MAX" denotes the maximum statistic based on "LEPAGE" and "WMD" for the one-sided alternative
### "YES" denotes Interim Analysis 
### "NO" denotes NO Interim Analysis 
### "WITH" denotes the statistic with weighted constant in two-stage
### "WITHOUT" denotes the statistic with unweighted constant in two-stage
### "BK" denotes the adaptive design of Bauer and Kohne (1994)
### "LW" denotes the adaptive design of Lehmacher and Wassmer (1999)


	if(Interim == "NO"){
	###############################
	##### No Interim Analysis #####
	###############################
	
	dname1 <- paste(deparse(substitute(xx)), "and", deparse(substitute(yy))) 

	if(method == "LEPAGE"){
		method.name1 <- c("One-sided Lepage statistic (No Interim Analysis)")
		N0.IA.OLS1 = round(ONE.first.LEPAGE.Stat(xx, yy)$T.first.WAB, 4) ### Statistic
		names(N0.IA.OLS1) <- "STATISTIC"
		N0.IA.P.Val.OLS1 = round(ONE.first.LEPAGE.Stat(xx, yy)$pvalue.first.WAB, 4) ### P-value
		names(N0.IA.P.Val.OLS1) <- "P-value of STATISTIC"
		res1 <- list(method = method.name1, data.name = dname1, statistic = N0.IA.OLS1, p.value = N0.IA.P.Val.OLS1)
	} else if(method == "WMD"){
		method.name1 <- c("One-sided Wilcoxon-Mood statistic (No Interim Analysis)")
		N0.IA.OLS2 = round(ONE.first.LEPAGE.Stat(xx, yy)$T.first.WMD, 4) ### Statistic
		names(N0.IA.OLS2) <- "STATISTIC"
		N0.IA.P.Val.OLS2 = round(ONE.first.LEPAGE.Stat(xx, yy)$pvalue.first.WMD, 4) ### P-value
		names(N0.IA.P.Val.OLS2) <- "P-value of STATISTIC"
		res1 <- list(method = method.name1, data.name = dname1, statistic = N0.IA.OLS2, p.value = N0.IA.P.Val.OLS2)
	} else if(method == "ADP"){
		method.name1 <- c("One-sided Adaptive statistic (No Interim Analysis)")
		N0.IA.ADP = round(ONE.first.LEPAGE.Stat(xx, yy)$T.first.ADP, 4) ### Statistic
		names(N0.IA.ADP) <- "STATISTIC"
		N0.IA.P.Val.ADP = round(ONE.first.LEPAGE.Stat(xx, yy)$pvalue.first.ADP, 4) ### P-value
		names(N0.IA.P.Val.ADP) <- "P-value of STATISTIC"
		res1 <- list(method = method.name1, data.name = dname1, statistic = N0.IA.ADP, p.value = N0.IA.P.Val.ADP)
	} else if(method == "MAX"){
		method.name1 <- c("One-sided Maximum statistic (No Interim Analysis)")
		N0.IA.MAX = round(ONE.first.LEPAGE.Stat(xx, yy)$T.first.MAX, 4) ### Statistic
		names(N0.IA.MAX) <- "STATISTIC"
		N0.IA.P.Val.MAX = round(ONE.first.LEPAGE.Stat(xx, yy)$pvalue.first.MAX, 4) ### P-value
		names(N0.IA.P.Val.MAX) <- "P-value of STATISTIC"
		res1 <- list(method = method.name1, data.name = dname1, statistic = N0.IA.MAX, p.value = N0.IA.P.Val.MAX)
	}

	} else if(Interim == "YES"){
	##############################
	##### Two Stage Analysis #####
	##############################

	dname2 <- paste(deparse(substitute(x1)), "and", deparse(substitute(y1)))
	dname3 <- paste(deparse(substitute(x2)), "and", deparse(substitute(y2)))

	weight.wab1 = ONE.first.LEPAGE.Stat(x1, y1)$Weight.WAB1
	weight.wab2 = ONE.first.LEPAGE.Stat(x1, y1)$Weight.WAB2
	weight.wmd1 = ONE.first.LEPAGE.Stat(x1, y1)$Weight.WMD1
	weight.wmd2 = ONE.first.LEPAGE.Stat(x1, y1)$Weight.WMD2


	if(design =="BK"){

		if(method == "LEPAGE"){
			FIRST.OLS1 = round(ONE.first.LEPAGE.Stat(x1, y1)$T.first.WAB, 4)
			names(FIRST.OLS1) <- "STATISTIC"
			pval1.wab = round(ONE.first.LEPAGE.Stat(x1, y1)$pvalue.first.WAB, 4)

			if(pval1.wab > 0.5){
				method.name1 <- c("One-sided Lepage Statistic based on Bauer and Kohne (1994) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.OLS1, p.value = pval1.wab)
			} else if(pval1.wab < p.val.form.BK(significance.level)$cr.prob.alpha1){
				method.name1 <- c("One-sided Lepage statistic based on Bauer and Kohne (1994) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.OLS1, p.value = pval1.wab)
			} else {
				if(weighted == "WITHOUT"){	
					method.name1 <- c("One-sided Unweighted Lepage Statistic based on Bauer and Kohne (1994) (Second Stage)")
					stat3.wab = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, "OLS1.WAB")$T.second.WAB, 4)
					pval3.wab = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, "OLS1.WAB")$pvalue.second.WAB, 4)
					pvalunweight.wab = round(pval1.wab*pval3.wab, 4)
					names(pvalunweight.wab) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pvalunweight.wab, p.value = 1-pchisq(-2*log(pvalunweight.wab), 4))
				} else if(weighted == "WITH"){
					method.name1 <- c("One-sided Weighted Lepage Statistic based on Bauer and Kohne (1994) (Second Stage)")
					stat2.wab = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, "OLS1.WAB")$T.second.WAB, 4)
					pval2.wab = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, "OLS1.WAB")$pvalue.second.WAB, 4)
					pval.wab = round(pval1.wab*pval2.wab, 4)
					names(pval.wab) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pval.wab, p.value = 1-pchisq(-2*log(pval.wab), 4))
				}
			}
		} else if(method == "WMD"){
			FIRST.OLS2 = round(ONE.first.LEPAGE.Stat(x1, y1)$T.first.WMD, 4)
			names(FIRST.OLS2) <- "STATISTIC"
			pval1.wmd = round(ONE.first.LEPAGE.Stat(x1, y1)$pvalue.first.WMD, 4)

			if(pval1.wmd > 0.5){
				method.name1 <- c("One-sided Wilcoxon-Mood Statistic based on Bauer and Kohne (1994) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.OLS2, p.value = pval1.wmd)
			} else if(pval1.wmd < p.val.form.BK(significance.level)$cr.prob.alpha1){
				method.name1 <- c("One-sided Wilcoxon-Mood statistic based on Bauer and Kohne (1994) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.OLS2, p.value = pval1.wmd)
			} else {
				if(weighted == "WITHOUT"){	
					method.name1 <- c("One-sided Unweighted Wilcoxon-Mood Statistic based on Bauer and Kohne (1994) (Second Stage)")
					stat3.wmd = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, "OLS2.WMD")$T.second.WMD, 4)
					pval3.wmd = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, "OLS2.WMD")$pvalue.second.WMD, 4)
					pvalunweight.wmd = round(pval1.wmd*pval3.wmd, 4)
					names(pvalunweight.wmd) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pvalunweight.wmd, p.value = 1-pchisq(-2*log(pvalunweight.wmd), 4))
				} else if(weighted == "WITH"){	
					method.name1 <- c("One-sided Weighted Wilcoxon-Mood Statistic based on Bauer and Kohne (1994) (Second Stage)")
					stat2.wmd = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, "OLS2.WMD")$T.second.WMD, 4)
					pval2.wmd = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, "OLS2.WMD")$pvalue.second.WMD, 4)
					pval.wmd = round(pval1.wmd*pval2.wmd, 4)
					names(pval.wmd) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pval.wmd, p.value = 1-pchisq(-2*log(pval.wmd), 4))
				}
			}
		} else if(method == "ADP"){
			FIRST.ADP = round(ONE.first.LEPAGE.Stat(x1, y1)$T.first.ADP, 4)
			names(FIRST.ADP) <- "STATISTIC"
			pval1.adp = round(ONE.first.LEPAGE.Stat(x1, y1)$pvalue.first.ADP, 4)

			if(Q2(c(x1, y1)) > 7){
				ADAPT = "OLS1.WAB"
			} else {
				ADAPT = "OLS2.WMD"
			}

			if(pval1.adp > 0.5){
				method.name1 <- c("One-sided Adaptive Statistic based on Bauer and Kohne (1994) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.ADP, p.value = pval1.adp)
			} else if(pval1.adp < p.val.form.BK(significance.level)$cr.prob.alpha1){
				method.name1 <- c("One-sided Adaptive statistic based on Bauer and Kohne (1994) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.ADP, p.value = pval1.adp)
			} else {
				if(weighted == "WITHOUT"){	
					method.name1 <- c("One-sided Unweighted Adaptive Statistic based on Bauer and Kohne (1994) (Second Stage)")
					stat3.adp = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, ADAPT)$T.second.ADP, 4)
					pval3.adp = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, ADAPT)$pvalue.second.ADP, 4)
					pvalunweight.adp = round(pval1.adp*pval3.adp, 4)
					names(pvalunweight.adp) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pvalunweight.adp, p.value = 1-pchisq(-2*log(pvalunweight.adp), 4))
				} else if(weighted == "WITH"){	
					method.name1 <- c("One-sided Weighted Adaptive Statistic based on Bauer and Kohne (1994) (Second Stage)")
					stat2.adp = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, ADAPT)$T.second.ADP, 4)
					pval2.adp = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, ADAPT)$pvalue.second.ADP, 4)
					pval.adp = round(pval1.adp*pval2.adp, 4)
					names(pval.adp) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pval.adp, p.value = 1-pchisq(-2*log(pval.adp), 4))
				}
			}
		} else if(method == "MAX"){
			FIRST.MAX = round(ONE.first.LEPAGE.Stat(x1, y1)$T.first.MAX, 4)
			names(FIRST.MAX) <- "STATISTIC"
			pval1.max = round(ONE.first.LEPAGE.Stat(x1, y1)$pvalue.first.MAX, 4)

			if(ONE.first.LEPAGE.Stat(x1, y1)$T.first.WAB > ONE.first.LEPAGE.Stat(x1, y1)$T.first.WMD){
				MAXIMUM = "OLS1.WAB"
			} else {
				MAXIMUM = "OLS2.WMD"	
			}

			if(pval1.max > 0.5){
				method.name1 <- c("One-sided Maximum Statistic based on Bauer and Kohne (1994) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.MAX, p.value = pval1.max)
			} else if(pval1.max < p.val.form.BK(significance.level)$cr.prob.alpha1){
				method.name1 <- c("One-sided Maximum statistic based on Bauer and Kohne (1994) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.MAX, p.value = pval1.max)
			} else {
				if(weighted == "WITHOUT"){	
					method.name1 <- c("One-sided Unweighted Maximum Statistic based on Bauer and Kohne (1994) (Second Stage)")
					stat3.max = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, MAXIMUM)$T.second.MAX, 4)
					pval3.max = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, MAXIMUM)$pvalue.second.MAX, 4)
					pvalunweight.max = round(pval1.max*pval3.max, 4)
					names(pvalunweight.max) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pvalunweight.max, p.value = 1-pchisq(-2*log(pvalunweight.max), 4))
				} else if(weighted == "WITH"){	
					method.name1 <- c("One-sided Weighted Maximum Statistic based on Bauer and Kohne (1994) (Second Stage)")
					stat2.max = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, MAXIMUM)$T.second.MAX, 4)
					pval2.max = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, MAXIMUM)$pvalue.second.MAX, 4)
					pval.max = round(pval1.max*pval2.max, 4)
					names(pval.max) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pval.max, p.value = 1-pchisq(-2*log(pval.max), 4))
				}
			}
		}

	} else if(design == "LW"){

		if(method == "LEPAGE"){
			FIRST.OLS1 = round(ONE.first.LEPAGE.Stat(x1, y1)$T.first.WAB, 4)
			names(FIRST.OLS1) <- "STATISTIC"
			pval1.wab = round(ONE.first.LEPAGE.Stat(x1, y1)$pvalue.first.WAB, 4)

			if(pval1.wab > 0.5){
				method.name1 <- c("One-sided Lepage Statistic based on Lehmacher and Wassmer (1999) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.OLS1, p.value = pval1.wab)
			} else if(pval1.wab < p.val.form.LW(significance.level)$cr.prob.alpha1){
				method.name1 <- c("One-sided Lepage statistic based on Lehmacher and Wassmer (1999) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.OLS1, p.value = pval1.wab)
			} else {
				if(weighted == "WITHOUT"){	
					method.name1 <- c("One-sided Unweighted Lepage Statistic based on Lehmacher and Wassmer (1999) (Second Stage)")
					stat3.wab = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, "OLS1.WAB")$T.second.WAB, 4)
					pval3.wab = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, "OLS1.WAB")$pvalue.second.WAB, 4)
					pvalunweight.wab = round(1 - pnorm(1/sqrt(2)*qnorm(1 - pval1.wab) + 1/sqrt(2)*qnorm(1 - pval3.wab)), 4)
					names(pvalunweight.wab) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pvalunweight.wab, p.value = pvalunweight.wab)
				} else if(weighted == "WITH"){
					method.name1 <- c("One-sided Weighted Lepage Statistic based on Lehmacher and Wassmer (1999) (Second Stage)")
					stat2.wab = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, "OLS1.WAB")$T.second.WAB, 4)
					pval2.wab = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, "OLS1.WAB")$pvalue.second.WAB, 4)
					pval.wab = round(1 - pnorm(1/sqrt(2)*qnorm(1 - pval1.wab) + 1/sqrt(2)*qnorm(1 - pval2.wab)), 4)
					names(pval.wab) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pval.wab, p.value = pval.wab)
				}
			}
		} else if(method == "WMD"){
			FIRST.OLS2 = round(ONE.first.LEPAGE.Stat(x1, y1)$T.first.WMD, 4)
			names(FIRST.OLS2) <- "STATISTIC"
			pval1.wmd = round(ONE.first.LEPAGE.Stat(x1, y1)$pvalue.first.WMD, 4)

			if(pval1.wmd > 0.5){
				method.name1 <- c("One-sided Wilcoxon-Mood Statistic based on Lehmacher and Wassmer (1999) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.OLS2, p.value = pval1.wmd)
			} else if(pval1.wmd < p.val.form.LW(significance.level)$cr.prob.alpha1){
				method.name1 <- c("One-sided Wilcoxon-Mood statistic based on Lehmacher and Wassmer (1999) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.OLS2, p.value = pval1.wmd)
			} else {
				if(weighted == "WITHOUT"){	
					method.name1 <- c("One-sided Unweighted Wilcoxon-Mood Statistic based on Lehmacher and Wassmer (1999) (Second Stage)")
					stat3.wmd = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, "OLS2.WMD")$T.second.WMD, 4)
					pval3.wmd = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, "OLS2.WMD")$pvalue.second.WMD, 4)
					pvalunweight.wmd = round(1 - pnorm(1/sqrt(2)*qnorm(1 - pval1.wmd) + 1/sqrt(2)*qnorm(1 - pval3.wmd)), 4)
					names(pvalunweight.wmd) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pvalunweight.wmd, p.value = pvalunweight.wmd)
				} else if(weighted == "WITH"){	
					method.name1 <- c("One-sided Weighted Wilcoxon-Mood Statistic based on Lehmacher and Wassmer (1999) (Second Stage)")
					stat2.wmd = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, "OLS2.WMD")$T.second.WMD, 4)
					pval2.wmd = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, "OLS2.WMD")$pvalue.second.WMD, 4)
					pval.wmd = round(1 - pnorm(1/sqrt(2)*qnorm(1 - pval1.wmd) + 1/sqrt(2)*qnorm(1 - pval2.wmd)), 4)
					names(pval.wmd) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pval.wmd, p.value = pval.wmd)
				}
			}
		} else if(method == "ADP"){
			FIRST.ADP = round(ONE.first.LEPAGE.Stat(x1, y1)$T.first.ADP, 4)
			names(FIRST.ADP) <- "STATISTIC"
			pval1.adp = round(ONE.first.LEPAGE.Stat(x1, y1)$pvalue.first.ADP, 4)

			if(Q2(c(x1, y1)) > 7){
				ADAPT = "OLS1.WAB"
			} else {
				ADAPT = "OLS2.WMD"
			}

			if(pval1.adp > 0.5){
				method.name1 <- c("One-sided Adaptive Statistic based on Lehmacher and Wassmer (1999) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.ADP, p.value = pval1.adp)
			} else if(pval1.adp < p.val.form.LW(significance.level)$cr.prob.alpha1){
				method.name1 <- c("One-sided Adaptive statistic based on Lehmacher and Wassmer (1999) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.ADP, p.value = pval1.adp)
			} else {
				if(weighted == "WITHOUT"){	
					method.name1 <- c("One-sided Unweighted Adaptive Statistic based on Lehmacher and Wassmer (1999) (Second Stage)")
					stat3.adp = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, ADAPT)$T.second.ADP, 4)
					pval3.adp = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, ADAPT)$pvalue.second.ADP, 4)
					pvalunweight.adp = round(1 - pnorm(1/sqrt(2)*qnorm(1 - pval1.adp) + 1/sqrt(2)*qnorm(1 - pval3.adp)), 4)
					names(pvalunweight.adp) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pvalunweight.adp, p.value = pvalunweight.adp)
				} else if(weighted == "WITH"){	
					method.name1 <- c("One-sided Weighted Adaptive Statistic based on Lehmacher and Wassmer (1999) (Second Stage)")
					stat2.adp = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, ADAPT)$T.second.ADP, 4)
					pval2.adp = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, ADAPT)$pvalue.second.ADP, 4)
					pval.adp = round(1 - pnorm(1/sqrt(2)*qnorm(1 - pval1.adp) + 1/sqrt(2)*qnorm(1 - pval2.adp)), 4)
					names(pval.adp) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pval.adp, p.value = pval.adp)
				}
			}
		} else if(method == "MAX"){
			FIRST.MAX = round(ONE.first.LEPAGE.Stat(x1, y1)$T.first.MAX, 4)
			names(FIRST.MAX) <- "STATISTIC"
			pval1.max = round(ONE.first.LEPAGE.Stat(x1, y1)$pvalue.first.MAX, 4)

			if(ONE.first.LEPAGE.Stat(x1, y1)$T.first.WAB > ONE.first.LEPAGE.Stat(x1, y1)$T.first.WMD){
				MAXIMUM = "OLS1.WAB"
			} else {
				MAXIMUM = "OLS2.WMD"	
			}

			if(pval1.max > 0.5){
				method.name1 <- c("One-sided Maximum Statistic based on Lehmacher and Wassmer (1999) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.MAX, p.value = pval1.max)
			} else if(pval1.max < p.val.form.LW(significance.level)$cr.prob.alpha1){
				method.name1 <- c("One-sided Maximum statistic based on Lehmacher and Wassmer (1999) (First Stage)")
				res1 <- list(method = method.name1, data.name = dname2, statistic = FIRST.MAX, p.value = pval1.max)
			} else {
				if(weighted == "WITHOUT"){	
					method.name1 <- c("One-sided Unweighted Maximum Statistic based on Lehmacher and Wassmer (1999) (Second Stage)")
					stat3.max = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, MAXIMUM)$T.second.MAX, 4)
					pval3.max = round(ONE.second.LEPAGE.Stat(x2, y2, 1, 1, 1, 1, MAXIMUM)$pvalue.second.MAX, 4)
					pvalunweight.max = round(1 - pnorm(1/sqrt(2)*qnorm(1 - pval1.max) + 1/sqrt(2)*qnorm(1 - pval3.max)), 4)
					names(pvalunweight.max) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pvalunweight.max, p.value = pvalunweight.max)
				} else if(weighted == "WITH"){	
					method.name1 <- c("One-sided Weighted Maximum Statistic based on Lehmacher and Wassmer (1999) (Second Stage)")
					stat2.max = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, MAXIMUM)$T.second.MAX, 4)
					pval2.max = round(ONE.second.LEPAGE.Stat(x2, y2, weight.wab1, weight.wab2, weight.wmd1, weight.wmd2, MAXIMUM)$pvalue.second.MAX, 4)
					pval.max = round(1 - pnorm(1/sqrt(2)*qnorm(1 - pval1.max) + 1/sqrt(2)*qnorm(1 - pval2.max)), 4)
					names(pval.max) <- "STATISTIC"
					res1 <- list(method = method.name1, data.name = c(dname2, " AND ", dname3), statistic = pval.max, p.value = pval.max)
				}
			}
		}

	}
	}

	class(res1) <- "htest"
	return(res1)
}




####### EXAMPLE 1 #######
x1 = gbsg[gbsg$meno == "0" & gbsg$hormon =="0", "pgr"] ## 
x1 = x1[!is.na(x1)]
y1 = gbsg[gbsg$meno == "0" & gbsg$hormon =="1", "pgr"] ## First stage of first sample
y1 = y1[!is.na(y1)]


x2 = gbsg[gbsg$meno == "1" & gbsg$hormon =="0", "pgr" ] ## Second stage of second sample
x2 = x2[!is.na(x2)]
y2 = gbsg[gbsg$meno == "1" & gbsg$hormon =="1", "pgr" ] ## Second stage of first sample
y2 = y2[!is.na(y2)]


xx = c(x1,x2) ## Total of two stages of First sample
yy = c(y1,y2) ## Total of two stages of Second sample

###############



####### EXAMPLE 2 #######
ALLDATA = polyps
datatest = ALLDATA

x1 = datatest[datatest$treatment == "sulindac" & datatest$sex == "male", "number3m"] ## First stage of second sample
x1 = x1[!is.na(x1)]

y1 = datatest[datatest$treatment == "placebo" & datatest$sex == "male", "number3m"] ## First stage of first sample
y1 = y1[!is.na(y1)]


x2 = polyps[polyps$treatment == "sulindac" & polyps$sex == "male", "number12m"] ## Second stage of second sample
x2 = x2[!is.na(x2)]
y2 = polyps[polyps$treatment == "placebo" & polyps$sex == "male", "number12m"] ## Second stage of first sample
y2 = y2[!is.na(y2)]

xx = c(x1,x2) ## Total of two stages of First sample
yy = c(y1,y2) ## Total of two stages of Second sample
#########################


significance.level = 0.05


One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "LEPAGE", "NO", "WITHOUT", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "LEPAGE", "NO", "WITH", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "LEPAGE", "YES", "WITHOUT", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "LEPAGE", "YES", "WITH", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "LEPAGE", "YES", "WITHOUT", "LW")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "LEPAGE", "YES", "WITH", "LW")


One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "WMD", "NO", "WITHOUT", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "WMD", "NO", "WITH", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "WMD", "YES", "WITHOUT", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "WMD", "YES", "WITH", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "WMD", "YES", "WITHOUT", "LW")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "WMD", "YES", "WITH", "LW")


One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "ADP", "NO", "WITHOUT", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "ADP", "NO", "WITH", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "ADP", "YES", "WITHOUT", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "ADP", "YES", "WITH", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "ADP", "YES", "WITHOUT", "LW")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "ADP", "YES", "WITH", "LW")


One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "MAX", "NO", "WITHOUT", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "MAX", "NO", "WITH", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "MAX", "YES", "WITHOUT", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "MAX", "YES", "WITH", "BK")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "MAX", "YES", "WITHOUT", "LW")
One.Sided.Lepage.type.statistics(x1, y1, x2, y2, xx, yy, significance.level, "MAX", "YES", "WITH", "LW")

