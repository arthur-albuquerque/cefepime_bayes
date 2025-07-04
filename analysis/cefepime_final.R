if (!require("pacman")) install.packages("pacman")

pacman::p_load(here,
               rmeta,
               meta,
               bayesmeta,
               brms)       

##META -ANALYSIS
data= read.csv(here("data/cef_dataonly.csv"), header=T)
data
colnames(data)

##exclude sponsor studies
data2 <- subset(data, Comparator != "Control")

#exclude NA
data3 <- subset(data, N_cef != "NA")#withsponsor
data4 <- subset(data2, N_cef != "NA")#withoutsponsor


#c = meta.MH(lowdose$N_cef, lowdose$N_com, lowdose$death_cef, lowdose$deaths_com, names=lowdose$Study.ID)
#summary(c)

###BAYES META

# Apply continuity correction (0.5 to each cell)
#a <- data4$death_cef + 0.5
#b <- data4$N_cef - data4$death_cef + 0.5
#c <- data4$deaths_com + 0.5
#d <- data4$N_com - data4$deaths_com + 0.5

# Calculate logRR using corrected values
#data4$logRR <- log((a * d) / (b * c))

# Calculate standard error
#data4$SE <- sqrt(1/a + 1/b + 1/c + 1/d)

df <- escalc(measure="RR", ai=death_cef, n1i = N_cef, 
                   ci=deaths_com, n2i=N_com, 
                   slab=Study.ID, data=data4)

df$sei = sqrt(df$vi)

# Use the Turner et al. prior for a binary outcome with log odds ratio
TP <- TurnerEtAlPrior("all-cause mortality", "pharma", "pharma")
print(TP)

result <- bayesmeta(y = df$yi,
                    sigma = df$sei,
                    tau.prior=TP$dprior,
                    labels = df$Study.ID)

result
result$summary

# Extracting the necessary values for log RR (mu)
logRR_mean <- result$summary["mean", "mu"]
logRR_median <- result$summary["median", "mu"]
logRR_CI_lower <- result$summary["95% lower", "mu"]
logRR_CI_upper <- result$summary["95% upper", "mu"]

# Transform to RR scale by exponentiating
RR_mean <- exp(logRR_mean)
RR_median <- exp(logRR_median)
RR_CI_lower <- exp(logRR_CI_lower)
RR_CI_upper <- exp(logRR_CI_upper)

# Display the results
cat("Pooled RR (mean):", round(RR_mean, 2), "\n")
cat("Pooled RR (median):", round(RR_median, 2), "\n")
cat("95% CrI for RR: [", round(RR_CI_lower, 2), ",", round(RR_CI_upper, 2), "]\n")

1 - result$pposterior(mu=0) #posterior probability of mu>0
# 95% credible intervals for the effect mu:
result$post.interval(mu.level=0.95)


##META-REG
df <- escalc(measure="RR", ai=death_cef, n1i = N_cef, ci=deaths_com, n2i=N_com, slab=Study.ID, data=data4)
febneut <- (df$Intended.clinical.indication=="Febrile neutropenia")
pnuemonia <- (df$Intended.clinical.indication=="Pneumonia")
meningitis <- (df$Intended.clinical.indication=="Meningitis")
severebact <- (df$Intended.clinical.indication=="Severe bacterial infections")
uti <- (df$Intended.clinical.indication=="UTI")
other <- (df$Intended.clinical.indication=="Mixed")
X <- cbind("Febrile neutropenia"=as.numeric(febneut),
           "Pneumonia"=as.numeric(pnuemonia),
           "Meningitis"=as.numeric(meningitis),
           "UTI"=as.numeric(uti),
           "Mixed"=as.numeric(other),
           "Severe bacterial infections"=as.numeric(severebact))

bmr01 <- bmr(y=df$yi, sigma=sqrt(df$vi), X=X, labels = df$Study.ID)

carb <- (df$Comparator.group=="Carbapenem")
ceft <- (df$Comparator.group=="Cefotaxime/Ceftriaxone")
ceftaz <- (df$Comparator.group=="Ceftazidime")
tazo <- (df$Comparator.group=="Piperacillin-tazobactam")
X <- cbind("Carbapenem"=as.numeric(carb),
           "Cefotaxime/Ceftriaxone"=as.numeric(ceft),
           "Ceftazidime"=as.numeric(ceftaz),
           "Piperacillin-tazobactam"=as.numeric(tazo))
bmr02 <- bmr(y=df$yi, sigma=sqrt(df$vi), X=X, labels = df$Study.ID)

adults <- (df$Patient.population=="Adults")
kids <- (df$Patient.population=="Pediatrics")
X <- cbind("Adults"=as.numeric(adults), 
           "Pediatrics"=as.numeric(kids))
bmr03 <- bmr(y=df$yi, sigma=sqrt(df$vi), X=X, labels = df$Study.ID)

cefdose <- subset(df, Cefepime_q12 != "NA")
lowdose <- (cefdose$Cefepime_q12=="Low")
highdose <- (cefdose$Cefepime_q12=="High")
X <- cbind("Low"=as.numeric(lowdose), 
           "High"=as.numeric(highdose))
bmr04 <- bmr(y=cefdose$yi, sigma=sqrt(cefdose$vi), X=X, labels = cefdose$Study.ID)
