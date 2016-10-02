#----------------- S:H RATIOS in Ofav Control vs NH4

getwd()

# All libraries and sources
    source("STEPoneFunction.R")# R.Cunning steponeR function
    library(plyr)
    library(reshape2)
    library(ggplot2)



##----------------------------------------------------------##
# 1. qPCR RATIOS 
##----------------------------------------------------------##

# +++++++++
# Libraries                           
# source("STEPoneFunction.R")
# +++++++++

# # # Get the raw data for Ofav R.Cunning steponeR function----------


# Get list of plate files to read in
  Ofav.plates <- list.files(path="data", pattern="csv", full.names=T)
  Ofav.plates


# # # Run stepone function to get Ratios
Ofav.Out <- steponeR(files=Ofav.plates, target.ratios=c("A.Ofav", "B.Ofav","C.Ofav", "D.Ofav"), 
                     fluor.norm=list(A=-0.064, B=4.197, C=3.798, D=0, Ofav=7.416),
                     copy.number=list(A=1, B=1,C=50, D=3, Ofav=14),
                     ploidy=list(A=1, B=1, C=1, D=1, Ofav=2),
                     extract=list(A=0.813, B=0.813, C=0.813, D=0.813, Ofav=0.982))

# Target ratio results
Ofav<-Ofav.Out$result

##----------------------------------------------------------##
# DATA CLEANING part A

## 1. Check and remove NTC wells
ntc <- Ofav[which(Ofav$Sample.Name=="NTC"), ]
if(any(!is.na(ntc$CT))) warning ("Template detected in NTC: interpret data with caution")
Ofav <- droplevels(Ofav[!rownames(Ofav) %in% rownames(ntc), ])

## 2. Check and remove + Control wells
Positive <- Ofav[which(Ofav$Sample.Name=="Control"), ]
Ofav <- droplevels(Ofav[!rownames(Ofav) %in% rownames(Positive), ])

##----------------------------------------------------------##

# Calculate total S/H ratio, clade spacific ratios and Clade Proportions

# 1. If Clade only detected in one technical replicate, set its ratio to NA (becomes zero)
    Ofav$A.Ofav[which(Ofav$A.reps==1)] <- NA
    Ofav$B.Ofav[which(Ofav$B.reps==1)] <- NA
    Ofav$C.Ofav[which(Ofav$C.reps==1)] <- NA
    Ofav$D.Ofav[which(Ofav$D.reps==1)] <- NA

# 2. Rename cols and make NA=0
    colnames(Ofav)[which(colnames(Ofav) %in% c("A.Ofav", "B.Ofav","C.Ofav", "D.Ofav" ))] <- c("A.SH", "B.SH","C.SH", "D.SH")  

    Ofav$A.SH[is.na(Ofav$A.SH)] <- 0
    Ofav$B.SH[is.na(Ofav$B.SH)] <- 0
    Ofav$C.SH[is.na(Ofav$C.SH)] <- 0
    Ofav$D.SH[is.na(Ofav$D.SH)] <- 0

# 3. Get the ratios and log 10
    # Total ratio
    Ofav$tot.SH <- Ofav$A.SH + Ofav$B.SH + Ofav$C.SH + Ofav$D.SH    
    Ofav$logTot.SH <- log10(Ofav$tot.SH )  # Calculate log10 SH ratio
    
    # A LOG10 ratio
    Ofav$logA.SH <- log10(Ofav$A.SH)
    # B LOG10 ratio
    Ofav$logB.SH <- log10(Ofav$B.SH)
    # C LOG10 ratio
    Ofav$logC.SH <- log10(Ofav$C.SH)
    # D LOG10 ratio
    Ofav$logD.SH <- log10(Ofav$D.SH)  
  

# 4. Clade Proportion
      # D Proportion
      Ofav$D.Prp<-(Ofav$D.SH/Ofav$tot.SH)
      hist(Ofav$D.Prp)

      Ofav$logTot.SH[which(Ofav$tot.SH==0)] <- NA
      Ofav$logA.SH[which(Ofav$A.SH==0)] <- NA
      Ofav$logB.SH[which(Ofav$B.SH==0)] <- NA
      Ofav$logC.SH[which(Ofav$C.SH==0)] <- NA
      Ofav$logD.SH[which(Ofav$D.SH==0)] <- NA


##----------------------------------------------------------##
# DATA CLEANING part B

## 3. Check and remove NoHSratio samples
    NoHSratio <- Ofav[which(Ofav$tot.SH==0), ]
    Ofav <- droplevels(Ofav[!rownames(Ofav) %in% rownames(NoHSratio), ])
    
    
## 4. Chose bw samples ran more than once
    
  # Create a unique ID for each sample.run
    
    Ofav$ID<-paste(Ofav$Sample.Name,Ofav$File.Name, sep = ".")
    
  # Check the number of times that a sigle sample was ran  
    ReRunA <- Ofav[duplicated(Ofav$Sample.Name),]
    
    n_RunA <- data.frame(table(Ofav$Sample.Name))
    colnames(n_RunA)<-c("Sample.Name","RanA")
    Ofav<-join(Ofav, n_RunA, type = "left")
    
    DuplicatesA <- Ofav[(Ofav$RanA>1),]
    write.csv(DuplicatesA, file = 'DuplicatesA.csv')
    # Look and the duplicates and manueally select best runs. 
    # Save the ones you are removing as ToRem1.csv
    
  # Remove bad replicates
    ToRem1<-read.csv("ToRem1.csv") 
    # 09/01/16
    Ofav<-Ofav[!(Ofav$ID %in% ToRem1$ID),]
    
  # Check for replicates again--should have none
    
    n_RunB <- data.frame(table(Ofav$Sample.Name))
    ReRunB <- Ofav[duplicated(Ofav$Sample.Name),]
    colnames(n_RunB)<-c("Sample.Name","RanB")
    Ofav<-join(Ofav, n_RunB, type = "left")
    
    
  # List of dupplicated samples, should have 0 rows now
    DuplicatesB <- Ofav[(Ofav$RanB>1),]
    # write.csv(DuplicatesC, file = 'DuplicatesC.csv')
    
  # Check the data!
    summary(Ofav)
    
  # Export data if want local backup
  # write.csv(Ofav, file = 'OfavSH.csv')
    
    
#----------------------------------------------------------#
# DATA ANALYSIS
#----------------------------------------------------------#
    
  # 1.Labeling and Factors
    
    ## Import treatment info from last Lauren's file 
    Treatments<-read.csv("TGFB-Final_AP.csv", header = TRUE)
    
    ## Merge SH ratios and Treatments 
    
    Ofav.data<-join(Ofav, Treatments, by = "Sample.Name", 
                    type = "full", match = "all")
    
  # 2. Exploratory graphs 

    # Libraries   
    # library(ggplot2)
    # +++++++++
    
    HistoSH<-qplot(logTot.SH, data=Ofav.data, binwidth=0.15)  
    
    HistoSH + 
      facet_wrap(Treatment~Treatment.2) + 
      geom_density()
    
    hist_Tre <- ggplot(Ofav.data, aes(x=logTot.SH, fill=Treatment))
    
    hist_Tre + 
      geom_density(alpha = 0.2) + facet_grid(Treatment.2~.)
    
    # Ggplot shows by default median instead mean vlues
    logSHTreatment <- ggplot(Ofav.data, aes(factor(Treatment), logTot.SH)) + 
    logSHTreatment +
      geom_boxplot(aes(fill=factor(Treatment))) + 
      geom_jitter(width = 0.2) +
      facet_grid ((~Treatment.2))
    
    # Plot mean +- SD
    
    # Create a function
        min.mean.sd.max <- function(x) {
          r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
          names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
          r
        }
        
    # ggplot code using mean and SD
      PlotMeanAndSD <- ggplot(aes(y = logTot.SH, x = factor(Treatment)), data = Ofav.data) +
                  facet_grid ((~Treatment.2))
      
      PlotMeanAndSD <- PlotMeanAndSD + 
         stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") +
         geom_jitter(position=position_jitter(width=.2), size=2) + 
         xlab("Treatments") + ylab("Log10(SH cell ratio)")
  
  # 3. ANOVAS
    
    # Effectof nutrients alone:
      
    OnwWayANOVA<-aov(logTot.SH~Treatment, data=Ofav.data)
        summary(OnwWayANOVA)

    # Effects of nutrients * treatment # 2
    
    TwoWayANOVA<-aov(logTot.SH ~ Treatment*Treatment.2 , data=Ofav.data)
        summary(TwoWayANOVA)
        drop1(TwoWayANOVA,~.,test="F")

    