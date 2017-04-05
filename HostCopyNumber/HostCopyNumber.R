getwd()
library(plyr)

# 1. Amplification efficiency 

  # Ofav Actin


    Actine<-read.csv("Act.csv", header = T )
    Actine$Log<-log10(Actine$Dilution)
    
    LMact <- lm(Actine$Log~Actine$C_)
    summary(LMact)
    
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    #   (Intercept)  16.9898     0.3033   56.02 7.98e-14 ***
    #   Actine$Log   -3.0851     0.1002  -30.80 3.06e-11 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 0.5927 on 10 degrees of freedom
    # Multiple R-squared:  0.9896,	Adjusted R-squared:  0.9885 
    # F-statistic: 948.4 on 1 and 10 DF,  p-value: 3.057e-11
    
    
    plot(Actine$C_, Actine$Log) 
    abline(lm(Actine$Log~Actine$C_))
    
    Act.Efficiency<-((10^(-1/-3.0851))-1)*100
    Act.Efficiency  # 110.9301%
    
    
  # PAX Actin
    
    PaxC<-read.csv("PAX.csv", header = T )
    PaxC$Log<-log10(PaxC$Dilution)
    
    LMpax <- lm(PaxC$C_~PaxC$Log)
    summary(LMpax)
    
    
      # Residuals:
      #   Min       1Q   Median       3Q      Max 
      # -0.37119 -0.21937 -0.08436  0.11497  0.69272 
      # 
      # Coefficients:
      #   Estimate Std. Error t value Pr(>|t|)    
      # (Intercept)  20.8465     0.2205   94.55 9.43e-11 ***
      #   PaxC$Log     -3.2138     0.1179  -27.27 1.61e-07 ***
      #   ---
      #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 0.3727 on 6 degrees of freedom
      # (4 observations deleted due to missingness)
      # Multiple R-squared:  0.992,	Adjusted R-squared:  0.9907 
      # F-statistic: 743.6 on 1 and 6 DF,  p-value: 1.607e-07
    
    
    plot(PaxC$C_, PaxC$Log)
    abline(lm(PaxC$Log~PaxC$C_ ))
    
    Act.Efficiency<-((10^(-1/-3.2138))-1)*100
    Act.Efficiency # 104.719%
    
    PRim<-read.csv("Primers.csv", header = T )
    PRim$Log<-log10(PRim$Dilution)
    PRim<-na.omit(PRim)
    PRim$Primer <- factor (as.character(PRim$Primer))
    
    
    equation = function(PRim) {
      mod = lm(PRim$C_~PRim$Log , data=PRim)
      mod_sum = summary(mod)
      formula = sprintf("y= %.3f %+.3f*x", coef(mod)[1], coef(mod)[2])
      r = mod_sum$r.squared
      r2 = sprintf("r2= %.3f", r)
      x  = cor.test(~x + y)
      r0 = sprintf("r= %.3f", x[4])
      p1 = pf(mod_sum$fstatistic[1],mod_sum$fstatistic[2],mod_sum$fstatistic[3],lower.tail=F)
      p =sprintf("p = %.3f", p1)
      n0 = length(mod_sum$residual)
      n1 = sprintf("N = %.f", n0)
      data.frame(formula=formula, r=r0,r2=r2, p=p,n=n1, stringsAsFactors=FALSE)
    }
    
    equation_end = ddply(PRim, c("PRim"), equation) 
    summary(PRim)

        library(ggplot2)
ggplot(PRim, aes(PRim$C_, PRim$Log, colour=factor(Primer))) +
    geom_point() +
    geom_smooth(method=lm) + labs(x = "Ct", y = "Log10 (dilution)", colour = "Primer")
    annotate("text", c(18,18), c(-4,-5), label=equation_end$formula)

# 2. Calculate the host copy number

  
  # All libraries and sources
    source("STEPoneFunction.R")# R.Cunning steponeR function
    library(ggplot2)
  
  # Get list of plate files to read in
  host.plates <- list.files(path="data", pattern="csv$", full.names=T)
  host.plates
  
  
  # # # Run stepone function to get Ratios
  Host.Out <- steponeR(files=host.plates, target.ratios=("Actin_Of.OF_SC"), 
                       fluor.norm=NULL,
                       copy.number=NULL,
                       ploidy=list(OF_SC=2, Actin_Of=2),
                       extract=list(OF_SC=0.982, Actin_Of=0.982))
  
  # Target ratio results  
      Host<-Host.Out$result
      
  # 1. Check and remove NTC wells, clean data
      
      ntc <- Host[which(Host$Sample.Name=="0-NTC"), ]
      if(any(!is.na(ntc$CT))) warning("Template detected in NTC: interpret data with caution")
      Host <- droplevels(Host[!rownames(Host) %in% rownames(ntc), ])
      
      StDe5 <- Host[which((Host$OfavAct.CT.sd)>0.5|(Host$PAX_C2.CT.sd>0.5)), ]
      Host <- droplevels(Host[!rownames(Host) %in% rownames(StDe5), ])
  
  # 2. Separate spp and colonies
  
      Colonies <- rbind.fill(lapply(strsplit(as.character(Host$Sample.Name), split="-"), 
                                             function(X) data.frame(t(X))))
      colnames(Colonies) <- c("Colony", "Core")
      Host <- cbind(Colonies, Host[,-1])
  ##########
      Dilution <- rbind.fill(lapply(strsplit(as.character(Host$Colony), split="_"), 
                                  function(X) data.frame(t(X))))
      colnames(Dilution) <- c("Colony", "Dilution")
      Host <- cbind(Dilution, Host[,-2])
      
      summary(Host$Actin_Of.OF_SC)
      
      # RATIO EXPLORATION
      
     # par(mfrow=c(1,2))
      plot(x=Host$Colony, y=Host$Actin_Of.OF_SC)
     # plot(x=Host$Dilution, y=Host$OfavAct.PAX_C2)
      
      par(mfrow=c(2,2))
      
      # Mean CTs and sd bu Colony
      plot(x=Host$Colony, y=Host$OF_SC.CT.mean)
      plot(x=Host$Colony, y=Host$Actin_Of.CT.mean)
      plot(x=Host$Colony, y=Host$OF_SC.CT.sd)
      plot(x=Host$Colony, y=Host$Actin_Of.CT.sd)
      
  
      # # Mean CTs and sd by dilution
        # plot(x=Host$Dilution, y=Host$OfavAct.CT.mean)
        # plot(x=Host$Dilution, y=Host$PAX_C2.CT.mean)
        # plot(x=Host$Dilution, y=Host$OfavAct.CT.sd)
        # plot(x=Host$Dilution, y=Host$PAX_C2.CT.sd)
      # High variation bw and inside colonies. Dulition does not seem to affect the number
      
      par(mfrow=c(1,1))
      hist(Host$Actin_Of.OF_SC, breaks = 20)
      
    # plot(Host$OfavAct.CT.sd, Host$OfavAct.PAX_C2, col=Host$Colony)
    # plot(Host$PAX_C2.CT.sd, Host$OfavAct.PAX_C2, col=Host$Colony)
      
     
      ggplot(Host, aes(x=Actin_Of.OF_SC, colour=Colony)) +
      geom_density(alpha=0.25) +
      facet_wrap("Colony")
      
      
      # Preserves marginal densities
      ggplot(Host, aes(OfavAct.PAX_C2, ..count.., fill = Colony)) +
        geom_density(position = "stack")
      
      ggplot(Host, aes(OfavAct.PAX_C2, fill = Colony, colour = Colony)) +
        geom_density(alpha = 0.1)
      
    ggplot(Host, aes((factor(Colony)), OfavAct.PAX_C2 )) + geom_boxplot() + geom_point()  + theme(axis.title = element_text(size = 13), 
      axis.text = element_text(size = 12), 
      axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size = 12), 
      panel.background = element_rect(fill = NA)) +labs(title = "Mean copy number by parent colony", 
      x = "Colony", y = "OfavAct : PAX_C2 ")
   
    summary(Host$Actin_Of.OF_SC)
    ddply(Host, c("Colony"), summarize, 
          outVal =c(mean(Actin_Of.OF_SC), sd(Actin_Of.OF_SC)))
    
   AO<- aov(Actin_Of.OF_SC ~ Colony, data = Host)
    summary(AO)
    
    # library(TukeyC)
    # av<-with(Host, aov(OfavAct.PAX_C2 ~ Colony, data=Host))
    # summary(av)
    # 
    # Tukey <- with(Host, TukeyC(x=av, which= "x"))
    # 
    # summary(Tukey)
    TK<-TukeyHSD(AO, data=Host, conf.level = 0.99)
    plot(TK)
    
    
    require(graphics)
    
    summary(AoV <- aov(OfavAct.PAX_C2 ~ Colony, data = Host))
    TukeyHSD(AoV, ordered = TRUE)
    plot(TukeyHSD(AoV))
    