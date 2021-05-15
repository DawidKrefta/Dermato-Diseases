#install.packages("tseries")
#install.packages("lmtest")
library(readxl)
library(lmtest)
library(tseries)
library(plyr) # ldply function
library(Kendall) #SeasonalMannKendall, pwmk
library(modifiedmk)
library(forecast) #tbats
library(stats) #symnum

Derma <- read_excel("Derma.xlsx", col_types = c("date", 
                                                 "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
#----------- checking & deleting columns which contains 0.1 ------
  contains_NA<-as.logical(grepl("0.1",Derma)) #only to save which column has 0.1

  #creating vector Na_column with column names to delete
  Na_column<-c()
  for (i in 1:ncol(Derma)) {
    if(contains_NA[i]){Na_column<-append(Na_column,colnames(Derma)[i])}
  }
  Na_column
  NA_df<-Derma[Na_column]
  Derma[Derma==0.1]<-NA
  Derma = Derma[,colSums(is.na(Derma)) == 0]
  

#----------- Calculating slope [RSV/Year] using lm --------------------------
    #creating data frame to store results without Data
  df <- setNames(data.frame(matrix(ncol = ncol(Derma))),names(Derma)[-1])
    #loop
  rounding<-2
  for(i in names(Derma)[-1]){
      #1 - slope is given for 1 second - multiplying to get 1 year
    df[[1,i]] <- lm(get(i) ~ Data, Derma)$coefficients[2]*60*60*24*365
      #rounding
    df[[1,i]] <- round(as.numeric(df[[1,i]]),rounding)
      #2 - p-value of the slope
    df[[2,i]] <- summary(lm(get(i) ~ Data, Derma))$coefficients[,4][2]
      #3 - p-value stars
    df[[3,i]] <- symnum(as.numeric(df[[2,i]]), corr = FALSE, na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
      #rounding and format because round(14.0001,2)= 14 and we loose zeros
    df[[2,i]] <- format(round(as.numeric(df[[2,i]]),rounding),nsmall = rounding)
  }
#----------- Calculating SeasonalMannKendall tau --------------------------  
  for(i in names(Derma)[-1]){
    # Seasonal Man Kendall tau
    df[[4,i]]<-SeasonalMannKendall(ts(Derma[,i], start=c(2010, 1), end=c(2020, 12), frequency=12))$tau
    # tau p-value
    df[[5,i]] <- format(SeasonalMannKendall(ts(Derma[,i], start=c(2010, 1), end=c(2020, 12), frequency=12))$sl, nsmall = 4)
    df[[6,i]] <- symnum(as.numeric(df[[5,i]]), corr = FALSE, na = FALSE, cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
    df[[4,i]] <- format(round(as.numeric(df[[4,i]]),rounding), nsmall=rounding)
    df[[5,i]] <- format(round(as.numeric(df[[5,i]]),rounding), nsmall=rounding)
  }
row.names(df) <- c("[RSV/Year](lm)","[RSV/Year] - p-value","[RSV/Year] stars","Mann-Kendall tau","tau p-value","tau p-value stars")
  
#----------- TBATS?? --------------------------  
  for(i in names(Derma)[-1]){
    df[[7,i]] <- tbats(ts(Derma[,i], start=c(2010, 1), end=c(2020, 12), frequency=12),use.box.cox = TRUE)$seasonal.periods
  }
  #ts.plot(Derma$czyrak)
tbats(ts(Derma$alergia, start=c(2010, 1), end=c(2020, 12), frequency=12),use.box.cox = TRUE)$seasonal.periods
tbats(ts(Derma$czyrak, start=c(2010, 1), end=c(2020, 12), frequency=12),use.box.cox = TRUE)$seasonal.periods
  

#----------- Checking stationarity ---------------------------------
  ht <- lapply(Derma, adf.test, alternative="stationary", k=0)
  pvalues = t(ldply(ht, function(x){ x$p.value }))
  pvalues[2,]
  
#---------- Adding stationarity to df --------
df<-rbind(df,pvalues[2,]) 



