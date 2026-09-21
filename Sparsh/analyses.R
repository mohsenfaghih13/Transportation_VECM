#analysis

pkgs <- c(
  "readxl","flextable","officer","dplyr","lubridate","tsDyn","urca","vars",
  "tibble","tidyr","tseries","forecast","TSstudio","ggplot2","zoo","corrplot","FinTS",
  "strucchange","sandwich","lmtest"
)
need <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(need)) install.packages(need)
invisible(lapply(pkgs, library, character.only = TRUE))


library(urca)

#start writing everything to a log file
log_file <- paste0("analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
sink(log_file, split = TRUE)

raw_df <- read.csv("C:/Users/sparshs/Research Projects/TC-IC Research Note/Transportation_VECM/cleaned_master_timeseries.csv")
raw_df<-raw_df[-1,]
#Log transforming all price series prior to estimation

#Tried with and without deflation - results don't really change tbh

#With deflation
raw_df$LnTL<-log(raw_df$Trucking_LD_Truckload/raw_df$PPI_All_Commodities)
raw_df$LnLTL<-log(raw_df$Trucking_LD_LTL/raw_df$PPI_All_Commodities)
raw_df$LnLocal<-log(raw_df$Trucking_Local/raw_df$PPI_All_Commodities)
raw_df$LnAF<-log(raw_df$Airfreight_Scheduled/raw_df$PPI_All_Commodities)
raw_df$LnR<-log(raw_df$Rail_Transportation/raw_df$PPI_All_Commodities)
raw_df$LnIC<-log(raw_df$Total_Inventories/raw_df$PPI_All_Commodities)
raw_df$LnSea<-log(raw_df$Deep_Sea_Freight/raw_df$PPI_All_Commodities)
raw_df$LnIF<-log(raw_df$Inland_Water_Freight/raw_df$PPI_All_Commodities)
raw_df$LnTotal<-log(raw_df$TC_Aggregate_WPU30/raw_df$PPI_All_Commodities)
#Without deflation
# raw_df$LnTL<-log(raw_df$Trucking_LD_Truckload)
# raw_df$LnLTL<-log(raw_df$Trucking_LD_LTL)
# raw_df$LnLocal<-log(raw_df$Trucking_Local)
# raw_df$LnAF<-log(raw_df$Airfreight_Scheduled)
# raw_df$LnR<-log(raw_df$Rail_Transportation)
# raw_df$LnIC<-log(raw_df$Total_Inventories)
# raw_df$LnSea<-log(raw_df$Deep_Sea_Freight)


#ADF Testing====
#All series should be I(1): non-stationary in levels, stationary in first differences

#Testing in levels
adf_TL<-ur.df(raw_df$Trucking_LD_Truckload, type = "drift",lags=4)
adf_LTL<-ur.df(raw_df$Trucking_LD_LTL, type = "drift",lags=4)
adf_Local<-ur.df(raw_df$Trucking_Local, type = "drift",lags=4)
adf_AF<-ur.df(raw_df$Airfreight_Scheduled, type = "drift",lags=4)
adf_R<-ur.df(raw_df$Rail_Transportation, type = "drift",lags=4)
adf_IC<-ur.df(raw_df$Total_Inventories, type = "drift",lags=4)
adf_Sea<-ur.df(raw_df$Deep_Sea_Freight, type = "drift",lags=4)
adf_IF<-ur.df(raw_df$Inland_Water_Freight, type = "drift",lags=4)
summary(adf_TL)
summary(adf_LTL)
summary(adf_Local)
summary(adf_AF)
summary(adf_R)
summary(adf_IC)
summary(adf_Sea)
summary(adf_IF)


diff_TL<-100*diff(raw_df$LnTL)
diff_LTL<-100*diff(raw_df$LnLTL)
diff_Local<-100*diff(raw_df$LnLocal)
diff_AF<-100*diff(raw_df$LnAF)
diff_R<-100*diff(raw_df$LnR)
diff_IC<-100*diff(raw_df$LnIC)
diff_Sea<-100*diff(raw_df$LnSea)
diff_IF<-100*diff(raw_df$LnIF)

#Testing in first differences
adf_diff_TL<-ur.df(diff_TL, type = "drift",lags=4)
adf_diff_LTL<-ur.df(diff_LTL, type = "drift",lags=4)
adf_diff_Local<-ur.df(diff_Local, type = "drift",lags=4)
adf_diff_AF<-ur.df(diff_AF, type = "drift",lags=4)
adf_diff_R<-ur.df(diff_R, type = "drift",lags=4)
adf_diff_IC<-ur.df(diff_IC, type = "drift",lags=4)
adf_diff_Sea<-ur.df(diff_Sea, type = "drift",lags=4)
adf_diff_IF<-ur.df(diff_IF, type = "drift",lags=4)

summary(adf_diff_TL)
summary(adf_diff_LTL)
summary(adf_diff_Local)
summary(adf_diff_AF)
summary(adf_diff_R)
summary(adf_diff_IC)
summary(adf_diff_Sea)
summary(adf_diff_IF)

#Everything is I(1), good to move forward with analysis

#Cointegration ====
#TL
TL_df<-data.frame(
  TL = raw_df$LnTL,
  IC = raw_df$LnIC
)
lag_selection <- VARselect(TL_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_TL<-ca.jo(TL_df,type="trace",ecdet="trend",K=5)
summary(johansen_test_TL) 
johansen_test_TL<-ca.jo(TL_df,type="trace",ecdet="trend",K=2)
summary(johansen_test_TL) 
#No Relationship
#VECM_TL <- VECM(TL_df, lag=10, r=1,
#                     estim = "ML",
#                     LRinclude = "const"
#)
#summary(VECM_TL)
#LTL
LTL_df<-data.frame(
  LTL = raw_df$LnLTL,
  IC = raw_df$LnIC
)
lag_selection <- VARselect(LTL_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_LTL<-ca.jo(LTL_df,type="trace",ecdet="trend", K=5)
summary(johansen_test_LTL)
johansen_test_LTL<-ca.jo(LTL_df,type="trace",ecdet="trend", K=3)
summary(johansen_test_LTL) 

VECM_LTL <- VECM(LTL_df, lag=2,
                 r=1, estim = "ML",
                 LRinclude = "trend")
summary(VECM_LTL)
VECM_LTL_alt <- VECM(LTL_df, lag=0,
                     r=1, estim = "ML",
                     LRinclude = "trend")
summary(VECM_LTL_alt)
#LTL Weak Exogeneity
alrtest(johansen_test_LTL, A = matrix(c(0,1), nrow=2), r=1)
#IC Weak Exogeneity
alrtest(johansen_test_LTL, A = matrix(c(1,0), nrow=2), r=1)
fit <- cajorls(johansen_test_LTL, r=1)
Box.test(residuals(fit$rlm)[,"LTL.d"], lag=12, type="Ljung-Box")
Box.test(residuals(fit$rlm)[,"IC.d"], lag=12, type="Ljung-Box")
#Local
Local_df<-data.frame(
  Local = raw_df$LnLocal,
  IC = raw_df$LnIC
)
lag_selection <- VARselect(Local_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_Local<-ca.jo(Local_df,type="trace",ecdet="trend", K=6)
summary(johansen_test_Local)
johansen_test_Local<-ca.jo(Local_df,type="trace",ecdet="trend", K=2)
summary(johansen_test_Local)

#R
R_df<-data.frame(
  R=raw_df$LnR,
  IC=raw_df$LnIC
)
lag_selection <- VARselect(R_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_R<-ca.jo(R_df,type="trace",ecdet="trend", K=3)
summary(johansen_test_R)

#no cointegration
#AF
AF_df<-data.frame(
  AF = raw_df$LnAF,
  IC = raw_df$LnIC
)
lag_selection <- VARselect(AF_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_AF<-ca.jo(AF_df,type="trace",ecdet="trend", K=6)
summary(johansen_test_AF)
johansen_test_AF<-ca.jo(AF_df,type="trace",ecdet="trend", K=2)
summary(johansen_test_AF)
#Sea
Sea_df<-data.frame(
  Sea = raw_df$LnSea,
  IC = raw_df$LnIC
)
lag_selection <- VARselect(Sea_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_Sea<-ca.jo(Sea_df,type="trace",ecdet="trend", K=6)
summary(johansen_test_Sea)
#IF
IF_df<-data.frame(
  IF = raw_df$LnIF,
  IC = raw_df$LnIC
)
lag_selection <- VARselect(IF_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_IF<-ca.jo(IF_df,type="trace",ecdet="trend", K=2)
summary(johansen_test_IF)
#nothing
#Full System
Full_df<-data.frame(
  TL=raw_df$LnTL,
  LTL=raw_df$LnLTL,
  Local=raw_df$LnLocal,
  R=raw_df$LnR,
  AF=raw_df$LnAF,
  Sea=raw_df$LnSea,
  IC=raw_df$LnIC
)
lag_selection <- VARselect(Full_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_Full<-ca.jo(Full_df,type="trace",ecdet="trend", K=1)
summary(johansen_test_Full)
VECM_Full <- VECM(Full_df, lag=3,
                  r=1, estim = "ML",
                  LRinclude = "trend")
summary(VECM_Full)
#Total System
total_df<-data.frame(
  Total=raw_df$LnTotal,
  IC=raw_df$LnIC
)
total_df <- na.omit(total_df)
lag_selection <- VARselect(total_df, lag.max = 10, type = "trend")
print(lag_selection)
johansen_test_total<-ca.jo(total_df,type="trace",ecdet="trend", K=5)
summary(johansen_test_total)

#Dummies ====
#robustness only
#ca.jo critical values aren't adjusted for dumvar, so these are indicative
raw_df$Date<-as.Date(paste0(raw_df$Year,"-",sprintf("%02d",raw_df$Month),"-01"))
raw_df$pandemic<-as.numeric(raw_df$Date >= as.Date("2020-01-01") & raw_df$Date <= as.Date("2022-12-01"))
raw_df$crisis<-as.numeric(raw_df$Date >= as.Date("2008-01-01") & raw_df$Date <= as.Date("2009-12-01"))

#drop any dummy with no variation
build_dumvar<-function(df, cols){
  m<-as.matrix(df[,cols,drop=FALSE])
  keep<-apply(m,2,function(x) length(unique(x))>1)
  if(!any(keep)) return(NULL)
  m[,keep,drop=FALSE]
}

dum_pan<-build_dumvar(raw_df,"pandemic")
dum_both<-build_dumvar(raw_df,c("pandemic","crisis"))

#TL with dummies
johansen_test_TL_pan<-ca.jo(TL_df,type="trace",ecdet="trend",K=5,dumvar=dum_pan)
summary(johansen_test_TL_pan)
johansen_test_TL_both<-ca.jo(TL_df,type="trace",ecdet="trend",K=5,dumvar=dum_both)
summary(johansen_test_TL_both)

#LTL with dummies
johansen_test_LTL_pan<-ca.jo(LTL_df,type="trace",ecdet="trend",K=3,dumvar=dum_pan)
summary(johansen_test_LTL_pan)
johansen_test_LTL_both<-ca.jo(LTL_df,type="trace",ecdet="trend",K=3,dumvar=dum_both)
summary(johansen_test_LTL_both)

#Local with dummies
johansen_test_Local_pan<-ca.jo(Local_df,type="trace",ecdet="trend",K=6,dumvar=dum_pan)
summary(johansen_test_Local_pan)
johansen_test_Local_both<-ca.jo(Local_df,type="trace",ecdet="trend",K=6,dumvar=dum_both)
summary(johansen_test_Local_both)

#R with dummies
johansen_test_R_pan<-ca.jo(R_df,type="trace",ecdet="trend",K=3,dumvar=dum_pan)
summary(johansen_test_R_pan)
johansen_test_R_both<-ca.jo(R_df,type="trace",ecdet="trend",K=3,dumvar=dum_both)
summary(johansen_test_R_both)

#AF with dummies
johansen_test_AF_pan<-ca.jo(AF_df,type="trace",ecdet="trend",K=6,dumvar=dum_pan)
summary(johansen_test_AF_pan)
johansen_test_AF_both<-ca.jo(AF_df,type="trace",ecdet="trend",K=6,dumvar=dum_both)
summary(johansen_test_AF_both)

#Sea with dummies
johansen_test_Sea_pan<-ca.jo(Sea_df,type="trace",ecdet="trend",K=6,dumvar=dum_pan)
summary(johansen_test_Sea_pan)
johansen_test_Sea_both<-ca.jo(Sea_df,type="trace",ecdet="trend",K=6,dumvar=dum_both)
summary(johansen_test_Sea_both)

#IF with dummies
johansen_test_IF_pan<-ca.jo(IF_df,type="trace",ecdet="trend",K=2,dumvar=dum_pan)
summary(johansen_test_IF_pan)
johansen_test_IF_both<-ca.jo(IF_df,type="trace",ecdet="trend",K=2,dumvar=dum_both)
summary(johansen_test_IF_both)

#Full system with dummies
johansen_test_Full_pan<-ca.jo(Full_df,type="trace",ecdet="trend",K=3,dumvar=dum_pan)
summary(johansen_test_Full_pan)
johansen_test_Full_both<-ca.jo(Full_df,type="trace",ecdet="trend",K=3,dumvar=dum_both)
summary(johansen_test_Full_both)

#Total with dummies
#total_df got na.omit'd so the dummy rows need to match it
total_keep<-!is.na(raw_df$LnTotal) & !is.na(raw_df$LnIC)
dum_pan_total<-build_dumvar(raw_df[total_keep,],"pandemic")
dum_both_total<-build_dumvar(raw_df[total_keep,],c("pandemic","crisis"))
johansen_test_total_pan<-ca.jo(total_df,type="trace",ecdet="trend",K=5,dumvar=dum_pan_total)
summary(johansen_test_total_pan)
johansen_test_total_both<-ca.jo(total_df,type="trace",ecdet="trend",K=5,dumvar=dum_both_total)
summary(johansen_test_total_both)


#Only dummy cointegrations to run
#LTL, pandemic + crisis (also check pandemic-only, both clear)
VECM_LTL_pan <- VECM(LTL_df, lag=2, r=1, estim="ML", LRinclude="trend", exogen=dum_pan)
summary(VECM_LTL_pan)
VECM_LTL_both <- VECM(LTL_df, lag=2, r=1, estim="ML", LRinclude="trend", exogen=dum_both)
summary(VECM_LTL_both)

#Local, pandemic + crisis only (pandemic-alone fell short at 24.96)
VECM_Local_both <- VECM(Local_df, lag=5, r=1, estim="ML", LRinclude="trend", exogen=dum_both)
summary(VECM_Local_both)

#IF, pandemic and pandemic + crisis (K=2 confirmed for both)
VECM_IF_pan <- VECM(IF_df, lag=1, r=1, estim="ML", LRinclude="trend", exogen=dum_pan)
summary(VECM_IF_pan)
VECM_IF_both <- VECM(IF_df, lag=1, r=1, estim="ML", LRinclude="trend", exogen=dum_both)
summary(VECM_IF_both)


#stop logging
sink()