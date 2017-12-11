RelDiffWithYesterday <- function(Data, minData=0.01,na.rm=F){
# RelDiff <- RelDiffWithYesterday(Preise,MinPreis);
# Transformation  von einem Daten-Vektor mit der Zeit in Relative Differenzen
# INPUT
# Data 							Vektor oder matrix tageweise in den Spalten. [1:n,1:d] 
# OPTIONAL
# minData       Data = max(Data,MinData), default MinData = 0.01
# na.rm					Set NaNs with zero or not?
# 
# OUTPUT
# RelDiff(1:n,1:d)  
 
# MT 2016
  if(is.matrix(Data)){
		nDays =dim(Data)[1]
		nQuotations=dim(Data)[1]
	}else{
		nQuotations = 1
		nDays <-length(Data)
	}


# Gestern = [1 1:(AnzTage-1)]';
	yesterday = c(1,1:(nDays-1)) # index von Gestern

# GestrigeWert = Wert(Gestern,:);
 if (nQuotations>1){
	 yData = Data[yesterday,]
 }else{
   yData = Data[yesterday];
 } # end  if (nQuotations>1)

# 
# Nenner  = (Preise+GestrigePreise);
	deno <- Data+yData

# NennerZeroInd = find(Nenner<MinPreis);
# Nenner(NennerZeroInd) = nan;  % underflow verhindern
if(na.rm){
		ind=which(deno<minData,arr.ind=T)
}else{
		deno[deno<minData]<-NaN;
}
# 
# RelDiff = 2* (Preise - GestrigePreise)./Nenner;

	relDiff <- 2*(Data - yData) / deno
if(na.rm){
	if(length(ind)>0)
		relDiff[ind]=0
}

 	return(relDiff)

 }

