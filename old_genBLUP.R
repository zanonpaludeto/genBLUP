genBLUP <- function(data, varResp, treatment = c("Prog","Clone"), plotType = c("LP","STP"), dominance = TRUE, 
                    fixed = c("Rep","Proc"), random = c("Rep","Proc","Control"), 
                    genPar_digits, otimizeSelection = FALSE, maxIndProgeny = NULL, 
                    maxProgenyBlock = NULL, excludeControl = NULL, excludeCod = NULL, directory = NULL){
  
  # loading packages --------------------------------------------------------
  
  pacman::p_load(future,sommer,stringr,tidyverse,ggroups)
  
  # stops and warnings ------------------------------------------------------
  
  if(!treatment %in% c("Prog","Clone")){
    stop("ERROR: The argument 'treatment' must be 'Prog' or 'Clone'")
  }  
  
  if(!plotType %in% c("LP","STP")){
    stop("ERROR: The argument 'plotType' must be 'LP' or 'STP'")
  } 
  
  if(!varResp %in% names(data)){
    stop("ERROR: varResp non-existent, please check your data")
  }
  
  if(is.null(directory)){
    warning("Directory is 'NULL', so no outputs were created")
  }
  
  # exploratory analysis ----------------------------------------------------
  
  resp <- data[,which(names(data)==varResp)]
  data$resp <- resp
  
  expAnalysis <- function(data){
    NAs <- sum(is.na(data))
    data <- data[!is.na(data)]
    exp <- round(c(mean(data),sd(data),min(data),max(data),length(data)), genPar_digits)
    exp <- round(c(exp,(exp[2]/exp[1])*100), genPar_digits)
    if (any(NAs)){
      exp <- c(exp,NAs)} else {exp <- c(exp,0)}
    exp <- round(c(exp,((exp[5]/(exp[5]+exp[7]))*100)), genPar_digits)
    names(exp) <- c("Mean", "sd", "Min", "Max", "Vivas", "CV%", "NAs","Sob%")
    return(exp)
  }
  respMeans <- expAnalysis(resp)
  
  groupMeans <- function(data, group){
    sub_group <- subset(data, select=group)
    names <- c(names(sub_group),"Mean", "sd", "Min", "Max", "Vivas", "CV%", "NAs","Sob%")
    names(sub_group) <- "grp"
    grpMatrix <- as.matrix(aggregate(resp, by=list(sub_group$grp), FUN=expAnalysis))
    grpMeans <- setNames(as.data.frame(grpMatrix),names)
  }
  treatMeans <- groupMeans(data,treatment)
  repMeans <- groupMeans(data,"Rep")
  
  if(any(random=="Proc")||any(fixed=="Proc")){
    procMeans <- groupMeans(data,"Proc")
  }
  
  # factorizing variables and treeCheck -------------------------------------
  
  if(plotType=="LP"){
    if(length(random)>=1){
      fct <- c(treatment,"Parc",fixed,random) 
      data[fct] <- lapply(data[fct], factor)
    }
    if(is.null(random)){
      fct <- c(treatment,"Parc",fixed) 
      data[fct] <- lapply(data[fct], factor)
    }}
  
  if(plotType=="STP"){
    
    treeCheck <- function(data){
      
      arvList <- list()
      
      for(i in levels(factor(data$Rep))){
        if(treatment=="Clone"){
          arvList[i] <- as.data.frame(ave(data[data$Rep==i,]$Clone==data[data$Rep==i,]$Clone, 
                                          data[data$Rep==i,]$Clone, FUN=cumsum))
        }else{
          arvList[i] <- as.data.frame(ave(data[data$Rep==i,]$Prog==data[data$Rep==i,]$Prog, 
                                          data[data$Rep==i,]$Prog, FUN=cumsum))
        }}
      
      Arv <- unlist(arvList)
      data$Arv <- Arv
      return(data)
    }
    
    data <- treeCheck(data)
    
    if(length(random)>=1){
      fct <- c(treatment,fixed,random) 
      data[fct] <- lapply(data[fct], factor)
    }
    
    if(is.null(random)){
      fct <- c(treatment,fixed) 
      data[fct] <- lapply(data[fct], factor)
    }}
  
  # pedigree matrix ---------------------------------------------------------
  
  if(plotType=="LP"){  
    data$Ind <- paste0(data$Rep,sep = ".",data[,names(data) %in% treatment],sep=".",data$Parc,sep=".",data$Arv)
  }
  if(plotType=="STP"){
    data$Ind <- paste0(data$Rep,sep = ".",data[,names(data) %in% treatment],sep=".",data$Arv)  
  }
  
  if(treatment=="Prog"){
    
    makeA_op <- function(data){
      
      ped <- data.frame(data,sire=0) %>% setNames(c("Ind","dam","sire")) %>% relocate(dam, .after = sire) %>% arrange(dam)
      if(any(duplicated(ped[,1])==TRUE)){
        stop("ERROR: There is duplicated individuals, please check your data")
      }
      exc <- length(levels(factor(ped$dam)))
      A <- data.frame(levels(factor(ped$dam)),0,0) %>% setNames(c("Ind","dam","sire")) %>% rbind(.,ped) %>% 
        ggroups::buildA(.) %>% .[-1:-exc,-1:-exc]
      return(A)
    }
    
    A <- makeA_op(data[,c("Ind","Prog")])
    data <- data[match(rownames(A),data$Ind),]
    
    resp <- data[,which(names(data)==varResp)]
    data$resp <- resp
    
    if(identical(data$Ind,rownames(A))==FALSE){
      stop("ERROR: Pedigree and data orders doesn't match, something strange happened")
      
    }
    
    # dominance matrix --------------------------------------------------------
    
    if(dominance==TRUE){
      
      D_matrix <- data.frame(data[,c("Ind","Prog")],sire=0) %>% setNames(c("Ind","dam","sire")) %>% relocate(dam, .after = sire) %>% 
        arrange(dam) %>% ggroups::buildD(.,A)
      if(any(duplicated(rownames(D_matrix))==TRUE)){
        stop("ERROR: There is duplicated individuals, please check your data")
      }
      D_matrix <- D_matrix[!is.na(data$resp),!is.na(data$resp)]
    }
    
    A <- A[!is.na(data$resp),!is.na(data$resp)]
    
  }
  
  # statistical modeling ----------------------------------------------------
  
  
  future::plan(multisession)
  
  if(plotType=="LP"){
    
    if(length(fixed)==1&length(random)==0){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F),c("Ind","Treat","Parc","resp","Fixed1"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1,
                                  random= ~ Treat + Parc,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1,
                                 random= ~ vs(Ind,Gu=A) + Parc,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + (1|Parc), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + (1|Treat), data=df)
    }
    
    if(length(fixed)==1&length(random)==1){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","Parc","resp","Fixed1","Random1"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1,
                                  random= ~ Treat + Parc + Random1,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1,
                                 random= ~ vs(Ind,Gu=A) + Parc + Random1,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc + Random1,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc) + (1|Random1), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + (1|Parc) + (1|Random1), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Random1), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc), data=df)
    }
    
    if(length(fixed)==1&length(random)==2){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","Parc","resp","Fixed1","Random1","Random2"))
      
      if(any(random=="Proc")){
        if(df$Random2==data$Proc){
          df$Random2 <- df$Random1
          df$Random1 <- data$Proc
        }
      }
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1,
                                  random= ~ Treat + Parc + Random1 + Random2,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1,
                                 random= ~ vs(Ind,Gu=A) + Parc + Random1 + Random2,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc + Random1 + Random2,
                                   rcov=~units,
                                   data=df)
          
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc) + (1|Random1) + (1|Random2), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + (1|Parc) + (1|Random1) + (1|Random2), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Random1) + (1|Random2), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc)+ (1|Random2), data=df)
      mSig4 <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc) + (1|Random1), data=df)
    }
    
    if(length(fixed)==2&length(random)==0){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","Parc","resp","Fixed1","Fixed2"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                  random= ~ Treat + Parc,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                 random= ~ vs(Ind,Gu=A) + Parc,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Parc), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat), data=df)
    }
    
    if(length(fixed)==2&length(random)==1){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","Parc","resp","Fixed1","Fixed2","Random1"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                  random= ~ Treat + Parc + Random1,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                 random= ~ vs(Ind,Gu=A) + Parc + Random1,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc + Random1,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc) + (1|Random1), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Parc) + (1|Random1), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Random1), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc), data=df)
    }
    
    if(length(fixed)==2&length(random)==2){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","Parc","resp","Fixed1","Fixed2","Random1","Random2"))
      
      if(any(random=="Proc")){
        if(df$Random2==data$Proc){
          df$Random2 <- df$Random1
          df$Random1 <- data$Proc
        }
      }
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                  random= ~ Treat + Parc + Random1 + Random2,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                 random= ~ vs(Ind,Gu=A) + Parc + Random1 + Random2,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc + Random1 + Random2,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc) + (1|Random1) + (1|Random2), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Parc) + (1|Random1) + (1|Random2), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Random1) + (1|Random2), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc) + (1|Random2), data=df)
      mSig4 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc) + (1|Random1), data=df)
    }
    
    if(length(fixed)==3&length(random)==0){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","Parc","resp","Fixed1","Fixed2","Fixed3"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random= ~ Treat + Parc,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                 random= ~ vs(Ind,Gu=A) + Parc,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Parc), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat), data=df)
    }
    
    if(length(fixed)==3&length(random)==1){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","Parc","resp","Fixed1","Fixed2","Fixed3","Random1"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random= ~ Treat + Parc + Random1,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                 random= ~ vs(Ind,Gu=A) + Parc + Random1,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc + Random1,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc) + (1|Random1), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Parc) + (1|Random1), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Random1), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc), data=df)
    }
    
    if(length(fixed)==3&length(random)==2){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$Parc,data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","Parc","resp","Fixed1","Fixed2","Fixed3","Random1","Random2"))
      
      if(any(random=="Proc")){
        if(df$Random2==data$Proc){
          df$Random2 <- df$Random1
          df$Random1 <- data$Proc
        }
      }
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random= ~ Treat + Parc + Random1 + Random2,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                 random= ~ vs(Ind,Gu=A) + Parc + Random1 + Random2,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                   random= ~ vs(Ind,Gu=D_matrix) + Parc + Random1 + Random2,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc) + (1|Random1) + (1|Random2), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Parc) + (1|Random1) + (1|Random2), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Random1) + (1|Random2), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc) + (1|Random2), data=df)
      mSig4 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc) + (1|Random1), data=df)
    }
  }
  
  ###################### Single-tree plot
  
  if(plotType=="STP"){
    
    if(length(fixed)==1&length(random)==0){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1,
                                  random= ~ Treat,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1,
                                 random= ~ vs(Ind,Gu=A),
                                 rcov=~units,
                                 data=df)
        
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1,
                                   random= ~ vs(Ind,Gu=D_matrix),
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat), data=df)
      mSig1 <- lm(resp ~ Fixed1, data=df)
    }
    
    if(length(fixed)==1&length(random)==1){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1","Random1"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1,
                                  random= ~ Treat + Random1,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1,
                                 random= ~ vs(Ind,Gu=A) + Random1,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1,
                                   random= ~ vs(Ind,Gu=D_matrix) + Random1,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Random1), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + (1|Random1), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + (1|Treat), data=df)
    }
    
    if(length(fixed)==1&length(random)==2){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1","Random1","Random2"))
      
      if(any(random=="Proc")){
        if(df$Random2==data$Proc){
          df$Random2 <- df$Random1
          df$Random1 <- data$Proc
        }
      }
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1,
                                  random= ~ Treat + Random1 + Random2,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1,
                                 random= ~ vs(Ind,Gu=A) + Random1 + Random2,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1,
                                   random= ~ vs(Ind,Gu=D_matrix) + Random1 + Random2,
                                   rcov=~units,
                                   data=df)
          
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Random1) + (1|Random2), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + (1|Random1) + (1|Random2), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Random2), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Random1), data=df)
    }
    
    if(length(fixed)==2&length(random)==0){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1","Fixed2"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                  random= ~ Treat,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                 random= ~ vs(Ind,Gu=A),
                                 rcov=~units,
                                 data=df)
        
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                   random= ~ vs(Ind,Gu=D_matrix),
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat), data=df)
      mSig1 <- lm(resp ~ Fixed1 + Fixed2, data=df)
    }
    
    if(length(fixed)==2&length(random)==1){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1","Fixed2","Random1"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                  random= ~ Treat + Random1,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                 random= ~ vs(Ind,Gu=A) + Random1,
                                 rcov=~units,
                                 data=df)
        
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1, Fixed2,
                                   random= ~ vs(Ind,Gu=D_matrix) + Random1,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Random1), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Random1), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat), data=df)
    }
    
    if(length(fixed)==2&length(random)==2){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1","Fixed2","Random1","Random2"))
      
      if(any(random=="Proc")){
        if(df$Random2==data$Proc){
          df$Random2 <- df$Random1
          df$Random1 <- data$Proc
        }
      }
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                  random= ~ Treat + Random1 + Random2,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                 random= ~ vs(Ind,Gu=A) + Random1 + Random2,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2,
                                   random= ~ vs(Ind,Gu=D_matrix) + Random1 + Random2,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Random1) + (1|Random2), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Random1) + (1|Random2), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Random2), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Random1), data=df)
    }
    
    if(length(fixed)==3&length(random)==0){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1","Fixed2","Fixed3"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random= ~ Treat,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                 random= ~ vs(Ind,Gu=A),
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                   random= ~ vs(Ind,Gu=D_matrix),
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat), data=df)
      mSig1 <- lm(resp ~ Fixed1 + Fixed2, data=df)
    }
    
    if(length(fixed)==3&length(random)==1){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1","Fixed2","Fixed3","Random1"))
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random= ~ Treat + Random1,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                 random= ~ vs(Ind,Gu=A) + Random1,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                   random= ~ vs(Ind,Gu=D_matrix) + Random1,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Random1), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Random1), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat), data=df)
    }
    
    if(length(fixed)==3&length(random)==2){
      
      df <- setNames(data.frame(data$Ind,data[,names(data) %in% treatment],data$resp,data[,names(data) %in% fixed],data[,names(data) %in% random],
                                stringsAsFactors = F), c("Ind","Treat","resp","Fixed1","Fixed2","Fixed3","Random1","Random2"))
      
      if(any(random=="Proc")){
        if(df$Random2==data$Proc){
          df$Random2 <- df$Random1
          df$Random1 <- data$Proc
        }
      }
      
      Progmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random= ~ Treat + Random1 + Random2,
                                  rcov=~units,
                                  data=df)
      if(treatment=="Prog"){
        Amodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                 random= ~ vs(Ind,Gu=A) + Random1 + Random2,
                                 rcov=~units,
                                 data=df)
        if(dominance==TRUE){
          Dmodel %<-% sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
                                   random= ~ vs(Ind,Gu=D_matrix) + Random1 + Random2,
                                   rcov=~units,
                                   data=df)
        }}
      mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Random1) + (1|Random2), data=df)
      mSig1 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Random1) + (1|Random2), data=df)
      mSig2 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Random2), data=df)
      mSig3 <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Random1), data=df)
    }
  }
  
  if(treatment=="Prog"){
    
    if(dominance==TRUE){
      model.list <- lapply(c("Amodel","Progmodel","Dmodel"),get,envir=sys.frame(sys.parent(0)))
    }else{
      model.list <- lapply(c("Amodel","Progmodel"),get,envir=sys.frame(sys.parent(0))) 
    }
    
    mAdd <- model.list[[1]]
    mProg <- model.list[[2]]
    if(dominance==TRUE){
      mDom <- model.list[[3]]
    }
    
    # accuracy by PEV ---------------------------------------------------------
    
    dfacc <- setNames(data.frame(sqrt(1-(diag(mProg$PevU$`Treat`$resp)/mProg$sigmaVector[1]))),"Accuracy") %>% 
      rownames_to_column(.,var="Prog") %>% mutate(Prog = str_replace_all(Prog, 'Treat', ""))
    accProg <- mean(dfacc$Accuracy, na.rm=T)
    
    
    dfaccInd <- setNames(data.frame(sqrt(1-(diag(mAdd$PevU$`u:Ind`$resp)/mAdd$sigmaVector[1]))),"Accuracy") %>% 
      rownames_to_column(.,var="Ind")
    
    if(plotType=="LP"){ 
      dfaccInd <- separate(data=dfaccInd, Ind, c("Block","Progeny","Plot","Tree"), sep="\\.")
    }
    if(plotType=="STP"){
      dfaccInd <- separate(data=dfaccInd, Ind, c("Block","Progeny","Tree"), sep="\\.") 
    }
    
    accInd <- mean(dfaccInd$Accuracy, na.rm=T)
    
    # Genetic Parameters
    
    vA <- mAdd$sigmaVector[1]
    SEvA <- sqrt(diag(mAdd$sigmaSE))[1]
    
    # Average Mean dn CVgi  
    Mean <- mean(df$resp,na.rm=T)
    CVgi <- sqrt(vA)/Mean*100
    
    # Linear plot
    if(plotType=="LP"){
      nArv <- length(levels(factor(data$Arv)))
      nRep <- length(levels(factor(data$Rep)))
      
      if(length(random)==0){
        vParc <- mAdd$sigmaVector[2]
        vE <- mAdd$sigmaVector[3]
        SEvParc <- sqrt(diag(mAdd$sigmaSE))[2]
        SEvE <- sqrt(diag(mAdd$sigmaSE))[3]
        vPhen <- vA + vParc + vE
        c2Parc <- vParc/vPhen
        h2aSE <- sommer::vpredict(mAdd, h2a ~ (V1) / (V1+V2+V3))
        h2dSE <- sommer::vpredict(mAdd, h2d ~ (0.75*V1) / (0.75*V1+V3))
        h2m <- (0.25*vA) / (0.25*vA+(vParc/nRep)+vE/(nRep*nArv))
        CVe = (sqrt((0.75*vA+vE)/nArv+vParc))/Mean*100
        
        # genPar
        genParNames<- c("vA","vParc","vE","vPhen","h2a","h2d","h2m","c2Parc","accProg","accInd","CVgi%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2aSE$Estimate,h2dSE$Estimate,h2m,c2Parc,accProg,accInd,CVgi,CVe,Mean), 
                                   SE=c(SEvA,SEvParc,SEvE,NA,h2aSE$SE,h2dSE$SE,matrix(NA,nrow=7,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }
      if(length(random)==1){
        vParc <- mAdd$sigmaVector[2]
        vRdm1 <- mAdd$sigmaVector[3]
        vE <- mAdd$sigmaVector[4]
        SEvParc <- sqrt(diag(mAdd$sigmaSE))[2]
        SEvRdm1 <- sqrt(diag(mAdd$sigmaSE))[3]
        SEvE <- sqrt(diag(mAdd$sigmaSE))[4]
        vPhen <- vA + vParc + vRdm1 + vE
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        h2aSE <- sommer::vpredict(mAdd, h2a ~ (V1) / (V1+V2+V3+V4))
        h2dSE <- sommer::vpredict(mAdd, h2d ~ (0.75*V1) / (0.75*V1+V4))
        CVe = (sqrt((0.75*vA+vE)/nArv+vParc))/Mean*100
        
        # genPar
        genParNames <- c("vA","vParc","vRdm1","vE","vPhen","h2a","h2d","c2Parc","c2Rdm1","accProg","accInd","CVgi%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vE,vPhen,h2aSE$Estimate,h2dSE$Estimate,c2Parc,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                   SE=c(SEvA,SEvParc,SEvRdm1,SEvE,NA,h2aSE$SE,h2dSE$SE,matrix(NA,nrow=7,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }
      if(length(random)==2){
        vParc <- mAdd$sigmaVector[2]
        vRdm1 <- mAdd$sigmaVector[3]
        vRdm2 <- mAdd$sigmaVector[4]
        vE <- mAdd$sigmaVector[5]
        SEvParc <- sqrt(diag(mAdd$sigmaSE))[2]
        SEvRdm1 <- sqrt(diag(mAdd$sigmaSE))[3]
        SEvRdm1 <- sqrt(diag(mAdd$sigmaSE))[4]
        SEvE <- sqrt(diag(mAdd$sigmaSE))[5]
        vPhen <- vA + vParc + vRdm1 + vRdm2 + vE
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        h2aSE <- sommer::vpredict(mAdd, h2a ~ (V1) / (V1+V2+V3+V4+V5))
        h2dSE <- sommer::vpredict(mAdd, h2d ~ (0.75*V1) / (0.75*V1+V5))
        CVe = (sqrt((0.75*vA+vE)/nArv+vParc))/Mean*100
        
        # genPar
        genParNames <- c("vA","vParc","vRdm1","vRdm2","vE","vPhen","h2a","h2d","c2Parc","c2Rdm1","c2Rdm2","accProg","accInd","CVgi%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2aSE$Estimate,h2dSE$Estimate,c2Parc,c2Rdm1,c2Rdm2,accProg,accInd,CVgi,CVe,Mean), 
                                   SE=c(SEvA,SEvParc,SEvRdm1,SEvRdm2,SEvE,NA,h2aSE$SE,h2dSE$SE,matrix(NA,nrow=8,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }}
    # Single tree plot
    if(plotType=="STP"){
      nRep <- length(levels(factor(data$Rep)))
      
      
      if(length(random)==0){
        vE <- mAdd$sigmaVector[2]
        SEvE <- sqrt(diag(mAdd$sigmaSE))[2]
        vPhen <- vA + vE
        h2aSE <- sommer::vpredict(mAdd, h2a ~ (V1) / (V1+V2))
        h2m <- (0.25*vA) / (0.25*vA+vE/(nRep))
        CVe = sqrt(0.75*vA+vE)/Mean*100
        
        # genPar
        genParNames <- c("vA","vE","vPhen","h2a","h2m","accProg","accInd","CVgi%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vA,vE,vPhen,h2aSE$Estimate,h2m,accProg,accInd,CVgi,CVe,Mean), 
                                   SE=c(SEvA,SEvE,NA,h2aSE$SE,matrix(NA,nrow=6,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }
      if(length(random)==1){
        vRdm1 <- mAdd$sigmaVector[2]
        vE <- mAdd$sigmaVector[3]
        SEvRdm1 <- sqrt(diag(mAdd$sigmaSE))[2]
        SEvE <- sqrt(diag(mAdd$sigmaSE))[3]
        vPhen <- vA + vRdm1 + vE
        c2Rdm1 <- vRdm1/vPhen
        h2aSE <- sommer::vpredict(mAdd, h2a ~ (V1) / (V1+V2+V3))
        CVe = sqrt(0.75*vA+vE)/Mean*100
        
        # genPar
        genParNames <- c("vA","vRdm1","vE","vPhen","h2a","c2Rdm1","accProg","accInd","CVgi%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2aSE$Estimate,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                   SE=c(SEvA,SEvRdm1,SEvE,NA,h2aSE$SE,matrix(NA,nrow=6,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }
      if(length(random)==2){
        vRdm1 <- mAdd$sigmaVector[2]
        vRdm2 <- mAdd$sigmaVector[3]
        vE <- mAdd$sigmaVector[4]
        SEvRdm1 <- sqrt(diag(mAdd$sigmaSE))[2]
        SEvRdm2 <- sqrt(diag(mAdd$sigmaSE))[3]
        SEvE <- sqrt(diag(mAdd$sigmaSE))[4]
        vPhen <- vA + vRdm1 + vRdm2 + vE
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        h2aSE <- sommer::vpredict(mAdd, h2a ~ (V1) / (V1+V2+V3+V4))
        CVe = sqrt(0.75*vA+vE)/Mean*100
        
        # genPar
        genParNames <- c("vA","vRdm1","vRdm2","vE","vPhen","h2a","c2Rdm1","c2Rdm2","accProg","accInd","CVgi%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2aSE$Estimate,c2Rdm1,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                   SE=c(SEvA,SEvRdm1,SEvE,NA,h2aSE$SE,matrix(NA,nrow=7,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }}
    
    # individual_and_progeny_BLUP ---------------------------------------------
    
    require(tidyr)
    
    # Creating Cod if it doesn't exists
    if(!("Cod" %in% colnames(data))){
      data$Cod = " "
    }
    
    # Progeny BLUP and overlapping generations
    progBLUP <- mProg$U$Treat$resp %>% data.frame() %>% setNames("BLUP") %>% rownames_to_column(.,"Progeny") %>% 
      mutate(Progeny = str_replace_all(Progeny, 'Treat', ""), BLUP = BLUP*2) %>% arrange(desc(BLUP))
    
        # Individual BLUP
    indBLUP <- mAdd$U$`u:Ind`$resp %>% data.frame %>% setNames("BLUP") %>% rownames_to_column(.,"Ind")
    
    # Phenotypic and cod data
    phenoCodData <- data %>% select(Ind,Cod,resp) %>% rename(.,f = resp)
    
    if(any(random=="Proc")||any(fixed=="Proc")){
      phenoCodData <- setNames(subset(data, select = c("Ind","Cod","Proc","resp")),c("Ind","Cod","Proc","f"))
    }
    
    dfblup <- merge(indBLUP,phenoCodData,by="Ind")
    
    if(dominance==TRUE){
      domBLUP <- mDom$U$`u:Ind`$resp %>% data.frame %>% setNames("d") %>% rownames_to_column(.,"Ind")
      dfblup <- merge(dfblup,domBLUP,by="Ind")
      dfblup$g <- dfblup$BLUP+dfblup$d
    }
    
    if(plotType=="LP"){ 
      dfblupInd <- separate(data=dfblup, col = Ind, into = c("Block","Progeny","Plot","Tree"), sep="\\.")
      overBLUP <- progBLUP %>% mutate(Block = 0, Plot = 0, Tree = 0, Cod = "Parent") %>% relocate(Block,Progeny,Plot,Tree,BLUP,Cod)
    }
    if(plotType=="STP"){
      dfblupInd <- separate(data=dfblup, col = Ind, into = c("Block","Progeny","Tree"), sep="\\.")
      overBLUP <- progBLUP %>% mutate(Block = 0, Tree = 0, Cod = "Parent") %>% relocate(Block,Progeny,Tree,BLUP,Cod)
    }
    
    # Provenance BLUP
    if(any(random=="Proc")){
      procBLUP <- mAdd$U$Random1$resp %>% data.frame() %>% setNames("procBLUP") %>% rownames_to_column(.,"Provenance") %>% 
        mutate(Provenance = str_replace_all(Provenance, 'Random1', "")) %>% arrange(desc(procBLUP))
      
      dfblupInd <- dfblupInd %>%  merge(.,procBLUP, by.x="Proc",by.y="Provenance")
      dfblupInd$BLUP <- dfblupInd$BLUP+dfblupInd$procBLUP
      dfblupInd <- dfblupInd[,-which(names(dfblupInd)%in%"procBLUP")]
    }
    
    # Organizing overlapping generation BLUP
    if(dominance){
      overBLUP <- dfblupInd %>% select (-c(d,f,g)) %>% rbind(.,overBLUP) %>% arrange(desc(BLUP))
    }else{
      overBLUP <- dfblupInd %>% select (-f) %>% rbind(.,overBLUP) %>% arrange(desc(BLUP))
    }
    
    # BLUP + mean and arrange
    
    dfblupInd$"u+a" <- dfblupInd$BLUP+Mean
    blupOrd <- arrange(dfblupInd,desc(dfblupInd$BLUP))
    
    # Absolute and percentual genetic gain
    gain <- matrix(nrow=nrow(blupOrd), ncol=1)
    
    for(i in 1:nrow(blupOrd)){
      gain[i,] <- mean(blupOrd$BLUP[1:i])
    }
    blupOrd$Gain <- gain
    blupOrd$"Gain%" <- (blupOrd$Gain/Mean)*100
    
    # effective_number_size_and_genetic_gain ----------------------------------
    
    dfNe <- matrix(nrow=nrow(blupOrd), ncol=1)
    
    for(i in 1:nrow(blupOrd)){
      varkf <- as.vector(table(blupOrd$Progeny[1:i]))
      dfNe[i,] <- var(varkf)
    }
    blupOrd$Np <- cumsum(!duplicated(blupOrd$Progeny))
    blupOrd$seq <- seq(1:nrow(blupOrd))
    blupOrd$Kf <- blupOrd$seq/blupOrd$Np
    dfNe[is.na(dfNe)] <- 0
    blupOrd$varkf <- dfNe
    blupOrd$Ne <- (4*blupOrd$Np*blupOrd$Kf)/(3+blupOrd$Kf+(blupOrd$varkf/blupOrd$Kf))
    blupOrd <- blupOrd[,!names(blupOrd) %in% c("Np","seq","Kf","varkf")]
    
    
    if(dominance==TRUE){
      if(plotType=="LP"&(any(random=="Proc")||any(fixed=="Proc"))){
        order <- c("Proc","Block","Progeny","Plot","Tree","f","BLUP","u+a","Gain","Gain%","Ne","d","g","Cod")
      }
      if(plotType=="LP"&((any(random=="Proc")||any(fixed=="Proc")))==FALSE){
        order <- c("Block","Progeny","Plot","Tree","f","BLUP","u+a","Gain","Gain%","Ne","d","g","Cod")
      }
      if(plotType=="STP"&(any(random=="Proc")||any(fixed=="Proc"))){
        order <- c("Proc","Block","Progeny","Tree","f","BLUP","u+a","Gain","Gain%","Ne","d","g","Cod")
      }
      if(plotType=="STP"&((any(random=="Proc")||any(fixed=="Proc")))==FALSE){
        order <- c("Block","Progeny","Tree","f","BLUP","u+a","Gain","Gain%","Ne","d","g","Cod")
      }
    }else{
      if(plotType=="LP"&(any(random=="Proc")||any(fixed=="Proc"))){
        order <- c("Proc","Block","Progeny","Plot","Tree","f","BLUP","u+a","Gain","Gain%","Ne","Cod")
      }
      if(plotType=="LP"&((any(random=="Proc")||any(fixed=="Proc")))==FALSE){
        order <- c("Block","Progeny","Plot","Tree","f","BLUP","u+a","Gain","Gain%","Ne","Cod")
      }
      if(plotType=="STP"&(any(random=="Proc")||any(fixed=="Proc"))){
        order <- c("Proc","Block","Progeny","Tree","f","BLUP","u+a","Gain","Gain%","Ne","Cod")
      }
      if(plotType=="STP"&((any(random=="Proc")||any(fixed=="Proc")))==FALSE){
        order <- c("Block","Progeny","Tree","f","BLUP","u+a","Gain","Gain%","Ne","Cod")
      }}
    blupOrdered <- blupOrd[,order]
    rownames(blupOrdered) <- NULL
    
    if(otimizeSelection==TRUE){
      
      rankBLUP <- blupOrdered
      if(!is.null(excludeControl)){
        rankBLUP <- rankBLUP[!rankBLUP$Progeny %in% excludeControl,]
      }
      if(!is.null(excludeCod)){
        rankBLUP <- rankBLUP[!rankBLUP$Cod %in% excludeCod,]
      }
      
      #Filtering using the maximum number of individuals per progeny admitted in the whole otimization
      rankBLUP$Np_csum <- ave(rankBLUP$Progeny==rankBLUP$Progeny, 
                              rankBLUP$Progeny, FUN=cumsum) 
      rankBLUP <- rankBLUP[!rankBLUP$Np_csum > maxIndProgeny,] 
      
      #Filtering using the maximum number of progenies admitted at the same block repetition
      rankBLUP$R_P <- paste0(rankBLUP$Block,rankBLUP$Progeny)
      rankBLUP$R_Pcsum <- ave(rankBLUP$R_P==rankBLUP$R_P, 
                              rankBLUP$R_P, FUN=cumsum) 
      rankBLUP <- rankBLUP[!rankBLUP$R_Pcsum > maxProgenyBlock,]
      
      rankBLUP <- rankBLUP[, !names(rankBLUP) %in% c("Np_csum","R_P","R_Pcsum")] 
      
      #Genetic Gain
      gain <- matrix(nrow=nrow(rankBLUP), ncol=1)
      
      for(i in 1:nrow(rankBLUP)){
        gain[i,] <- mean(rankBLUP$BLUP[1:i])
      }
      rankBLUP$Gain <- gain
      rankBLUP$"Gain%" <- (rankBLUP$Gain/mean(resp, na.rm=T))*100
      
      #Effective number size
      dfNe <- matrix(nrow=nrow(rankBLUP), ncol=1)
      
      for(i in 1:nrow(rankBLUP)){
        varkf <- as.vector(table(rankBLUP$Progeny[1:i]))
        dfNe[i,] <- var(varkf)
      }
      rankBLUP$Np <- cumsum(!duplicated(rankBLUP$Progeny))
      rankBLUP$seq <- seq(1:nrow(rankBLUP))
      rankBLUP$Kf <- rankBLUP$seq/rankBLUP$Np
      dfNe[is.na(dfNe)] <- 0
      rankBLUP$varkf <- dfNe
      rankBLUP$Ne <- (4*rankBLUP$Np*rankBLUP$Kf)/(3+rankBLUP$Kf+(rankBLUP$varkf/rankBLUP$Kf))
      otimizedBLUP <- rankBLUP[,!names(rankBLUP) %in% c("Np","seq","Kf","varkf")]
      rownames(otimizedBLUP) <- NULL
    }
    
    
    # final_output_list -------------------------------------------------------
    
    genParBLUP <- list()
    
    # exploratory_analysis -------------------------------------------------
    genParBLUP$expAnalysis$respMeans <- respMeans
    genParBLUP$expAnalysis$progMeans <- treatMeans
    genParBLUP$expAnalysis$repMeans <- repMeans
    
    # model_output_and_significance_tests -------------------------------------------------
    genParBLUP$Model$mmerModel <- mAdd
    genParBLUP$Model$mSig$fixedSig <- anova(mSig)
    genParBLUP$Model$mSig$randomSig <- lmerTest::ranova(mSig, reduce.terms = F)
    
    # genetic_parameters_accuracy_and_Blup -------------------------------------------------
    genParBLUP$genPar <- genPar
    genParBLUP$accuracyPEV$progPEV <- dfacc
    genParBLUP$accuracyPEV$indPEV <- dfaccInd
    genParBLUP$BLUP$progBLUP <- progBLUP
    genParBLUP$BLUP$indBLUP <- blupOrdered
    genParBLUP$BLUP$overBLUP <- overBLUP
    
    if(any(random=="Proc")){
      genParBLUP$BLUP$procBLUP <- dfblupProc
      genParBLUP$expAnalysis$procMeans <- procMeans
    }
    
    if(otimizeSelection==TRUE){
      genParBLUP$BLUP$otimizedBLUP <- otimizedBLUP
    }
  }
  
  if(treatment=="Clone"){
    
    mProg <- lapply("Progmodel",get,envir=sys.frame(sys.parent(0)))[[1]] 
    
    # accuracyPEV
    
    dfacc <- setNames(data.frame(sqrt(1-(diag(mProg$PevU$`Treat`$resp)/mProg$sigmaVector[1]))),"Accuracy") %>% 
      rownames_to_column(.,var="Clone") %>% mutate(Clone = str_replace_all(Clone, 'Treat', ""))
    accClone <- mean(dfacc$Accuracy, na.rm=T)
    
    # Genetic Parameters
    vG <- mProg$sigmaVector[1]
    SEvG <- sqrt(diag(mProg$sigmaSE))[1]
    
    # Average Mean and CVg  
    Mean <- mean(df$resp,na.rm=T)
    CVg <- sqrt(vG)/Mean*100
    
    # Linear plot
    if(plotType=="LP"){
      nArv <- length(levels(factor(data$Arv)))
      nRep <- length(levels(factor(data$Rep)))
      
      if(length(random)==0){
        
        vParc <- mProg$sigmaVector[2]
        vE <- mProg$sigmaVector[3]
        SEvParc <- sqrt(diag(mProg$sigmaSE))[2]
        SEvE <- sqrt(diag(mProg$sigmaSE))[3]
        vPhen <- vG + vParc + vE
        c2Parc <- vParc/vPhen
        h2GSE <- sommer::vpredict(mProg, h2G ~ (V1) / (V1+V2+V3))
        h2mc <- (vG) / (vG+(vParc/nRep)+vE/(nRep*nArv))
        CVe = (sqrt((3*vG+vE)/nArv+vParc))/Mean*100
        
        # genPar
        genParNames <- c("vG","vParc","vE","vPhen","h2G","h2mc","c2Parc","accClone","CVg%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vG,vParc,vE,vPhen,h2GSE$Estimate,h2mc,c2Parc,accClone,CVg,CVe,Mean), 
                                   SE=c(SEvG,SEvParc,SEvE,NA,h2GSE$SE,matrix(NA,nrow=6,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }
      if(length(random)==1){
        
        vParc <- mProg$sigmaVector[2]
        vRdm1 <- mProg$sigmaVector[3]
        vE <- mProg$sigmaVector[4]
        SEvParc <- sqrt(diag(mProg$sigmaSE))[2]
        SEvRdm1 <- sqrt(diag(mProg$sigmaSE))[3]
        SEvE <- sqrt(diag(mProg$sigmaSE))[4]
        vPhen <- vG + vParc + vRdm1 + vE
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        h2GSE <- sommer::vpredict(mProg, h2g ~ (V1) / (V1+V2+V3+V4))
        CVe = (sqrt((3*vG+vE)/nArv+vParc))/Mean*100
        
        # genPar
        genParNames <- c("vG","vParc","vRdm1","vE","vPhen","h2G","c2Parc","c2Rdm1","accClone","CVg%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2GSE$Estimate,c2Rdm1,c2Parc,accClone,CVg,CVe,Mean), 
                                   SE=c(SEvG,SEvParc,SEvRdm1,SEvE,NA,h2GSE$SE,matrix(NA,nrow=6,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "  
        
      }
      if(length(random)==1){
        
        vParc <- mProg$sigmaVector[2]
        vRdm1 <- mProg$sigmaVector[3]
        vRdm2 <- mProg$sigmaVector[4]
        vE <- mProg$sigmaVector[5]
        SEvParc <- sqrt(diag(mProg$sigmaSE))[2]
        SEvRdm1 <- sqrt(diag(mProg$sigmaSE))[3]
        SEvRdm2 <- sqrt(diag(mProg$sigmaSE))[4]
        SEvE <- sqrt(diag(mProg$sigmaSE))[5]
        vPhen <- vG + vParc + vRdm1 + vRdm2 + vE
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm1 <- vRdm2/vPhen
        h2gSE <- sommer::vpredict(mProg, h2g ~ (V1) / (V1+V2+V3+V4+V4))
        CVe = (sqrt((3*vG+vE)/nArv+vParc))/Mean*100
        
        # genPar
        genParNames <- c("vG","vParc","vRdm1","vRdm2","vE","vPhen","h2G","c2Parc","c2Rdm1","c2Rdm2","accClone","CVg%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2GSE$Estimate,c2Rdm1,c2Rdm2,c2Parc,accClone,CVg,CVe,Mean), 
                                   SE=c(SEvG,SEvParc,SEvE,SEvRdm1,SEvRdm2,NA,h2GSE$SE,matrix(NA,nrow=7,ncol=1)), 
                                   row.names = genParNames,genPar_digits))
        genPar[is.na(genPar)] <- " "  
        
      }}
    
    if(plotType=="STP"){
      nRep <- length(levels(factor(data$Rep)))
      
      if(length(random)==0){
        vE <- mProg$sigmaVector[2]
        SEvE <- sqrt(diag(mProg$sigmaSE))[2]
        vPhen <- vG + vE
        h2GSE <- sommer::vpredict(mProg, h2a ~ (V1) / (V1+V2))
        h2mc <- (vG) / (vG+(vE/nRep))
        CVe = sqrt(3*vG+vE)/Mean*100
        
        # genPar
        genParNames <- c("vG","vE","vPhen","h2G","h2mc","accClone","CVg%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vG,vE,vPhen,h2GSE$Estimate,h2mc,accClone,CVg,CVe,Mean), 
                                   SE=c(SEvG,SEvE,NA,h2GSE$SE,matrix(NA,nrow=5,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }
      if(length(random)==1){
        vRdm1 <- mProg$sigmaVector[2]
        vE <- mProg$sigmaVector[3]
        vPhen <- vG + vRdm1 + vE
        SEvRdm1 <- sqrt(diag(mProg$sigmaSE))[2]
        SEvE <- sqrt(diag(mProg$sigmaSE))[3]
        c2Rdm1 <- vRdm1/vPhen
        h2GSE <- sommer::vpredict(mProg, h2a ~ (V1) / (V1+V2+V3))
        CVe = sqrt(3*vG+vE)/Mean*100
        
        # genPar
        genParNames <- c("vG","vRdm1","vE","vPhen","h2G","c2Rdm1","accClone","CVg%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vG,vRdm1,vE,vPhen,h2GSE$Estimate,c2Rdm1,accClone,CVg,CVe,Mean), 
                                   SE=c(SEvG,SEvRdm1,SEvE,NA,h2GSE$SE,matrix(NA,nrow=5,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }  
      if(length(random)==2){
        vRdm1 <- mProg$sigmaVector[2]
        vRdm2 <- mProg$sigmaVector[3]
        vE <- mProg$sigmaVector[4]
        SEvRdm1 <- sqrt(diag(mProg$sigmaSE))[2]
        SEvRdm2 <- sqrt(diag(mProg$sigmaSE))[3]
        SEvE <- sqrt(diag(mProg$sigmaSE))[4]
        vPhen <- vG + vRdm1 + vRdm2 + vE
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        h2GSE <- sommer::vpredict(mProg, h2a ~ (V1) / (V1+V2+V3+V4))
        CVe = sqrt(3*vG+vE)/Mean*100
        
        # genPar
        genParNames <- c("vG","vRdm1","vRdm2","vE","vPhen","h2G","c2Rdm1","c2Rdm2","accClone","CVg%","CVe%","Mean")
        genPar <- round(data.frame(Estimates=c(vG,vRdm1,vRdm2,vE,vPhen,h2GSE$Estimate,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean), 
                                   SE=c(SEvA,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE$SE,matrix(NA,nrow=6,ncol=1)), 
                                   row.names = genParNames),genPar_digits)
        genPar[is.na(genPar)] <- " "
      }
      
    }
    
    # BLUP_Individual
    
    # Criando Cod caso nao existaa
    if(!("Cod" %in% colnames(data))){
      data$Cod = " "
    }
    
    # Clone BLUP
    cloneBLUP <- mProg$U$Treat$resp %>% data.frame() %>% setNames("BLUP") %>% rownames_to_column(.,"Clone") %>% 
      mutate(Clone = str_replace_all(Clone, 'Treat', "")) %>% arrange(desc(BLUP))
    
    # Provenance BLUP
    if(any(random=="Proc")){
      procBLUP <- mAdd$U$Random1$resp %>% data.frame() %>% setNames("procBLUP") %>% rownames_to_column(.,"Provenance") %>% 
        mutate(Provenance = str_replace_all(Provenance, 'Random1', "")) %>% arrange(desc(procBLUP))
    }
    
    # Create final output list
    genParBLUP <- list()
    
    # Exploratory analysis
    genParBLUP$expAnalysis$cloneMeans <- treatMeans
    genParBLUP$expAnalysis$repMeans <- repMeans
    
    # Model output
    genParBLUP$Model$mmerModel <- mProg
    genParBLUP$Model$mSig$fixedSig <- anova(mSig)
    
    genParBLUP$Model$mSig$randomSig <- lmerTest::ranova(mSig, reduce.terms = F)
    
    # Genetic parameters
    genParBLUP$genPar <- genPar
    
    # Clone accuracy
    genParBLUP$accuracyPEV$clonePEV <- dfacc
    # BLUP
    genParBLUP$BLUP$cloneBLUP <- cloneBLUP
    
    if(any(random=="Proc")){
      genParBLUP$BLUP$procBLUP <- procBLUP
      genParBLUP$expAnalysis$procMeans <- procMeans
    }
  }
  
  # output ------------------------------------------------------------------
  
  if(!is.null(directory)){
    
    options(max.print=999999)
    
    dir_0 <- getwd()
    
    if(!is.null(directory)){
      if(dir.exists(file.path(getwd(), directory))==FALSE){  ## Criar diret?rio, caso ele n?o exista e como diret?rio de trabalho
        dir.create(directory)
        setwd(directory)
      }else{setwd(directory)}
    }
    
    genParFile <- paste0(varResp,"_genPar",".txt") # o objeto var foi criado anteriormente
    
    if (file.exists(genParFile)){
      file.remove(genParFile)
    }
    
    # Building string model to output
    
    if(is.null(random)){
      if(plotType=="LP"){
        strMod <- paste0(varResp," = ","f",fixed," + rGen + rParc")}
      else{
        strMod <- paste0(varResp," = ","f",fixed," + rGen")
      }
    }else{
      rnd = glue::glue_collapse(random, " + r", last = " + r")
      if(plotType=="LP"){
        strMod <- paste0(varResp," = ","f",fixed," + rGen + rParc + ","r",rnd)}
      else{
        strMod <- paste0(varResp," = ","f",fixed," + rGen + ","r",rnd)
      }
    }
    
    # sink output
    
    sink(genParFile, append=TRUE, type = "output")
    
    cat("\n----------------------------------------------------------------------------------------------\n")
    cat("                                       |genParBLUP Analysis| \n"                                   )
    cat("\n----------------------------------------------------------------------------------------------\n")
    cat("Date and Time: ")
    cat(format(Sys.time()))
    cat("\n")
    cat("Model: ")
    cat(strMod)
    cat("\n")
    cat("----------------------------------------------------------------------------------------------\n\n")
    cat("---------------------------------------- Exploratory Analysis --------------------------------\n\n")
    print(respMeans)
    cat("\n")
    cat("------------------------------------------ Block repetition ----------------------------------\n\n")
    print(repMeans)
    cat("\n")
    cat("--------------------------------------------- Treatment --------------------------------------\n\n")
    print(treatMeans)
    cat("\n")
    if(any(random=="Proc")||any(fixed=="Proc")){
      cat("-------------------------------------------- Provenance --------------------------------------\n\n")
      print(procMeans)
      cat("\n")
    }
    cat("----------------------------------------- Genetic Parameters ---------------------------------\n\n")
    print(genPar)
    cat("\n")
    cat("\n")
    cat("---------------------------------------------- |BLUP|  ---------------------------------------\n\n")
    cat("\n")
    if(any(random=="Proc")){
      cat("------------------------------------------ Provenance BLUP -----------------------------------\n\n")
      print(dfblupProc)
      cat("\n")
    }
    cat("------------------------------------------- Treatment BLUP -----------------------------------\n\n")
    cat("\n")
    if(treatment=="Clone"){
      print(cloneBLUP)}
    if(treatment=="Prog"){
      print(progBLUP)
      cat("\n")
      cat("------------------------------------------ Individual BLUP -----------------------------------\n\n")
      print(blupOrdered)
      cat("\n")
      cat("---------------------------------- Overlapping Generations Selection -----------------------------\n\n")
      print(overBLUP)
    }
    
    sink()
    
    # Accuracy output #
    
    accFile <- paste0(varResp,"_acc",".txt")
    
    if (file.exists(accFile)){
      file.remove(accFile)
    }
    
    sink(accFile, append=TRUE, type = "output")
    
    cat("\n----------------------------------------------------------------------------------------------\n")
    cat("                                             |Accuracy| \n"                                        )
    cat("----------------------------------------------------------------------------------------------\n\n")
    cat("----------------------------------------- Treatment accuracy ---------------------------------\n\n")
    print(dfacc)
    if(treatment=="Prog"){
      cat("\n")
      cat("----------------------------------------- Individual accuracy --------------------------------\n\n")
      print(dfaccInd)
    }
    
    sink()
    
    # Significance output #
    
    sigFile <- paste0(varResp,"_sig",".txt")
    
    if (file.exists(sigFile)){
      file.remove(sigFile)
    }
    
    sink(sigFile, append=TRUE, type = "output")
    
    cat("\n----------------------------------------------------------------------------------------------\n")
    cat("                                     |Significance of effects| \n"                                 )
    cat("----------------------------------------------------------------------------------------------\n\n")
    cat("------------------------------------------- Fixed effects ------------------------------------\n\n")
    print(anova(mSig))
    cat("----------------------------------------------------------------------------------------------\n\n")
    cat("------------------------------------------ Random effects ------------------------------------\n\n")
    print(lmerTest::ranova(mSig))
    
    sink()
    
    # Otimized BLUP #
    
    if(otimizeSelection==TRUE){
      
      blupOtimized <- paste0(varResp,"_otimizedBLUP",".txt")
      
      if(file.exists(blupOtimized)){
        file.remove(blupOtimized)
      }
      
      sink(blupOtimized, append=TRUE, type = "output")
      
      cat("\n----------------------------------------------------------------------------------------------\n")
      cat("                                           |Otimized BLUP| \n"                                     )
      cat("----------------------------------------------------------------------------------------------\n\n")
      cat("--------------------------------------------- Arguments --------------------------------------\n\n")
      cat("Max. Individuals per Progeny:\n\n")
      print(maxIndProgeny) 
      cat("Max. Progenies per Block:\n\n")
      print(maxProgenyBlock) 
      cat("Controls Excluded:\n\n")
      print(excludeControl)
      cat("----------------------------------------------------------------------------------------------\n\n")
      cat("---------------------------------------- Individual Selection --------------------------------\n\n")
      print(otimizedBLUP)
      
      sink()
      
    }
    
    setwd(dir_0)
    
  }
  class(genParBLUP) <- "genList"
  invisible(genParBLUP)
}