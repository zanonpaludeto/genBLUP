genBLUP <- function(data, varResp, treatment, plotType, fixed = "Rep", random = NULL, method = "ai", 
                    excludeControl = NULL, genPar_digits, codPerc = NULL, optimizeSelection = FALSE, maxIndProgeny = NULL, 
                    maxProgenyBlock = NULL, excludeCod = NULL, directory = NULL){
  
  # loading packages --------------------------------------------------------
  
  pacman::p_load(lme4,stringr,tidyverse,ggroups,dplyr,msm)
  
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
  
  if(!is.numeric(data[,varResp])){
    stop("ERROR: varResp is not numeric, please check your data")
  }
  
  if(treatment=="Clone"&optimizeSelection==T){
    stop("ERROR: You can't optimize individual selection in a clonal analysis")
  }
  
  if(treatment=="Prog"&optimizeSelection==T&!is.numeric(maxIndProgeny)|!is.numeric(maxProgenyBlock)){
    stop("ERROR: If you want to optimize selection, please choose numeric parameters for maxIndProgeny and maxProgenyBlock")
  }
  if(plotType=="LP"&!"Parc"%in%random){
    stop("ERROR: If you're trying to adjust an LP model, you must specify the plot effect as 'Parc' in the random argument")
  }
  
  if(method=="ai"||is.null(method)){
    cat("AI-REML algorithm was selected\n")
  }else{
    cat("EM-REML algorithm was selected\n")
  }
  
  if(!is.null(excludeControl)&!any(data[,treatment]%in%excludeControl)){
    stop(paste0("ERROR: The specified controls in excludeControl argument are not in the ",treatment," column of your dataset"))
  }
  
  # exploratory analysis ----------------------------------------------------
  
  # preparing data
  
  if(!is.null(excludeControl)){
    data <- data %>% filter(!.[,treatment]%in%excludeControl)
    cat(paste(c("Controls:",excludeControl, "where removed using excludeControl argument\n"))) 
  }
  
  # Creating Cod if it doesn't exist
  if(!("Cod" %in% colnames(data))){
    data$Cod = " "
  }
  if(all(c("Cod1","Cod2")%in%colnames(data))){
    data$Cod <- paste0(data$Cod1,";",data$Cod2)
  }
  
  data <- data %>% arrange(get(treatment))
  
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
  
  if(!is.null(codPerc)){
    countTreat <- data %>% count(get(treatment)) %>% setNames(c(treatment,"nObs"))
    suppressMessages(codData <- separate_rows(data,Cod) %>% count(!!as.name(treatment),Cod) %>% filter(Cod!=".") %>% pivot_wider(., names_from = Cod,values_from = n) %>% 
                       left_join(countTreat,.) %>% relocate(nObs, .after = last_col()))
    
    for(i in codPerc){
      if(i %in% names(codData)){
        percCod <- setNames(data.frame(codData[,i]/codData$nObs*100),paste0(i,"%"))
        codData <- cbind(codData,percCod)
      }
    }
    if(any(codPerc%in%names(codData))){
      treatPerc <- codData %>% dplyr::select(as.name(treatment),(which(names(.) == "nObs") + 1):last_col())
      treatMeans <- cbind(treatMeans,treatPerc[,-1])
      treatMeans[is.na(treatMeans)] <- 0
    }
  }
  
  if(any(random=="Proc")||any(fixed=="Proc")){
    procMeans <- groupMeans(data,"Proc")
  }
  
  # factorizing_variables_and_treeCheck -------------------------------------
  
  # factorizing variables
  fct <- c(treatment,fixed,random) 
  data[fct] <- lapply(data[fct], factor)
  
  # treeCheck
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
  }
  
  # pedigree matrix ---------------------------------------------------------
  
  if(treatment=="Prog"){
    
    #creating Ind column
    if(plotType=="LP"){  
      data$Ind <- paste0(data$Rep,sep = ".",data[,names(data) %in% treatment],sep=".",data$Parc,sep=".",data$Arv)
    }else{data$Ind <- paste0(data$Rep,sep = ".",data[,names(data) %in% treatment],sep=".",data$Arv)  
    }
    
    #counting nProg (dams)
    uniqueProg <- length(na.omit(unique(data$Prog)))
    
    if("GD"%in%names(data)&"sire"%in%names(data)){
      
      cat("Grandparents and sires were found and considered in the pedigree\n")
      
      uniqueGD <- length(na.omit(unique(data$GD)))
      nGD <- data.frame(GD = unique(na.omit(data$GD)), nGD = seq(1,uniqueGD))
      uniqueSire <- length(na.omit(unique(data$sire)))
      nSire <- data.frame(sire = unique(na.omit(data$sire)), nSire = seq(uniqueGD+1,uniqueGD+uniqueSire))
      
      ancPed <- data.frame(idNum=c(nGD$nGD,nSire$nSire),sire=0,dam=0)
      
      nProg <- data.frame(Prog = sort(unique(data$Prog)), nProg = as.numeric(seq(uniqueGD+uniqueSire+1,uniqueProg+uniqueGD+uniqueSire)))
    }
    if(!"GD"%in%names(data)&"sire"%in%names(data)){
      
      cat("Sires were found and considered in the pedigree\n")
      
      uniqueSire <- length(na.omit(unique(data$sire)))
      nSire <- data.frame(sire = unique(na.omit(data$sire)), nSire = seq(1,uniqueSire))
      
      ancPed <- data.frame(idNum=nSire$nSire,sire=0,dam=0)
      
      nProg <- data.frame(Prog = sort(unique(data$Prog)), nProg = as.numeric(seq(uniqueSire+1,uniqueProg+uniqueSire)))
    }
    if("GD"%in%names(data)&!"sire"%in%names(data)){
      
      cat("Grandparents were found and considered in the pedigree\n")
      
      uniqueGD <- length(na.omit(unique(data$GD)))
      nGD <- data.frame(GD = unique(na.omit(data$GD)), nGD = seq(1,uniqueGD))
      
      ancPed <- data.frame(idNum=nGD$nGD,sire=0,dam=0)
      
      nProg <- data.frame(Prog = sort(unique(data$Prog)), nProg = as.numeric(seq(uniqueGD+1,uniqueProg+uniqueGD)))
    }
    if(!"GD"%in%names(data)&!"sire"%in%names(data)){
      nProg <- data.frame(Prog = sort(unique(data$Prog)), nProg = as.numeric(seq(1,uniqueProg)))
    }
    
    nInd <- seq(1001,1000+nrow(data))
    data <- data %>% mutate(idNum=as.numeric(nInd))
    data <- dplyr::left_join(data,nProg,by="Prog")
    
    breedR.createPed <- function(pedData){
      
      if("sire"%in%names(data)){
        data <- dplyr::left_join(data,nSire,by="sire")
        ped <- data.frame(pedData,sire=data$nSire) %>% setNames(c("idNum","dam","sire")) %>% relocate(dam, .after = sire)
      }else{
        ped <- data.frame(pedData,sire=0) %>% setNames(c("idNum","dam","sire")) %>% relocate(dam, .after = sire)
      }
      
      #pedigree with grand dam information
      if(any(names(data)=="GD")){
        
        dfGD <- data[c("nProg","GD")] %>% group_by(nProg) %>% distinct() %>% arrange(nProg)
        if(nrow(dfGD)>nrow(nProg)){
          stop("ERROR: There is progenies with more than one grandparent assisgned in the dataset, please check your data")
        }
        ped <- dfGD %>% left_join(.,nGD, by="GD") %>% dplyr::select(-GD) %>% mutate(sire=0, .after = nProg) %>% 
          setNames(c("idNum","sire","dam")) %>% bind_rows(ancPed,.,ped)
      }
      
      if(any(duplicated(ped[,1])==TRUE)){
        stop("ERROR: There is duplicated individuals, please check your data")
      }
      
      ped[is.na(ped)] <- 0
      return(ped)
    }
    
    ped <- breedR.createPed(data[,c("idNum","nProg")])
    
    # dominance matrix --------------------------------------------------------
    
    # if(dominance==TRUE){
    # 
    #   A <- buildA(ped)
    #   D <- ggroups::buildD(ped,A)
    #     if(any(duplicated(rownames(D))==TRUE)){
    #       stop("ERROR: There is duplicated individuals, please check your data")
    #     }}
  }
  
  # statistical modelling ---------------------------------------------------
  
  #creating all possible object formulas
  fixef <- as.formula(paste0("get(varResp) ~ 1 +", paste(fixed, collapse=" + ")))
  
  # if(plotType=="LP"){
  #   if(is.null(random)&treatment=="Prog"){
  #     ranef <- as.formula("~Parc")
  #   }
  #   if(!is.null(random)&treatment=="Prog"){
  #     ranef <- as.formula(paste0("~Parc + ", paste(random, collapse=" + ")))
  #   }
  #   if(is.null(random)&treatment=="Clone"){
  #     ranef <- as.formula(paste0("~Parc + ", treatment))
  #   }
  #   if(!is.null(random)&treatment=="Clone"){
  #     ranef <- as.formula(paste0("~Parc +", treatment, "+", paste(random, collapse=" + ")))
  #   }
  # }else{
  #   if(is.null(random)&treatment=="Prog"){
  #     ranef <- NULL
  #   }
  #   if(!is.null(random)&treatment=="Prog"){
  #     ranef <- as.formula(paste0("~", paste(random, collapse=" + ")))
  #   }
  #   if(is.null(random)&treatment=="Clone"){
  #     ranef <- as.formula(paste0("~", treatment))
  #   }
  #   if(!is.null(random)&treatment=="Clone"){
  #     ranef <- as.formula(paste0("~", treatment, "+", paste(random, collapse=" + ")))
  #   }
  # }
  if(!is.null(random)&"Parc"%in%random){
    random <- c(random[which(random=="Parc")],random[which(random!="Parc")])
  }
  if(treatment=="Prog"&!is.null(random)){
    ranef <- as.formula(paste0("~", paste(random, collapse=" + ")))
  }
  if(treatment=="Clone"&!is.null(random)){
    ranef <- as.formula(paste0("~", treatment, "+", paste(random, collapse=" + ")))
  }
  if(treatment=="Clone"&is.null(random)){
    ranef <- as.formula(paste0("~", treatment))
  }
  if(treatment=="Prog"&is.null(random)){
    ranef <- NULL
  }
  
  #suppressing messages and warnings
  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
  defaultW <- getOption("warn") 
  options(warn = -1)
  convergence = T
  
  while(convergence){
    suppressMessages(
      if(treatment=="Prog"){
        tryCatch(mAdd <- breedR::remlf90(fixef,
                                         random = ranef,
                                         genetic = list(model = "add_animal",
                                                        pedigree = ped,id = "idNum"),
                                         data = data,method = method,
                                         progsf90.options = 'EM-REML 10'), error = function(e) {convergence <<- TRUE})
      }else{
        tryCatch(mClone <- breedR::remlf90(fixef,
                                           random = ranef,
                                           data = data,method = method,
                                           progsf90.options = 'EM-REML 10'), error = function(e) {convergence <<- TRUE})
      })
    
    if(is.null(random)){
      lmerModel <- as.formula(paste0("get(varResp) ~ ",paste(fixed, collapse=" + ")," + ", 
                                     paste0("(1|", treatment, ")"))) 
      
    }else{
      lmerModel <- as.formula(paste0("get(varResp) ~ ",paste(fixed, collapse=" + ")," + ", 
                                     paste0("(1|", treatment, ") + ", paste("(1|",random, collapse=") + "),")")))
    }
    
    mSig <- lme4::lmer(lmerModel, data=data, control = control)
    
    if(!exists("mAdd")&!exists("mClone")){
      cat("AI-REML algorithm got an error, switching to EM-REML\n")
      method = "em"
    }else{
      convergence <- ifelse(treatment=="Prog",any(is.na(mAdd$var[,1])),any(is.na(mClone$var[,1])))
      if(convergence){
        cat("AI-REML algorithm did not converge, switching to EM-REML\n")
        method = "em"
      }
    }
  }
  options(warn = defaultW)
  
  if(treatment=="Prog"){
    
    # reliability_and_individual_blup -----------------------------------------
    
    if("GD"%in%names(data)){
      r2GD <- mAdd$ranef$genetic[[1]] %>% 
        mutate(r2=1-(s.e.)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
        filter(row_number() %in% nGD$nGD) %>% add_column(GD = nGD$GD, .before = "value") %>% mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
        rename(a = value)
    }
    if("sire"%in%names(data)){
      r2Sire <- mAdd$ranef$genetic[[1]] %>% 
        mutate(r2=1-(s.e.)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
        filter(row_number() %in% nSire$nSire) %>% add_column(Sire = nSire$sire, .before = "value") %>% mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
        rename(a = value)
    }
    r2Prog <- mAdd$ranef$genetic[[1]] %>% 
      mutate(r2=1-(s.e.)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
      filter(row_number() %in% nProg$nProg) %>% add_column(Prog = nProg$Prog, .before = "value") %>% mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
      rename(a = value)
    
    r2Ind <- mAdd$ranef$genetic[[1]] %>% 
      mutate(r2=1-s.e.^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
      filter(row_number() > max(nProg$nProg)) %>% add_column(Ind = data$Ind, .before = "value") %>% filter(!is.na(resp)) %>% rename(a = value)
    
    r2Ind_df <- r2Ind %>% 
      {if (plotType=="LP") separate(.,Ind,c("Block","Progeny","Plot","Tree"),  sep="\\.") else separate(.,Ind,c("Block","Progeny","Tree"), sep="\\.")}
    
    accInd <- mean(sqrt(1-((r2Ind$s.e.)^2)/mAdd$var["genetic",1]), na.rm=T)
    accProg <- mean(sqrt(1-(((r2Prog$s.e./2))^2)/(mAdd$var["genetic",1]/4)), na.rm=T)
    
    
    # genetic_parameters ------------------------------------------------------
    
    vA <- mAdd$var["genetic",1]
    vE <- mAdd$var["Residual",1]
    Mean <- mean(data$resp,na.rm=T)
    nRep <- length(unique(data$Rep))
    CVgi <- sqrt(vA)/Mean*100
    
    # Linear plot
    if(plotType=="LP"){
      nArv <- length(unique(data$Arv))
      vParc <- mAdd$var[which(random=="Parc"),1]
      h2d <- (0.75*vA) / (0.75*vA+vE)
      CVe = (sqrt((0.75*vA+vE)/nArv+vParc))/Mean*100
      
      if(length(random)==1){
        vPhen <- vA + vParc + vE
        c2Parc <- vParc/vPhen
        h2a <- vA/vPhen
        h2m <- (0.25*vA) / (0.25*vA+(vParc/nRep)+vE/(nRep*nArv))
        
        genParNames <- c("vA","vParc","vE","vPhen","h2a","h2m","h2d",
                         "c2Parc","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var[random,2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2a,h2m,h2d,c2Parc,accProg,accInd,CVgi,CVe,Mean), 
                                     SE=c(SEvA,SEvParc,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2a,h2m,h2d,c2Parc,accProg,accInd,CVgi,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==2){
        vRdm1 <- mAdd$var[random[2],1]
        vPhen <- vA + vParc + vRdm1 + vE
        h2a <- vA/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        genParNames <- c("vA","vParc",paste0("v",random[2]),"vE","vPhen","h2a","h2d",
                         "c2Parc",paste0("c2",random[2]),"accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var["Parc",2]
          SEvRdm1 <- mAdd$var[2,2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                     SE=c(SEvA,SEvParc,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,accProg,accInd,CVgi,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==3){
        vRdm1 <- mAdd$var[random[2],1]
        vRdm2 <- mAdd$var[random[3],1]
        vPhen <- vA + vParc + vRdm1 + vRdm2  + vE
        h2a <- vA/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        
        genParNames <- c("vA","vParc",paste0("v",random[2]),paste0("v",random[3]),"vE","vPhen","h2a","h2d",
                         "c2Parc",paste0("c2",random[2]),paste0("c2",random[3]),"accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var["Parc",2]
          SEvRdm1 <- mAdd$var[random[2],2]
          SEvRdm2 <- mAdd$var[random[3],2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2Rdm2,accProg,accInd,CVgi,CVe,Mean), 
                                     SE=c(SEvA,SEvParc,SEvRdm1,SEvRdm2,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2Rdm2,accProg,accInd,CVgi,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
    }
    # Single tree plot
    if(plotType=="STP"){
      CVe = sqrt(0.75*vA+vE)/Mean*100
      
      if(length(random)==0){
        vPhen <- vA + vE
        h2a <- vA/vPhen
        h2m <- (0.25*vA) / (0.25*vA+vE/(nRep))
        
        genParNames <- c("vA","vE","vPhen","h2a","h2m","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          genPar <- round(data.frame(Estimates=c(vA,vE,vPhen,h2a,h2m,accProg,accInd,CVgi,CVe,Mean), 
                                     SE=c(SEvA,SEvE,NA,h2aSE,matrix(NA,nrow=6,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vA,vE,vPhen,h2a,h2m,accProg,accInd,CVgi,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
      }
      if(length(random)==1){
        vRdm1 <- mAdd$var[random,1]
        vPhen <- vA + vRdm1 + vE
        h2a <- vA/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        genParNames <- c("vA",paste0("v",random),"vE","vPhen","h2a","c2Rdm1","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvRdm1 <- mAdd$var[random,2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                     SE=c(SEvA,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=6,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,accProg,accInd,CVgi,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
      }
      if(length(random)==2){
        vRdm1 <- mAdd$var[random[1],1]
        vRdm2 <- mAdd$var[random[2],1]
        vPhen <- vA + vRdm1 + vRdm2 + vE
        h2a <- vA/vPhen
        
        genParNames <- c("vA",paste0("v",random[1]),paste0("v",random[2]),"vE","vPhen","h2a",
                         "c2Rdm1","c2Rdm2","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvRdm1 <- mAdd$var[random[1],2]
          SEvRdm2 <- mAdd$var[random[2],2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                     SE=c(SEvA,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=7,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
          
        }else{
          genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2aSE$Estimate,c2Rdm1,c2Rdm1,accProg,accInd,CVgi,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }}
    }
    genPar[is.na(genPar)] <- " "
    
    # BLUP_dataframes ---------------------------------------------------------
    
    # Progeny BLUP
    progBLUP <- r2Prog %>% dplyr::select(Prog,a) %>% rename(Progeny=Prog) %>% arrange(desc(a))
    
    # Individual_BLUP_genetic_gain_and_effectve_population_size --------------- 
    indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
      {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
      {if (plotType=="LP") separate(.,Ind,c("Block","Progeny","Plot","Tree"),  sep="\\.") else separate(.,Ind,c("Block","Progeny","Tree"), sep="\\.")} %>% 
      mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a))
    
    # Provenance BLUP
    if(any(random=="Proc")){
      procBLUP <- mAdd$ranef$Proc[[1]] %>% rownames_to_column(.,"Proc") %>% rename(BLUP=value) %>% arrange(desc(BLUP))
      indBLUP <- indBLUP %>% left_join(.,procBLUP[,1:2],by="Proc") %>% mutate(a=a+BLUP) %>% mutate("u+a"=Mean+a) %>% dplyr::select(.,-last_col()) %>% 
        arrange(desc(a))
      progBLUP <- progBLUP %>% left_join(.,data[,c("Prog","Proc")], by= c("Progeny"="Prog")) %>% unique() %>% left_join(.,procBLUP[,1:2], by="Proc") %>% 
        mutate(a=a+BLUP) %>% dplyr::select(-c(Proc,BLUP)) %>% arrange(desc(a))
    }
    
    # Absolute and percentual genetic gain
    gain <- matrix(nrow=nrow(indBLUP), ncol=1)
    for(i in 1:nrow(indBLUP)){
      gain[i,] <- mean(indBLUP$a[1:i])
    }
    indBLUP <- indBLUP %>% mutate(Gain=gain) %>% mutate("Gain%"=gain/Mean*100)
    
    # Effective population size
    dfNe <- matrix(nrow=nrow(indBLUP), ncol=1)
    
    for(i in 1:nrow(indBLUP)){
      varkf <- as.vector(table(indBLUP$Progeny[1:i]))
      dfNe[i,] <- var(varkf)
    }
    dfNe[is.na(dfNe)] <- 0
    
    indBLUP <- indBLUP %>% mutate(Np = cumsum(!duplicated(indBLUP$Progeny))) %>% mutate(seq = seq(1:nrow(indBLUP))) %>% mutate(Kf = seq/Np) %>% 
      mutate(varkf = dfNe) %>% mutate(Ne = (4*Np*Kf)/(3+Kf+(varkf/Kf))) %>% dplyr::select(.,-c("Np","seq","Kf","varkf")) %>% 
      relocate(Cod, .after = last_col())
    
    #if(dominance==TRUE){
    #domBLUP <- mDom$U$`u:Ind`$resp %>% data.frame %>% setNames("d") %>% rownames_to_column(.,"Ind")
    #dfblup <- merge(dfblup,domBLUP,by="Ind")
    #dfblup$g <- dfblup$BLUP+dfblup$d
    #}
    
    # Overlapping generations
    overBLUP <- progBLUP %>% {if (plotType=="LP") mutate(.,Block = 0, Plot = 0, Tree = 0, Cod = "Parent") %>% 
        relocate(Block,Progeny,Plot,Tree,a,Cod) else mutate(.,Block = 0, Tree = 0, Cod = "Parent") %>% relocate(Block,Progeny,Tree,a,Cod)} %>% 
      rbind(.,indBLUP[,colnames(.)]) %>% arrange(desc(a))
    
    # optimize_selection -------------------------------------------------------
    
    if(optimizeSelection==TRUE){
      
      rankBLUP <- indBLUP
      
      if(!is.null(excludeCod)){
        rankBLUP <- indBLUP[!indBLUP$Cod %in% excludeCod,]
      }
      
      # Filtering using the maximum number of individuals per progeny, then the max number of progenies addmited in the same block
      rankBLUP <- rankBLUP %>% mutate(Np_csum = ave(Progeny==Progeny, Progeny, FUN=cumsum)) %>% filter(Np_csum <= maxIndProgeny) %>% 
        mutate(RP = paste0(Block,Progeny)) %>% mutate(RP_csum = ave(RP==RP,RP, FUN=cumsum)) %>% 
        filter(RP_csum <= maxProgenyBlock) %>% dplyr::select(.,-c("Np_csum","RP","RP_csum"))
      
      # recalculating effective population size and genetic gain
      
      dfNe <- matrix(nrow=nrow(rankBLUP), ncol=1)
      
      for(i in 1:nrow(rankBLUP)){
        varkf <- as.vector(table(rankBLUP$Progeny[1:i]))
        dfNe[i,] <- var(varkf)
      }
      dfNe[is.na(dfNe)] <- 0
      
      optimizedBLUP <- rankBLUP %>% mutate(Np = cumsum(!duplicated(rankBLUP$Progeny))) %>% mutate(seq = seq(1:nrow(rankBLUP))) %>% mutate(Kf = seq/Np) %>% 
        mutate(varkf = dfNe) %>% mutate(Ne = (4*Np*Kf)/(3+Kf+(varkf/Kf))) %>% dplyr::select(.,-c("Np","seq","Kf","varkf")) %>% 
        relocate(Cod, .after = last_col())
    }
    
    
    # final_output_list -------------------------------------------------------
    
    genParBLUP <- list()
    
    # exploratory_analysis -------------------------------------------------
    genParBLUP$expAnalysis$respMeans <- respMeans; genParBLUP$expAnalysis$progMeans <- treatMeans; genParBLUP$expAnalysis$repMeans <- repMeans
    
    # model_output_and_significance_tests -------------------------------------------------
    genParBLUP$Model$remlMethod <- method; genParBLUP$Model$remlModel <- mAdd; genParBLUP$Model$mSig$fixedSig <- anova(mSig)
    genParBLUP$Model$mSig$randomSig <- lmerTest::ranova(mSig, reduce.terms = F)
    
    # genetic_parameters_accuracy_and_Blup -------------------------------------------------
    genParBLUP$genPar <- genPar; genParBLUP$blupAccuracy$progAccuracy <- r2Prog; genParBLUP$blupAccuracy$indAccuracy <- r2Ind_df 
    genParBLUP$BLUP$progBLUP <- progBLUP; genParBLUP$BLUP$indBLUP <- indBLUP; genParBLUP$BLUP$overBLUP <- overBLUP
    
    if(any(random=="Proc")){
      genParBLUP$BLUP$procBLUP <- procBLUP
    }
    if(any(fixed=="Proc")){
      genParBLUP$expAnalysis$procMeans <- procMeans
    }
    
    if(optimizeSelection==TRUE){
      genParBLUP$BLUP$optimizedBLUP <- optimizedBLUP
    }
  }
  if(treatment=="Clone"){
    
    # reliability_and_individual_blup -----------------------------------------
    
    r2Clone <- mClone$ranef$Clone[[1]] %>% 
      mutate(r2=1-(s.e./2)^2/(diag(diag(length(mClone$ranef$Clone[[1]][,1])))*as.data.frame(mClone$var)["Clone",1])) %>% 
      rownames_to_column("Clone") %>% rename(g = value)
    
    accClone <- mean(sqrt(1-((r2Clone$s.e.)^2)/(mClone$var["Clone",1])), na.rm=T)
    
    # Genetic Parameters
    vG <- mClone$var["Clone",1]
    vE <- mClone$var["Residual",1]
    Mean <- mean(data$resp,na.rm=T)
    nRep <- length(unique(data$Rep))
    CVg <- sqrt(vG)/Mean*100
    
    if(plotType=="LP"){
      nArv <- length(unique(data$Arv))
      vParc <- mClone$var["Parc",1]
      CVe = (sqrt((3*vG+vE)/nArv+vParc))/Mean*100
      
      if(length(random)==1){
        vPhen <- vG + vParc + vE
        h2G <- vG/vPhen
        h2mc <- (vG) / (vG+(vParc/nRep)+vE/(nRep*nArv))
        c2Parc <- vParc/vPhen
        
        genParNames <- c("vG","vParc","vE","vPhen","h2G","h2mc","c2Parc","accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          SEvParc <- mClone$var["Parc",2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3),
                                    c(vG,vParc,vE),
                                    mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vParc,vE,vPhen,h2G,h2mc,c2Parc,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvParc,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vParc,vE,vPhen,h2G,h2mc,c2Parc,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==2){
        vRdm1 <- mClone$var[random[2],1]
        vPhen <- vG + vParc + vRdm1 + vE
        h2G <- vG/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        genParNames <- c("vG","vParc",paste0("v",random[2]),"vE","vPhen","h2G","c2Parc",
                         paste0("c2",random[2]),"accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          SEvParc <- mClone$var["Parc",2]
          SEvRdm1 <- mClone$var[random[2],2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4),
                                    c(vG,vParc,vRdm1,vE),
                                    mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2G,c2Parc,c2Rdm1,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvParc,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2G,c2Parc,c2Rdm1,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==3){
        rdmIndex <- which(random!="Parc")
        vRdm1 <- mClone$var[random[2],1]
        vRdm2 <- mClone$var[random[3],1]
        vPhen <- vG + vParc + vRdm1 + vRdm2 + vE
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        h2G <- vG/vPhen
        
        genParNames <- c("vG","vParc",paste0("v",random[2]),paste0("v",random[3]),"vE","vPhen","h2G",
                         "c2Parc",paste0("c2",random[2]),paste0("c2",random[3]),"accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          SEvParc <- mClone$var["Parc",2]
          SEvRdm1 <- mClone$var[random[2],2]
          SEvRdm2 <- mClone$var[random[3],2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4+x5),
                                    c(vG,vParc,vRdm1,vRdm2,vE),
                                    mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2Parc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvParc,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=7,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2Parc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }}
    
    if(plotType=="STP"){
      CVe = sqrt(3*vG+vE)/Mean*100
      
      if(length(random)==0){
        vPhen <- vG + vE
        h2G <- vG/vPhen
        h2mc <- (vG) / (vG+(vE/nRep))
        
        genParNames <- c("vG","vE","vPhen","h2G","h2mc","accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- msm::deltamethod(~ x1/(x1+x2),
                                    c(vG,vE),
                                    mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==1){
        vRdm1 <- mClone$var[random,1]
        vPhen <- vG + vRdm1 + vE
        h2G <- vG/vPhen
        h2mc <- (vG) / (vG+(vE/nRep))
        c2Rdm1 <- vRdm1/vPhen
        
        genParNames <- c("vG",paste0("v",random),"vE","vPhen","h2G","h2mc",paste0("c2",random),"accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          SEvRdm1 <- mClone$var[random,2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3),
                                    c(vG,vRdm1,vE),
                                    mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vRdm1,vE,vPhen,h2G,h2mc,c2Rdm1,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vRdm1,vE,vPhen,h2G,h2mc,c2Rdm1,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }  
      if(length(random)==2){
        vRdm1 <- mClone$var[random[1],1]
        vRdm2 <- mClone$var[random[2],1]
        vPhen <- vG + vRdm1 + vRdm2 + vE
        h2G <- vG/vPhen
        h2mc <- (vG) / (vG+(vE/nRep))
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        
        genParNames <- c("vG",paste0("v",random[1]),paste0("v",random[1]),"vE","vPhen","h2G","h2mc",
                         paste0("c2",random[1]),paste0("c2",random[2]),"accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          SEvRdm1 <- mClone$var[random[1],2]
          SEvRdm2 <- mClone$var[random[2],2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4),
                                    c(vG,vRdm1,vE),
                                    mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
      }
    }
    genPar[is.na(genPar)] <- " "
    
    # BLUP_dataframes ---------------------------------------------------------
    
    # Clone BLUP
    cloneBLUP <- r2Clone %>% dplyr::select(Clone,g) %>% arrange(desc(g))
    
    # Provenance BLUP
    if(any(random=="Proc")){
      procBLUP <- mAdd$ranef$Random1[[1]] %>% rownames_to_column(.,"Proc") %>% rename(BLUP=value) %>% arrange(desc(BLUP))
    }
    
    # final_output_list -------------------------------------------------------
    genParBLUP <- list()
    
    # exploratory_analysis -------------------------------------------------
    genParBLUP$expAnalysis$respMeans <- respMeans; genParBLUP$expAnalysis$cloneMeans <- treatMeans
    genParBLUP$expAnalysis$repMeans <- repMeans
    
    # model_output_and_significance_tests -------------------------------------------------
    genParBLUP$Model$remlMethod <- method; genParBLUP$Model$remlModel <- mClone; genParBLUP$Model$mSig$fixedSig <- anova(mSig)
    genParBLUP$Model$mSig$randomSig <- lmerTest::ranova(mSig, reduce.terms = F)
    
    # genetic_parameters_accuracy_and_Blup -------------------------------------------------
    genParBLUP$genPar <- genPar; genParBLUP$blupAccuracy$cloneAccuracy <- r2Clone;
    genParBLUP$BLUP$cloneBLUP <- cloneBLUP
    
    if(any(random=="Proc")){
      genParBLUP$BLUP$procBLUP <- procBLUP
    }
    if(any(fixed=="Proc")){
      genParBLUP$expAnalysis$procMeans <- procMeans
    }
  }
  
  # output ------------------------------------------------------------------
  
  if(!is.null(directory)){
    
    options(max.print=999999)
    
    dir_0 <- getwd()
    
    if(!is.null(directory)){
      if(dir.exists(file.path(getwd(), directory))==FALSE){
        dir.create(directory)
        setwd(directory)
      }else{setwd(directory)}
    }
    
    genParFile <- paste0(directory,"_",varResp,"_genPar",".txt")
    
    if (file.exists(genParFile)){
      file.remove(genParFile)
    }
    
    # Building string model to output
    # if(is.null(random)){
    #   if(plotType=="LP"){
    #     strMod <- paste0(varResp," = ","f",fixed," + rGen + rParc")}
    #   else{
    #     strMod <- paste0(varResp," = ","f",fixed," + rGen")
    #   }
    # }else{
    #   rnd = glue::glue_collapse(random, " + r", last = " + r")
    #   if(plotType=="LP"){
    #     strMod <- paste0(varResp," = ","f",fixed," + rGen + rParc + ","r",rnd)}
    #   else{
    #     strMod <- paste0(varResp," = ","f",fixed," + rGen + ","r",rnd)
    #   }
    # }
    
    if(!is.null(random)){
      strFix <- paste(paste0("f",fixed),collapse="+")
      strRnd <- paste(paste0("r",random),collapse="+")
      strMod <- paste0(varResp,"=",strFix,"+","r",treatment,"+",strRnd)
    }else{
      strFix <- paste(paste0("f",fixed),collapse="+")
      strMod <- paste0(varResp,"=",strFix,"+","r",treatment)
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
    cat("Controls Excluded:\n\n")
    print(excludeControl)
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
      cat("------------------------------------------ Provenance BLUP --------------------------------------\n\n")
      print(procBLUP)
      cat("\n")
    }
    cat("------------------------------------------ Treatment BLUP ------------------------------------\n\n")
    cat("\n")
    if(treatment=="Clone"){
      print(cloneBLUP)}
    if(treatment=="Prog"){
      print(progBLUP)
      cat("\n")
      cat("---------------------------------------- Individual BLUP -----------------------------------\n\n")
      print(indBLUP)
      cat("\n")
      cat("--------------------------------- Overlapping Generations Selection -----------------------------\n\n")
      print(overBLUP)
    }
    
    sink()
    
    # Accuracy output #
    
    accFile <- paste0(directory,"_",varResp,"_acc",".txt")
    
    if (file.exists(accFile)){
      file.remove(accFile)
    }
    
    sink(accFile, append=TRUE, type = "output")
    
    cat("\n----------------------------------------------------------------------------------------------\n")
    cat("                                             |Accuracy| \n"                                        )
    cat("----------------------------------------------------------------------------------------------\n\n")
    if("GD"%in%names(data)){
      cat("\n")
      cat("--------------------------------------- Grandmother accuracy ---------------------------------\n\n")
      print(r2GD)
    }
    if("sire"%in%names(data)){
      cat("\n")
      cat("----------------------------------------- Sires accuracy -------------------------------------\n\n")
      print(r2Sire)
    }
    cat("\n")
    cat("----------------------------------------- Treatment accuracy ---------------------------------\n\n")
    if(treatment=="Clone"){
      cat("\n")
      print(r2Clone)}
    if(treatment=="Prog"){
      print(r2Prog)
      cat("\n")
      cat("----------------------------------------- Individual accuracy --------------------------------\n\n")
      print(r2Ind_df)
    }
    
    sink()
    
    # Significance output #
    
    sigFile <- paste0(directory,"_",varResp,"_sig",".txt")
    
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
    
    # optimized BLUP #
    
    if(treatment=="Prog" & optimizeSelection==TRUE){
      
      blupoptimized <- paste0(directory,"_",varResp,"_optimizedBLUP",".txt")
      
      if(file.exists(blupoptimized)){
        file.remove(blupoptimized)
      }
      
      sink(blupoptimized, append=TRUE, type = "output")
      
      cat("\n----------------------------------------------------------------------------------------------\n")
      cat("                                           |optimized BLUP| \n"                                     )
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
      print(optimizedBLUP)
      
      sink()
      
    }
    
    setwd(dir_0)
    
  }
  
  if(is.null(directory)){
    cat("Directory argument is 'NULL', so no txt outputs were created\n")
  }else{
    cat("Outputs were created at", paste0(getwd(),"/",directory))
  }
  
  class(genParBLUP) <- "genList"
  invisible(genParBLUP)
}
