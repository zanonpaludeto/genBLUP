genBLUP <- function(data, varResp, treatment, plotType, fixed = "Rep", random = NULL, method = "ai", GxE = F, PxE = F,
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
  if(plotType=="LP"&is.null(random)){
    cat(paste0("You choose to adjust an LP model without specifying 'Parc' in the random argument, 
the function will automatically build Rep:",treatment," cross interactions and proceed the analysis \n \n"))
  }
  if(PxE&!GxE){
    stop("ERROR: Argument 'GxE' is FALSE. You need to run a GxE analysis to estimate a Provenance x Environment effect")
  }
  if(PxE&!"Proc" %in% random){
    stop("ERROR: To estimate a Provenance x Environment effect, you need to specify 'Proc' in the random argument")
  }
  if("Proc"%in%random & !"Proc"%in%names(data)){
    stop("ERROR: You cannot try to estimate provenance effect without a column names 'Proc' in your data")
  }
  
  if(optimizeSelection&GxE){
    warning("WARNING: You can't optimize selection in a GxE analysis, the argument was changed to FALSE")
    optimizeSelection=F
  }
  
  if(!is.null(excludeControl)&!any(data[,treatment]%in%excludeControl)){
    stop(paste0("ERROR: The specified controls in excludeControl argument are not in the ",treatment," column of your dataset"))
  }
  
  # messages ----------------------------------------------------------------
  
  if(GxE) cat("Genotype x Environment Analysis: \n") else cat("Individual site Analysis: \n")
  
  if(method=="ai"||is.null(method)){
    cat("AI-REML algorithm was selected\n")
  }else{
    cat("EM-REML algorithm was selected\n")
  }
  
  # exploratory analysis ----------------------------------------------------
  
  # preparing data
  
  if(!is.null(excludeControl)){
    data <- data %>% filter(!.[,treatment]%in%excludeControl)
    cat(paste(c("Controls:",excludeControl, "has been removed using excludeControl argument\n"))) 
  }
  
  # Creating Cod if it doesn't exist
  if(!("Cod" %in% colnames(data))){
    data$Cod = "."
  }
  if(all(c("Cod1","Cod2")%in%colnames(data))){
    data$Cod <- paste0(data$Cod1,";",data$Cod2)
  }
  
  # creating arv column if it doesn't exist
  if(!"Arv"%in%names(data)){
    data$Arv <- 1
  }
  
  if(GxE){
    #creating basic GxE interactions
    data <- data %>% mutate(EnvRep = factor(paste0(Env,"_x_",Rep)))
    if(treatment=="Prog"){
      data <- data %>% mutate(EnvProg = factor(paste0(Env,"_x_",Prog)))
      if(plotType=="LP") data <- data %>% mutate(Parc = factor(paste0(Env,"_x_",Rep,"_x_",Prog)))
    }else{
      data <- data %>% mutate(EnvClone = factor(paste0(Env,"_x_",Clone)))
      if(plotType=="LP") data <- data %>% mutate(Parc = factor(paste0(Env,"_x_",Rep,"_x_",Clone)))
    }
    if("Proc"%in%random&PxE){
      data <- data %>% mutate(EnvProc = factor(paste0(Env,"_x_",Proc)))
      random <- c(random,"EnvProc")  
    }
  }else{
    if(plotType=="LP"&is.null(random)){
      random = "Parc"
      if(treatment=="Prog") data <- data %>% mutate(Parc = factor(paste0(Rep,"_x_",Prog)))
      if(treatment=="LP") data <- data %>% mutate(Parc = factor(paste0(Rep,"_x_",Prog)))
    }
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
  suppressMessages(treatMeans <- groupMeans(data,treatment))
  suppressMessages(repMeans <- groupMeans(data,"Rep"))
  
  if(GxE){
    envMeans <- groupMeans(data,"Env")
    repMeans <- groupMeans(data,"EnvRep") %>% separate(.,EnvRep, into = c("Env","Rep"), sep="_x_")
  }
  
  if(any(random=="Proc")||any(fixed=="Proc")){
    suppressMessages(procMeans <- groupMeans(data,"Proc"))
  }
  
  if(!is.null(codPerc)&any(codPerc%in%data$Cod)){
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
  
  # factorizing_variables_and_treeCheck -------------------------------------
  
  # factorizing variables
  if(GxE) fct <- c("Env",treatment,fixed,random) else fct <- c(treatment,fixed,random)
  
  data[fct] <- lapply(data[fct], factor)
  
  # # treeCheck
  # if(plotType=="STP"){
  #   
  #   treeCheck <- function(data){
  #     
  #     arvList <- list()
  #     
  #     for(i in levels(factor(data$Rep))){
  #       if(treatment=="Clone"){
  #         arvList[i] <- as.data.frame(ave(data[data$Rep==i,]$Clone==data[data$Rep==i,]$Clone, 
  #                                         data[data$Rep==i,]$Clone, FUN=cumsum))
  #       }else{
  #         arvList[i] <- as.data.frame(ave(data[data$Rep==i,]$Prog==data[data$Rep==i,]$Prog, 
  #                                         data[data$Rep==i,]$Prog, FUN=cumsum))
  #       }}
  #     
  #     Arv <- unlist(arvList)
  #     data$Arv <- Arv
  #     return(data)
  #   }
  #   
  #   data <- treeCheck(data)
  # }
  
  # pedigree matrix ---------------------------------------------------------
  
  if(treatment=="Prog"){
    
    #creating Ind column
    data$seq <- paste0(rep("I"),seq(from=1, to=nrow(data))) #creating sequential index to differentiate individuals
    
    if(GxE){
      if("sire"%in%names(data)){
        if(plotType=="LP"){  
          data$Ind <- paste0(data$seq,"-x-",data$Env,"-x-",data$Rep,"-x-",data$sire,"-x-",data[,names(data) %in% treatment],"-x-",data$Parc,"-x-",data$Arv)
        }else{
          data$Ind <- paste0(data$seq,"-x-",data$Env,"-x-",data$Rep,"-x-",data$sire,"-x-",data[,names(data) %in% treatment],"-x-",data$Arv)  
        }
      }else{
        if(plotType=="LP"){  
          data$Ind <- paste0(data$seq,"-x-",data$Env,"-x-",data$Rep,"-x-",data[,names(data) %in% treatment],"-x-",data$Parc,"-x-",data$Arv)
        }else{
          data$Ind <- paste0(data$seq,"-x-",data$Env,"-x-",data$Rep,"-x-",data[,names(data) %in% treatment],"-x-",data$Arv)  
        }
      }
      
    }else{
      if("sire"%in%names(data)){
        if(plotType=="LP"){  
          data$Ind <- paste0(data$seq,"-x-",data$Rep,"-x-",data$sire,"-x-",data[,names(data) %in% treatment],"-x-",data$Parc,"-x-",data$Arv)
        }else{data$Ind <- paste0(data$seq,"-x-",data$Rep,"-x-",data$sire,"-x-",data[,names(data) %in% treatment],"-x-",data$Arv)  
        }
      }else{
        if(plotType=="LP"){  
          data$Ind <- paste0(data$seq,"-x-",data$Rep,"-x-",data[,names(data) %in% treatment],"-x-",data$Parc,"-x-",data$Arv)
        }else{data$Ind <- paste0(data$seq,"-x-",data$Rep,"-x-",data[,names(data) %in% treatment],"-x-",data$Arv)  
        }
      }}
    
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
  if(GxE){
    
    # fixed effects
    if(length(fixed)==1&fixed=="Rep"){
      fixef <- as.formula(paste0("get(varResp) ~ Env + EnvRep"))
    }else{
      fixef <- as.formula(paste0("get(varResp) ~ Env + EnvRep", paste("+",fixed[which(fixed!="Rep")], collapse=" + ")))
    }
    
    # random effects
    if(!is.null(random)&"Parc"%in%random){
      random <- c(random[which(random=="Parc")],random[which(random!="Parc")])
    }
    if(treatment=="Prog"&!is.null(random)){
      ranef <- as.formula(paste0("~ EnvProg","+", paste(random, collapse=" + ")))
    }
    if(treatment=="Clone"&!is.null(random)){
      ranef <- as.formula(paste0("~ Clone + EnvClone", "+", paste(random, collapse=" + ")))
    }
    if(treatment=="Clone"&is.null(random)){
      ranef <- as.formula("~ Clone + EnvClone")
    }
    if(treatment=="Prog"&is.null(random)){
      ranef <- as.formula("~ EnvProg")
    } 
  }else{
    
    #fixed effects
    if(!is.null(fixed)){
      fixef <- as.formula(paste0("get(varResp) ~ 1 +", paste(fixed, collapse=" + ")))
    }else{
      fixef <- as.formula(paste0("get(varResp) ~ 1"))
    }
    
    # random effects
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
  }
  
  # breedR modelling --------------------------------------------------------
  
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
    
    #significance of effects
    
    if(GxE){
      if(is.null(random)){
        lmerModel <- as.formula(paste0(format(fixef)," + ", 
                                       paste0("(1|", treatment, ")"," + ","(1|", paste0("Env",treatment), ")"))) 
        
      }else{
        lmerModel <- as.formula(paste0(paste0("get(varResp) ~ Env + EnvRep")," + ", 
                                       paste0("(1|", treatment, ")","+","(1|", paste0("Env",treatment), ") + ")
                                       , paste("(1|", random, collapse=") + "),")"))
      }
    }else{
      if(is.null(random)){
        lmerModel <- as.formula(paste0("get(varResp) ~ ",paste(fixed, collapse=" + ")," + ", 
                                       paste0("(1|", treatment, ")"))) 
        
      }else{
        lmerModel <- as.formula(paste0("get(varResp) ~ ",paste(fixed, collapse=" + ")," + ", 
                                       paste0("(1|", treatment, ") + ", paste("(1|",random, collapse=") + "),")")))
      }
    }
    suppressMessages(mSig <- lme4::lmer(lmerModel, data=data, control = control))
    
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
    
    r2Prog <- mAdd$ranef$genetic[[1]] %>% 
      mutate(r2=1-(s.e.)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
      filter(row_number() %in% nProg$nProg) %>% add_column(Prog = nProg$Prog, .before = "value") %>% mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
      rename(a = value)
    
    r2Ind <- mAdd$ranef$genetic[[1]] %>% 
      mutate(r2=1-s.e.^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
      filter(row_number() > max(nProg$nProg)) %>% add_column(Ind = data$Ind, .before = "value") %>% filter(!is.na(resp)) %>% rename(a = value)
    
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
      
      if(GxE){
        r2Ind_df <- r2Ind %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Env","Block","Sire","Progeny","Plot","Tree"),  sep="-x-") else 
            separate(.,Ind,c("Seq","Env","Block","Sire","Progeny","Tree"), sep="-x-")}
      }else{
        r2Ind_df <- r2Ind %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Block","Sire","Progeny","Plot","Tree"),  sep="-x-") else 
            separate(.,Ind,c("Seq","Block","Sire","Progeny","Tree"), sep="-x-")}
      }
      
    }else{
      
      if(GxE){
        r2Ind_df <- r2Ind %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Env","Block","Progeny","Plot","Tree"),  sep="-x-") else 
            separate(.,Ind,c("Seq","Env","Block","Progeny","Tree"), sep="-x-")}
        
      }else{
        r2Ind_df <- r2Ind %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Block","Progeny","Plot","Tree"),  sep="-x-") else 
            separate(.,Ind,c("Seq","Block","Progeny","Tree"), sep="-x-")}
      }}  
    
    
    accInd <- mean(sqrt(1-((r2Ind$s.e.)^2)/mAdd$var["genetic",1]), na.rm=T)
    accProg <- mean(sqrt(1-(((r2Prog$s.e./2))^2)/(mAdd$var["genetic",1]/4)), na.rm=T)
    
    
    # genetic_parameters ------------------------------------------------------
    
    vA <- mAdd$var["genetic",1]
    vE <- mAdd$var["Residual",1]
    Mean <- mean(data$resp,na.rm=T)
    nRep <- length(unique(data$Rep))
    CVgi <- sqrt(vA)/Mean*100
    
    if(GxE){
      vGxE <- mAdd$var["EnvProg",1]
      nEnv <- length(unique(data$Env))
      rgloc <- (vA/4)/((vA/4)+vGxE)
    }
    
    # Linear plot
    if(plotType=="LP"){
      nArv <- length(unique(data$Arv))
      vParc <- mAdd$var[which(random=="Parc"),1]
      h2d <- (0.75*vA) / (0.75*vA+vE)
      CVe = (sqrt((0.75*vA+vE)/nArv+vParc))/Mean*100
      
      if(length(random)==1){
        if(GxE){
          vPhen <- vA + vGxE + vParc + vE
          c2GxE <- vGxE/vPhen
          h2m <- (0.25*vA) / (0.25*vA+(vGxE/nEnv)+(vParc/nRep)+vE/(nRep*nArv))
          genParNames <- c("vA","vGxE","vParc","vE","vPhen","h2a","h2m","h2d",
                           "c2Parc","c2GxE","rgloc","accProg","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vParc + vE
          h2m <- (0.25*vA) / (0.25*vA+(vParc/nRep)+vE/(nRep*nArv))
          genParNames <- c("vA","vParc","vE","vPhen","h2a","h2m","h2d",
                           "c2Parc","accProg","accInd","CVgi%","CVe%","Mean")
        }
        c2Parc <- vParc/vPhen
        h2a <- vA/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          if(GxE) SEvGxE <- mAdd$var["EnvProg",2]
          SEvParc <- mAdd$var[random,2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vParc,vGxE,vE,vPhen,h2a,h2m,h2d,c2Parc,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvParc,SEvGxE,SEvE,NA,h2aSE,matrix(NA,nrow=10,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2a,h2m,h2d,c2Parc,accProg,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvParc,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vParc,vGxE,vE,vPhen,h2a,h2m,h2d,c2Parc,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2a,h2m,h2d,c2Parc,accProg,accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==2){
        vRdm1 <- mAdd$var[random[2],1]
        if(GxE){
          vPhen <- vA + vParc + vRdm1 + vGxE + vE
          c2GxE <- vGxE/vPhen
          genParNames <- c("vA","vParc",paste0("v",random[2]),"vGxE","vE","vPhen","h2a","h2d",
                           "c2Parc",paste0("c2",random[2]),"c2GxE","rgloc","accProg","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vParc + vRdm1 + vE
          genParNames <- c("vA","vParc",paste0("v",random[2]),"vE","vPhen","h2a","h2d",
                           "c2Parc",paste0("c2",random[2]),"accProg","accInd","CVgi%","CVe%","Mean")
        }
        h2a <- vA/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var["Parc",2]
          SEvRdm1 <- mAdd$var[2,2]
          if(GxE) SEvGxE <- mAdd$var["EnvProg",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vGxE,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvParc,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvParc,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vGxE,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,accProg,accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==3){
        vRdm1 <- mAdd$var[random[2],1]
        vRdm2 <- mAdd$var[random[3],1]
        if(GxE){
          vPhen <- vA + vParc + vRdm1 + vRdm2 + vGxE + vE
          genParNames <- c("vA","vParc",paste0("v",random[2]),paste0("v",random[3]),"vGxE","vE","vPhen","h2a","h2d",
                           "c2Parc",paste0("c2",random[2]),paste0("c2",random[3]),"c2GxE","rgloc","accProg","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vParc + vRdm1 + vRdm2  + vE
          genParNames <- c("vA","vParc",paste0("v",random[2]),paste0("v",random[3]),"vE","vPhen","h2a","h2d",
                           "c2Parc",paste0("c2",random[2]),paste0("c2",random[3]),"accProg","accInd","CVgi%","CVe%","Mean")
        }
        h2a <- vA/vPhen
        c2GxE <- vGxE/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var["Parc",2]
          SEvRdm1 <- mAdd$var[random[2],2]
          SEvRdm2 <- mAdd$var[random[3],2]
          if(GxE) SEvGxE <- mAdd$var["EnvProg",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vGxE,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2Rdm2,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean), 
                                   SE=c(SEvA,SEvParc,SEvRdm1,SEvRdm2,SEvE,NA,h2aSE,matrix(NA,nrow=11,ncol=1)), 
                                   row.names = genParNames),genPar_digits) else 
                                     genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2Rdm2,accProg,accInd,CVgi,CVe,Mean), 
                                                                SE=c(SEvA,SEvParc,SEvRdm1,SEvRdm2,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)), 
                                                                row.names = genParNames),genPar_digits)
        }else{
          if(GxE) round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vGxE,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2Rdm2,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean),
                                   row.names = genParNames),genPar_digits) else 
                                     genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2Rdm2,accProg,accInd,CVgi,CVe,Mean),
                                                                row.names = genParNames),genPar_digits)
        }
      }
    }
    # Single tree plot
    if(plotType=="STP"){
      CVe = sqrt(0.75*vA+vE)/Mean*100
      
      if(length(random)==0){
        if(GxE){
          vPhen <- vA + vGxE + vE
          c2GxE <- vGxE/vPhen
          h2m <- (0.25*vA) / (0.25*vA+(vGxE/nEnv)+vE/(nRep))
          genParNames <- c("vA","vGxE","vE","vPhen","h2a","h2d","h2m","accProg","accInd","c2GxE","rgloc","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vE
          genParNames <- c("vA","vGxE","vE","vPhen","h2a","h2d","h2m","accProg","accInd","CVgi%","CVe%","Mean")
        }
        h2a <- vA/vPhen
        h2d <- (0.75*vA) / (0.75*vA+vE)
        h2m <- (0.25*vA) / (0.25*vA+vE/(nRep))
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          if(GxE) SEvGxE <- mAdd$var["EnvProg",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vGxE,vE,vPhen,h2a,h2d,h2m,accProg,accInd,c2GxE,rgloc,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvGxE,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vE,vPhen,h2a,h2d,h2m,accProg,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvE,NA,h2aSE,matrix(NA,nrow=7,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vGxE,vE,vPhen,h2a,h2d,h2m,accProg,accInd,c2GxE,rgloc,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vE,vPhen,h2a,h2d,h2m,accProg,accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
      }
      if(length(random)==1){
        vRdm1 <- mAdd$var[random,1]
        if(GxE){
          vPhen <- vA + vGxE + vRdm1 + vE
          c2GxE <- vGxE/vPhen
          genParNames <- c("vA",paste0("v",random),"vGxE","vE","vPhen","h2a",paste0("c2",random),
                           "c2GxE","rgloc","accProg","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vRdm1 + vE
          genParNames <- c("vA",paste0("v",random),"vE","vPhen","h2a","c2Rdm1","accProg","accInd","CVgi%","CVe%","Mean")
        }
        
        h2a <- vA/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvRdm1 <- mAdd$var[random,2]
          if(GxE) SEvGxE <- mAdd$var["EnvProg",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vRdm1,vGxE,vE,vPhen,h2a,c2Rdm1,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvRdm1,SEvGxE,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=6,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vRdm1,vGxE,vE,vPhen,h2a,c2Rdm1,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,accProg,accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
      }
      if(length(random)==2){
        vRdm1 <- mAdd$var[random[1],1]
        vRdm2 <- mAdd$var[random[2],1]
        if(GxE){
          vPhen <- vA + vRdm1 + vRdm2 + vGxE + vE
          c2GxE <- vGxE/vPhen
          genParNames <- c("vA",paste0("v",random[1]),paste0("v",random[2]),"vGxE","vE","vPhen","h2a",
                           paste0("c2",random[1]),paste0("c2",random[2]),"c2GxE","rgloc","accProg","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vRdm1 + vRdm2 + vE
          genParNames <- c("vA",paste0("v",random[1]),paste0("v",random[2]),"vE","vPhen","h2a",
                           "c2Rdm1","c2Rdm2","accProg","accInd","CVgi%","CVe%","Mean")
        }
        
        h2a <- vA/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvRdm1 <- mAdd$var[random[1],2]
          SEvRdm2 <- mAdd$var[random[2],2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vRdm1,vRdm2,vGxE,vE,vPhen,h2a,c2Rdm1,c2Rdm2,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvRdm1,SEvRdm2,SEvGxE,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,c2Rdm1,accProg,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=7,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
          
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vRdm1,vRdm2,vGxE,vE,vPhen,h2a,c2Rdm1,c2Rdm2,c2GxE,rgloc,accProg,accInd,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vRdm1,vE,vPhen,h2aSE$Estimate,c2Rdm1,c2Rdm1,accProg,accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }}
    }
    genPar[is.na(genPar)] <- " "
    
    # BLUP_dataframes ---------------------------------------------------------
    
    # Progeny BLUP
    progBLUP <- r2Prog %>% dplyr::select(Prog,a) %>% rename(Progeny=Prog) %>% arrange(desc(a))
    
    # Individual_BLUP_genetic_gain_and_effectve_population_size --------------- 
    
    if("sire"%in%names(data)){
      
      if(GxE){
        indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
          {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Env","Block","Sire","Progeny","Plot","Tree"),  sep="-x-") else separate(.,Ind,c("Seq","Env","Block","Sire","Progeny","Tree"), sep="-x-")} %>% 
          mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a)) %>% dplyr::select(-Seq)
      }else{
        indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
          {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Block","Sire","Progeny","Plot","Tree"),  sep="-x-") else separate(.,Ind,c("Seq","Block","Sire","Progeny","Tree"), sep="-x-")} %>% 
          mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a)) %>% dplyr::select(-Seq)
      }
    }else{
      
      if(GxE){
        indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
          {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Env","Block","Progeny","Plot","Tree"),  sep="-x-") else separate(.,Ind,c("Seq","Env","Block","Progeny","Tree"), sep="-x-")} %>% 
          mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a)) %>% dplyr::select(-Seq)
      }else{
        indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
          {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Block","Progeny","Plot","Tree"),  sep="-x-") else separate(.,Ind,c("Seq","Block","Progeny","Tree"), sep="-x-")} %>% 
          mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a)) %>% dplyr::select(-Seq)
      }
    }
    
    # Provenance BLUP
    if(any(random=="Proc")){
      procBLUP <- mAdd$ranef$Proc[[1]] %>% rownames_to_column(.,"Proc") %>% rename(BLUP=value) %>% arrange(desc(BLUP))
      indBLUP <- indBLUP %>% left_join(.,procBLUP[,1:2],by="Proc") %>% mutate(a=a+BLUP) %>% mutate("u+a"=Mean+a) %>% dplyr::select(.,-last_col()) %>% 
        arrange(desc(a))
      progBLUP <- progBLUP %>% left_join(.,data[,c("Prog","Proc")], by= c("Progeny"="Prog")) %>% unique() %>% left_join(.,procBLUP[,1:2], by="Proc") %>% 
        mutate(a=a+BLUP) %>% dplyr::select(-c(Proc,BLUP)) %>% arrange(desc(a))
      if(PxE){
        procgeBLUP <- mAdd$ranef$EnvProc %>% as.data.frame() %>% rownames_to_column() %>% separate(.,rowname,c("Env","Proc"),sep = "_x_") %>% 
          left_join(.,procBLUP[,1:2],by="Proc") %>% mutate("p+pe"=value+BLUP) %>% dplyr::select(-c(s.e.,BLUP,value)) %>% arrange(desc(.[,3])) %>% 
          split(.,f=~Env)
      }
    }
    
    if(GxE){
      
      # Progeny BLUP (g) plus ge
      if(GxE) geBLUP <- mAdd$ranef$EnvProg %>% as.data.frame() %>% rownames_to_column() %>% separate(.,rowname,c("Env","Progeny"),sep = "_x_") %>% 
          left_join(.,progBLUP,by="Progeny") %>% left_join(.,envMeans[,1:2], by="Env") %>% mutate(a=a/2, Mean=as.numeric(Mean),"g+ge"=a+value,"g+ge+u"=a+value+Mean) %>% 
          dplyr::select(-c(s.e.,a,value,Mean)) %>% arrange(desc(.[,3])) %>% split(.,f=~Env)
      
      # BLUP indexes ------------------------------------------------------------
      
      # MHVG
      n <- table(plyr::ldply(geBLUP, data.frame)[,-1]$Progeny) %>% as.data.frame
      dataMHVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% group_by(Progeny) %>% summarise(MHVG=sum(1/g.ge.u)) %>% left_join(.,n,by=c("Progeny"="Var1")) %>% 
        mutate(MHVG=Freq/MHVG) %>% dplyr::select(-Freq)
      
      # PRVG
      dataPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% mutate(Mean=as.numeric(Mean),prvgIndex=.[,4]/Mean) %>% 
        dplyr::select(-5) %>% group_by(Progeny) %>% summarise(PRVG=mean(prvgIndex)) %>% mutate(PRVG=PRVG*Mean)
      
      # MHPRVG
      dataMHPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% mutate(Mean=as.numeric(Mean),prvgIndex=1/(.[,4]/Mean)) %>% 
        group_by(Progeny) %>% summarise(MHPRVG=sum(prvgIndex)) %>% left_join(.,n,by=c("Progeny"="Var1")) %>% mutate(MHPRVG=Freq/MHPRVG*Mean) %>% 
        dplyr::select(-Freq)
      
      indexBLUP <- dataMHVG %>% mutate(r_MHVG=order(order(MHVG, decreasing = T))) %>% left_join(.,dataPRVG,by="Progeny") %>% 
        mutate(r_PRVG=order(order(PRVG, decreasing = T))) %>% left_join(.,dataMHPRVG,by="Progeny") %>% 
        mutate(r_MHPRVG=order(order(MHPRVG, decreasing = T))) %>% arrange(-MHPRVG) %>% as.data.frame()
    }
    
    
    # Genetic gain and Effective Population Size ------------------------------
    
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
    
    # Overlapping generations
    if(GxE){
      overBLUP <- progBLUP %>% {if (plotType=="LP") mutate(.,Env = 0, Block = 0, Plot = 0, Tree = 0, Cod = "Parent") %>% 
          relocate(Env,Block,Progeny,Plot,Tree,a,Cod) else mutate(.,Env = 0, Block = 0, Tree = 0, Cod = "Parent") %>% relocate(Env,Block,Progeny,Tree,a,Cod)} %>% 
        rbind(.,indBLUP[,colnames(.)]) %>% arrange(desc(a))
    }else{
      overBLUP <- progBLUP %>% {if (plotType=="LP") mutate(.,Block = 0, Plot = 0, Tree = 0, Cod = "Parent") %>% 
          relocate(Block,Progeny,Plot,Tree,a,Cod) else mutate(.,Block = 0, Tree = 0, Cod = "Parent") %>% relocate(Block,Progeny,Tree,a,Cod)} %>% 
        rbind(.,indBLUP[,colnames(.)]) %>% arrange(desc(a))
    }
    
    
    # optimize_selection -------------------------------------------------------
    
    if(!GxE){
      
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
      }}
    
    # final_output_list -------------------------------------------------------
    
    genParBLUP <- list()
    
    # exploratory_analysis -------------------------------------------------
    genParBLUP$expAnalysis$respMeans <- respMeans; if(GxE) genParBLUP$expAnalysis$envMeans <- envMeans; 
    if(any(random=="Proc")||any(fixed=="Proc")) genParBLUP$expAnalysis$procMeans <- procMeans <- procMeans; 
    genParBLUP$expAnalysis$repMeans <- repMeans; genParBLUP$expAnalysis$progMeans <- treatMeans; 
    
    # model_output_and_significance_tests -------------------------------------------------
    genParBLUP$Model$remlMethod <- method; genParBLUP$Model$remlModel <- mAdd; genParBLUP$Model$mSig$fixedSig <- anova(mSig)
    genParBLUP$Model$mSig$randomSig <- suppressMessages(lmerTest::ranova(mSig, reduce.terms = F))
    
    # genetic_parameters_accuracy_and_Blup -------------------------------------------------
    genParBLUP$genPar <- genPar; genParBLUP$blupAccuracy$progAccuracy <- r2Prog; genParBLUP$blupAccuracy$indAccuracy <- r2Ind_df 
    genParBLUP$BLUP$progBLUP <- progBLUP; if(GxE) genParBLUP$BLUP$geBLUP <- geBLUP; if(GxE) genParBLUP$BLUP$blupIndex <- indexBLUP; 
    if("Proc"%in%random) genParBLUP$BLUP$procBLUP <- procBLUP; if(PxE) genParBLUP$BLUP$procgeBLUP <- procgeBLUP;
    genParBLUP$BLUP$indBLUP <- indBLUP; genParBLUP$BLUP$overBLUP <- overBLUP
    
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
    
    if(GxE){
      vGxE <- mClone$var["EnvClone",1]
      nEnv <- length(unique(data$Env))
      rgloc <- vG/(vG+vGxE)
    }
    
    if(plotType=="LP"){
      nArv <- length(unique(data$Arv))
      vParc <- mClone$var["Parc",1]
      CVe = (sqrt((3*vG+vE)/nArv+vParc))/Mean*100
      
      if(length(random)==1){
        if(GxE){
          vPhen <- vG + vGxE + vParc + vE
          c2GxE <- vGxE/vPhen
          h2mc <- (vG) / (vG+(vGxE/nEnv)+(vParc/nRep)+vE/(nRep*nArv))
          genParNames <- c("vG","vGxE","vParc","vE","vPhen","h2G","h2mc","c2GxE","c2Parc","rgloc","accClone","CVg%","CVe%","Mean")
          
        }else{
          vPhen <- vG + vParc + vE
          h2mc <- (vG) / (vG+(vParc/nRep)+vE/(nRep*nArv))
          genParNames <- c("vG","vParc","vE","vPhen","h2G","h2mc","c2Parc","accClone","CVg%","CVe%","Mean")
        }
        c2Parc <- vParc/vPhen
        h2G <- vG/vPhen
        
        genParNames <- c("vG","vParc","vE","vPhen","h2G","h2mc","c2Parc","accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          if(GxE){
            SEvGxE <- mClone$var["EnvClone",2]
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4),
                                      c(vG,vGxE,vParc,vE),
                                      mClone$reml$invAI)
          }else{
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3),
                                      c(vG,vParc,vE),
                                      mClone$reml$invAI)
          } 
          SEvParc <- mClone$var["Parc",2]
          SEvE <- mClone$var["Residual",2]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vE,vPhen,h2G,h2mc,c2GxE,c2Parc,rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvParc,SEvE,NA,h2GSE,matrix(NA,nrow=8,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vParc,vE,vPhen,h2G,h2mc,c2Parc,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvParc,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vE,vPhen,h2G,h2mc,c2GxE,c2Parc,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vParc,vE,vPhen,h2G,h2mc,c2Parc,accClone,CVg,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==2){
        rdmIndex <- which(random!="Parc")
        if(GxE){
          vPhen <- vG + vGxE + vParc + vRdm1 + vE
          c2GxE <- vGxE/vPhen
          genParNames <- c("vG","vGxE","vParc",paste0("v",random[rdmIndex[1]]),"vE","vPhen","h2G","c2GxE","c2Parc",
                           paste0("c2",random[rdmIndex[1]]),"rgloc","accClone","CVg%","CVe%","Mean")
        }else{
          vPhen <- vG + vParc + vRdm1 + vE
          genParNames <- c("vG","vParc",paste0("v",random[rdmIndex[1]]),"vE","vPhen","h2G","c2Parc",
                           paste0("c2",random[rdmIndex[1]]),"accClone","CVg%","CVe%","Mean")
        }
        vRdm1 <- mClone$var[random[2],1]
        vPhen <- vG + vParc + vRdm1 + vE
        h2G <- vG/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          SEvParc <- mClone$var["Parc",2]
          SEvRdm1 <- mClone$var[random[rdmIndex[1]],2]
          if(GxE){
            SEvGxE <- mClone$var["EnvProg",2]
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4+x5),
                                      c(vG,vGxE,vParc,vRdm1,vE),
                                      mClone$reml$invAI) 
          }else{
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4),
                                      c(vG,vParc,vRdm1,vE),
                                      mClone$reml$invAI) 
          } 
          SEvE <- mClone$var["Residual",2]
          
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vRdm1,vE,vPhen,h2G,c2GxE,c2Parc,c2Rdm1,rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvParc,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=8,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2G,c2Parc,c2Rdm1,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvParc,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vRdm1,vE,vPhen,h2G,c2GxE,c2Parc,c2Rdm1,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2G,c2Parc,c2Rdm1,accClone,CVg,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==3){
        rdmIndex <- which(random!="Parc")
        vRdm1 <- mClone$var[random[rdmIndex[1]],1]
        vRdm2 <- mClone$var[random[rdmIndex[2]],1]
        if(GxE){
          vPhen <- vG + vGxE + vParc + vRdm1 + vRdm2 + vE
          c2GxE <- vGxE/vPhen
          genParNames <- c("vG","vGxE","vParc",paste0("v",random[rdmIndex[1]]),paste0("v",random[rdmIndex[2]]),"vE","vPhen","h2G","c2GxE",
                           "c2Parc",paste0("c2",random[rdmIndex[1]]),paste0("c2",random[rdmIndex[2]]),"rgloc","accClone","CVg%","CVe%","Mean")
        }else{
          vPhen <- vG + vParc + vRdm1 + vRdm2 + vE
          genParNames <- c("vG","vParc",paste0("v",random[rdmIndex[1]]),paste0("v",random[rdmIndex[2]]),"vE","vPhen","h2G",
                           "c2Parc",paste0("c2",random[rdmIndex[1]]),paste0("c2",random[rdmIndex[2]]),"rgloc","accClone","CVg%","CVe%","Mean")
        }
        
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        h2G <- vG/vPhen
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          SEvParc <- mClone$var["Parc",2]
          if(GxE){
            SEvGxE <- mClone$var["EnvClone",2]
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4+x5+x6),
                                      c(vG,vGxE,vParc,vRdm1,vRdm2,vE),
                                      mClone$reml$invAI) 
          }else{
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4+x5),
                                      c(vG,vParc,vRdm1,vRdm2,vE),
                                      mClone$reml$invAI) 
          }
          SEvRdm1 <- mClone$var[random[2],2]
          SEvRdm2 <- mClone$var[random[3],2]
          SEvE <- mClone$var["Residual",2]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2GxE,c2Parc,c2Rdm1,c2Rdm2,rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvParc,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE,matrix(NA,nrow=9,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2Parc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvParc,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=7,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2GxE,c2Parc,c2Rdm1,c2Rdm2,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2Parc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
      }}
    
    if(plotType=="STP"){
      CVe = sqrt(3*vG+vE)/Mean*100
      
      if(length(random)==0){
        if(GxE){
          vPhen <- vG + vGxE + vE
          c2GxE <- vGxE/vPhen
          h2mc <- (vG) / (vG+(vGxE/nEnv)+(vE/nRep))
          genParNames <- c("vG","vGxE","vE","vPhen","h2G","h2mc","c2GxE","rgloc","accClone","CVg%","CVe%","Mean")
        }else{
          vPhen <- vG + vE
          h2mc <- (vG) / (vG+(vE/nRep))
          genParNames <- c("vG","vE","vPhen","h2G","h2mc","accClone","CVg%","CVe%","Mean")
        }
        h2G <- vG/vPhen
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          if(GxE){
            SEvGxE <- mClone$var["EnvClone",2]
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3),
                                      c(vG,vGxE,vE),
                                      mClone$reml$invAI)
          }else{
            h2GSE <- msm::deltamethod(~ x1/(x1+x2),
                                      c(vG,vE),
                                      mClone$reml$invAI)
          } 
          SEvE <- mClone$var["Residual",2]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vE,vPhen,h2G,h2mc,c2GxE,rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvE,NA,h2GSE,matrix(NA,nrow=7,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vE,vPhen,h2G,h2mc,c2GxE,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==1){
        vRdm1 <- mClone$var[random,1]
        if(GxE){
          vPhen <- vG + vGxE + vRdm1 + vE
          c2GxE <- vGxE/vPhen
          h2mc <- (vG) / (vG+(vGxE/nEnv)+(vE/nRep))
          genParNames <- c("vG","vGxE",paste0("v",random),"vE","vPhen","h2G","h2mc","c2GxE",paste0("c2",random),"rgloc","accClone","CVg%","CVe%","Mean")
        }else{
          vPhen <- vG + vRdm1 + vE
          h2mc <- (vG) / (vG+(vE/nRep))
          genParNames <- c("vG",paste0("v",random),"vE","vPhen","h2G","h2mc",paste0("c2",random),"accClone","CVg%","CVe%","Mean")
        }
        
        h2G <- vG/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          if(GxE){
            SEvGxE <- mClone$var["EnvClone",2]
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4),
                                      c(vG,vGxE,vRdm1,vE),
                                      mClone$reml$invAI)
          }else{
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3),
                                      c(vG,vRdm1,vE),
                                      mClone$reml$invAI)
          }
          SEvRdm1 <- mClone$var[random,2]
          SEvE <- mClone$var["Residual",2]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vRdm1,vE,vPhen,h2G,h2mc,c2GxE,c2Rdm1,rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=8,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vRdm1,vE,vPhen,h2G,h2mc,c2Rdm1,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vRdm1,vE,vPhen,h2G,h2mc,c2GxE,c2Rdm1,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vRdm1,vE,vPhen,h2G,h2mc,c2Rdm1,accClone,CVg,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }  
      if(length(random)==2){
        vRdm1 <- mClone$var[random[1],1]
        vRdm2 <- mClone$var[random[2],1]
        if(GxE){
          vPhen <- vG + vGxE + vRdm1 + vRdm2 + vE
          c2GxE <- vGxE/vPhen
          h2mc <- (vG) / (vG+(vGxE/nEnv)+(vE/nRep))
          genParNames <- c("vG","vGxE",paste0("v",random[1]),paste0("v",random[2]),"vE","vPhen","h2G","h2mc","c2GxE",
                           paste0("c2",random[1]),paste0("c2",random[2]),"rgloc","accClone","CVg%","CVe%","Mean")
        }else{
          vPhen <- vG + vRdm1 + vRdm2 + vE
          h2mc <- (vG) / (vG+(vE/nRep))
          genParNames <- c("vG",paste0("v",random[1]),paste0("v",random[1]),"vE","vPhen","h2G","h2mc",
                           paste0("c2",random[1]),paste0("c2",random[2]),"accClone","CVg%","CVe%","Mean")
        }
        h2G <- vG/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        
        if(method=="ai"){
          SEvG <- mClone$var["Clone",2]
          if(GxE){
            SEvGxE <- mClone$var["EnvClone",2]
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4+x5),
                                      c(vG,vGxE,vRdm1,vRdm2,vE),
                                      mClone$reml$invAI) 
          }else{
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4),
                                      c(vG,vRdm1,vRdm2,vE),
                                      mClone$reml$invAI) 
          } 
          SEvRdm1 <- mClone$var[random[1],2]
          SEvRdm2 <- mClone$var[random[2],2]
          SEvE <- mClone$var["Residual",2]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2GxE,c2Rdm1,c2Rdm2,rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE,matrix(NA,nrow=9,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vG,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2GxE,c2Rdm1,c2Rdm2,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
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
    if("Proc"%in%random){
      procBLUP <- mClone$ranef$Proc %>% as.data.frame() %>% rownames_to_column(.,"Proc") %>% rename(BLUP=value) %>% arrange(desc(BLUP))
      cloneBLUP <- cloneBLUP %>% left_join(.,data[,c("Clone","Proc")], by= "Clone") %>% unique() %>% left_join(.,procBLUP[,1:2], by="Proc") %>% 
        mutate(g=g+BLUP) %>% dplyr::select(-c(Proc,BLUP)) %>% arrange(desc(g))
      if(PxE) procgeBLUP <- mClone$ranef$EnvProc %>% as.data.frame() %>% rownames_to_column() %>% separate(.,rowname,c("Env","Proc"),sep = "_x_") %>% 
        left_join(.,procBLUP[,1:2],by="Proc") %>% mutate("p+pe"=value+BLUP) %>% dplyr::select(-c(s.e.,BLUP,value)) %>% arrange(desc(.[,3])) %>% 
        split(.,f=~Env)
    }
    
    if(GxE){
      # Clone BLUP (g) plus ge
      geBLUP <- mClone$ranef$EnvClone[[1]] %>% rownames_to_column() %>% separate(.,rowname,c("Env","Clone"),sep = "_x_") %>% 
        left_join(.,cloneBLUP,by="Clone") %>% left_join(.,envMeans[,1:2], by="Env") %>% mutate(Mean=as.numeric(Mean),"g+ge"=g+value,"g+ge+u"=g+value+Mean) %>% 
        dplyr::select(-c(s.e.,g,value,Mean)) %>% arrange(desc(.[,3])) %>% split(.,f=~Env)
      
      # BLUP indexes ------------------------------------------------------------
      
      # MHVG
      n <- table(plyr::ldply(geBLUP, data.frame)[,-1]$Clone) %>% as.data.frame
      dataMHVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% group_by(Clone) %>% summarise(MHVG=sum(1/g.ge.u)) %>% left_join(.,n,by=c("Clone"="Var1")) %>% 
        mutate(MHVG=Freq/MHVG) %>% dplyr::select(-Freq)
      
      # PRVG
      dataPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% 
        mutate(Mean=as.numeric(Mean),prvgIndex=.[,4]/Mean) %>% dplyr::select(-5) %>% group_by(Clone) %>% 
        summarise(PRVG=mean(prvgIndex)) %>% mutate(PRVG=PRVG*Mean)
      
      # MHPRVG
      dataMHPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% 
        mutate(Mean=as.numeric(Mean),prvgIndex=1/(.[,4]/Mean)) %>% group_by(Clone) %>% summarise(MHPRVG=sum(prvgIndex)) %>% 
        left_join(.,n,by=c("Clone"="Var1")) %>% mutate(MHPRVG=Freq/MHPRVG*Mean) %>% 
        dplyr::select(-Freq)
      
      indexBLUP <- dataMHVG %>% mutate(r_MHVG=order(order(MHVG, decreasing = T))) %>% left_join(.,dataPRVG,by="Clone") %>% 
        mutate(r_PRVG=order(order(PRVG, decreasing = T))) %>% left_join(.,dataMHPRVG,by="Clone") %>% 
        mutate(r_MHPRVG=order(order(MHPRVG, decreasing = T))) %>% arrange(-MHPRVG) %>% as.data.frame()
    }
    
    # final_output_list -------------------------------------------------------
    genParBLUP <- list()
    
    # exploratory_analysis -------------------------------------------------
    genParBLUP$expAnalysis$respMeans <- respMeans; if(GxE) genParBLUP$expAnalysis$envMeans <- envMeans; 
    if(any(random=="Proc")||any(fixed=="Proc")) genParBLUP$expAnalysis$procMeans <- procMeans; 
    genParBLUP$expAnalysis$repMeans <- repMeans; genParBLUP$expAnalysis$cloneMeans <- treatMeans
    
    # model_output_and_significance_tests -------------------------------------------------
    genParBLUP$Model$remlMethod <- method; genParBLUP$Model$remlModel <- mClone; genParBLUP$Model$mSig$fixedSig <- anova(mSig)
    genParBLUP$Model$mSig$randomSig <- suppressMessages(lmerTest::ranova(mSig, reduce.terms = F))
    
    # genetic_parameters_accuracy_and_Blup -------------------------------------------------
    genParBLUP$genPar <- genPar; genParBLUP$blupAccuracy$cloneAccuracy <- r2Clone; genParBLUP$BLUP$cloneBLUP <- cloneBLUP; 
    if(GxE) genParBLUP$BLUP$geBLUP <- geBLUP; if(GxE) genParBLUP$BLUP$blupIndex <- indexBLUP
    if("Proc"%in%random) genParBLUP$BLUP$procBLUP <- procBLUP; if(PxE) genParBLUP$BLUP$procgeBLUP <- procgeBLUP;
    
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
    
    if(GxE){
      if("Rep"%in%fixed){
        fixed <- c("Env","EnvRep",fixed[which(fixed!="Rep")])
      }
      
      if(!is.null(random)){
        strFix <- paste(paste0("f",fixed),collapse=" + ")
        strRnd <- paste(paste0("r",random),collapse=" + ")
        if(treatment=="Prog"){
          strMod <- paste0(varResp," = ",strFix," + ","r",treatment," + rEnvProg"," + ",strRnd)
        }else{
          strMod <- paste0(varResp," = ",strFix," + ","r",treatment," + rEnvClone"," + ",strRnd)
        }
      }else{
        strFix <- paste(paste0("f",fixed),collapse=" + ")
        if(treatment=="Prog"){
          strMod <- paste0(varResp," = ",strFix," + r",treatment," + r","EnvProg")
        }else{
          strMod <- paste0(varResp," = ",strFix," + r",treatment," + r","EnvClone")
        }
      }
    }else{
      if(!is.null(random)){
        strFix <- paste(paste0("f",fixed),collapse="+")
        strRnd <- paste(paste0("r",random),collapse="+")
        strMod <- paste0(varResp,"=",strFix,"+","r",treatment,"+",strRnd)
      }else{
        strFix <- paste(paste0("f",fixed),collapse="+")
        strMod <- paste0(varResp,"=",strFix,"+","r",treatment)
      }
    }
    
    # sink output
    
    sink(genParFile, append=TRUE, type = "output")
    
    cat("\n----------------------------------------------------------------------------------------------\n")
    cat("                                        |genBLUP Analysis| \n"                                     )
    if(GxE){cat("                                 Genotype x Environment Analysis \n"                       )
      cat("Sites: \n")
      print(levels(data$Env))
    }else{
      cat("                                       Individual Analysis \n"                                   )
    }
    cat("\n")
    cat("Date and Time: ")
    cat(format(Sys.time()))
    cat("\n")
    cat("Model: ")
    cat(strMod)
    cat("\n")
    cat("Controls Excluded: ")
    print(excludeControl)
    cat("----------------------------------------------------------------------------------------------\n\n")
    cat("---------------------------------------- Exploratory Analysis --------------------------------\n\n")
    print(respMeans)
    cat("\n")
    if(GxE){
      cat("------------------------------------------ Environmental  ------------------------------------\n\n")
      print(envMeans)
      cat("\n")
    }
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
      if(PxE){
        cat("-------------------------------- Provenance BLUP - Each Environments ----------------------------\n\n")
        print(procgeBLUP)
        cat("\n")
      }
    }
    cat("------------------------------------------ Treatment BLUP ------------------------------------\n\n")
    cat("\n")
    if(treatment=="Clone"){
      print(cloneBLUP)}
    if(treatment=="Prog"){
      print(progBLUP)
    }
    if(GxE){
      cat("\n")
      cat("----------------------------------- Treatment - Each environment -----------------------------\n\n")
      print(geBLUP)
      cat("\n")
      cat("----------------------------------------- BLUP Indexes -------------------------------------\n\n")
      print(indexBLUP)
    }
    if(treatment=="Prog"){
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
    print(suppressMessages(lmerTest::ranova(mSig)))
    
    sink()
    
    # optimized BLUP #
    
    if(!GxE & treatment=="Prog" & optimizeSelection==TRUE){
      
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