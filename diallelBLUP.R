diallelBLUP <- function(data, varResp, plotType=NULL, envCol=NULL, repCol=NULL, plotCol=NULL, damCol=NULL, sireCol=NULL, crossCol=NULL,
                        provCol=NULL, fixed = "Rep", random = NULL, method = "ai", GxE = F, PxE = F, excludeControl = NULL, 
                        gp_digits, codPerc = NULL, optimizeSelection = FALSE, maxIndProgeny = NULL, maxProgenyBlock = NULL, 
                        excludeCod = NULL, directory = NULL){
  
  # loading packages --------------------------------------------------------
  
  pacman::p_load(crayon,lme4,stringr,tidyverse,ggroups,dplyr,msm,breedR)
  
  # stops and warnings ------------------------------------------------------
  
  if(!plotType %in% c("LP","STP")){
    stop("ERROR: The argument 'plotType' must be 'LP' or 'STP'")
  } 
  
  if(!varResp %in% names(data)){
    stop("ERROR: varResp non-existent, please check your data")
  }
  
  if(!is.numeric(data[,varResp])){
    stop("ERROR: varResp is not numeric, please check your data")
  }
  
  if(plotType=="LP"&!"Plot"%in%random){
    cat(paste0(bold("WARNING:"), " The plotType argument was set to ", bold(cyan("LP")), " but nothing was especified in ", 
               bold('plotCol'), " argument. The function will build it using rep:cross interactions \n"))
  }
  if(PxE&!GxE){
    stop("ERROR: Argument 'GxE' is FALSE. You need to run a GxE analysis to estimate a Provenance x Environment effect")
  }
  
  #certifying that rep is in fixed or random effects if specified in repCol
  if(!is.null(repCol)){
    if(!repCol%in%c(fixed,random)){
      stop("ERROR: If you specify block replicates column at repCol argument, you need to specify it at fixed or random effects")
    }}
  
  #certifying that provenance is in fixed or random effects if specified in repCol
  if(!is.null(provCol)){
    if(!provCol%in%c(fixed,random)){
      stop("ERROR: If you specify provenance effect at provCol argument, you need to specify it at fixed or random effects")
    }}
  
  if(!is.null(excludeControl)&!any(data[,crossCol]%in%excludeControl)){
    stop(paste0("ERROR: The specified controls in excludeControl argument are not in the ",crossCol," column of your dataset"))
  }
  
  # creating_columns --------------------------------------------------------
  
  #basic columns
  data$sire <- data[,sireCol]
  data$dam <- data[,damCol]
  data$family <- data[,crossCol]
  
  #creating Env col if envCol is specified
  if(GxE) data$Env <- data[,envCol]
  
  #changing repCol to "Rep" and adding it in fixed or random (user specified)
  if(!is.null(repCol)){
    data$Rep <- data[,repCol]
    if(repCol%in%fixed){
      fixed <- c("Rep",fixed)
      if(any(duplicated(fixed))){
        fixed <- fixed[-which(duplicated(fixed))]
      }else{
        fixed <- fixed[-which(fixed==repCol)]
      }}
    if(repCol%in%random){
      random <- c("Rep",random)
      if(any(duplicated(random))){
        random <- random[-which(duplicated(random))]
      }else{
        random <- random[-which(random==repCol)]
      }}
  }
  
  #changing provCol to "prov" and adding it in fixed or random (user specified)
  if(!is.null(provCol)){
    data$prov <- data[,provCol]
    if(provCol%in%fixed){
      fixed <- c("prov",fixed)
      if(any(duplicated(fixed))){
        fixed <- fixed[-which(duplicated(fixed))]
      }else{
        fixed <- fixed[-which(fixed==provCol)]
      }}
    if(provCol%in%random){
      random <- c("prov",random)
      if(any(duplicated(random))){
        random <- random[-which(duplicated(random))]
      }else{
        random <- random[-which(random==provCol)]
      }}
  }
  
  # messages ----------------------------------------------------------------
  
  if(GxE) cat("diallelBLUP - Genotype x Environment Analysis: \n") else cat("diallelBLUP - Individual site Analysis: \n")
  
  if(method=="ai"||is.null(method)){
    cat("AI-REML algorithm was selected\n")
  }else{
    cat("EM-REML algorithm was selected\n")
  }
  # rearranging data --------------------------------------------------------
  
  # checking if any crosses have two or more dam and/or sires
  checkParent <- data %>% mutate(damSire=paste0(dam,"_x_",sire)) %>% dplyr::select(family,damSire) %>% 
    group_by(family) %>% unique() %>% filter(n()>1) %>% .$family %>% unique() %>% as.character()
  
  if(!is_empty(checkParent)){
    stop(paste0("The cross ", checkParent, " have two or more sires and/or dams, please check your data \n"))
  }
  
  #excluding controls in the family (cross) column
  if(!is.null(excludeControl)){
    data <- data %>% filter(!.$family%in%excludeControl)
    cat(paste(c("Controls:", excludeControl , "has been removed using excludeControl argument\n"))) 
  }
  
  # Creating Cod and other interactions if it doesn't exist
  if(!("Cod" %in% colnames(data))){
    data$Cod = " "
  }
  if(all(c("Cod1","Cod2")%in%colnames(data))){
    data$Cod <- paste0(data$Cod1,";",data$Cod2)
  }
  if(plotType=="LP"&!"Plot"%in%random){
    data$Plot <- paste0(data$Rep,"_x_",data$family)
  }
  
  if(GxE){
    #creating basic GxE interactions
    data <- data %>% mutate(EnvRep = factor(paste0(Env,"_x_",Rep)))
    data <- data %>% mutate(EnvFam = factor(paste0(Env,"_x_",family)))
    if(plotType=="LP") data <- data %>% mutate(Plot = factor(paste0(Env,"_x_",Rep,"_x_",family)))
    
    if("prov"%in%random&PxE){
      data <- data %>% mutate(EnvProv = factor(paste0(Env,"_x_",prov)))
      random <- c(random,"EnvProv") 
    }
  }
  
  data <- data %>% arrange("family")
  
  resp <- data[,which(names(data)==varResp)]
  data$resp <- resp
  
  # exploratory analysis ----------------------------------------------------  
  
  #create vector to round expAnalysis results
  roundVector <- c(rep(gp_digits,4),0,1,0,1)
  
  expAnalysis <- function(data){
    NAs <- sum(is.na(data))
    data <- data[!is.na(data)]
    exp <- round(c(mean(data),sd(data),min(data),max(data),length(data)), gp_digits)
    exp <- round(c(exp,(exp[2]/exp[1])*100), gp_digits)
    if (any(NAs)){
      exp <- c(exp,NAs)} else {exp <- c(exp,0)}
    exp <- round(c(exp,((exp[5]/(exp[5]+exp[7]))*100)), gp_digits)
    exp <- mapply(function(x,y){as.character(round(x,y))}, exp, roundVector)
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
  suppressMessages(treatMeans <- groupMeans(data,"family"))
  suppressMessages(repMeans <- groupMeans(data,"Rep"))
  
  if(GxE){
    envMeans <- groupMeans(data,"Env")
    repMeans <- groupMeans(data,"EnvRep") %>% separate(.,EnvRep, into = c("Env","Rep"), sep="_x_")
  }
  
  if(!is.null(provCol)){
    suppressMessages(provMeans <- groupMeans(data,"prov"))
  }
  
  if(!is.null(codPerc)){
    countTreat <- data %>% count(family) %>% setNames(c("family","nObs"))
    suppressMessages(codData <- separate_rows(data,Cod) %>% count(family,Cod) %>% filter(Cod!=".") %>% pivot_wider(., names_from = Cod,values_from = n) %>% 
                       left_join(countTreat,.) %>% relocate(nObs, .after = last_col()))
    nCod <- ncol(codData)+1
    
    for(i in codPerc){
      if(i %in% names(codData)){
        percCod <- setNames(data.frame(codData[,i]/codData$nObs*100),paste0(i,"%"))
        codData <- cbind(codData,percCod)
      }
    }
    
    codData[,nCod:ncol(codData)] <- round(codData[,nCod:ncol(codData)],2)
    
    if(any(codPerc%in%names(codData))){
      treatPerc <- codData %>% dplyr::select("family",(which(names(.) == "nObs") + 1):last_col())
      treatMeans <- cbind(treatMeans,treatPerc[,-1])
      treatMeans[is.na(treatMeans)] <- 0
    }
  }
  
  # factorizing_variables_and_treeCheck -------------------------------------
  
  # factorizing variables
  if(GxE) fct <- c("Env","sire","dam","family",fixed,random) else fct <- c("sire","dam","family",fixed,random)
  
  data[fct] <- lapply(data[fct], factor)
  
  # pedigree matrix ---------------------------------------------------------
  
  #creating Ind column
  
  if(GxE){
    if(plotType=="LP"){  
      data$Ind <- paste0(data$Env,"-x-",data$Rep,"-x-",data$sire,"-x-",data$dam,"-x-",data$family,"-x-",data$Plot,"-x-",data$Arv)
    }else{
      data$Ind <- paste0(data$Env,"-x-",data$Rep,"-x-",data$sire,"-x-",data$dam,"-x-",data$family,"-x-",data$Arv)  
    }
  }else{
    if(plotType=="LP"){  
      data$Ind <- paste0(data$Rep,"-x-",data$sire,"-x-",data$dam,"-x-",data$family,"-x-",data$Plot,"-x-",data$Arv)
    }else{data$Ind <- paste0(data$Rep,"-x-",data$sire,"-x-",data$dam,"-x-",data$family,"-x-",data$Arv)  
    }}
  
  if(any(duplicated(data$Ind))){
    stop("There are duplicated individuals in your data, please check it")
  }
  
  # pedigree ----------------------------------------------------------------
  
  #creating unique combinations for each individual
  data <- data %>% rownames_to_column("nId") %>% 
    mutate(nId=as.numeric(nId)+1000)
  
  #ancestors pedigree (n)
  nPed <- unlist(data[,c("sire","dam")], use.names = F) %>% unique() %>% 
    data.frame(parentCode=.) %>% drop_na() %>% mutate(nId=rownames(.)) %>% 
    relocate(nId,.before=parentCode)
  
  data <- data %>% left_join(nPed,by = c("sire"="parentCode")) %>% rename("nSire"="nId.y") %>% 
    left_join(nPed,by = c("dam"="parentCode")) %>% rename("nDam"="nId","nId"="nId.x") %>% 
    mutate(nIdd=nId)
  
  #building numeric pedigree
  ped <- nPed %>% mutate(sire=NA,dam=NA) %>% dplyr::select(-parentCode) %>%  
    rbind(.,data.frame(nId=data$nId,sire=as.numeric(data$nSire),dam=as.numeric(data$nDam))) %>% 
    mutate(nId=as.numeric(nId))
  
  # uniSire <- unique(data$sire) %>% .[!is.na(.)]
  # uniDam <- unique(data$dam) %>% .[!is.na(.)]
  # nInd <- seq(1001,1000+nrow(data))
  # 
  # uniPar <- data.frame(ID=unique(c(uniSire,uniDam))) %>% mutate(nPar=seq(1,length(ID)))
  # 
  # if("grandmother"%in%names(data)){
  #   
  #   #cat("Grandmothers were found and considered in the pedigree\n")
  #   
  #   #uniqueGD <- unique(data$grandmother) %>% .[!is.na(.)]
  #   #nGD <- data.frame(GD = unique(na.omit(data$GD)), nGD = seq(1,uniqueGD))
  #   
  #   #ancPed <- data.frame(idNum=nGD$nGD,sire=0,dam=0)
  #   
  #   #nProg <- data.frame(family = sort(unique(data$family)), nProg = as.numeric(seq(uniqueGD+1,uniqueProg+uniqueGD)))
  # }else{
  #   uniPar <- data.frame(ID=unique(c(uniSire,uniDam))) %>% mutate(nPar=seq(1,length(ID)))
  #   
  #   data <- data %>% mutate(idNum=as.numeric(nInd)) %>% dplyr::left_join(.,uniPar,by=c("sire"="ID")) %>% rename(.,"nSire"="nPar") %>% 
  #     dplyr::left_join(.,uniPar,by=c("dam"="ID")) %>% rename(.,"nDam"="nPar") %>% mutate(iddNum=idNum)
  #   
  #   ped <- data.frame(ind=uniPar$nPar,sire=0,dam=0) %>% rbind(data.frame(ind=data$idNum,sire=data$nSire,dam=data$nDam)) %>% 
  #     mutate_all(~replace(., is.na(.), 0))
  # }
  
  # dominance matrix --------------------------------------------------------
  
  D <- ped %>% replace(is.na(.), 0) %>% AGHmatrix::Amatrix(data = ., ploidy = 2, dominance = T)
  D <- D[-c(1:nrow(nPed)), -c(1:nrow(nPed))]
  D <- as(D, "dgCMatrix")
  
  Z <- as(Matrix::sparse.model.matrix(~ 0 + nIdd, data %>% 
                                        mutate(nIdd = as.factor(nIdd))), "indMatrix")
  
  # statistical modelling ---------------------------------------------------
  
  #creating all possible object formulas
  if(GxE){
    
    #changing Rep to EnvRep
    if("Rep"%in%fixed){
      fixed[which(fixed=="Rep")] <- "EnvRep"
    }
    if("Rep"%in%random){
      random[which(random=="Rep")] <- "EnvRep"
    }
    
    #fixed effects
    if(!is.null(fixed)){
      fixef <- as.formula(paste0("resp ~ 1 +", paste(fixed, collapse=" + ")))
    }else{
      fixef <- as.formula(paste0("resp ~ 1"))
    }
    
    # random effects
    if(!is.null(random)&"Plot"%in%random){
      random <- c(random[which(random=="Plot")],random[which(random!="Plot")])
    }
    if(!is.null(random)){
      ranef <- as.formula(paste0("~ EnvFam","+", paste(random, collapse=" + ")))
      ranef2 <- as.formula(paste0("~ overlay(dam,sire) + family + EnvFam","+", paste(random, collapse=" + ")))
    }else{
      ranef <- as.formula("~ EnvFam")
      ranef2 <- as.formula(paste0("~ overlay(dam,sire) + family + EnvFam"))
    }
  }else{
    #fixed effects
    if(!is.null(fixed)){
      fixef <- as.formula(paste0("resp ~ 1 +", paste(fixed, collapse=" + ")))
    }else{
      fixef <- as.formula(paste0("resp ~ 1"))
    }
    
    # random effects
    if(!is.null(random)&"Plot"%in%random){
      random <- c(random[which(random=="Plot")],random[which(random!="Plot")])
    }
    if(!is.null(random)){
      ranef <- as.formula(paste0("~", paste(random, collapse=" + ")))
      ranef2 <- as.formula(paste0("~ overlay(dam,sire) + family","+", paste(random, collapse=" + ")))
    }else{
      ranef <- NULL
      ranef2 <- as.formula(paste0("~ overlay(dam,sire) + family"))
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
      tryCatch(mAdd <- remlf90(fixef,
                               random = ranef,
                               genetic = list(model = "add_animal",
                                              pedigree = ped,id = "nId"),
                               generic = list(dominance = list(Z, D)),
                               data = data,method = method,
                               progsf90.options = 'EM-REML 10'), error = function(e) {convergence <<- TRUE}))
    if(!exists("mAdd")){
      cat("AI-REML algorithm got an error, switching to EM-REML\n")
      method = "em"
    }else{
      convergence <- ifelse(any(is.na(mAdd$var[,1])),TRUE,FALSE)
      if(convergence){
        cat("AI-REML algorithm did not converge, switching to EM-REML\n")
        method = "em"
      }}
    suppressMessages(tryCatch(mCruz <- sommer::mmer(fixef,
                                                    random= ranef2,
                                                    rcov=~units,
                                                    data=data, verbose=F, date.warning=F), 
                              error = function(e) {convergence <<- TRUE}))
    
    #significance of effects
    
    if(GxE){
      if(is.null(random)){
        lmerModel <- as.formula(paste0(format(fixef)," + ", 
                                       paste0("(1|family) + (1|EnvFam)"))) 
        
      }else{
        lmerModel <- as.formula(paste0(format(fixef)," + ", 
                                       paste0("(1|family) + (1|EnvFam) +"),
                                       paste("(1|", random, collapse=") + "),")"))
      }
    }else{
      if(is.null(random)){
        lmerModel <- as.formula(paste0(format(fixef)," + ", 
                                       paste0("(1|family)"))) 
        
      }else{
        lmerModel <- as.formula(paste0(format(fixef)," + ", 
                                       paste0("(1|family) +"),
                                       paste("(1|", random, collapse=") + "),")"))
      }
    }
    suppressMessages(mSig <- lme4::lmer(lmerModel, data=data, control = control))
    
  }
  options(warn = defaultW)
  
  # reliability_accuracy_and_individual_blup -----------------------------------------
  
  # specific_combining_ability_SCA
  r2SCA <- data.frame(a=mCruz$U$family$resp,r2=1-(diag(mCruz$PevU$`family`$resp)/mCruz$sigmaVector["family.resp-resp"])) %>% 
    rownames_to_column("Cross") %>% mutate(Cross=gsub("family","",Cross),s.e.=sqrt(diag(mCruz$PevU$`family`$resp)),s.e.=as.numeric(s.e.)) %>% 
    relocate(s.e.,.before = r2) %>% arrange(desc(a)) %>% rename(SCA = a)
  accSCA <- mean(sqrt(r2SCA$r2), na.rm=T)
  
  # general_combining_ability_GCA
  r2GCA <- mAdd$ranef$genetic[[1]] %>% 
    mutate(r2=1-s.e.^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
    filter(row_number() <= max(nrow(nPed))) %>% add_column(Parent = nPed$parentCode, .before = "value") %>% rename(GCA = value) %>% 
    arrange(desc(GCA))
  accGCA <- mean(sqrt(r2GCA$r2), na.rm=T)
  
  #individual_blup
  r2Ind <- mAdd$ranef$genetic[[1]] %>% 
    mutate(r2=1-s.e.^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
    filter(row_number() > max(nrow(nPed))) %>% add_column(Ind = data$Ind, .before = "value") %>% filter(!is.na(resp)) %>% rename(a = value)
  accInd <- mean(sqrt(r2Ind$r2), na.rm=T)
  
  if("grandmother"%in%names(data)){
    # r2GD <- mAdd$ranef$genetic[[1]] %>% 
    #   mutate(r2=1-(s.e.)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
    #   filter(row_number() %in% nGD$nGD) %>% add_column(GD = nGD$GD, .before = "value") %>% mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
    #   rename(a = value)
  }
  
  # genetic_parameters ------------------------------------------------------
  
  #function to build a named list and facilitate genPar organization
  named.list <- function(...) {
    l <- list(...)
    exprs <- lapply(substitute(list(...))[-1], deparse)
    names(l) <- exprs
    l
  } 
  
  # extracting vc's and calculating genetic parameters ----------------------
  
  vA <- mAdd$var["genetic",1]
  vD <- mAdd$var["dominance",1]
  vFam <-  as.numeric(mCruz$sigmaVector["family.resp-resp"])
  vE <- mAdd$var["Residual",1]
  vGCA <- mCruz$sigmaVector["overlay(dam, sire).resp-resp"]
  vPhen <- mAdd$var %>% as.data.frame() %>% dplyr::select('Estimated variances') %>% sum()
  Mean <- mean(data$resp,na.rm=T)
  
  if(!is.null(repCol)) nRep <- length(unique(data[,repCol]))
  
  tmp_gpl_1 <- named.list(vA,vD,vFam,vE) #named list number one (common variance components)
  
  #calculating basic genetic parameters
  h2a <- vA/vPhen
  h2G <- (vA+vD)/vPhen
  h2d <- (0.75*vA) / (0.75*vA+vE)
  c2D <- vD/vPhen
  bakersRatio <- as.numeric(2*vGCA/(2*vGCA+vFam)) #Baker's ratio
  CVgi <- sqrt(vA)/Mean*100
  
  if(plotType=="STP"){
    h2m <- (vA+vD) / (vA+vD+vE/(nRep))
    CVe = sqrt((0.75*vA+vE))/Mean*100
  }else{
    nArv <- length(unique(data$Arv))
    h2m <- (vA+vD) / (vA+vD+(vPlot/nRep)+vE/(nRep*nArv))
    CVe = (sqrt((0.75*vA+vE)/nArv+vPlot))/Mean*100
  }
  
  tmp_gpl_2 <- named.list(h2a,h2G,h2m,h2d) #named list number two (heritabilities)
  tmp_gpl_3 <- named.list(bakersRatio,c2D) #named list number three (heritabilities)
  
  #additional variance components
  if(length(rownames(mAdd$var))>3){
    
    if(length(rownames(mAdd$var))==4){
      dummyVar <- matrix(data=c(0,0),nrow=1,ncol=2) %>% `rownames<-`("dummyVar")
      mAdd_var <- rbind(mAdd$var,dummyVar)
      addVc <- as.data.frame(mAdd_var[setdiff(rownames(mAdd_var),c("genetic","dominance","Residual")),]) %>% 
        slice_head()
    }else{
      addVc <- as.data.frame(mAdd$var[setdiff(rownames(mAdd$var),c("genetic","dominance","Residual")),])
    }
    
    addVc_list <- list() #list to store additional variance components
    addc2_list <- list() #list to store additional c2
    
    for(i in 1:nrow(addVc)){
      assign(paste0("v",rownames(addVc)[i]), addVc[i,1]) #extracting additional variance components
      assign(paste0("c2",rownames(addVc)[i]),get(paste0("v",rownames(addVc)[i]))/vPhen) #calculating c2 for every additional vc
      
      addVc_i <- paste0("v",rownames(addVc)[i])
      addc2_i <- paste0("c2",rownames(addVc)[i])
      
      addVc_list[[addVc_i]] <- get(addVc_i)
      addc2_list[[addc2_i]] <- get(addc2_i)
    }
    
    tmp_gpl_1 <- c(tmp_gpl_1,addVc_list,named.list(vPhen))
    tmp_gpl_3 <- c(tmp_gpl_3,addc2_list)
  }else{
    tmp_gpl_1 <- c(tmp_gpl_1,named.list(vPhen))
  }
  
  #if GxE
  if(GxE){
    nEnv <- length(unique(data[,envCol]))
    rgloc <- (vA+vD)/((vA+vD)+vEnvFam)
    nRep <- length(unique(data[,"EnvRep"]))
    
    tmp_gpl_3 <- c(tmp_gpl_3,named.list(rgloc)) #adding rgloc to tmp_gpl_3 list
    
    if(plotType=="STP"){
      h2m <- (vA+vD) / (vA+vD+(vEnvFam/nEnv)+vE/(nRep*nEnv))
    }else{
      h2m <- (vA+vD) / (vA+vD+(vEnvFam/nEnv)+(vPlot/nRep)+vE/(nRep*nEnv*nArv))
    }
  }
  
  tmp_gp_list <- c(tmp_gpl_1,tmp_gpl_2,tmp_gpl_3,
                   named.list(accSCA,accGCA,accInd,CVgi,CVe,Mean))
  
  genPar <- data.frame(Estimate=unlist(tmp_gp_list)) %>% round(.,gp_digits)
  
  if(method=="ai"){
    #common standard errors
    vec_se <- c(mAdd$var["genetic",2],mAdd$var["dominance",2],sqrt(diag(mCruz$sigmaSE))[[2]],
                mAdd$var["Residual",2])
    
    #standard error from additional variance components
    if(length(rownames(mAdd$var))>3){
      for(i in 1:nrow(addVc)){
        vec_se <- c(vec_se,addVc[i,2])
      }
    }
    #build vector of standard error and fill missing values with NA
    vec_se <- c(vec_se,NA,mAdd$funvars["sample sd",1])
    vec_se <- c(vec_se,rep(NA,nrow(genPar)-length(vec_se))) %>% round(.,gp_digits)
    
    genPar$SE <- vec_se 
    genPar[is.na(genPar)] <- " "
  }
  
  # BLUP_dataframes ---------------------------------------------------------
  
  # Individual_BLUP_genetic_gain_and_effectve_population_size --------------- 
  
  if(GxE){
    indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
      {if ("prov"%in%random||"prov"%in%fixed) dplyr::select(.,prov,Ind,resp,a,Cod) else dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>%
      mutate("u+a" = Mean+a, d=mAdd$ranef$dominance[[1]][,1], g=a+d) %>% 
      {if (plotType=="LP") separate(.,Ind,c("Env","Block","Sire","Dam","Cross","Plot","Tree"),  sep="-x-") else 
        separate(.,Ind,c("Env","Block","Sire","Dam","Cross","Tree"), sep="-x-")} %>% drop_na() %>% arrange(desc(a)) %>% relocate(Cod,.after = last_col())
  }else{
    indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
      {if ("prov"%in%random||"prov"%in%fixed) dplyr::select(.,prov,Ind,resp,a,Cod) else dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>%
      mutate("u+a" = Mean+a, d=mAdd$ranef$dominance[[1]][,1], g=a+d) %>% 
      {if (plotType=="LP") separate(.,Ind,c("Block","Sire","Dam","Cross","Plot","Tree"),  sep="-x-") else 
        separate(.,Ind,c("Block","Sire","Dam","Cross","Tree"), sep="-x-")} %>% drop_na() %>% arrange(desc(a)) %>% relocate(Cod,.after = last_col())
  }
  
  # Provenance BLUP
  if("prov"%in%random){
    provBLUP <- mAdd$ranef$prov[[1]] %>% rownames_to_column(.,"prov") %>% rename(BLUP=value) %>% arrange(desc(BLUP))
    indBLUP <- indBLUP %>% left_join(.,provBLUP[,1:2],by="prov") %>% mutate(a=a+BLUP,"u+a"=Mean+a, g=a+d) %>% dplyr::select(.,-last_col()) %>% 
      arrange(desc(a))
    # parentBLUP <- parentBLUP %>% left_join(.,data[,c("family","prov")], by= c("Parent"="family")) %>% unique() %>% left_join(.,procBLUP[,1:2], by="Proc") %>% 
    #   mutate(a=a+BLUP) %>% dplyr::select(-c(Proc,BLUP)) %>% arrange(desc(a))
    if(PxE){
      procgeBLUP <- mAdd$ranef$EnvProv %>% as.data.frame() %>% rownames_to_column() %>% separate(.,rowname,c("Env","Prov"),sep = "_x_") %>% 
        left_join(.,provBLUP[,1:2],by="Prov") %>% mutate("p+pe"=value+BLUP) %>% dplyr::select(-c(s.e.,BLUP,value)) %>% arrange(desc(.[,3])) %>% 
        split(.,f=~Env)
    }
  }
  
  if(GxE){
    
    # Cross BLUP plus ge
    geBLUP <- mCruz$U$EnvFam %>% as.data.frame() %>% rownames_to_column() %>% 
      separate(.,rowname,c("Env","Cross"),sep = "_x_") %>% 
      mutate(Env=gsub('EnvFam','',.$Env)) %>% rename(g=resp) %>% left_join(.,r2SCA[,1:2],by="Cross") %>% 
      left_join(.,envMeans[,1:2], by="Env") %>% mutate(Mean=as.numeric(Mean),"g+ge"=SCA+g,"g+ge+u"=SCA+g+Mean) %>% 
      dplyr::select(-c(SCA,g,Mean)) %>% arrange(desc(.[,3])) %>% split(.,f=~Env)
    
    # BLUP indexes ------------------------------------------------------------
    
    # MHVG
    n <- table(plyr::ldply(geBLUP, data.frame)[,-1]$Cross) %>% as.data.frame
    dataMHVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% group_by(Cross) %>% summarise(MHVG=sum(1/g.ge.u)) %>% 
      left_join(.,n,by=c("Cross"="Var1")) %>% mutate(MHVG=Freq/MHVG) %>% dplyr::select(-Freq)
    
    # PRVG
    dataPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% 
      mutate(Mean=as.numeric(Mean),prvgIndex=.[,4]/Mean) %>% 
      dplyr::select(-5) %>% group_by(Cross) %>% summarise(PRVG=mean(prvgIndex)) %>% mutate(PRVG=PRVG*Mean)
    
    # MHPRVG
    dataMHPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% 
      mutate(Mean=as.numeric(Mean),prvgIndex=1/(.[,4]/Mean)) %>% group_by(Cross) %>% 
      summarise(MHPRVG=sum(prvgIndex)) %>% left_join(.,n,by=c("Cross"="Var1")) %>% 
      mutate(MHPRVG=Freq/MHPRVG*Mean) %>% dplyr::select(-Freq)
    
    indexBLUP <- dataMHVG %>% mutate(r_MHVG=order(order(MHVG, decreasing = T))) %>% left_join(.,dataPRVG,by="Cross") %>% 
      mutate(r_PRVG=order(order(PRVG, decreasing = T))) %>% left_join(.,dataMHPRVG,by="Cross") %>% 
      mutate(r_MHPRVG=order(order(MHPRVG, decreasing = T))) %>% arrange(-MHPRVG) %>% as.data.frame()
  }
  
  
  # genetic gain and overlapping generation ---------------------------------
  
  # Absolute and percentual genetic gain (individual blup)
  gain <- matrix(nrow=nrow(indBLUP), ncol=1)
  for(i in 1:nrow(indBLUP)){
    gain[i,] <- mean(indBLUP$a[1:i])
  }
  indBLUP <- indBLUP %>% mutate(Gain=gain) %>% mutate("Gain%"=gain/Mean*100) %>% relocate(Cod,.after=last_col())
  
  # Overlapping generations
  overBLUP <- indBLUP %>% mutate(Parent=0) %>% relocate(Parent, .after = Block)
  if(GxE){
    overBLUP <- r2GCA %>% dplyr::select(-c(s.e.,r2)) %>% rename("a"="GCA") %>% {if (plotType=="LP") 
      mutate(.,Env = 0, Block = 0, Sire = 0, Dam = 0, Cross =0, Plot = 0, Tree = 0, Cod = "Parent") %>% 
        relocate(Env,Block,Parent,Plot,Sire,Dam,Cross,Tree,a,Cod) else 
          mutate(.,Env = 0, Block = 0, Sire = 0, Dam = 0, Cross =0, Tree = 0, Cod = "Parent") %>% 
        relocate(Env,Block,Parent,Sire,Dam,Cross,Tree,a,Cod)} %>% rbind(.,overBLUP[,colnames(.)]) %>% arrange(desc(a))
  }else{
    overBLUP <- r2GCA %>% dplyr::select(-c(s.e.,r2)) %>% rename("a"="GCA") %>% {if (plotType=="LP") 
      mutate(.,Block = 0, Sire = 0, Dam = 0, Cross = 0, Plot = 0, Tree = 0, Cod = "Parent") %>% 
        relocate(Block,Parent,Sire,Dam,Cross,Plot,Tree,a,Cod) else 
          mutate(.,Block = 0, Sire = 0, Dam = 0, Cross =0, Tree = 0, Cod = "Parent") %>% 
        relocate(Block,Parent,Sire,Dam,Cross,Tree,a,Cod)} %>% rbind(.,overBLUP[,colnames(.)]) %>% arrange(desc(a))
  }
  # Absolute and percentual genetic gain (overlapping)
  gain <- matrix(nrow=nrow(overBLUP), ncol=1)
  for(i in 1:nrow(overBLUP)){
    gain[i,] <- mean(overBLUP$a[1:i])
  }
  overBLUP <- overBLUP %>% mutate(Gain=gain) %>% mutate("Gain%"=gain/Mean*100) %>% relocate(Cod,.after=last_col())
  
  # final_output_list -------------------------------------------------------
  
  genParBLUP <- list()
  
  # exploratory_analysis -------------------------------------------------
  genParBLUP$expAnalysis$respMeans <- respMeans; if(GxE) genParBLUP$expAnalysis$envMeans <- envMeans; 
  if("prov"%in%c(fixed,random)) genParBLUP$expAnalysis$provMeans <- provMeans; 
  genParBLUP$expAnalysis$repMeans <- repMeans; genParBLUP$expAnalysis$crossMeans <- treatMeans; 
  
  # model_output_and_significance_tests -------------------------------------------------
  genParBLUP$Model$remlMethod <- method; genParBLUP$Model$remlModel <- mAdd; genParBLUP$Model$mSig$fixedSig <- anova(mSig)
  genParBLUP$Model$mSig$randomSig <- suppressMessages(lmerTest::ranova(mSig, reduce.terms = F))
  
  # genetic_parameters_accuracy_and_Blup -------------------------------------------------
  genParBLUP$genPar <- genPar; genParBLUP$blupAccuracy$indAccuracy <- r2Ind;
  genParBLUP$BLUP$GCA <- r2GCA; genParBLUP$BLUP$SCA <- r2SCA; if(GxE) genParBLUP$BLUP$geBLUP <- geBLUP; 
  if(GxE) genParBLUP$BLUP$blupIndex <- indexBLUP 
  if("prov"%in%c(fixed,random)) genParBLUP$BLUP$provBLUP <- provBLUP; if(PxE) genParBLUP$BLUP$provgeBLUP <- provgeBLUP;
  genParBLUP$BLUP$indBLUP <- indBLUP; genParBLUP$BLUP$overBLUP <- overBLUP
  
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
    
    #creating genParFile
    genParFile <- paste0(directory,"_",varResp,"_genPar",".txt")
    
    if (file.exists(genParFile)){
      file.remove(genParFile)
    }
    
    #building model string to show in genParFile
    if(!is.null(fixed)){
      strFix <- paste(paste0("f",fixed),collapse=" + ")
    }else{
      strFix <- "f1"
    }
    rVars <- c("Zu","Zd",all.vars(ranef))
    strRnd <- paste(paste0("r",rVars),collapse=" + ")
    
    strMod <- paste0(varResp," = ",strFix," + ",strRnd)
    
    
    sink(genParFile, append=TRUE, type = "output")
    
    cat("\n----------------------------------------------------------------------------------------------\n")
    cat("                                   |diallelBLUP Analysis| \n"                                      )
    if(GxE){cat("                             Genotype x Environment Analysis \n"                           )
      cat("Sites: \n")
      print(levels(data$Env))
    }else{
      cat("                                     Individual Analysis \n"                                     )
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
    cat("-------------------------------------- Exploratory Analysis ----------------------------------\n\n")
    print(respMeans)
    cat("\n")
    if(GxE){
      cat("---------------------------------------- Environmental ---------------------------------------\n\n")
      print(envMeans)
      cat("\n")
    }
    cat("---------------------------------------- Block repetition ------------------------------------\n\n")
    print(repMeans)
    cat("\n")
    cat("--------------------------------------------- Parent -----------------------------------------\n\n")
    print(treatMeans)
    cat("\n")
    if("prov"%in%c(fixed,random)){
      cat("------------------------------------------ Provenance -----------------------------------------\n\n")
      print(provMeans)
      cat("\n")
    }
    cat("--------------------------------------- Genetic Parameters ----------------------------------\n\n")
    print(genPar)
    cat("\n")
    cat("\n")
    cat("-------------------------------------------- |BLUP|  ----------------------------------------\n\n")
    cat("\n")
    if("prov"%in%random){
      cat("-------------------------------------- Provenance BLUP -----------------------------------------\n\n")
      print(provBLUP)
      cat("\n")
      if(PxE){
        cat("------------------------------- Provenance BLUP - Each Environments -----------------------------\n\n")
        print(provgeBLUP)
        cat("\n")
      }
    }
    cat("---------------------------------------- Parent (GCA) BLUP -----------------------------------\n\n")
    cat("\n")
    print(r2GCA)
    cat("\n")
    cat("--------------------------------------------- SCA BLUP ---------------------------------------\n\n")
    print(r2SCA)
    
    if(GxE){
      cat("\n")
      cat("---------------------------------- Treatment - Each environment -------------------------------\n\n")
      print(geBLUP)
      cat("\n")
      cat("---------------------------------------- BLUP Indexes -----------------------------------------\n\n")
      print(indexBLUP)
    }
    cat("\n")
    cat("----------------------------------------- Individual BLUP ---------------------------------------\n\n")
    print(indBLUP)
    cat("\n")
    cat("--------------------------------- Overlapping Generations Selection -----------------------------\n\n")
    print(overBLUP)
    
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
    cat("\n")
    cat("----------------------------------------- Individual accuracy --------------------------------\n\n")
    print(r2Ind)
    
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
    if(!is.null(fixed)){
      cat("------------------------------------------- Fixed effects ------------------------------------\n\n")
      print(anova(mSig))
    }else{
      cat("There is no fixed effects in our model.")
    }
    cat("----------------------------------------------------------------------------------------------\n\n")
    cat("------------------------------------------ Random effects ------------------------------------\n\n")
    print(suppressMessages(lmerTest::ranova(mSig)))
    
    sink()
    
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