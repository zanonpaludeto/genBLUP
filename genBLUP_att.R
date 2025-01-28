genBLUP <- function(data, varResp,envCol = NULL, treeCol = NULL, plotCol = NULL, repCol, matGenCol, matGen, provCol = NULL, 
                    plotType, fixed = "Rep",codCol = NULL, random = NULL, method = "ai", GxE = F, PxE = F, 
                    excludeControl = NULL, genPar_digits = 6, plotRanking = NULL,
                    codPerc = NULL, optimizeSelection = FALSE, maxIndProgeny = NULL, 
                    maxProgenyBlock = NULL, excludeCod = NULL, directory = NULL){
  
  # loading packages --------------------------------------------------------
  
  pacman::p_load(lme4,stringr,tidyverse,ggroups,dplyr,msm,lmerTest)
  
   # stops and warnings ------------------------------------------------------
  
  if(!matGen %in% c("family","clone")){
    stop("ERROR: The argument 'matGen' must be 'family' or 'clone'")
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
  
  if(matGen=="clone"&optimizeSelection==T){
    stop("ERROR: You can't optimize individual selection in a clonal analysis")
  }
  
  if(!is.null(repCol)){
    if(!repCol%in%fixed&!repCol%in%random){
      stop("ERROR: If you specify the block replicates column at repCol argument, you need 
           to specify it at fixed or random effects")
    }}
  
  if(!is.null(provCol)){
    if(!provCol%in%fixed&!provCol%in%random){
      stop("ERROR: If you specify a provenance effect column at provCol argument, you need 
           to specify it at fixed or random effects")
    }}
  
  if(optimizeSelection==T&(!is.numeric(maxIndProgeny)|!is.numeric(maxProgenyBlock))){
    stop("ERROR: If you want to optimize selection, please choose numeric parameters for 
         maxIndProgeny and maxProgenyBlock")
  }
  if(PxE&!GxE){
    stop("ERROR: Argument 'GxE' is FALSE. You need to run a GxE analysis to estimate a 
         Provenance x Environment effect")
  }
  if(PxE&!"Proc" %in% random){
    stop("ERROR: To estimate a Provenance x Environment effect, you need to specify 
         'Proc' in the random argument")
  }
  if("Proc"%in%random & !"Proc"%in%names(data)){
    stop("ERROR: You cannot try to estimate provenance effect without a column names 
         'Proc' in your data")
  }
  
  if(optimizeSelection&GxE){
    warning("WARNING: You can't optimize selection in a GxE analysis, the argument was 
            changed to FALSE")
    optimizeSelection=F
  }
  
  if(!is.null(excludeControl)&!any(data[,matGenCol]%in%excludeControl)){
    stop(paste0("ERROR: The specified controls in excludeControl argument are not in the ",matGenCol," column of your dataset"))
  }
  
  # creating_columns --------------------------------------------------------
  
  # matGen column
  if(matGen=="family"){
    data$family <- data[,matGenCol]
  }else{
    data$clone <- data[,matGenCol]
  }  
  
  # block replicates column (fixed)
  if(!is.null(repCol)){
    if(repCol%in%fixed){
      data$Rep <- data[,repCol]
      fixed <- c("Rep",fixed)
      if(any(duplicated(fixed))){
        fixed <- fixed[-which(duplicated(fixed))]
      }else{
        fixed <- fixed[-which(fixed==repCol)]
      }
    }}
  
  # block replicates column (random)
  if(!is.null(repCol)){
    if(repCol%in%random){
      data$Rep <- data[,repCol]
      random <- c("Rep",random)
      if(any(duplicated(random))){
        random <- random[-which(duplicated(random))]
      }else{
        random <- random[-which(random==repCol)]
      }
    }}
  
  #creating Env col if envCol is specified
  if(GxE) data$Env <- data[,envCol]
  
  #creating tree column if plotType is LP and treeCol is specified
  if(plotType=="LP") data$Arv <- data[,treeCol]
  
  #creating "Parc" column if plotCol is specified and plotType is LP
  if(plotType=="LP"&!is.null(plotCol)){
    data$Parc <- data[,plotCol]
    random <- c(random,"Parc")
  }
  
  #creating provenance column if 'provCol' is specified
  if(!is.null(provCol)){
    data$Proc <- data[,provCol]
    
    #correcting fixed effects if provCol is in fixed
    if(provCol%in%fixed){
      fixed <- fixed[-which(fixed==provCol)]
      fixed <- c(fixed,"Proc")
    }
    #correcting random effects if provCol is in random
    if(provCol%in%random){
      random <- random[-which(random==provCol)]
      random <- c(random,"Proc")
    }
  }
  
  #creating Cod column
  if(is.null(codCol)){
    data$Cod = " "
  }else{
    if(length(codCol)>1){
      data$Cod <- apply(data[,codCol], 1,paste, collapse = ";")
    }
    if(length(codCol)==1){
      data$Cod = data[,codCol]
    }
  }
  
  # messages ----------------------------------------------------------------
  
  if(GxE) cat("Genotype x Environment Analysis: \n") else cat("Individual site Analysis: \n")
  
  if(method=="ai"||is.null(method)){
    cat("AI-REML algorithm was selected\n")
  }else{
    cat("EM-REML algorithm was selected\n")
  }
  
  if(plotType=="LP"&is.null(plotCol)){
    cat(paste0("You choose to adjust an LP model without specifying 'plotCol' argument, 
               the function will automatically build Rep:", matGenCol," cross interactions 
               and proceed the analysis \n \n"))
    
    random <- c(random,"Parc")
  }
  
  if(!is.null(excludeControl)){
    data <- data %>% filter(!.[,matGen]%in%excludeControl)
    cat(paste(c("Controls:",excludeControl, "has been removed using excludeControl argument\n"))) 
  }
  
  # exploratory analysis ----------------------------------------------------
  
  # preparing data
  
  if(GxE){
    #creating basic GxE interactions
    if("Rep"%in%fixed|"Rep"%in%random){
      data <- data %>% mutate(EnvRep = factor(paste0(Env,"_x_",Rep)))
    }
    if(matGen=="family"){
      data <- data %>% mutate(EnvFam = factor(paste0(Env,"_x_",family)))
      if("Rep"%in%fixed){
        if(plotType=="LP") data <- data %>% mutate(Parc = factor(paste0(Env,"_x_",Rep,"_x_",family)))
      }
    }else{
      data <- data %>% mutate(EnvClone = factor(paste0(Env,"_x_",clone)))
      if("Rep"%in%fixed){
        if(plotType=="LP") data <- data %>% mutate(Parc = factor(paste0(Env,"_x_",Rep,"_x_",clone)))
      }
    }
    if("Proc"%in%random&PxE){
      data <- data %>% mutate(EnvProc = factor(paste0(Env,"_x_",Proc)))
      random <- c(random,"EnvProc")
    }
  }else{
    if(plotType=="LP"&"Parc"%in%random&is.null(plotCol)){
      if("Rep"%in%fixed){
        if(matGen=="family") data <- data %>% mutate(Parc = factor(paste0(Rep,"_x_",family)))
        if(matGen=="clone") data <- data %>% mutate(Parc = factor(paste0(Rep,"_x_",clone)))
      }}}
  
  data <- data %>% arrange(get(matGenCol))
  
  resp <- data[,which(names(data)==varResp)]
  data$resp <- resp
  
  #create vector to round expAnalysis results
  roundVector <- c(rep(genPar_digits,4),0,1,0,1)
  
  expAnalysis <- function(data){
    NAs <- sum(is.na(data))
    data <- data[!is.na(data)]
    exp <- round(c(mean(data),sd(data),min(data),max(data),length(data)), genPar_digits)
    exp <- round(c(exp,(exp[2]/exp[1])*100), genPar_digits)
    if (any(NAs)){
      exp <- c(exp,NAs)} else {exp <- c(exp,0)}
    exp <- round(c(exp,((exp[5]/(exp[5]+exp[7]))*100)), genPar_digits)
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
  suppressMessages(treatMeans <- groupMeans(data,matGenCol))
  
  suppressMessages(if("Rep"%in%fixed|"Rep"%in%random) repMeans <- groupMeans(data,"Rep"))
  
  if(GxE){
    envMeans <- groupMeans(data,"Env")
    if("Rep"%in%fixed|"Rep"%in%random){
      repMeans <- groupMeans(data,"EnvRep") %>% separate(.,EnvRep, into = c("Env","Rep"), sep="_x_")
    }
  }
  
  if(any(random=="Proc")||any(fixed=="Proc")){
    suppressMessages(procMeans <- groupMeans(data,"Proc"))
  }
  
  if(!is.null(codPerc)){
    countTreat <- data %>% dplyr::count(across(!!matGen)) %>% setNames(c(matGen,"nObs"))
    suppressMessages(codData <- separate_rows(data,Cod) %>% dplyr::count(across(!!matGen),Cod) %>% filter(Cod!=".") %>% 
                       pivot_wider(., names_from = Cod,values_from = n) %>% 
                       left_join(countTreat,.) %>% relocate(nObs, .after = last_col()))
    
    for(i in codPerc){
      if(i %in% names(codData)){
        percCod <- setNames(data.frame(codData[,i]/codData$nObs*100),paste0(i,"%")) %>% round(.,2)
        codData <- cbind(codData,percCod)
      }
    }
    if(any(codPerc%in%names(codData))){
      treatPerc <- codData %>% dplyr::select(as.name(matGen),(which(names(.) == "nObs") + 1):last_col())
      treatMeans <- cbind(treatMeans,treatPerc[,-1])
      treatMeans[is.na(treatMeans)] <- 0
    }
  }
  
  # factorizing_variables_and_treeCheck -------------------------------------
  
  # factorizing variables
  if(GxE) fct <- c("Env",matGen,fixed,random) else fct <- c(matGen,fixed,random)
  
  data[fct] <- lapply(data[fct], factor)
  
  # # treeCheck
  # if(plotType=="STP"){
  #   
  #   treeCheck <- function(data){
  #     
  #     arvList <- list()
  #     
  #     for(i in levels(factor(data$Rep))){
  #       if(matGen=="clone"){
  #         arvList[i] <- as.data.frame(ave(data[data$Rep==i,]$Clone==data[data$Rep==i,]$Clone, 
  #                                         data[data$Rep==i,]$Clone, FUN=cumsum))
  #       }else{
  #         arvList[i] <- as.data.frame(ave(data[data$Rep==i,]$family==data[data$Rep==i,]$family, 
  #                                         data[data$Rep==i,]$family, FUN=cumsum))
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
  
  if(matGen=="family"){
    
    #creating Ind column
    data$seq <- paste0(rep("I"),seq(from=1, to=nrow(data))) #creating sequential index to differentiate individuals
    
    if(GxE){
      if("sire"%in%names(data)){
        if(plotType=="LP"){  
          data$Ind <- paste(data$seq,data$Env,data$Rep,data$sire,data[,names(data) %in% matGen],data$Parc,data$Arv, 
                            sep = "-x-")
        }else{
          data$Ind <- paste(data$seq,data$Env,data$Rep,data$sire,data[,names(data) %in% matGen],data$Arv, sep = "-x-")  
        }
      }else{
        if(plotType=="LP"){  
          data$Ind <- paste(data$seq,data$Env,data$Rep,data[,names(data) %in% matGen],data$Parc,data$Arv, sep = "-x-")
        }else{
          data$Ind <- paste(data$seq,data$Env,data$Rep,data[,names(data) %in% matGen],data$Arv, sep = "-x-")  
        }
      }
      
    }else{
      if("sire"%in%names(data)){
        if(plotType=="LP"){  
          data$Ind <- paste(data$seq,data$Rep,data$sire,data[,names(data) %in% matGen],data$Parc,data$Arv, sep = "-x-")
        }else{data$Ind <- paste(data$seq,data$Rep,data$sire,data[,names(data) %in% matGen],data$Arv, sep = "-x-")  
        }
      }else{
        if(plotType=="LP"){  
          data$Ind <- paste(data$seq,data$Rep,data[,names(data) %in% matGen],data$Parc,data$Arv, sep = "-x-")
        }else{data$Ind <- paste(data$seq,data$Rep,data[,names(data) %in% matGen],data$Arv, sep = "-x-")  
        }
      }}
    
    #counting nFam (dams)
    uniqueFam <- length(na.omit(unique(data$family)))
    
    if("GD"%in%names(data)&"sire"%in%names(data)){
      
      cat("Grandparents and sires were found and considered in the pedigree\n")
      
      uniqueGD <- length(na.omit(unique(data$GD)))
      nGD <- data.frame(GD = unique(na.omit(data$GD)), nGD = seq(1,uniqueGD))
      uniqueSire <- length(na.omit(unique(data$sire)))
      nSire <- data.frame(sire = unique(na.omit(data$sire)), 
                          nSire = seq(uniqueGD+1,uniqueGD+uniqueSire))
      
      ancPed <- data.frame(idNum=c(nGD$nGD,nSire$nSire),sire=0,dam=0)
      
      nFamily <- data.frame(family = sort(unique(data$family)), 
                            nFamily = as.numeric(seq(uniqueGD+uniqueSire+1,uniqueFam+uniqueGD+uniqueSire)))
    }
    if(!"GD"%in%names(data)&"sire"%in%names(data)){
      
      cat("Sires were found and considered in the pedigree\n")
      
      uniqueSire <- length(na.omit(unique(data$sire)))
      nSire <- data.frame(sire = unique(na.omit(data$sire)), nSire = seq(1,uniqueSire))
      
      ancPed <- data.frame(idNum=nSire$nSire,sire=0,dam=0)
      
      nFamily <- data.frame(family = sort(unique(data$family)), 
                            nFamily = as.numeric(seq(uniqueSire+1,uniqueFam+uniqueSire)))
    }
    if("GD"%in%names(data)&!"sire"%in%names(data)){
      
      cat("Grandparents were found and considered in the pedigree\n")
      
      uniqueGD <- length(na.omit(unique(data$GD)))
      nGD <- data.frame(GD = unique(na.omit(data$GD)), nGD = seq(1,uniqueGD))
      
      ancPed <- data.frame(idNum=nGD$nGD,sire=0,dam=0)
      
      nFamily <- data.frame(family = sort(unique(data$family)), 
                            nFamily = as.numeric(seq(uniqueGD+1,uniqueFam+uniqueGD)))
    }
    if(!"GD"%in%names(data)&!"sire"%in%names(data)){
      nFamily <- data.frame(family = sort(unique(data$family)), 
                            nFamily = as.numeric(seq(1,uniqueFam)))
    }
    
    nInd <- seq(1001,1000+nrow(data))
    data <- data %>% mutate(idNum=as.numeric(nInd))
    data <- dplyr::left_join(data,nFamily,by="family")
    
    breedR.createPed <- function(pedData){
      
      if("sire"%in%names(data)){
        data <- dplyr::left_join(data,nSire,by="sire")
        ped <- data.frame(pedData,sire=data$nSire) %>% setNames(c("idNum","dam","sire")) %>% 
          relocate(dam, .after = sire)
      }else{
        ped <- data.frame(pedData,sire=0) %>% setNames(c("idNum","dam","sire")) %>% 
          relocate(dam, .after = sire)
      }
      
      #pedigree with grand dam information
      if(any(names(data)=="GD")){
        
        dfGD <- data[c("nFamily","GD")] %>% group_by(nFamily) %>% distinct() %>% 
          arrange(nFamily)
        if(nrow(dfGD)>nrow(nFamily)){
          stop("ERROR: There is families with more than one grandparent assisgned 
               in the dataset, please check your data")
        }
        ped <- dfGD %>% left_join(.,nGD, by="GD") %>% dplyr::select(-GD) %>% 
          mutate(sire=0, .after = nFamily) %>% setNames(c("idNum","sire","dam")) %>% 
          bind_rows(ancPed,.,ped)
      }
      
      if(any(duplicated(ped[,1])==TRUE)){
        stop("ERROR: There is duplicated individuals, please check your data")
      }
      
      ped[is.na(ped)] <- 0
      return(ped)
    }
    
    ped <- breedR.createPed(data[,c("idNum","nFamily")])
    
  }
  
  # statistical modelling ---------------------------------------------------
  
  #creating all possible object formulas
  base_terms <- c()
  if(GxE){
    
    # fixed effects
    if(is.null(fixed)){
      fixef <- as.formula(paste0("get(varResp) ~ Env"))
    }
    if("Rep"%in%fixed){
      fixef <- as.formula(paste0("get(varResp) ~ Env + EnvRep"))
    }
    if(length(fixed)>1){
      fixef <- as.formula(paste0("get(varResp) ~ Env + EnvRep", 
                                 paste("+",fixed[which(fixed!="Rep")], collapse=" + ")))
    }
    # random effects
    if(!is.null(random)&"Parc"%in%random){
      random <- c(random[which(random=="Parc")],random[which(random!="Parc")])
    }
    # Adicionar EnvRep se "Rep" estiver em random
    if(!is.null(random) & "Rep" %in% random) {
      base_terms <- c(base_terms, "EnvRep")
      
      random[which(random=="Rep")] <- "EnvRep"
    }
    # Adicionar termos espec?ficos de matGen
    if(matGen == "family") {
      base_terms <- c(base_terms, "EnvFam")
    }else if(matGen == "clone") {
      base_terms <- c(base_terms, "clone", "EnvClone")
    }
    # Adicionar os termos de random, se existirem
    if (!is.null(random)) {
      base_terms <- c(base_terms, random)
      if("Rep"%in%base_terms){
        base_terms <- base_terms[-which(base_terms=="Rep")]
      }
    }
    
    # Criar a f?rmula final
    ranef <- as.formula(paste("~", paste(unique(base_terms), collapse = " + ")))
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
    if(matGen=="family"&!is.null(random)){
      ranef <- as.formula(paste0("~", paste(random, collapse=" + ")))
    }
    if(matGen=="clone"&!is.null(random)){
      ranef <- as.formula(paste0("~", matGen, "+", paste(random, collapse=" + ")))
    }
    if(matGen=="clone"&is.null(random)){
      ranef <- as.formula(paste0("~", matGen))
    }
    if(matGen=="family"&is.null(random)){
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
      if(matGen=="family"){
        tryCatch(mAdd <- breedR::remlf90(fixef,
                                         random = ranef,
                                         genetic = list(model = "add_animal",
                                                        pedigree = ped,id = "idNum"),
                                         data = data,method = method,
                                         progsf90.options = 'EM-REML 10'), 
                 error = function(e) {convergence <<- TRUE})
      }else{
        tryCatch(mClone <- breedR::remlf90(fixef,
                                           random = ranef,
                                           data = data,method = method,
                                           progsf90.options = 'EM-REML 10'), 
                 error = function(e) {convergence <<- TRUE})
      })
    
    #transform ranef to lme4 formula
    lmer_ranef <- ranef %>% as.character() %>% {strsplit(.[2], " \\+ ")[[1]]} %>% 
      paste0("(1|", ., ")") %>% paste(., collapse = " + ")
    
    #lme4 formula
    lmerModel <- as.formula(paste0(format(fixef)," + ", lmer_ranef))
    
    suppressMessages(mSig <- lmerTest::lmer(lmerModel, data=data, control = control))
    
    if(!exists("mAdd")&!exists("mClone")){
      cat("AI-REML algorithm got an error, switching to EM-REML\n")
      method = "em"
    }else{
      convergence <- ifelse(matGen=="family",any(is.na(mAdd$var[,1])),any(is.na(mClone$var[,1])))
      if(convergence){
        cat("AI-REML algorithm did not converge, switching to EM-REML\n")
        method = "em"
      }
    }
  }
  options(warn = defaultW)
  
  if(matGen=="family"){
    
    # reliability_and_individual_blup -----------------------------------------
    
    r2Gen <- mAdd$ranef$genetic[[1]] %>% 
      mutate(r2=1-(s.e.)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
      filter(row_number() %in% nFamily$nFamily) %>% 
      add_column(family = nFamily$family, .before = "value") %>% 
      mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
      rename(a = value)
    
    r2Ind <- mAdd$ranef$genetic[[1]] %>% 
      mutate(r2=1-s.e.^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
      filter(row_number() > max(nFamily$nFamily)) %>% 
      add_column(Ind = data$Ind, .before = "value") %>% filter(!is.na(resp)) %>% rename(a = value)
    
    if("GD"%in%names(data)){
      r2GD <- mAdd$ranef$genetic[[1]] %>% 
        mutate(r2=1-(s.e.)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
        filter(row_number() %in% nGD$nGD) %>% 
        add_column(GD = nGD$GD, .before = "value") %>% mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
        rename(a = value)
    }
    if("sire"%in%names(data)){
      r2Sire <- mAdd$ranef$genetic[[1]] %>% 
        mutate(r2=1-(s.e.)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["genetic",1])) %>% 
        filter(row_number() %in% nSire$nSire) %>% 
        add_column(Sire = nSire$sire, .before = "value") %>% mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
        rename(a = value)
      
      if(GxE){
        r2Ind_df <- r2Ind %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Env","Block","Sire","Family","Plot","Tree"),  
                                        sep="-x-") else separate(.,Ind,c("Seq","Env","Block","Sire","Family","Tree"), 
                                                                 sep="-x-")}
      }else{
        r2Ind_df <- r2Ind %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Block","Sire","Family","Plot","Tree"),  
                                        sep="-x-") else separate(.,Ind,c("Seq","Block","Sire","Family","Tree"), 
                                                                 sep="-x-")}
      }
      
    }else{
      
      if(GxE){
        r2Ind_df <- r2Ind %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Env","Block","Family","Plot","Tree"),  
                                        sep="-x-") else separate(.,Ind,c("Seq","Env","Block","Family","Tree"), 
                                                                 sep="-x-")}
        
      }else{
        r2Ind_df <- r2Ind %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Block","Family","Plot","Tree"),  
                                        sep="-x-") else separate(.,Ind,c("Seq","Block","Family","Tree"), 
                                                                 sep="-x-")}
      }}  
    
    
    accInd <- mean(sqrt(1-((r2Ind$s.e.)^2)/mAdd$var["genetic",1]), na.rm=T)
    accFam <- mean(sqrt(1-(((r2Gen$s.e./2))^2)/(mAdd$var["genetic",1]/4)), na.rm=T)
    
    
    # genetic_parameters ------------------------------------------------------
    
    vA <- mAdd$var["genetic",1]
    vE <- mAdd$var["Residual",1]
    Mean <- mean(data$resp,na.rm=T)
    nRep <- if("Rep"%in%fixed) length(unique(data$Rep))
    CVgi <- sqrt(vA)/Mean*100
    
    if(GxE){
      vGxE <- mAdd$var["EnvFam",1]
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
          
          if("Rep"%in%fixed|"Rep"%in%random){
            h2m <- (0.25*vA) / (0.25*vA+(vGxE/nEnv)+(vParc/nRep*nEnv)+vE/(nRep*nArv))
          }else{
            h2m <- (0.25*vA) / (0.25*vA+(vGxE/nEnv)+(vParc/nEnv)+vE)
          }
          
          h2m <- (0.25*vA) / (0.25*vA+(vGxE/nEnv)+(vParc/nRep)+vE/(nRep*nArv))
          genParNames <- c("vA","vGxE","vParc","vE","vPhen","h2a","h2m","h2d",
                           "c2Parc","c2GxE","rgloc","accFam","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vParc + vE
          
          if("Rep"%in%fixed|"Rep"%in%random){
            h2m <- (0.25*vA) / (0.25*vA+(vParc/nRep)+vE/(nRep*nArv))
          }else{
            h2m <- "There is no block replication"
          }
          genParNames <- c("vA","vParc","vE","vPhen","h2a","h2m","h2d",
                           "c2Parc","accFam","accInd","CVgi%","CVe%","Mean")
        }
        c2Parc <- vParc/vPhen
        h2a <- vA/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          if(GxE) SEvGxE <- mAdd$var["EnvFam",2]
          SEvParc <- mAdd$var[random,2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vParc,vGxE,vE,vPhen,h2a,h2m,h2d,c2Parc,c2GxE,rgloc,accFam,accInd,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvParc,SEvGxE,SEvE,NA,h2aSE,matrix(NA,nrow=10,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2a,h2m,h2d,c2Parc,accFam,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvParc,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vParc,vGxE,vE,vPhen,h2a,h2m,h2d,c2Parc,c2GxE,rgloc,accFam,accInd,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2a,h2m,h2d,c2Parc,accFam,accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==2){
        vRdm1 <- mAdd$var[random[2],1]
        if(GxE){
          vPhen <- vA + vParc + vRdm1 + vGxE + vE
          c2GxE <- vGxE/vPhen
          genParNames <- c("vA","vParc",paste0("v",random[2]),"vGxE","vE","vPhen","h2a","h2d",
                           "c2Parc",paste0("c2",random[2]),"c2GxE","rgloc","accFam","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vParc + vRdm1 + vE
          genParNames <- c("vA","vParc",paste0("v",random[2]),"vE","vPhen","h2a","h2d",
                           "c2Parc",paste0("c2",random[2]),"accFam","accInd","CVgi%","CVe%","Mean")
        }
        h2a <- vA/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var["Parc",2]
          SEvRdm1 <- mAdd$var[2,2]
          if(GxE) SEvGxE <- mAdd$var["EnvFam",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vGxE,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2GxE,rgloc,accFam,accInd,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvParc,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,accFam,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvParc,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vGxE,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,c2GxE,rgloc,accFam,accInd,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,accFam,accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==3){
        vRdm1 <- mAdd$var[random[2],1]
        vRdm2 <- mAdd$var[random[3],1]
        if(GxE){
          vPhen <- vA + vParc + vRdm1 + vRdm2 + vGxE + vE
          genParNames <- c("vA","vParc",paste0("v",random[2]),paste0("v",random[3]),"vGxE","vE","vPhen","h2a","h2d",
                           "c2Parc",paste0("c2",random[2]),paste0("c2",random[3]),"c2GxE","rgloc","accFam","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vParc + vRdm1 + vRdm2  + vE
          genParNames <- c("vA","vParc",paste0("v",random[2]),paste0("v",random[3]),"vE","vPhen","h2a","h2d",
                           "c2Parc",paste0("c2",random[2]),paste0("c2",random[3]),"accFam","accInd","CVgi%","CVe%","Mean")
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
          if(GxE) SEvGxE <- mAdd$var["EnvFam",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vGxE,vE,vPhen,h2a,h2d,c2Parc,
                                               c2Rdm1,c2Rdm2,c2GxE,rgloc,accFam,accInd,CVgi,CVe,Mean), 
                                   SE=c(SEvA,SEvParc,SEvRdm1,SEvRdm2,SEvE,NA,h2aSE,matrix(NA,nrow=11,ncol=1)), 
                                   row.names = genParNames),genPar_digits) else 
                                     genPar <- round(data.frame(
                                       Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2a,h2d,c2Parc,
                                                   c2Rdm1,c2Rdm2,accFam,accInd,CVgi,CVe,Mean), 
                                                                SE=c(SEvA,SEvParc,SEvRdm1,SEvRdm2,
                                                                     SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)), 
                                                                row.names = genParNames),genPar_digits)
        }else{
          if(GxE) round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vGxE,vE,vPhen,h2a,h2d,c2Parc,
                                               c2Rdm1,c2Rdm2,c2GxE,rgloc,accFam,accInd,CVgi,CVe,Mean),
                                   row.names = genParNames),genPar_digits) else 
                                     genPar <- round(data.frame(
                                       Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2a,h2d,c2Parc,c2Rdm1,
                                                   c2Rdm2,accFam,accInd,CVgi,CVe,Mean),
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
          if("Rep"%in%fixed|"Rep"%in%random){
            h2m <- (0.25*vA) / (0.25*vA+(vGxE/nEnv)+vE/(nRep))
          }else{
            h2m <- (0.25*vA) / (0.25*vA+(vGxE/nEnv)+vE)
          }
          genParNames <- c("vA","vGxE","vE","vPhen","h2a","h2d","h2m","accFam","accInd","c2GxE",
                           "rgloc","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vE
          genParNames <- c("vA","vE","vPhen","h2a","h2d","h2m","accFam","accInd","CVgi%","CVe%","Mean")
        }
        h2a <- vA/vPhen
        h2d <- (0.75*vA) / (0.75*vA+vE)
        
        if("Rep"%in%fixed){
          h2m <- (0.25*vA) / (0.25*vA+vE/(nRep))
        }else{
          h2m <- "There is no block replicates"
        }
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          if(GxE) SEvGxE <- mAdd$var["EnvFam",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vGxE,vE,vPhen,h2a,h2d,h2m,
                                                         accFam,accInd,c2GxE,rgloc,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvGxE,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vA,vE,vPhen,h2a,h2d,h2m,accFam,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvE,NA,h2aSE,matrix(NA,nrow=7,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vGxE,vE,vPhen,h2a,h2d,h2m,accFam,accInd,c2GxE,
                                                         rgloc,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vA,vE,vPhen,h2a,h2d,h2m,accFam,accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
      }
      if(length(random)==1){
        vRdm1 <- mAdd$var[random,1]
        
        if(GxE){
          vPhen <- vA + vGxE + vRdm1 + vE
          c2GxE <- vGxE/vPhen
          genParNames <- c("vA",paste0("v",random),"vGxE","vE","vPhen","h2a",paste0("c2",random),
                           "c2GxE","rgloc","accFam","accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vRdm1 + vE
          genParNames <- c("vA",paste0("v",random),"vE","vPhen","h2a","c2Rdm1","accFam","accInd","CVgi%","CVe%","Mean")
        }
        
        h2a <- vA/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvRdm1 <- mAdd$var[random,2]
          if(GxE) SEvGxE <- mAdd$var["EnvFam",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vRdm1,vGxE,vE,vPhen,h2a,c2Rdm1,c2GxE,
                                                         rgloc,accFam,accInd,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvRdm1,SEvGxE,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,accFam,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvRdm1,SEvE,NA,
                                                                               h2aSE,matrix(NA,nrow=6,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vRdm1,vGxE,vE,vPhen,h2a,c2Rdm1,c2GxE,
                                                         rgloc,accFam,accInd,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,accFam,accInd,CVgi,CVe,Mean),
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
                           paste0("c2",random[1]),paste0("c2",random[2]),"c2GxE","rgloc","accFam",
                           "accInd","CVgi%","CVe%","Mean")
        }else{
          vPhen <- vA + vRdm1 + vRdm2 + vE
          genParNames <- c("vA",paste0("v",random[1]),paste0("v",random[2]),"vE","vPhen","h2a",
                           "c2Rdm1","c2Rdm2","accFam","accInd","CVgi%","CVe%","Mean")
        }
        
        h2a <- vA/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          if(GxE) SEvGxE <- mAdd$var["EnvFam",2]
          SEvRdm1 <- mAdd$var[random[1],2]
          SEvRdm2 <- mAdd$var[random[2],2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vRdm1,vRdm2,vGxE,vE,vPhen,h2a,c2Rdm1,c2Rdm2,c2GxE,rgloc,
                                                         accFam,accInd,CVgi,CVe,Mean), 
                                             SE=c(SEvA,SEvRdm1,SEvRdm2,SEvGxE,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vA,vRdm1,vE,vPhen,h2a,c2Rdm1,c2Rdm1,accFam,accInd,CVgi,CVe,Mean), 
                                                                          SE=c(SEvA,SEvRdm1,SEvE,NA,h2aSE,matrix(NA,nrow=7,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
          
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vA,vRdm1,vRdm2,vGxE,vE,vPhen,h2a,c2Rdm1,c2Rdm2,c2GxE,rgloc,
                                                         accFam,accInd,CVgi,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vA,vRdm1,vE,vPhen,h2aSE$Estimate,c2Rdm1,c2Rdm1,accFam,
                                                             accInd,CVgi,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        }
    }
    genPar[is.na(genPar)] <- " "
    
    # BLUP_dataframes ---------------------------------------------------------
    
    # Family BLUP
    progBLUP <- r2Gen %>% dplyr::select(family,a) %>% rename(Family=family) %>% arrange(desc(a))
    
    # Individual_BLUP_genetic_gain_and_effectve_population_size --------------- 
    
    if("sire"%in%names(data)){
      
      if(GxE){
        indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
          {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else 
            dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Env","Block","Sire","Family","Plot","Tree"),  sep="-x-") else 
            separate(.,Ind,c("Seq","Env","Block","Sire","Family","Tree"), sep="-x-")} %>% 
          mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a)) %>% dplyr::select(-Seq)
      }else{
        indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
          {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else 
            dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Block","Sire","Family","Plot","Tree"),  sep="-x-") else 
            separate(.,Ind,c("Seq","Block","Sire","Family","Tree"), sep="-x-")} %>% 
          mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a)) %>% dplyr::select(-Seq)
      }
    }else{
      
      if(GxE){
        indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
          {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else 
            dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Env","Block","Family","Plot","Tree"),  sep="-x-") else 
            separate(.,Ind,c("Seq","Env","Block","Family","Tree"), sep="-x-")} %>% 
          mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a)) %>% dplyr::select(-Seq)
      }else{
        indBLUP <- data %>% left_join(.,r2Ind[,1:2],by="Ind") %>% 
          {if (any(random=="Proc")||any(fixed=="Proc")) dplyr::select(.,Proc,Ind,resp,a,Cod) else 
            dplyr::select(.,Ind,resp,a,Cod)} %>% rename(f=resp) %>% 
          {if (plotType=="LP") separate(.,Ind,c("Seq","Block","Family","Plot","Tree"),  sep="-x-") else 
            separate(.,Ind,c("Seq","Block","Family","Tree"), sep="-x-")} %>% 
          mutate("u+a" = Mean+a) %>% drop_na() %>% arrange(desc(a)) %>% dplyr::select(-Seq)
      }
    }
    
    # Provenance BLUP
    if(any(random=="Proc")){
      procBLUP <- mAdd$ranef$Proc[[1]] %>% rownames_to_column(.,"Proc") %>% rename(BLUP=value) %>% arrange(desc(BLUP))
      indBLUP <- indBLUP %>% left_join(.,procBLUP[,1:2],by="Proc") %>% mutate(a=a+BLUP) %>% 
        mutate("u+a"=Mean+a) %>% dplyr::select(.,-last_col()) %>% 
        arrange(desc(a))
      progBLUP <- progBLUP %>% left_join(.,data[,c("family","Proc")], by= c("Family"="family")) %>% 
        unique() %>% left_join(.,procBLUP[,1:2], by="Proc") %>% 
        mutate(a=a+BLUP) %>% dplyr::select(-c(Proc,BLUP)) %>% arrange(desc(a))
      if(PxE){
        procgeBLUP <- mAdd$ranef$EnvProc %>% as.data.frame() %>% rownames_to_column() %>% 
          separate(.,rowname,c("Env","Proc"),sep = "_x_") %>% 
          left_join(.,procBLUP[,1:2],by="Proc") %>% mutate("p+pe"=value+BLUP) %>% 
          dplyr::select(-c(s.e.,BLUP,value)) %>% arrange(desc(.[,3])) %>% 
          split(.,f=~Env)
      }
    }
    
    if(GxE){
      
      # Family BLUP (g) plus ge
      if(GxE) geBLUP <- mAdd$ranef$EnvFam %>% as.data.frame() %>% rownames_to_column() %>% 
          separate(.,rowname,c("Env","Family"),sep = "_x_") %>% 
          left_join(.,progBLUP,by="Family") %>% left_join(.,envMeans[,1:2], by="Env") %>% 
          mutate(a=a/2, Mean=as.numeric(Mean),"g+ge"=a+value,"g+ge+u"=a+value+Mean) %>% 
          dplyr::select(-c(s.e.,a,value,Mean)) %>% arrange(desc(.[,3])) %>% split(.,f=~Env)
      
      # BLUP indexes ------------------------------------------------------------
      
      # MHVG
      n <- table(plyr::ldply(geBLUP, data.frame)[,-1]$Family) %>% as.data.frame
      dataMHVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% group_by(Family) %>% 
        summarise(MHVG=sum(1/g.ge.u)) %>% left_join(.,n,by=c("Family"="Var1")) %>% 
        mutate(MHVG=Freq/MHVG) %>% dplyr::select(-Freq)
      
      # PRVG
      dataPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% 
        mutate(Mean=as.numeric(Mean),prvgIndex=.[,4]/Mean) %>% 
        dplyr::select(-5) %>% group_by(Family) %>% summarise(PRVG=mean(prvgIndex)) %>% 
        mutate(PRVG=PRVG*Mean)
      
      # MHPRVG
      dataMHPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% 
        mutate(Mean=as.numeric(Mean),prvgIndex=1/(.[,4]/Mean)) %>% 
        group_by(Family) %>% summarise(MHPRVG=sum(prvgIndex)) %>% left_join(.,n,by=c("Family"="Var1")) %>% 
        mutate(MHPRVG=Freq/MHPRVG*Mean) %>% 
        dplyr::select(-Freq)
      
      indexBLUP <- dataMHVG %>% mutate(r_MHVG=order(order(MHVG, decreasing = T))) %>% 
        left_join(.,dataPRVG,by="Family") %>% 
        mutate(r_PRVG=order(order(PRVG, decreasing = T))) %>% left_join(.,dataMHPRVG,by="Family") %>% 
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
      varkf <- as.vector(table(indBLUP$Family[1:i]))
      dfNe[i,] <- var(varkf)
    }
    dfNe[is.na(dfNe)] <- 0
    
    indBLUP <- indBLUP %>% mutate(Np = cumsum(!duplicated(indBLUP$Family))) %>% 
      mutate(seq = seq(1:nrow(indBLUP))) %>% mutate(Kf = seq/Np) %>% 
      mutate(varkf = dfNe) %>% mutate(Ne = (4*Np*Kf)/(3+Kf+(varkf/Kf))) %>% 
      dplyr::select(.,-c("Np","seq","Kf","varkf")) %>% 
      relocate(Cod, .after = last_col())
    
    # Overlapping generations
    if(GxE){
      overBLUP <- progBLUP %>% {if (plotType=="LP") mutate(.,Env = 0, Block = 0, Plot = 0, Tree = 0, Cod = "Parent") %>% 
          relocate(Env,Block,Family,Plot,Tree,a,Cod) else mutate(.,Env = 0, Block = 0, Tree = 0, Cod = "Parent") %>% 
          relocate(Env,Block,Family,Tree,a,Cod)} %>% 
        rbind(.,indBLUP[,colnames(.)]) %>% arrange(desc(a))
    }else{
      overBLUP <- progBLUP %>% {if (plotType=="LP") mutate(.,Block = 0, Plot = 0, Tree = 0, Cod = "Parent") %>% 
          relocate(Block,Family,Plot,Tree,a,Cod) else mutate(.,Block = 0, Tree = 0, Cod = "Parent") %>% 
          relocate(Block,Family,Tree,a,Cod)} %>% 
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
        rankBLUP <- rankBLUP %>% mutate(Np_csum = ave(Family==Family, Family, FUN=cumsum)) %>% 
          filter(Np_csum <= maxIndProgeny) %>% 
          mutate(RP = paste0(Block,Family)) %>% mutate(RP_csum = ave(RP==RP,RP, FUN=cumsum)) %>% 
          filter(RP_csum <= maxProgenyBlock) %>% dplyr::select(.,-c("Np_csum","RP","RP_csum"))
        
        # recalculating effective population size and genetic gain
        
        dfNe <- matrix(nrow=nrow(rankBLUP), ncol=1)
        
        for(i in 1:nrow(rankBLUP)){
          varkf <- as.vector(table(rankBLUP$Family[1:i]))
          dfNe[i,] <- var(varkf)
        }
        dfNe[is.na(dfNe)] <- 0
        
        optimizedBLUP <- rankBLUP %>% mutate(Np = cumsum(!duplicated(rankBLUP$Family))) %>% 
          mutate(seq = seq(1:nrow(rankBLUP))) %>% mutate(Kf = seq/Np) %>% 
          mutate(varkf = dfNe) %>% mutate(Ne = (4*Np*Kf)/(3+Kf+(varkf/Kf))) %>% 
          dplyr::select(.,-c("Np","seq","Kf","varkf")) %>% 
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
    genParBLUP$genPar <- genPar; genParBLUP$blupAccuracy$progAccuracy <- r2Gen; 
    genParBLUP$blupAccuracy$indAccuracy <- r2Ind_df 
    genParBLUP$BLUP$progBLUP <- progBLUP; if(GxE) genParBLUP$BLUP$geBLUP <- geBLUP; if(GxE) 
      genParBLUP$BLUP$blupIndex <- indexBLUP; 
    if("Proc"%in%random) genParBLUP$BLUP$procBLUP <- procBLUP; if(PxE) 
      genParBLUP$BLUP$procgeBLUP <- procgeBLUP;
    genParBLUP$BLUP$indBLUP <- indBLUP; genParBLUP$BLUP$overBLUP <- overBLUP
    
    if(optimizeSelection==TRUE){
      genParBLUP$BLUP$optimizedBLUP <- optimizedBLUP
    }
  }
  if(matGen=="clone"){
    
    # reliability_and_individual_blup -----------------------------------------
    
    r2Gen <- mClone$ranef$clone[[1]] %>% 
      mutate(r2=1-(s.e./2)^2/(diag(diag(length(mClone$ranef$clone[[1]][,1])))*as.data.frame(mClone$var)["clone",1])) %>% 
      rownames_to_column("clone") %>% dplyr::rename(g = value)
    
    accClone <- mean(sqrt(1-((r2Gen$s.e.)^2)/(mClone$var["clone",1])), na.rm=T)
    
    # Genetic Parameters
    vG <- mClone$var["clone",1]
    vE <- mClone$var["Residual",1]
    Mean <- mean(data$resp,na.rm=T)
    if(GxE){
      nRep <- length(unique(data$EnvRep))
    }else{
      nRep <- length(unique(data$Rep))
    }
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
          c2Parc <- vParc/vPhen
          h2G <- vG/vPhen
          c2GxE <- vGxE/vPhen
          if("Rep"%in%fixed){
            h2mc <- (vG) / (vG+(vGxE/nEnv)+(vParc/nRep)+vE/(nRep*nArv))
          }else{
            h2mc <- (vG) / (vG+(vGxE/nEnv)+vE)
          }
          genParNames <- c("vG","vGxE","vParc","vE","vPhen","h2G","h2mc","c2GxE","c2Parc","rgloc","accClone","CVg%","CVe%","Mean")
          
        }else{
          vPhen <- vG + vParc + vE
          c2Parc <- vParc/vPhen
          h2G <- vG/vPhen
          if("Rep"%in%fixed){
            h2mc <- (vG) / (vG+(vParc/nRep)+vE/(nRep*nArv))
          }else{
            h2mc <- "There is no block replicates"
          }
          genParNames <- c("vG","vParc","vE","vPhen","h2G","h2mc","c2Parc","accClone","CVg%","CVe%","Mean")
        }

        if(method=="ai"){
          SEvG <- mClone$var["clone",2]
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
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vParc,vE,vPhen,h2G,h2mc,c2Parc,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvParc,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vE,vPhen,h2G,h2mc,c2GxE,c2Parc,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vParc,vE,vPhen,h2G,h2mc,c2Parc,accClone,CVg,CVe,Mean),
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
          SEvG <- mClone$var["clone",2]
          SEvParc <- mClone$var["Parc",2]
          SEvRdm1 <- mClone$var[random[rdmIndex[1]],2]
          if(GxE){
            SEvGxE <- mClone$var["EnvFam",2]
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4+x5),
                                      c(vG,vGxE,vParc,vRdm1,vE),
                                      mClone$reml$invAI) 
          }else{
            h2GSE <- msm::deltamethod(~ x1/(x1+x2+x3+x4),
                                      c(vG,vParc,vRdm1,vE),
                                      mClone$reml$invAI) 
          } 
          SEvE <- mClone$var["Residual",2]
          
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vRdm1,vE,vPhen,h2G,c2GxE,
                                                         c2Parc,c2Rdm1,rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvParc,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=8,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2G,c2Parc,
                                                             c2Rdm1,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvParc,SEvRdm1,
                                                                               SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vRdm1,vE,vPhen,h2G,c2GxE,
                                                         c2Parc,c2Rdm1,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2G,c2Parc,
                                                             c2Rdm1,accClone,CVg,CVe,Mean),
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
          genParNames <- c("vG","vGxE","vParc",paste0("v",random[rdmIndex[1]]),
                           paste0("v",random[rdmIndex[2]]),"vE","vPhen","h2G","c2GxE",
                           "c2Parc",paste0("c2",random[rdmIndex[1]]),
                           paste0("c2",random[rdmIndex[2]]),"rgloc","accClone","CVg%","CVe%","Mean")
        }else{
          vPhen <- vG + vParc + vRdm1 + vRdm2 + vE
          genParNames <- c("vG","vParc",paste0("v",random[rdmIndex[1]]),
                           paste0("v",random[rdmIndex[2]]),"vE","vPhen","h2G",
                           "c2Parc",paste0("c2",random[rdmIndex[1]]),
                           paste0("c2",random[rdmIndex[2]]),"rgloc","accClone","CVg%","CVe%","Mean")
        }
        
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        h2G <- vG/vPhen
        
        if(method=="ai"){
          SEvG <- mClone$var["clone",2]
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
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vRdm1,vRdm2,vE,vPhen,h2G,
                                                         c2GxE,c2Parc,c2Rdm1,c2Rdm2,rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvParc,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE,matrix(NA,nrow=9,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2Parc,
                                                             c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvParc,SEvRdm1,SEvE,
                                                                               NA,h2GSE,matrix(NA,nrow=7,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vParc,vRdm1,vRdm2,vE,vPhen,h2G,
                                                         c2GxE,c2Parc,c2Rdm1,c2Rdm2,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2G,
                                                             c2Parc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
      }}
    
    if(plotType=="STP"){
      CVe = sqrt(3*vG+vE)/Mean*100
      
      if(length(random)==0){
        if(GxE){
          vPhen <- vG + vGxE + vE
          c2GxE <- vGxE/vPhen
          if("Rep"%in%fixed){
            h2mc <- (vG) / (vG+(vGxE/nEnv)+(vE/nRep))
          }else{
            h2mc <- (vG) / (vG+(vGxE/nEnv)+vE)
          }
          
          genParNames <- c("vG","vGxE","vE","vPhen","h2G","h2mc","c2GxE","rgloc","accClone","CVg%","CVe%","Mean")
        }else{
          vPhen <- vG + vE
          h2mc <- (vG) / (vG+(vE/nRep))
          genParNames <- c("vG","vE","vPhen","h2G","h2mc","accClone","CVg%","CVe%","Mean")
        }
        h2G <- vG/vPhen
        
        if(method=="ai"){
          SEvG <- mClone$var["clone",2]
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
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vE,vPhen,h2G,h2mc,c2GxE,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==1){
        vRdm1 <- mClone$var[random,1]
        if(GxE){
          vPhen <- vG + vGxE + vRdm1 + vE
          c2GxE <- vGxE/vPhen
          h2mc <- (vG) / (vG+(vGxE/nEnv)+(vE/nRep))
          genParNames <- c("vG","vGxE",paste0("v",random),"vE","vPhen","h2G","h2mc","c2GxE",paste0("c2",random),
                           "rgloc","accClone","CVg%","CVe%","Mean")
        }else{
          vPhen <- vG + vRdm1 + vE
          h2mc <- (vG) / (vG+(vE/nRep))
          genParNames <- c("vG",paste0("v",random),"vE","vPhen","h2G","h2mc",paste0("c2",random),"accClone","CVg%","CVe%","Mean")
        }
        
        h2G <- vG/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        if(method=="ai"){
          SEvG <- mClone$var["clone",2]
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
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vRdm1,vE,vPhen,h2G,h2mc,c2Rdm1,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vRdm1,vE,vPhen,h2G,h2mc,c2GxE,c2Rdm1,rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vRdm1,vE,vPhen,h2G,h2mc,c2Rdm1,accClone,CVg,CVe,Mean),
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
          SEvG <- mClone$var["clone",2]
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
          
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2GxE,c2Rdm1,c2Rdm2,
                                                         rgloc,accClone,CVg,CVe,Mean), 
                                             SE=c(SEvG,SEvGxE,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE,matrix(NA,nrow=9,ncol=1)), 
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean), 
                                                                          SE=c(SEvG,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE,
                                                                               matrix(NA,nrow=7,ncol=1)), 
                                                                          row.names = genParNames),genPar_digits)
        }else{
          if(GxE) genPar <- round(data.frame(Estimates=c(vG,vGxE,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2GxE,c2Rdm1,c2Rdm2,
                                                         rgloc,accClone,CVg,CVe,Mean),
                                             row.names = genParNames),genPar_digits) else 
                                               genPar <- round(data.frame(
                                                 Estimates=c(vG,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,c2Rdm1,c2Rdm2,accClone,CVg,CVe,Mean),
                                                                          row.names = genParNames),genPar_digits)
        }
      }
    }
    genPar[is.na(genPar)] <- " "
    
    # BLUP_dataframes ---------------------------------------------------------
    
    # Clone BLUP
    cloneBLUP <- r2Gen %>% dplyr::select(clone,g) %>% arrange(desc(g))
    
    # Provenance BLUP
    if("Proc"%in%random){
      procBLUP <- mClone$ranef$Proc %>% as.data.frame() %>% rownames_to_column(.,"Proc") %>% 
        rename(BLUP=value) %>% arrange(desc(BLUP))
      cloneBLUP <- cloneBLUP %>% left_join(.,data[,c("clone","Proc")], by= "clone") %>% 
        unique() %>% left_join(.,procBLUP[,1:2], by="Proc") %>% 
        mutate(g=g+BLUP) %>% dplyr::select(-c(Proc,BLUP)) %>% arrange(desc(g))
      if(PxE) procgeBLUP <- mClone$ranef$EnvProc %>% as.data.frame() %>% rownames_to_column() %>% 
        separate(.,rowname,c("Env","Proc"),sep = "_x_") %>% 
        left_join(.,procBLUP[,1:2],by="Proc") %>% mutate("p+pe"=value+BLUP) %>% 
        dplyr::select(-c(s.e.,BLUP,value)) %>% arrange(desc(.[,3])) %>% 
        split(.,f=~Env)
    }
    
    if(GxE){
      # Clone BLUP (g) plus ge
      geBLUP <- mClone$ranef$EnvClone[[1]] %>% rownames_to_column() %>% 
        separate(.,rowname,c("Env","clone"),sep = "_x_") %>% 
        left_join(.,cloneBLUP,by="clone") %>% left_join(.,envMeans[,1:2], by="Env") %>% 
        mutate(Mean=as.numeric(Mean),"g+ge"=g+value,"g+ge+u"=g+value+Mean) %>% 
        dplyr::select(-c(s.e.,g,value,Mean)) %>% arrange(desc(.[,3])) %>% split(.,f=~Env)
      
      # BLUP indexes ------------------------------------------------------------
      
      # MHVG
      n <- table(plyr::ldply(geBLUP, data.frame)[,-1]$clone) %>% as.data.frame
      dataMHVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% group_by(clone) %>% 
        dplyr::summarise(MHVG=sum(1/g.ge.u)) %>% left_join(.,n,by=c("clone"="Var1")) %>% 
        mutate(MHVG=Freq/MHVG) %>% dplyr::select(-Freq)
      
      # PRVG
      dataPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% 
        mutate(Mean=as.numeric(Mean),prvgIndex=.[,4]/Mean) %>% dplyr::select(-5) %>% group_by(clone) %>% 
        dplyr::summarise(PRVG=mean(prvgIndex)) %>% mutate(PRVG=PRVG*Mean)
      
      # MHPRVG
      dataMHPRVG <- plyr::ldply(geBLUP, data.frame)[,-1] %>% left_join(.,envMeans[,c(1:2)],by="Env") %>% 
        mutate(Mean=as.numeric(Mean),prvgIndex=1/(.[,4]/Mean)) %>% group_by(clone) %>% 
        dplyr::summarise(MHPRVG=sum(prvgIndex)) %>% 
        left_join(.,n,by=c("clone"="Var1")) %>% mutate(MHPRVG=Freq/MHPRVG*Mean) %>% 
        dplyr::select(-Freq)
      
      indexBLUP <- dataMHVG %>% mutate(r_MHVG=order(order(MHVG, decreasing = T))) %>% left_join(.,dataPRVG,by="clone") %>% 
        mutate(r_PRVG=order(order(PRVG, decreasing = T))) %>% left_join(.,dataMHPRVG,by="clone") %>% 
        mutate(r_MHPRVG=order(order(MHPRVG, decreasing = T))) %>% arrange(-MHPRVG) %>% as.data.frame()
    }
    
    # final_output_list -------------------------------------------------------
    genParBLUP <- list()
    
    # exploratory_analysis -------------------------------------------------
    genParBLUP$expAnalysis$respMeans <- respMeans; if(GxE) genParBLUP$expAnalysis$envMeans <- envMeans; 
    if(any(random=="Proc")||any(fixed=="Proc")) genParBLUP$expAnalysis$procMeans <- procMeans; 
    if("Rep"%in%fixed) genParBLUP$expAnalysis$repMeans <- repMeans; genParBLUP$expAnalysis$cloneMeans <- treatMeans
    
    # model_output_and_significance_tests -------------------------------------------------
    genParBLUP$Model$remlMethod <- method; genParBLUP$Model$remlModel <- mClone; genParBLUP$Model$mSig$fixedSig <- anova(mSig)
    genParBLUP$Model$mSig$randomSig <- suppressMessages(lmerTest::ranova(mSig, reduce.terms = F))
    
    # genetic_parameters_accuracy_and_Blup -------------------------------------------------
    genParBLUP$genPar <- genPar; genParBLUP$blupAccuracy$cloneAccuracy <- r2Gen; genParBLUP$BLUP$cloneBLUP <- cloneBLUP; 
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
    
    # plot_outputs ----------------------------------------------------------
    
    if(plotRanking & !is.null(directory)){
      
      # BLUP family or clone plot
      pDf <- r2Gen %>% rename(gen=1, y=2) %>%  arrange(y) %>% 
        mutate(gen = factor(gen, levels = gen), ci=s.e.*1.76,)
      
      pDf$check <- ifelse(pDf$y > 0, "Above", "Below")
      
      scale <- 12-log(length(levels(pDf$gen)),2)
      
      rank_plot <- ggplot(pDf, aes(x = .data$gen, y = .data$y, color = .data$check)) + 
        coord_flip() + 
        xlab(matGen) + 
        ylab(paste0("BLUP - ", varResp)) + theme_minimal() +
        theme(axis.text.y = element_text(size = scale + 1),
              axis.title = element_text(size = scale + 5, face = "bold"),
              legend.position = c(0.95, 0.05),
              legend.justification = c(1, 0), 
              legend.title = element_blank(), 
              legend.background = element_blank() 
        ) +
        geom_hline(yintercept = 0, lwd = I(7/12), colour = I(grDevices::hsv(0/12, 7/12, 7/12)), 
                   alpha = I(5/12)) +
        geom_errorbar(aes(ymin = !!sym("y") - !!sym("ci"), ymax = !!sym("y") + !!sym("ci"), 
                          color = !!sym("check")), width = 0) +
        geom_point(size = 3) +
        scale_color_manual(values = c("Above" = "darkgreen", "Below" = "darkred"))
      
      ggsave(filename = paste0(directory, "_rank_plot_", varResp, ".pdf"),
             plot=rank_plot, width=8, height=10, units="in")
      
      # if GxE ---------------------------------------------------------------
      
      if(GxE){
        
        pDf_ge <- plyr::ldply(geBLUP, data.frame) %>% 
          rename(gen=3) %>% select(Env, gen, `g.ge`)
        
        env_sel <- table(pDf_ge$Env) %>% as.data.frame() %>% 
          filter(Freq == max(Freq)) %>% slice(1) %>% pull(Var1) %>% 
          as.character()
        
        gen_category <- pDf_ge %>% filter(Env==env_sel) %>% 
          arrange(desc(`g.ge`)) %>%
          mutate(rank = row_number(),
                 hjust = 1.1,
                 category_env_sel = case_when(
                   rank <= quantile(rank, 1/3) ~ "Melhores",
                   rank > quantile(rank, 1/3) & rank <= quantile(rank, 2/3) ~ "Medianos",
                   rank > quantile(rank, 1/3) ~ "Piores", 
                   TRUE ~ NA_character_
                 )) %>% dplyr::select(gen, category_env_sel) # edita aqui se for manter assim mesmo
        
        ranked_data <- pDf_ge %>%
          group_by(Env) %>%
          arrange(desc(`g.ge`)) %>%
          mutate(rank = row_number(),
                 hjust = case_when(
                   Env == env_sel ~ 1.1,
                   TRUE ~ -0.1)) %>%
          ungroup() %>% left_join(.,gen_category, by = "gen")
        
        unique_env <- unique(ranked_data$Env) %>% .[.!=env_sel] %>% 
          as.character() %>% sort()
        
        ranked_data$hjust[which(ranked_data$Env==last(unique_env))] <- -0.1
        
        ranked_data$Env <- factor(ranked_data$Env, levels = c(env_sel, unique_env))
        
        rank_env_plot <- ggplot(ranked_data, aes(x = Env, y = -rank, group = gen, color = category_env_sel)) +  # -rank para inverter
          geom_line(linewidth = 1) +      # Linhas conectando os grupos
          geom_point(size = 3) +     # Pontos em cada grupo
          geom_text(aes(label = gen, hjust = hjust), size = 4) + # Ajustar rótulos
          scale_color_manual(
            values = c("Melhores" = "darkgreen", "Medianos" = "#dab600", "Piores" = "darkred")
          ) + 
          theme_minimal() +
          theme(
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(size=12,face="bold")
          )
        
        ggsave(filename = paste0(directory,"_rank_env_plot_", varResp, ".pdf"),
               plot=rank_env_plot, width=8, height=10, units="in")
      }
    }
    
    
    # txt output ------------------------------------------------------------

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
        if(matGen=="family"){
          strMod <- paste0(varResp," = ",strFix," + ","r",matGen," + rEnvFam"," + ",strRnd)
        }else{
          strMod <- paste0(varResp," = ",strFix," + ","r",matGen," + rEnvClone"," + ",strRnd)
        }
      }else{
        strFix <- paste(paste0("f",fixed),collapse=" + ")
        if(matGen=="family"){
          strMod <- paste0(varResp," = ",strFix," + r",matGen," + r","EnvFam")
        }else{
          strMod <- paste0(varResp," = ",strFix," + r",matGen," + r","EnvClone")
        }
      }
    }else{
      if(!is.null(random)){
        strFix <- paste(paste0("f",fixed),collapse="+")
        strRnd <- paste(paste0("r",random),collapse="+")
        strMod <- paste0(varResp,"=",strFix,"+","r",matGen,"+",strRnd)
      }else{
        strFix <- paste(paste0("f",fixed),collapse="+")
        strMod <- paste0(varResp,"=",strFix,"+","r",matGen)
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
    if("Rep"%in%fixed) print(repMeans)
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
    if(matGen=="clone"){
      print(cloneBLUP)}
    if(matGen=="family"){
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
    if(matGen=="family"){
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
    if(matGen=="clone"){
      cat("\n")
      print(r2Gen)}
    if(matGen=="family"){
      print(r2Gen)
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
    
    if(!GxE & matGen=="family" & optimizeSelection==TRUE){
      
      blupoptimized <- paste0(directory,"_",varResp,"_optimizedBLUP",".txt")
      
      if(file.exists(blupoptimized)){
        file.remove(blupoptimized)
      }
      
      sink(blupoptimized, append=TRUE, type = "output")
      
      cat("\n----------------------------------------------------------------------------------------------\n")
      cat("                                           |optimized BLUP| \n"                                     )
      cat("----------------------------------------------------------------------------------------------\n\n")
      cat("--------------------------------------------- Arguments --------------------------------------\n\n")
      cat("Max. Individuals per Family:\n\n")
      print(maxIndProgeny) 
      cat("Max. Families per Block:\n\n")
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