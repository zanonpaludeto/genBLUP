genBLUP <- function(data, varResp, treatment, plotType, fixed = "Rep", random = NULL, method = "ai", 
                    genPar_digits, otimizeSelection = FALSE, maxIndProgeny = NULL, 
                    maxProgenyBlock = NULL, excludeControl = NULL, excludeCod = NULL, directory = NULL){
  
  # loading packages --------------------------------------------------------
  
  pacman::p_load(future,lme4,stringr,tidyverse,ggroups,dplyr)
  
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
    warning("Directory is 'NULL', so no txt outputs were created")
  }
  if(treatment=="Clone"&otimizeSelection==T){
    stop("ERROR: You can't otimize selection in a clonal analysis")
  }  
  if(method=="ai"||is.null(method)){
    cat("AI-REML algorithm was selected\n")
  }else{
    cat("EM-REML algorith was selected\n")
  }
  
  # exploratory analysis ----------------------------------------------------
  
  data <- data[order(data[,which(colnames(data)==treatment)]),]
  
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
  
  # factorizing_variables_and_treeCheck -------------------------------------
  
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
  
  # preparing data
  
  # Creating Cod if it doesn't exist
  if(!("Cod" %in% colnames(data))){
    data$Cod = " "
  }
  
  #creating Ind column
  if(plotType=="LP"){  
    data$Ind <- paste0(data$Rep,sep = ".",data[,names(data) %in% treatment],sep=".",data$Parc,sep=".",data$Arv)
  }else{data$Ind <- paste0(data$Rep,sep = ".",data[,names(data) %in% treatment],sep=".",data$Arv)  
  }
  
  data <- data %>% mutate(idNum=1)
  
  if(treatment=="Prog"){
    
    #counting nProg (dams)
    nProg <- data.frame(Prog = sort(unique(data$Prog)), nProg = as.numeric(seq(1,length(unique(data$Prog)))))
    nInd <- seq(dim(nProg)[1]+1,dim(data)[1]+(dim(nProg)[1]))
    data <- data %>% mutate(idNum=as.numeric(nInd))
    data <- dplyr::left_join(data,nProg,by="Prog")
    
    breedR.createPed <- function(pedData){
      
      ped <- data.frame(pedData,sire=0) %>% setNames(c("idNum","dam","sire")) %>% relocate(dam, .after = sire)
      
      if(any(duplicated(ped[,1])==TRUE)){
        stop("ERROR: There is duplicated individuals, please check your data")
      }
      ped <- data.frame(as.numeric(levels(factor(ped$dam))),0,0) %>% setNames(c("idNum","sire","dam")) %>% rbind(.,ped)
      return(ped)
    }
    
    ped <- breedR.createPed(data[,c("idNum","nProg")])
    
    # dominance matrix --------------------------------------------------------
    
    # if(dominance==TRUE){
    #   
    #   A <- buildA(ped)
    #   D <- ggroups::buildD(ped,A)
    #   if(any(duplicated(rownames(D))==TRUE)){
    #     stop("ERROR: There is duplicated individuals, please check your data")
    #   }}
  }
  
  
  # statistical modelling ---------------------------------------------------
  
  # linear_plot -------------------------------------------------------------
  
  #suppressing messages and warnings
  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
  defaultW <- getOption("warn") 
  options(warn = -1)
  convergence = F
  
  while(!convergence){
    
    suppressMessages(if(plotType=="LP"){
      
      if(length(fixed)==1&length(random)==0){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed) %>% rename("Fixed1"=fixed,"Treat"=3)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1,
                                  random = ~ Parc,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1,
                                    random = ~ Treat + Parc,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        #if(dominance==TRUE){
        #mDom <- sommer::mmer(resp ~ Fixed1,
        #random= ~ vs(idNum,Gu=D) + Parc,
        #rcov=~units,
        #data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc), data=df, control = control)
      }
      if(length(fixed)==1&length(random)==1){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed,random) %>% rename("Fixed1"=fixed,"Treat"=3,"Random1"=random)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1,
                                  random = ~ Parc + Random1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1,
                                    random = ~ Treat + Parc + Random1,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1,
        #                          random= ~ vs(idNum,Gu=D) + Parc + Random1,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc) + (1|Random1), data=df, control = control)
      }
      
      if(length(fixed)==1&length(random)==2){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed,random) %>% rename("Fixed1"=fixed,"Treat"=3,"Random"=random)
        
        if(any(random=="Proc")){
          if(df$Random2==data$Proc){
            df$Random2 <- df$Random1
            df$Random1 <- data$Proc
          }
        }
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1,
                                  random = ~ Parc + Random1 + Random2,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1,
                                    random = ~ Treat + Parc + Random1 + Random2,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1,
        #                          random= ~ vs(idNum,Gu=D) + Parc + Random1 + Random2,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        #   
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Parc) + (1|Random1) + (1|Random2), data=df, control = control)
      }
      
      if(length(fixed)==2&length(random)==0){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed) %>% rename("Fixed"=fixed,"Treat"=3)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                  random = ~ Parc,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                    random = ~ Treat + Parc,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2,
        #                          random= ~ vs(idNum,Gu=D) + Parc,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc), data=df, control = control)
      }
      
      if(length(fixed)==2&length(random)==1){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed,random) %>% rename("Fixed"=fixed,"Treat"=3,"Random1"=random)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                  random = ~ Parc + Random1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                    random = ~ Treat + Parc + Random1,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2,
        #                          random= ~ vs(idNum,Gu=D) + Parc + Random1,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc) + (1|Random1), data=df, control = control)
      }
      
      if(length(fixed)==2&length(random)==2){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed,random) %>% rename("Fixed"=fixed,"Treat"=3,"Random"=random)
        
        if(any(random=="Proc")){
          if(df$Random2==data$Proc){
            df$Random2 <- df$Random1
            df$Random1 <- data$Proc
          }
        }
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                  random = ~ Parc + Random1 + Random2,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                    random = ~ Treat + Parc + Random1 + Random2,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2,
        #                          random= ~ vs(idNum,Gu=D) + Parc + Random1 + Random2,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Parc) + (1|Random1) + (1|Random2), data=df, control = control)
      }
      
      if(length(fixed)==3&length(random)==0){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed) %>% rename("Fixed"=fixed,"Treat"=3)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random = ~ Parc,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                    random = ~ Treat + Parc,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
        #                          random= ~ vs(idNum,Gu=D) + Parc,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc), data=df, control = control)
      }
      
      if(length(fixed)==3&length(random)==1){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed,random) %>% rename("Fixed"=fixed,"Treat"=3,"Random1"=random)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random = ~ Parc + Random1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                    random = ~ Treat + Parc + Random1,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
        #                          random= ~ vs(idNum,Gu=D) + Parc + Random1,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc) + (1|Random1), data=df, control = control)
      }
      
      if(length(fixed)==3&length(random)==2){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"Parc","resp",fixed,random) %>% rename("Fixed"=fixed,"Treat"=3,"Random"=random)
        
        if(any(random=="Proc")){
          if(df$Random2==data$Proc){
            df$Random2 <- df$Random1
            df$Random1 <- data$Proc
          }
        }
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random = ~ Parc + Random1 + Random2,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                    random = ~ Treat + Parc + Random1 + Random2,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
        #                          random= ~ vs(idNum,Gu=D) + Parc + Random1 + Random2,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Parc) + (1|Random1) + (1|Random2), data=df, control = control)
      }
    })
    
    # single_tree_plot --------------------------------------------------------
    
    suppressMessages(if(plotType=="STP"){
      
      if(length(fixed)==1&length(random)==0){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed) %>% rename("Fixed1"=fixed,"Treat"=3)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1,
                                  random = ~ 1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1,
                                    random = ~ Treat,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1,
        #                          random= ~ vs(idNum,Gu=D),
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat), data=df, control = control)
      }
      
      if(length(fixed)==1&length(random)==1){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed,random) %>% rename("Fixed1"=fixed,"Treat"=3,"Random1"=random)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1,
                                  random = ~ Random1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1,
                                    random = ~ Treat + Random1,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1,
        #                          random= ~ vs(idNum,Gu=D) + Random1,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Random1), data=df, control = control)
      }
      
      if(length(fixed)==1&length(random)==2){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed,random) %>% rename("Fixed1"=fixed,"Treat"=3,"Random"=random)
        
        if(any(random=="Proc")){
          if(df$Random2==data$Proc){
            df$Random2 <- df$Random1
            df$Random1 <- data$Proc
          }
        }
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1,
                                  random = ~ Random1 + Random2,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1,
                                    random = ~ Treat + Random1 + Random2,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1,
        #                          random= ~ vs(idNum,Gu=D) + Random1 + Random2,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        #   
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + (1|Treat) + (1|Random1) + (1|Random2), data=df, control = control)
      }
      
      if(length(fixed)==2&length(random)==0){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed) %>% rename("Fixed"=fixed,"Treat"=3)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                  random = ~ 1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                    random = ~ Treat,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2,
        #                          random= ~ vs(idNum,Gu=D),
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat), data=df, control = control)
      }
      
      if(length(fixed)==2&length(random)==1){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed,random) %>% rename("Fixed"=fixed,"Treat"=3,"Random1"=random)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                  random = ~ Random1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                    random = ~ Treat + Random1,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')}      
        
        
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2,
        #                          random= ~ vs(idNum,Gu=D) + Random1,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Random1), data=df, control = control)
      }
      
      if(length(fixed)==2&length(random)==2){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed,random) %>% rename("Fixed"=fixed,"Treat"=3,"Random"=random)
        
        if(any(random=="Proc")){
          if(df$Random2==data$Proc){
            df$Random2 <- df$Random1
            df$Random1 <- data$Proc
          }
        }
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                  random = ~ Random1 + Random2,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2,
                                    random = ~ Treat + Random1 + Random2,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2,
        #                          random= ~ vs(idNum,Gu=D) + Random1 + Random2,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + (1|Treat) + (1|Random1) + (1|Random2), data=df, control = control)
      }
      
      if(length(fixed)==3&length(random)==0){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed) %>% rename("Fixed"=fixed,"Treat"=3)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random = ~ 1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                    random = ~ Treat,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
        #                          random= ~ vs(idNum,Gu=D),
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat), data=df, control = control)
      }
      
      if(length(fixed)==3&length(random)==1){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed,random) %>% rename("Fixed"=fixed,"Treat"=3,"Random1"=random)
        
        if(treatment=="Prog"){
          mAdd <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                  random = ~ Random1,
                                  genetic = list(model = "add_animal",
                                                 pedigree = ped,id = "idNum"),
                                  data = df,method = method,
                                  progsf90.options = 'EM-REML 10')
          
        }else{
          mClone <- breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                    random = ~ Treat + Random1,
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
        #                          random= ~ vs(idNum,Gu=D) + Random1,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Random1), data=df, control = control)
      }
      
      if(length(fixed)==3&length(random)==2){
        
        df <- data %>% dplyr::select("Ind","idNum",treatment,"resp",fixed,random) %>% rename("Fixed"=fixed,"Treat"=3,"Random"=random)
        
        if(any(random=="Proc")){
          if(df$Random2==data$Proc){
            df$Random2 <- df$Random1
            df$Random1 <- data$Proc
          }
        }
        if(treatment=="Prog"){
          mAdd %<-% breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                    random = ~ Random1 + Random2,
                                    genetic = list(model = "add_animal",
                                                   pedigree = ped,id = "idNum"),
                                    data = df,method = method,
                                    progsf90.options = 'EM-REML 10')
          
        }else{
          mClone %<-% breedR::remlf90(fixed = resp ~ Fixed1 + Fixed2 + Fixed3,
                                      random = ~ Treat + Random1 + Random2,
                                      data = df,method = method,
                                      progsf90.options = 'EM-REML 10')
        }
        # if(dominance==TRUE){
        #   mDom <- sommer::mmer(resp ~ Fixed1 + Fixed2 + Fixed3,
        #                          random= ~ vs(idNum,Gu=D) + Random1 + Random2,
        #                          rcov=~units,
        #                          data=df %>% mutate(idNum = as.factor(idNum)), verbose = F)
        # }
        mSig <- lme4::lmer(resp ~ Fixed1 + Fixed2 + Fixed3 + (1|Treat) + (1|Random1) + (1|Random2), data=df, control = control)
      }
    })
    
    convergence <- ifelse(treatment=="Prog",any(!is.na(mAdd$var[,1])),any(!is.na(mClone$var[,1])))
    if(convergence==F){
      cat("AI-REML algorithm did not converge, switching to EM-REML\n")
      method = "em"
    }}
  options(warn = defaultW)
  
  if(treatment=="Prog"){
    
    # reliability_and_individual_blup -----------------------------------------
    
    r2Prog <- mAdd$ranef$genetic[[1]] %>% 
      mutate(r2=1-(s.e./2)^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["Residual",1])) %>% 
      filter(row_number() < nrow(nProg)+1) %>% add_column(Prog = nProg$Prog, .before = "value") %>% mutate(value=value) %>% mutate(s.e.=s.e.) %>% 
      rename(a = value)
    
    r2Ind <- mAdd$ranef$genetic[[1]] %>% 
      mutate(r2=1-s.e.^2/(diag(diag(length(mAdd$ranef$genetic[[1]][,1])))*as.data.frame(mAdd$var)["Residual",1])) %>% 
      filter(row_number() > nrow(nProg)) %>% add_column(Ind = df$Ind, .before = "value") %>% filter(!is.na(resp)) %>% rename(a = value)
    r2Ind_df <- r2Ind %>% 
      {if (plotType=="LP") separate(.,Ind,c("Block","Progeny","Plot","Tree"),  sep="\\.") else separate(.,Ind,c("Block","Progeny","Tree"), sep="\\.")}
    
    accInd <- mean(sqrt(1-((r2Ind$s.e.)^2)/mAdd$var["genetic",1]), na.rm=T)
    accProg <- mean(sqrt(1-(((r2Prog$s.e./2))^2)/(mAdd$var["genetic",1]/4)), na.rm=T)
    
    
    # genetic_parameters ------------------------------------------------------
    
    vA <- mAdd$var["genetic",1]
    vE <- mAdd$var["Residual",1]
    Mean <- mean(df$resp,na.rm=T)
    nRep <- length(unique(data$Rep))
    CVgi <- sqrt(vA)/Mean*100
    
    # Linear plot
    if(plotType=="LP"){
      nArv <- length(unique(data$Arv))
      vParc <- mAdd$var["Parc",1]
      h2d <- (0.75*vA) / (0.75*vA+vE)
      CVe = (sqrt((0.75*vA+vE)/nArv+vParc))/Mean*100
      
      if(length(random)==0){
        vPhen <- vA + vParc + vE
        c2Parc <- vParc/vPhen
        h2a <- vA/vPhen
        h2m <- (0.25*vA) / (0.25*vA+(vParc/nRep)+vE/(nRep*nArv))
        
        genParNames<- c("vA","vParc","vE","vPhen","h2a","h2d","h2m","c2Parc","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var["Parc",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2a,h2d,h2m,c2Parc,accProg,accInd,CVgi,CVe,Mean), 
                                     SE=c(SEvA,SEvParc,SEvE,NA,h2aSE,matrix(NA,nrow=8,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vA,vParc,vE,vPhen,h2a,h2d,h2m,c2Parc,accProg,accInd,CVgi,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==1){
        vRdm1 <- mAdd$var["Random1",1]
        vPhen <- vA + vParc + vRdm1 + vE
        h2a <- vA/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        genParNames <- c("vA","vParc","vRdm1","vE","vPhen","h2a","h2d","c2Parc","c2Rdm1","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var["Parc",2]
          SEvRdm1 <- mAdd$var["Random1",2]
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
      if(length(random)==2){
        vRdm1 <- mAdd$var["Random1",1]
        vRdm2 <- mAdd$var["Random2",1]
        vPhen <- vA + vParc + vRdm1 + vRdm2 + vE
        h2a <- vA/vPhen
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        c2Rdm2 <- vRdm2/vPhen
        
        genParNames <- c("vA","vParc","vRdm1","vRdm2","vE","vPhen","h2a","h2d","c2Parc","c2Rdm1","c2Rdm2","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvParc <- mAdd$var["Parc",2]
          SEvRdm1 <- mAdd$var["Random1",2]
          SEvRdm2 <- mAdd$var["Random2",2]
          SEvE <- mAdd$var["Residual",2]
          h2aSE <- mAdd$funvars["sample sd",1]
          
          genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2aSE,h2dSE,c2Parc,c2Rdm1,c2Rdm2,accProg,accInd,CVgi,CVe,Mean), 
                                     SE=c(SEvA,SEvParc,SEvRdm1,SEvRdm2,SEvE,NA,h2aSE,matrix(NA,nrow=9,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vA,vParc,vRdm1,vRdm2,vE,vPhen,h2aSE,h2dSE,c2Parc,c2Rdm1,c2Rdm2,accProg,accInd,CVgi,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }}
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
        vRdm1 <- mAdd$var["Random1",1]
        vPhen <- vA + vRdm1 + vE
        h2a <- vA/vPhen
        c2Rdm1 <- vRdm1/vPhen
        
        genParNames <- c("vA","vRdm1","vE","vPhen","h2a","c2Rdm1","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvA <- mAdd$var["genetic",2]
          SEvRdm1 <- mAdd$var["Random1",2]
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
        vRdm1 <- mAdd$var["Random1",1]
        vRdm2 <- mAdd$var["Random2",1]
        vPhen <- vA + vRdm1 + vRdm2 + vE
        h2a <- vA/vPhen
        
        genParNames <- c("vA","vRdm1","vRdm2","vE","vPhen","h2a","c2Rdm1","c2Rdm2","accProg","accInd","CVgi%","CVe%","Mean")
        
        if(method=="ai"){
          h2a <- mAdd$funvars["mean",1]
          SEvRdm1 <- mAdd$var["Random1",2]
          SEvRdm2 <- mAdd$var["Random2",2]
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
    
    # Provenance BLUP
    if(any(random=="Proc")){
      procBLUP <- mAdd$ranef$Random1[[1]] %>% rownames_to_column(.,"Proc") %>% rename(BLUP=value) %>% arrange(desc(BLUP))
      indBLUP <- indBLUP %>% left_join(.,procBLUP[,1:2],by="Proc") %>% mutate(a=a+BLUP) %>% dplyr::select(.,-last_col())
    }
    
    # otimize_selection -------------------------------------------------------
    
    if(otimizeSelection==TRUE){
      
      rankBLUP <- indBLUP
      
      if(!is.null(excludeControl)){
        rankBLUP <- indBLUP[!indBLUP$Progeny %in% excludeControl,]
      }
      if(!is.null(excludeCod)){
        rankBLUP <- indBLUP[!indBLUP$Cod %in% excludeCod,]
      }
      
      # Filtering using the maximum number of individuals per progeny, then the max number of progenies addmited in the same block
      rankBLUP <- rankBLUP %>% mutate(Np_csum = ave(Progeny==Progeny, Progeny, FUN=cumsum)) %>% filter(Np_csum <= maxIndProgeny) %>% 
        mutate(RP = paste0(Block,Progeny)) %>% mutate(RP_csum = ave(RP==RP,RP, FUN=cumsum)) %>% 
        filter(RP_csum <= maxProgenyBlock) %>% select(.,-c("Np_csum","RP","RP_csum"))
      
      # recalculating effective population size and genetic gain
      
      dfNe <- matrix(nrow=nrow(rankBLUP), ncol=1)
      
      for(i in 1:nrow(rankBLUP)){
        varkf <- as.vector(table(rankBLUP$Progeny[1:i]))
        dfNe[i,] <- var(varkf)
      }
      dfNe[is.na(dfNe)] <- 0
      
      otimizedBLUP <- rankBLUP %>% mutate(Np = cumsum(!duplicated(rankBLUP$Progeny))) %>% mutate(seq = seq(1:nrow(rankBLUP))) %>% mutate(Kf = seq/Np) %>% 
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
    
    if(otimizeSelection==TRUE){
      genParBLUP$BLUP$otimizedBLUP <- otimizedBLUP
    }
  }
  if(treatment=="Clone"){
    
    # reliability_and_individual_blup -----------------------------------------
    
    r2Clone <- mClone$ranef$Treat[[1]] %>% 
      mutate(r2=1-(s.e./2)^2/(diag(diag(length(mClone$ranef$Treat[[1]][,1])))*as.data.frame(mClone$var)["Residual",1])) %>% 
      rownames_to_column("Clone") %>% rename(g = value)
    
    accClone <- mean(sqrt(1-(((r2Clone$s.e./2))^2)/(mClone$var["Treat",1]/4)), na.rm=T)
    
    # Genetic Parameters
    vG <- mClone$var["Treat",1]
    vE <- mClone$var["Residual",1]
    Mean <- mean(df$resp,na.rm=T)
    nRep <- length(unique(data$Rep))
    CVg <- sqrt(vG)/Mean*100
    
    if(plotType=="LP"){
      nArv <- length(unique(data$Arv))
      vParc <- mClone$var["Parc",1]
      CVe = (sqrt((3*vG+vE)/nArv+vParc))/Mean*100
      
      if(length(random)==0){
        vPhen <- vG + vParc + vE
        c2Parc <- vParc/vPhen
        h2G <- vG/vPhen
        h2mc <- (vG) / (vG+(vParc/nRep)+vE/(nRep*nArv))
        
        genParNames <- c("vG","vParc","vE","vPhen","h2G","h2mc","c2Parc","accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Treat",2]
          SEvParc <- mClone$var["Parc",2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- deltamethod(~ x1/(x1+x2+x3),
                               c(vG,vParc,vE),
                               mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vParc,vE,vPhen,h2GSE,h2mc,c2Parc,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvParc,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vParc,vE,vPhen,h2G,h2mc,c2Parc,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==1){
        vRdm1 <- mClone$var["Random1",1]
        vPhen <- vG + vParc + vRdm1 + vE
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        h2G <- vG/vPhen
        
        genParNames <- c("vG","vParc","vRdm1","vE","vPhen","h2G","c2Parc","accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Treat",2]
          SEvParc <- mClone$var["Parc",2]
          SEvRdm1 <- mClone$var["Random1",2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- deltamethod(~ x1/(x1+x2+x3+x4),
                               c(vG,vParc,vRdm1,vE),
                               mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2G,c2Rdm1,c2Parc,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvParc,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vE,vPhen,h2G,c2Rdm1,c2Parc,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }
      if(length(random)==2){
        vRdm1 <- mClone$var["Random1",1]
        vRdm2 <- mClone$var["Random2",1]
        vPhen <- vG + vParc + vRdm1 + vRdm2 + vE
        c2Parc <- vParc/vPhen
        c2Rdm1 <- vRdm1/vPhen
        h2G <- vG/vPhen
        
        genParNames <- c("vG","vParc","vRdm1","vRdm2","vE","vPhen","h2G","c2Parc","accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Treat",2]
          SEvParc <- mClone$var["Parc",2]
          SEvRdm1 <- mClone$var["Random1",2]
          SEvRdm2 <- mClone$var["Random2",2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- deltamethod(~ x1/(x1+x2+x3+x4+x5),
                               c(vG,vParc,vRdm1,vRdm2,vE),
                               mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2Rdm1,c2Parc,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvParc,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=6,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vParc,vRdm1,vRdm2,vE,vPhen,h2G,c2Rdm1,c2Parc,accClone,CVg,CVe,Mean),
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
          SEvG <- mClone$var["Treat",2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- deltamethod(~ x1/(x1+x2),
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
        vRdm1 <- mClone$var["Random1",1]
        vPhen <- vG + vRdm1 + vE
        h2G <- vG/vPhen
        h2mc <- (vG) / (vG+(vE/nRep))
        
        genParNames <- c("vG","vRdm1","vE","vPhen","h2G","h2mc","accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Treat",2]
          SEvRdm1 <- mClone$var["Random1",2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- deltamethod(~ x1/(x1+x2+x3),
                               c(vG,vRdm1,vE),
                               mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vRdm1,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvRdm1,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vRdm1,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean),
                                     row.names = genParNames),genPar_digits)
        }
        
      }  
      if(length(random)==2){
        vRdm1 <- mClone$var["Random1",1]
        vRdm2 <- mClone$var["Random2",1]
        vPhen <- vG + vRdm1 + vRdm2 + vE
        h2G <- vG/vPhen
        h2mc <- (vG) / (vG+(vE/nRep))
        
        genParNames <- c("vG","vRdm1","vRdm2","vE","vPhen","h2G","h2mc","accClone","CVg%","CVe%","Mean")
        
        if(method=="ai"){
          SEvG <- mClone$var["Treat",2]
          SEvRdm1 <- mClone$var["Random1",2]
          SEvRdm2 <- mClone$var["Random2",2]
          SEvE <- mClone$var["Residual",2]
          h2GSE <- deltamethod(~ x1/(x1+x2+x3+x4),
                               c(vG,vRdm1,vE),
                               mClone$reml$invAI) 
          
          genPar <- round(data.frame(Estimates=c(vG,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean), 
                                     SE=c(SEvG,SEvRdm1,SEvRdm2,SEvE,NA,h2GSE,matrix(NA,nrow=5,ncol=1)), 
                                     row.names = genParNames),genPar_digits)
        }else{
          genPar <- round(data.frame(Estimates=c(vG,vRdm1,vRdm2,vE,vPhen,h2G,h2mc,accClone,CVg,CVe,Mean),
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
    
    genParFile <- paste0(varResp,"_genPar",".txt")
    
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
    
    accFile <- paste0(varResp,"_acc",".txt")
    
    if (file.exists(accFile)){
      file.remove(accFile)
    }
    
    sink(accFile, append=TRUE, type = "output")
    
    cat("\n----------------------------------------------------------------------------------------------\n")
    cat("                                             |Accuracy| \n"                                        )
    cat("----------------------------------------------------------------------------------------------\n\n")
    cat("----------------------------------------- Treatment accuracy ---------------------------------\n\n")
    if(treatment=="Clone"){
      print(r2Clone)}
    if(treatment=="Prog"){
      print(r2Prog)
      cat("\n")
      cat("----------------------------------------- Individual accuracy --------------------------------\n\n")
      print(r2Ind_df)
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
    
    if(treatment=="Prog" & otimizeSelection==TRUE){
      
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