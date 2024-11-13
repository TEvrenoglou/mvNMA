structure <- function(data,sm){
  
  event_combined <- cbind.data.frame(data[,grep("event",names(data))])
  
  #data_new <- cbind.data.frame(studlab,treat,n,event_combined)
  
  data$treat <- factor(data$treat)
  
  #data_new$treat1 <- relevel(data_new$treat1,ref = reference.group)
  
  ## create labels 
  
  labs <- character2numeric_new(data$treat)
  
  labtreat <- c(labs$treat)
  
  names(labtreat) <- labs$label
  
  ### Match treatments with their labels
  data$label <- NA
  
  for(i in 1:nrow(data)){
    
    for(j in 1:nrow(labs)){
      
      if(data$treat[i]==labs$treat[j]){
        
        data$label[i]=labs$label[j]
        
      }
    }
  }
  
  ## create a numeric studlab for each study
  
  u <- unique(data$studlab)
  
  A <- as.data.frame(table(data$studlab))
  
  A <- A[match(u, A$Var1),]
  
  row.names(A) <- 1:nrow(A)
  
  studlab_numeric <- as.numeric(rep(row.names(A),A$Freq))
  
  data$studlab_numeric <- studlab_numeric
  
  #data_new$treat1 <- NULL
  
  labels <- data$label
  
  ## calculate treatment effects for each outcome
  
  p <- list()
  
  r <- list()
  
  r1 <- list()
  
  p1 <- list()
  
  for(j in 1:ncol(event_combined)){
    
    r[[j]] <- cbind.data.frame(data$studlab,
                               data$studlab_numeric,
                               data$treat,
                               data$n,
                               event_combined[,j],
                               labels)  
    
    names(r[[j]]) <- c("studlab","studlab_numeric","treat","n","event","label")
  
  
    p[[j]] <- pairwise(studlab = studlab, event =event , n=n, data = r[[j]],treat = treat,sm=sm,keep.all.comparisons = T)
    
    p[[j]] <- multi_arm(p[[j]])
    
    p[[j]]$outcome <- j
    
     }
  
  data_final <- list.rbind(p)
  
  data_final<- data_final %>% 
    group_by(studlab_numeric) %>% 
    arrange(desc(label1),.by_group = T)
  
  out <- list("new_data" = data_final,
              "data" = data,
              "labtreat" = labtreat,
              "n_outcomes" = length(p)
              #"reference.group" = reference.group
              )
  
  return(out)
}



character2numeric_new <- function(treat){
  
  treats <- cbind.data.frame(levels(treat))
  
  treats$numeric <- 1:nrow(treats)
  
  names(treats) <- c("treat","label")
  
  return(treats)
  
}


catch <- function(argname, matchcall, data, encl)
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)


jags_data <- function(dat){
  
  ## number of studies  
  
  Ns <- length(unique(dat$data$studlab))    
  
  ## labtreat 
  
  labtreat <- dat$labtreat
  
  ## number of outcomes
  
  n_outcomes <- dat$n_outcomes
  
  ## arms per study
  
  arms <- as.data.frame(table(dat$data$studlab))
  
  names(arms) <- c("study","narms")
  
  max_arms <- max(arms$narms)
  
  ## number of 2 arms studies
  
  two_arm <- length(which(arms$narms==2))
  
  # new_data <- character2numeric(data)
  
  treat_data <- create_T(dat$new_data,max_arms=max_arms)
  
  # jags_data <- long2jags(data = dat$data,
  #                        event = event1,
  #                        treat = treat,
  #                        reference.group = reference.group,
  #                        n= n,
  #                        studlab = studlab
  # )

  
  # treat_data <- jags_data$T
  # 
  # for(i in 1:nrow(treat_data)){
  #   
  #   for(j in 1:ncol(treat_data)){
  #     
  #     treat_data[i,j] <- ifelse(is.na(treat_data[i,j]),0,treat_data[i,j])  
  #     
  #   }
  #   
  # }
  
  
  ## extract vector with treatment effects
  
  y <- dat$new_data$TE
  
  var <- list()
  
  names_vec <- c()
  
  for(i in 1:n_outcomes){
    
    var[[i]] <- dat$new_data[dat$new_data$outcome==i,]$seTE^2
    
    var[[i]] <- ifelse(is.na(var[[i]]),10000,var[[i]])
    
    names_vec[i] <- paste("var",i,sep = "")
    
  }
  
  var_f <- list.cbind(var) 
  
  var_f <- as.data.frame(var_f)
  
  names(var_f) <- names_vec
  
  # attach(var_f)
  
  dat_f <- list("y"=y,
                "var" = var_f,
                "T" = treat_data,
                "Ns" = Ns,
                "N2h" = two_arm,
                "na" = arms,
                "labtreat" = labtreat
  )
  
  return(dat_f)
  
}

multi_arm <- function(dat){
  
  u <- unique(dat$studlab)    
  
  r <- list()
  
  t <- list()
  
  E <- list()
  
  for(i in 1:length(u)){
    
    r[[i]] <- dat %>% 
      filter(studlab == u[i])
    
    if(nrow(r[[i]])>2){
      
      t[[i]] <- as.data.frame(table(r[[i]]$treat2))
      
      E[[i]] <- t[[i]]$Var1[which(t[[i]]$Freq==max(t[[i]]$Freq))]
      
      
      r[[i]] <- r[[i]] %>% 
        filter(treat2 == E[[i]])
    }
    
  }
  
  harmonized <- list_rbind(r)
  
  
  return(harmonized)
}


create_T <- function(data1,max_arms){
  
  u <- unique(data1$studlab)  
  
  mat <- matrix(NA,nrow = length(u),ncol=max_arms)
  
  r <- list()
  
  treats <- list()
  
  for (i in 1:length(u)){
    
    r[[i]] <- data1 %>% 
      filter(studlab==u[i])
    
    
    treats[[i]] <- c(unique(r[[i]]$label2),unique(r[[i]]$label1))
    
    #treats[[i]] <- ifelse(length(treats[[i]])==max_arms,treats[[i]],c(treats[[i]],NA))
    
    mat[i,][1:length(treats[[i]])] <- treats[[i]]
    
    # replace NA's with 0
    #T[i,] <- ifelse(is.na(T[i,]),0,T[i,])
    
    
  }  
  
  return(mat)
  
}


