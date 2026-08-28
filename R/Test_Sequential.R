## Score 0: nothing
## Score 1: Classical
## Score 2: DTRW
## Score 3: LDM
## Score 4: YANG

######## 1- CLYD ########
DT_CLYD = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##Classical
  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision
  if(dec!="NO"){     return(score=1)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision }
  if(dec!="NO"){     return(score=3)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ##DTRW
  if(dec=="NO"){ dec=Test_DTRW_Indep(X, alpha=alpha)$decision }
  if(dec!="NO"){     return(score=2)}

  ## No model
  if(dec=="NO"){ return(score)}
}

######## 2 - CLDY ########
DT_CLDY = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##Classical
  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision }
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}


  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ## No model
  if(dec=="NO"){ return(score)}
}

######## 3 - CYLD ########
DT_CYLD = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##Classical
  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  }
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}
}

######## 4 - CYDL ########
DT_CYDL = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){

  score=0

  ##Classical
  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  }
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}


}

######## 5 - CDYL ########
DT_CDYL = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){

  score=0

  ##Classical
  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  }
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}
}

######## 6 - CDLY ########
DT_CDLY = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){

  score=0

  ##Classical
  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){
    dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  }
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 7- LCYD ########
DT_LCYD = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){

  score=0

  ## LDM
  dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 8- LCDY ########
DT_LCDY = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){

  score=0

  ## LDM
  dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 9- LYCD ########
DT_LYCD = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## LDM
  dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 10- LYDC ########
DT_LYDC = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){

  score=0

  ## LDM
  dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 11- LDYC ########
DT_LDYC = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## LDM
  dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){  dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ##Classical
  if(dec=="NO"){ dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 12- LDCY ########
DT_LDCY = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## LDM
  dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return(score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 13- YCLD ########
DT_YCLD = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## YANG
  dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision
  if(dec!="NO"){    return(score=4)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){dec=Test_DTRW_Indep(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 14- YCDL ########
DT_YCDL = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## YANG
  dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){dec=Test_DTRW_Indep(X, alpha=alpha)$decision }
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 15- YLCD ########
DT_YLCD = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## YANG
  dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision
  if(dec!="NO"){    return( score=4)}

  ## LDM
  dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){  dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ##DTRW
  if(dec=="NO"){ dec=Test_DTRW_Indep(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=2)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 16- YLDC ########
DT_YLDC = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## YANG
  dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision
  if(dec!="NO"){    return( score=4)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 17- YDLC ########
DT_YDLC = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## YANG
  dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision
  if(dec!="NO"){    return( score=4)}

  ##DTRW
  if(dec=="NO"){dec=Test_DTRW_Indep(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 18- YDCL ########
DT_YDCL = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ## YANG
  dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision
  if(dec!="NO"){    return( score=4)}

  ##DTRW
  if(dec=="NO"){
    dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  }
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 19- DCLY ########
DT_DCLY = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##DTRW
  dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision }
  if(dec!="NO"){    return( score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 20- DCYL ########
DT_DCYL = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##DTRW
  dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=2)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return( score=4)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}


  ## No model
  if(dec=="NO"){ return(score)}

}

######## 21- DYCL ########
DT_DYCL = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##DTRW
  dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision }
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision }
  if(dec!="NO"){    return( score=3)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 22- DYLC ########
DT_DYLC = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##DTRW
  dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=2)}

  ## YANG
  if(dec=="NO"){dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision }
  if(dec!="NO"){    return( score=4)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision }
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 23- DLYC ########
DT_DLYC = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##DTRW
  dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}

  ## YANG
  if(dec=="NO"){dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision}
  if(dec!="NO"){    return( score=4)}

  ##Classical
  if(dec=="NO"){dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## No model
  if(dec=="NO"){ return(score)}

}

######## 24- DLCY ########
DT_DLCY = function(X, alpha=0.05, RSq = 0.8, warmup = 2, K=NULL){
  score=0

  ##DTRW
  dec=Test_DTRW_Indep(X, alpha=alpha)$decision
  if(dec!="NO"){    return( score=2)}

  ## LDM
  if(dec=="NO"){ dec=Test_LDM_Regression(X, alpha=alpha, RSq = RSq)$decision}
  if(dec!="NO"){    return( score=3)}

  ##Classical
  if(dec=="NO"){ dec=Test_iid_BoxJenkins(X, alpha=alpha)$decision}
  if(dec!="NO"){    return( score=1)}

  ## YANG
  if(dec=="NO"){ dec=Test_YNM_Geom(X=X,alpha=alpha, warmup = warmup, K=K )$decision }
  if(dec!="NO"){    return( score=4)}

  ## No model
  if(dec=="NO"){ return(score)}

}

