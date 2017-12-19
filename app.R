#!/usr/bin/Rscript 
# Author: Henrik Zauber
# v 1.0000  

library(shiny)
library(data.table)
library(RSQLite)
library(shinyBS)
library(parallel)
library(optparse)

# Presettings--------

Use_parallel <- F
# ServerType <- "test"
ScoreThreshold <- 50
CurrentWD <- getwd()
includeMouse <- "true"
rescaledHydrophobicities <- "true"
AllGenes <- ""
AllSeq <- ""
AllSpec <- ""

TypeList <-list(PRM= c("Agilent TOF",
                       "Thermo Fusion",
                       "Thermo Q Exactive"
                       ),
                SRM=c("AB Sciex",
                      "Agilent",
                      "Bruker",
                      "Thermo Quantiva",
                      "Thermo TSQ",
                      "Waters",
                      "Waters Synapt transfer",
                      "Waters Synapt trap"
                      
                ))




# Get number of available Spectra and Genes
dbpath    <- "MagicTool.sqlite"
if((!exists("AllSeq")|0)&file.exists(dbpath)){ 
  db <- dbConnect(SQLite(),dbname = dbpath)
  InitialSequencesInfo <- dbGetQuery(db,paste("SELECT Sequence FROM Spectra WHERE Score >=",ScoreThreshold))
  InitialSequencesInfo2 <- dbGetQuery(db,"SELECT * FROM Sequences")
  
  InitialSequencesInfo2 <- InitialSequencesInfo2[InitialSequencesInfo2$MissingCleavages<=1&!(InitialSequencesInfo2$Gene_names== ""&InitialSequencesInfo2$Unique_identifiers==""),]
  
  InitialSequencesInfo <- InitialSequencesInfo[!is.na(match(InitialSequencesInfo$Sequence,InitialSequencesInfo2$Sequence)),]
  
  sp4 <- InitialSequencesInfo2$Gene_names[!is.na(match(InitialSequencesInfo2$Sequence,unique(InitialSequencesInfo)))]
  
  dbDisconnect(db)
  rm(db)
  AllSeq <- length(unique(unlist(InitialSequencesInfo)))
  AllSpec <- length(InitialSequencesInfo)
  
  sp4 <- InitialSequencesInfo2$Gene_names[!is.na(match(InitialSequencesInfo2$Sequence,unique(InitialSequencesInfo)))]
  humi <- unique(unlist(strsplit(InitialSequencesInfo2$Gene_names,";")))
  AllGenes <- length(unique(unlist(strsplit(sp4[sp4 != ""],";"))))
}

wdinit <- getwd()
mqpar <- ""
try(mqpar <- readLines("mqpar.xml",warn = F),silent = T)


# Functions --------
ExportWrapper <- function(SpecFilFil,input,Transitions,file,session,KmassShift= 8.0141988132,RmassShift = 10.008268600){
  
  
  if(input$SILAC != "none"){
    se <- unlist(unique(SpecFilFil$Sequence))
    # se[2] <- "KKRR"
    
    MassShift <- function(se,AA,MassShift,SRM){
      AAMatch <- gregexpr(AA,se)
      AAvec <- rep(0,length(se))
      AAvec[AAMatch!= "-1"] <- lengths(AAMatch[AAMatch!= "-1"])
      AAvec <- AAvec*MassShift
      AAposition <- sapply(AAMatch,paste,collapse = "#")
      return(cbind(AAvec,AAposition))
    }
    Kshift <-MassShift(se,"K",KmassShift)
    Rshift <-MassShift(se,"R",RmassShift)
    ReMatchVector <-match(SpecFilFil$Sequence,se)
    
    Kshift <- Kshift[ReMatchVector,]
    Rshift <- Rshift[ReMatchVector,]
    InsertShift <- as.numeric(SpecFilFil$mz)+(as.numeric(Kshift[,1])+as.numeric(Rshift[,1]))/as.numeric(SpecFilFil$Charge)
    SpecFilFilShifted <- SpecFilFil
    SpecFilFilShifted$mz <- InsertShift
    
    
    if(input$SRM == "SRM"){
      # only consider fragments, that carry a R or K
      SEL <- grep("[xyz]",SpecFilFil$Matches,invert = F)
      SpecFilFil <- SpecFilFil[SEL,]
      SpecFilFilShifted <- SpecFilFilShifted[SEL,]
      Kshift <- Kshift[SEL,]
      Rshift <- Rshift[SEL,]
      
      YMatches <- gsub("[+xyz]","",SpecFilFil$Matches)
      Kfactor <- DoubleFactor(Kshift)
      Rfactor <- DoubleFactor(Rshift)
      
      Kfactor <- Kfactor*KmassShift
      Rfactor <- Rfactor*RmassShift
      
      InsertShift <- as.numeric(SpecFilFil$Masses)+(as.numeric(Kfactor)+as.numeric(Rfactor))
      SpecFilFilShifted$Masses <- InsertShift
      SpecFilFilShifted$GeneSymbol <- paste(SpecFilFilShifted$GeneSymbol ,"heavy",sep = "_")
    }else{
      SpecFilFilShifted$GeneSymbol <- paste(SpecFilFilShifted$GeneSymbol ,"heavy",sep = "_")
      
    }
    SpecFilFilTemplate <- rbind(SpecFilFil,SpecFilFilShifted)
    
  }else{
    SpecFilFilTemplate <- SpecFilFil
  }
  
  SpecFilFil$Start[SpecFilFil$Start<0] <- 0
  SpecFilFil$Retention_time <- NULL
  
  
  MachineType   <- input$ExportTypes
  TemplateList <- "Error during Export."
  try(TemplateList <- ExportFormats(SpecFilFil = SpecFilFilTemplate,MachineType = input$ExportTypes,DwellTime = input$FixedDwellTime,Skyline_CE = 20,AbundanceBinnedDwellTimes = input$AbuBinDwell))
  Template <- TemplateList$Template
 
  rownames(Template) <- NULL
  
  tempdir = paste(tempdir(),basename(tempfile()),sep = "/")
  dir.create(tempdir,showWarnings = F)
  tempdir <- paste(tempdir,"InclusionList",sep = "/")
  
  dir.create(tempdir,showWarnings = F)
  try(file.copy("ReadMe.pdf",paste(tempdir,"ReadMe.pdf",sep = "/")),silent = T)
  export <- Progress$new(session)
  on.exit(export$close())
  
  export$set(message = 'Extracting Transitions',
             detail = 'This may take a while...')
  tab <- Transitions
  tab <- tab
  Templates <- Template
  wd <- getwd()
  setwd(tempdir)
  export$set(message = 'Extracting Transitions',
             detail = 'Done')

  dlist = list(SpecFilFil,tab,Template)
  export$set(message = 'Writing zip...',
             detail = 'busy')
  MSMSExport(tab,SpecFilFil,"./")
  for(i in 1:length(dlist)){
    sep = "\t"
    ending = ".txt"
    if(i == 3){
      sep = TemplateList$sep
      ending = paste(".",TemplateList$ending,sep = "")
    }
    write.table(x = as.matrix(dlist[[i]]), file = paste0(c("SpectraTable.txt","TransitionTable.txt",paste("InclusionList",ending,sep = ""))[i]), row.names = FALSE,sep = sep,quote = F)
    
    if(i == 1){
      try(temp <- data.frame(dlist[[i]]),silent = T)
      OutX <<- c()
      Psplit <- strsplit(as.character(temp$Proteins),";")
      FunOut <<- c()
      FastaOut <-sapply(  unique(unlist(strsplit(as.character(temp$Proteins),";")))
                          ,function(x){
                            # x <<-x 
                            
                            Sel <- sapply(Psplit,function(y){any(x == y)})
                            
                            FunOut <<- rbind(FunOut,cbind(x,paste(unique(temp$GeneSymbol[Sel]),collapse = ";"),as.character(temp$Sequence[Sel])))
                            return(NULL)
                          })
      
      AggP <-aggregate(FunOut[,c(1,2)],list(FunOut[,3]),function(x){paste(unique(x),collapse = ";")})
      
      AlmostaFasta <- (as.matrix(AggP))
      # stop()
      # stop()
      AlmostaFasta[,2]<- paste(paste(">",AlmostaFasta[,2],sep = ""),AlmostaFasta[,3],sep = "|")
      
      write.table(AlmostaFasta,paste0("Skyline_Insert_Peptides.txt"),quote = F,row.names = F,sep = "\t",col.names = F  )
    }
    
    
  }
  inputList <- as.list(input)
  inputList[lengths(inputList) ==0] <- "NA"
  inputList$RetInput <- NULL
  write(paste(names(inputList),sapply(inputList,paste,collapse = ";"),sep = ":\t"),"Picky_Parameters.txt")
  
  # save(inputList,file = "Picky_Settings")
  # zip(zipfile =  file, files = "../InclusionList",flags = "-9 -y -r -q")
  system(paste("zip -9 -y -r -q",file,paste(list.files(),collapse = " ")))
  
  export$set(message = 'Writing zip...',
             detail = 'Done')
  setwd(wd)
}

DoubleFactor <- function(CheckDouble){
  CheckDouble <<- CheckDouble
  Double <- grep("#",CheckDouble[,2])
  OutFac <- rep(0,dim(CheckDouble)[1])
  try(OutFac[as.numeric(CheckDouble[,2])>-1]<- 1,silent = T)
  if(length(Double) > 0){
    CheckDoubleTemp <- CheckDouble[Double,]
    CheckDoubleTempVec <- strsplit(CheckDoubleTemp[,2],"#")
    nse <- nchar(SpecFilFil$Sequence)
    DoubleFactorY <- sapply(1:length(Double),function(i){
      tempx <- as.numeric(CheckDoubleTempVec[[i]])
      tempM <- YMatches[Double][i]
      tempS <- nse[Double][i]
      tempMC <- tempS-as.numeric(tempM)
      Bi <- tempx > tempMC
      length(Bi[Bi])
    })
    
    OutFac[Double] <- DoubleFactorY
    
  }
  return(OutFac)
}

ExportFormats <- function(SpecFilFil,MachineType,DwellTime = 20,Skyline_CE = 20,AbundanceBinnedDwellTimes = T){
  if(AbundanceBinnedDwellTimes){
    DwellTime <- SpecFilFil$DwellTimes
  }else{
    DwellTime <- DwellTime
  }
  
  SpecFilFil$Picky_RET <- SpecFilFil$Start+(SpecFilFil$End-SpecFilFil$Start)/2
  
  if(MachineType == "Agilent TOF"){
    
    Template <- cbind("True",
                      SpecFilFil$mz,
                      SpecFilFil$Charge,
                      SpecFilFil$Picky_RET,
                      max(SpecFilFil$End-SpecFilFil$Start),
                      "Narrow (~1.3 m/z)",
                      Skyline_CE, # Check wether this could be otpimized by any equation...
                      "" # Acquisition Time (ms/spec), what is this?
    )
    colnames(Template) <- c("On","Prec. m/z","Z","Ret. Time (min)","Delta Ret. Time (min)","Iso. Width","Collision Energy","Acquisition Time (ms/spec)")
    # csv
    sepv = ","
    endv = "csv"
    headerv = T
  }
  
  if(MachineType == "Thermo Fusion"){
    Template <- cbind(SpecFilFil$mz,
                      SpecFilFil$Charge,
                      SpecFilFil$Start,
                      SpecFilFil$End,
                      SpecFilFil$CollisionEnergy
    )
    
    colnames(Template) <- c("m/z","z","t start (min)","t end (min)","CID Collision Energy (%)")
    #csv
    sepv = ","
    endv = "csv"
    headerv = T
  }
  
  if(MachineType == "Thermo Q Exactive"){
    Template <- cbind(SpecFilFil$mz,
                      "","","",
                      SpecFilFil$Charge,
                      "Positive",SpecFilFil$Start,SpecFilFil$End,
                      SpecFilFil$CollisionEnergy,
                      "NCE","",
                      paste(SpecFilFil$GeneSymbol,SpecFilFil$Proteins,SpecFilFil$Modified_sequence,SpecFilFil$SpecID,sep = "#"))
    colnames(Template) <- c("Mass [m/z]",
                            "Formula [M]","Formula type","Species",
                            "CS [z]",
                            "Polarity","Start [min]","End [min]",
                            "(N)CE",
                            "(N)CE type","MSX ID",
                            "Comment")
    #csv
    sepv = ","
    endv = "csv"
    headerv = T
  }
  
  # COntinue to check here!
  if(MachineType == "Waters Synapt transfer"){
    
    # aggregate(SpecFilFil$Masses)
    IDSel <- paste(SpecFilFil$SpecID,SpecFilFil$GeneSymbol)
    Template <- sapply(IDSel,function(x){
      x <- x
      tempx <- SpecFilFil[IDSel == x,]
      tempx <- tempx[order(tempx$Intensities_RAW,decreasing = T),]
      limit = 6
      if(dim(tempx)[1] > limit){
        tempx <- tempx[1:limit,]
      }
      if(dim(tempx)[1]< limit){
        INT  <- c(tempx$Intensities_RAW,rep(0,limit-dim(tempx)[1]))
      }
      tempxS <- tempx[1,]
      Template <- c("0",
                    tempxS$Start,
                    tempxS$End,
                    tempxS$mz,
                    INT,
                    4,
                    4,
                    20,
                    30,
                    30,
                    0,
                    0,
                    199,
                    paste(tempxS$GeneSymbol,tempxS$Proteins,tempxS$Modified_sequence,tempxS$SpecID,sep = "#")
                    
      )
      return(Template)
      
    })
    Template <- t(Template)
    rownames(Template) <- NULL
    colnames(Template) <- c("Function channel","Retention Time Start","Retention Time End","Set Mass","Mass Fragments 1","2","3","4","5","6","Trap CE Start","Trap CE End","Transfer CE Start","Transfer CE End","CV","EDCMass","DT Start","DT End","compound name")
     
    # template <- read.csv("Waters Synapt transfer_IL.mrm")
    sepv = ","
    endv = "mrm"
    headerv = T
    #mrm
  }

  if(MachineType == "Waters Synapt trap"){
    # template <- read.csv("Waters_Xevo_QTOF_IL.mrm")
    # aggregate(SpecFilFil$Masses)
    IDSel <- paste(SpecFilFil$SpecID,SpecFilFil$GeneSymbol)
    Template <- sapply(IDSel,function(x){
      x <- x
      tempx <- SpecFilFil[IDSel == x,]
      tempx <- tempx[order(tempx$Intensities_RAW,decreasing = T),]
      limit = 6
      if(dim(tempx)[1] > limit){
        tempx <- tempx[1:limit,]
      }
      if(dim(tempx)[1]< limit){
        INT  <- c(tempx$Intensities_RAW,rep(0,limit-dim(tempx)[1]))
      }
      tempxS <- tempx[1,]
      Template <- c("0",
                    tempxS$Start,
                    tempxS$End,
                    tempxS$mz,
                    INT,
                    20,
                    30,
                    2,
                    2,
                    30,
                    0,
                    0,
                    199,
                    paste(tempxS$GeneSymbol,tempxS$Proteins,tempxS$Modified_sequence,tempxS$SpecID,sep = "#")
                    
      )
      return(Template)
      
    })
    Template <- t(Template)
    rownames(Template) <- NULL
    colnames(Template) <- c("Function channel","Retention Time Start","Retention Time End","Set Mass","Mass Fragments 1","2","3","4","5","6","Trap CE Start","Trap CE End","Transfer CE Start","Transfer CE End","CV","EDCMass","DT Start","DT End","compound name")
    
    sepv = ","
    endv = "mrm"
    headerv = T
    #mrm
  }
  
  
  if(MachineType == "AB Sciex"){
    Template <- cbind(SpecFilFil$mz,
                      SpecFilFil$Masses,
                      DwellTime,#set to 20, dwell Time? 
                      paste(paste(">",SpecFilFil$Proteins,"|",SpecFilFil$GeneSymbol,sep = ""),SpecFilFil$Sequence,paste(paste("+",SpecFilFil$Charge,sep = ""),SpecFilFil$Matches,sep = ""),sep = "."),#comment
                      "",# no idea, something between 60 and 110 retention time?
                     Skyline_CE# no idea, something between 20 and 60# collision energy?
    )
    sepv = ","
    endv = "csv"
    headerv = T
    #csv, no header
  }
  
  if(MachineType == "Agilent"){
    Template <- cbind(paste(">",paste(SpecFilFil$Proteins,SpecFilFil$GeneSymbol,sep = "|"),sep = ""),
                      SpecFilFil$Sequence,
                      FALSE,# ISTD? WHat is this?
                      SpecFilFil$mz,#Precursor Ion
                      "Unit",
                      SpecFilFil$Masses,
                      "Unit",
                      DwellTime,#DwellTime
                      130,# Fragmentor?
                      Skyline_CE,# Collision Energy
                      4,#Cell Accelerator
                      SpecFilFil$Matches,
                      try(unlist(aggregate(SpecFilFil$Intensities_RAW,list(unlist(SpecFilFil$Proteins)),rank,ties.method = "random")$x)) #library Rank
                      # What is Fragmentor, ISTD and Cell Accelerator and Library Rank, where to put retention time  
                      # csv and header
                      
    )
    
    colnames(Template)    <-  c("Compound Group","Compound Name","ISTD?","Precursor Ion","MS1 Res","Product Ion","MS2 Res","Dwell","Fragmentor","Collision Energy","Cell Accelerator Voltage","Ion Name","Library Rank")
    sepv = ","
    endv = "csv"
    headerv = F              
    #csv, no header
  }
  
  if(MachineType == "Bruker"){
    # template <- read.csv("Bruker.csv",header = T,check.names = F)
    SeqID <-paste(SpecFilFil$Sequence,SpecFilFil$mz)
    SeqIDT <- table(SeqID)
    SeqIDT[match(SeqID,names(SeqIDT))]
    
    Template <- cbind(paste(SpecFilFil$Sequence,sapply(SpecFilFil$Charge,function(x){paste(rep("+",x),collapse = "")}),sep = ""),
                      SpecFilFil$Picky_RET,
                      max(SpecFilFil$End-SpecFilFil$Start),
                      "",#CAS Number
                      rank(SpecFilFil$Picky_RET,ties.method = "first"),#Retention Index?
                      "MRM",#Scan Type
                      "Positive",#polarity
                      100,# Scan Time in ms
                      "LCMS",
                      "",#source
                      "",#regulation
                      "",#Classification
                      paste(SpecFilFil$Proteins,SpecFilFil$GeneSymbol,sep = "|"),#comment
                      SeqIDT,#Transitions Count
                      SpecFilFil$mz,
                      "",
                      "",
                      SpecFilFil$Masses,
                      "",
                      "",
                      Skyline_CE, # Collision Energy
                      DwellTime,
                      1,
                      SpecFilFil$Masses,
                      0,
                      "",
                      "",
                      "",
                      "",
                      "",
                      "",
                      "",
                      "",
                      "",
                      "",
                      "",
                      ""
    )
    colnames(Template) <- c("Compound Name","Retention Time","Retention Time Window","CAS Number","Retention Index","Scan Type","Polarity","Scan Time (ms)","Separation Method","Source","Regulation","Classification","Comment","Transitions Count","Q1 First Mass","Q1 Last Mass","Q1 Resolution","Q3 First Mass","Q3 Last Mass","Q3 Resolution","Collision Energy","Dwell Time (ms)","Is Quantifier","Quantifier Ions","Is Qualifier","Qualifier Count","Qual Mass 1","Qual Ratio 1","Qual Mass 2","Qual Ratio 2","Qual Mass 3","Qual Ratio 3","Qual Mass 4","Qual Ratio 4","Qual Mass 5","Qual Ratio 5","GUID (Dont fill this Column)")
    # What means Scan Time
    sepv = ","
    endv = "csv"
    headerv = T
  }
  
  
  # maybe not shimadzu (: Quite Complicated)
  if(MachineType == "Shimadzu"){
    # template <- read.csv("Shimadzu.txt",sep = "\t",check.names = F)
    SeqID <-paste(SpecFilFil$Sequence,SpecFilFil$mz)
    SeqIDT <- table(SeqID)
    SeqIDT[match(SeqID,names(SeqIDT))]
    
    Template <- cbind(paste(SpecFilFil$Sequence,"_light",sep = ""),
                      SeqIDT,
                      SpecFilFil$mz,
                      SpecFilFil$Masses,
                      SpecFilFil$Picky_RET,
                      "",#Unit
                      0,0,0,0,
                      0,#SpecFilFil$Start,
                      0,#SpecFilFil$End,
                      0,0,0,0,0,0, 0, 0,
                      1,#Event
                      0,0,0,0,0,
                      1,#Correction Factor
                      1,#standard concentration Factor
                      SpecFilFil$Start,#Start Time
                      SpecFilFil$End,#End Time
                      15,#Acquisition Mode
                      0,
                      0.103,# Event Time?
                      paste(SpecFilFil$Proteins,SpecFilFil$GeneSymbol,sep = "|"),#comment
                      "","","",
                      3,#Target Pause Time
                      5,#Target Dwell Time
                      "","","",
                      Skyline_CE*-1,
                      "","",
                      SpecFilFil$mz,
                      SpecFilFil$Masses,
                      75,#Ref.(1) Ratio
                      0,
                      "",
                      
                      3,#Target Pause Time
                      5,#Target Dwell Time
                      "","","",
                      Skyline_CE*-1,
                      "","",
                      SpecFilFil$mz,
                      SpecFilFil$Masses,
                      75,#Ref.(1) Ratio
                      0,
                      "",
                      3,#Target Pause Time
                      5,#Target Dwell Time
                      "","","",
                      Skyline_CE*-1,
                      "","",
                      SpecFilFil$mz,
                      SpecFilFil$Masses,
                      75,#Ref.(1) Ratio
                      0,
                      "",
                      3,#Target Pause Time
                      5,#Target Dwell Time
                      "","","",
                      Skyline_CE*-1,
                      "","",
                      SpecFilFil$mz,
                      SpecFilFil$Masses,
                      75,#Ref.(1) Ratio
                      0,
                      "",
                      
                      3,#Target Pause Time
                      5,#Target Dwell Time
                      "","","",
                      Skyline_CE*-1,
                      "","",
                      SpecFilFil$mz,
                      SpecFilFil$Masses,
                      75,#Ref.(1) Ratio
                      0,
                      "",
                      3,#Target Pause Time
                      5,#Target Dwell Time
                      "","","",
                      Skyline_CE*-1,
                      "",""
                      
                      
                      
    )
    #Acquisition Mode
    sepv = "\t"
    endv = "txt"
    headerv = T
    
  }
  
  if(MachineType == "Thermo Quantiva"){
    Template <- cbind(paste(SpecFilFil$GeneSymbol,SpecFilFil$Proteins,SpecFilFil$Modified_sequence,SpecFilFil$Matches,sep = ""),
                      SpecFilFil$Picky_RET,
                      max(SpecFilFil$End-SpecFilFil$Start),
                      "Positive",
                      SpecFilFil$mz,
                      SpecFilFil$Masses,
                      SpecFilFil$CollisionEnergy
    )
    colnames(Template) <- c("Compound",
                            "Retention Time (min)","RT Window (min)","Polarity",
                            "Precursor (m/z)",
                            "Product (m/z)","Collision Energy (V)"
    )
    sepv = ","
    endv = "csv"
    headerv = T
  }
  
  if(MachineType == "Thermo TSQ"){
    Template <- cbind(paste(SpecFilFil$GeneSymbol,SpecFilFil$Proteins,SpecFilFil$Modified_sequence,SpecFilFil$Matches,sep = "_"),
                      SpecFilFil$Start,
                      SpecFilFil$End,
                      "Positive",
                      SpecFilFil$mz,
                      SpecFilFil$Masses,
                      SpecFilFil$CollisionEnergy,
                      DwellTime,
                      ""
    )
    colnames(Template) <- c("Compound",
                            "Start Time (min)","End Time (min)","Polarity",
                            "Precursor (m/z)",
                            "Product (m/z)","Collision Energy (V)","Dwell Time (ms)","RF Lens (V)"
    )
    sepv = ","
    endv = "csv"
    headerv = T
  }
  
  if(MachineType == "Waters"){
    Template <- cbind(paste(SpecFilFil$Proteins,SpecFilFil$GeneSymbol,sep = "|"),
                      paste(SpecFilFil$Sequence,SpecFilFil$Charge,sep = "."),
                      SpecFilFil$mz,
                      SpecFilFil$Picky_RET,
                      SpecFilFil$Masses,
                      Skyline_CE,
                      35,#cone voltage
                      SpecFilFil$Sequence,
                      SpecFilFil$Matches,
                      try(unlist(aggregate(SpecFilFil$Intensities_RAW,list(unlist(SpecFilFil$Proteins)),rank,ties.method = "random")$x)) #library Rank
                      
    )
    colnames(Template) <- c("protein.name","peptide.seq","precursor.mz","precursor.retT","product.m_z","collision_energy","cone_voltage","peptide_unmod.seq","ion_name","library_rank")
    sepv = ","
    endv = "csv"
    headerv = T
  }
  return(list(Template = Template,sep = sepv,ending = endv,header = headerv))
}

medwin <- 
  function(x,y,win = 8,fun = mean,type = "x",...){
    
    ord <- order(as.numeric(as.character(x)))
    x <- x[ord]
    y <- y[ord]
    ro <- ceiling(length(y)/win)
    
    resfun<- sapply(1:length(y),function(r){
      
      st <- r-win
      if(st < 1){st = 1}
      endt = (r+win)
      if(endt > length(y)){
        endt <- length(y)
      }
      yt <- y[st:endt]
      xt <- x[r]
      ytm <- fun(yt,...)
      yts <- sd(yt,na.rm = T)
      return(c(xt,ytm,yts))
    })
    resfun <- t(resfun)
    colnames(resfun) <- c("xmu","ymu","ystd")
    
    return(list(x=resfun[,1],y = resfun[,2],ysd=resfun[,3]))
  }

MSMSExport <- function(tr,IL,outputpath,name = "msms.txt",exportRaw  = T){
  ILunique <- IL
  ILunique$Matches <- NULL
  ILunique$Intensities <- NULL
  ILunique$Intensities_RAW <- NULL
  ILunique$ID <- NULL
  ILunique$Masses <- NULL
  IL <- unique(ILunique)

  wd<- getwd()
  setwd(outputpath)
  dir.create("Spectral_Library")
  setwd("Spectral_Library")
  
  sapply(list.files(),unlink)
  # unlink(name)
  write(mqpar,"mqpar.xml")
  write(paste(c("Raw File","Sequence","Modified Sequence","Mass","m/z","Charge","Masses","Matches","Intensities","Scan Number","Modifications","Retention Time","PEP","Score"),collapse = "\t"),name)
  write(paste(c("Raw File","Sequence","Modified Sequence","Mass","m/z","Charge","Masses","Matches","Intensities","Scan Number","Modifications","Retention Time","PEP","Score"),collapse = "\t"),RawName<- paste("RAW_Spectra_",name,sep = ""))
  
  it <- 1
  for(i in unique(tr$SpecID)){
    it <- it+1
    tempi <- tr[tr$SpecID == i,]
    tempIL <- IL[IL$SpecID == i,]
    tempiRaw <- tempi[!is.na(tempi$Intensities_RAW),]
    tempiRaw$Intensities <- tempiRaw$Intensities_RAW
    tempi <- tempi[!is.na(tempi$Intensities),]
    rountIT <- 0
    for(tempi in list(tempi,tempiRaw)){
      rountIT <- rountIT +1
      if(dim(tempi)[1] > 0){
        tempi <- tempi
        Masses <-paste(tempi$Masses,sep = ";",collapse = ";")
        Matches <- paste(tempi$Matches,collapse = ";")
        Intensities <- paste(tempi$Intensities,collapse = ";")
        Mass <- as.numeric(tempIL$mz)*as.numeric(tempIL$Charge)
        mz <- tempIL$mz
        Sequence <- tempIL$Sequence
        tempIL$Retention_time <- tempIL$Start+(tempIL$Start+tempIL$End)/2
        Modified.Sequence <- tempIL$Modified_sequence
        tempIL <- tempIL
        MassesTable = c("Picky",Sequence,Modified.Sequence,Mass,mz,tempIL$Charge,Masses,Matches,Intensities,round(as.numeric(tempIL$Retention_time)*100),"Unmodified",tempIL$Retention_time,"0","100")
        MassesTable <- MassesTable
        write(paste(MassesTable,collapse = "\t"),c(name,RawName)[rountIT],append = T)
      }
    }
  }
  
  gwd <- getwd()
  setwd(wd)
  return(gwd)
  
}
FragmentPlot <- function(msms,se){
  xions <- grep("^x",msms$Matches,value = T)
  yions <- grep("^y",msms$Matches,value = T)
  zions <- grep("^z",msms$Matches,value = T)
  
  aions <- grep("^a",msms$Matches,value = T)
  bions <- grep("^b",msms$Matches,value = T)
  cions <- grep("^c",msms$Matches,value = T)
  
  allL <- list(b = bions,x = xions,a = aions,y = yions,z = zions,c = cions)
  Ions <- lapply(allL,function(x){
    x <- gsub("[abcxyz]","",x)
    x <- strsplit(as.character(x),"[^0123456789]")
    x <- as.numeric(unique(sapply(x,function(x){x[length(x)]})))
    x <- unlist(x[!is.na(x)])
    if(length(x) == 0){x <- NA}else{x <- x}
  })
  
  Adds <- lapply(allL,function(x){
    x <- gsub("y|b|x|a","",x)
    x <- grep("[[:punct:]]",x,value = T)
    n <- strsplit(as.character(x),"[^0123456789]")
    n <- as.numeric(sapply(n,function(x){x[1]}))
    return(cbind(x,n))
  }
  )
  
  Adds <- lapply(Adds,function(x){
    x  <- apply(x,1,function(x){
      tx <- gsub(paste("^",x[2],sep = ""),"",x[1])
      return(c(tx,x[2]))
      
    })
    return(t(x))
    
  })
  Adds <- lapply(Adds,unique)
  # Adds <- unique(Adds)
  
  
  se <- gsub("_","",se)
  seSplit <- strsplit(se,"")
  
  TotalL <- nchar(as.character(se))
  Seq <- unlist(strsplit(as.character(se),""))
  Ycol <- rep("darkgrey", TotalL)
  Ycol[match(Ions$x,1: TotalL)] <- "dodgerblue3"
  Ycol[match(Ions$y,1: TotalL)] <- "dodgerblue3"
  Ycol[match(Ions$z,1: TotalL)] <- "dodgerblue3"
  
  Ycol <- rev(Ycol)
  Bcol <- rep("darkgrey", TotalL)
  Bcol[match(Ions$b,1: TotalL)] <- "red"
  Bcol[match(Ions$a,1: TotalL)] <- "red"
  Bcol[match(Ions$c,1: TotalL)] <- "red"
  
  par(mai = c(0,0,0,0))
  
  PlotL <- 21
  cexMod <- 0
  if(TotalL > 20){
    PlotL <- TotalL +1; cexMod <- (20/TotalL)
    if(cexMod < 0.5){
      cexMod <- 0.5
    }
    
  }
  
  aseries <- 1:TotalL
  
  aseries[match(Ions$c,1: TotalL)] <- paste("c",aseries[match(Ions$c,1: TotalL)],sep = "")
  aseries[match(Ions$b,1: TotalL)] <- paste("b",aseries[match(Ions$b,1: TotalL)],sep = "")
  
  aseries[match(Ions$a,1: TotalL)] <- paste("a",aseries[match(Ions$a,1: TotalL)],sep = "")
  
  xseries <- 1:TotalL
  xseries[match(Ions$z,1: TotalL)] <- paste("z",xseries[match(Ions$z,1: TotalL)],sep = "")
  xseries[match(Ions$y,1: TotalL)] <- paste("y",xseries[match(Ions$y,1: TotalL)],sep = "")
  
  xseries[match(Ions$x,1: TotalL)] <- paste("x",xseries[match(Ions$x,1: TotalL)],sep = "")
  
  
  for(i in 1: TotalL){
    # Sequence
    
    text(i+0.5,1.5,Seq[i],cex = 1.5-1.5*cexMod/0.7)
    # b ions
    lines(x<-c(i+0.2,i+0.8,i+0.8),y <- c(2,2,1),col = Bcol[i],lwd = 2)
    text(x<-c(i+0.5),y <- c(2),aseries[i],col = Bcol[i],pos = 3,offset = 0.2,cex = 1-1*cexMod*1.3)
    
    
    # b ions mods
    addB <- Adds$b
    if(length(addB) >0){
      if(any(addB[,2] == i)){
        texthu <- addB[addB[,2] == i,1]
        text(x<-c(i+0.5),y <- c(2),paste(texthu,collapse = "\n") ,col = "orange",pos = 3,offset = 1.1,cex = 0.7-0.5*cexMod,xpd = NA)
      }
      
    }
    
    # y ions
    lines(x<-c(i+0.1,i+0.1,i+0.7),y <- c(2,1,1),col = Ycol[i],lwd = 2)
    text(x<-c(i+0.5),y <- c(1),xseries[TotalL - i+1],col = Ycol[i],pos = 1,offset = 0.2,cex = 1-1*cexMod*1.3)

    # b ions mods
    addY <- Adds$y
    if(length(addY) >0){
      if(any(addY[,2] == (TotalL-i+1))){
        texthu <- addY[addY[,2] == (TotalL-i+1),1]
        text(x<-c(i+0.5),y <- c(1),paste(texthu,collapse = "\n") ,col = "orange",pos = 1,offset = 1.1,cex = 0.7-0.5*cexMod,xpd = NA)
      }
    }	
  }
}
SeSeq <- function(se){
  se <- gsub("_","",se)
  seSplit <- strsplit(se,"")
  
  TotalL <- nchar(as.character(se))

  PlotL <- 21
  cexMod <- 0
  if(TotalL > 20){
    PlotL <- TotalL +1
  }
  return(PlotL)
}

Wrapper <- function(WHAT,FROM,FILTERColumn= NULL,PATTERN,FUN = "in",MultiFun = "or"){
  p1 <-paste("select",WHAT,"from",FROM,sep = " ")
  PATTERN <- paste("'",PATTERN,"'",sep = "")
  
  if(length(FILTERColumn) >0){
    if(length(PATTERN) > 0){
      PATTERN <- PATTERN
      PATTERN <- paste("(",paste(PATTERN,collapse = ","),")")
    }
    
    p1 <- paste(p1,"where",FILTERColumn,FUN,PATTERN,sep = " ")
  }
}
FilterFun <- function(QueryTable,Fragmentation,MassAnalyzer,Charge){
  Frag <- grepl(paste(paste("^",Fragmentation,"$",sep = ""),collapse = "|"),QueryTable$Fragmentation)
  MaAn <- grepl(paste(paste("^",MassAnalyzer,"$",sep = ""),collapse = "|"),QueryTable$Mass_analyzer)
  Char <- grepl(paste(paste("^",Charge,"$",sep = ""),collapse = "|"),QueryTable$Charge)
  return(QueryTable[Frag&MaAn&Char,])
}
PeptideHydrophobicity <- function(Seq){
  AA <- c("W","F","L","I","M","V","Y","A","T","P","E","D","C","S","Q","G","N","R","H","K")
  Rc <- c(11,10.5,9.6,8.4,5.8,5,4,0.8,0.4,0.2,0,-0.5,-0.8,-0.8,-0.9,-0.9,-1.2,-1.3,-1.3,-1.9)
  RcNT <- c(-4,-7,-9,-8,-5.5,-5.5,-3,-1.5,5,4,7,9,4,5,1,5,5,8,4,4.6)
  rc <- cbind(AA,Rc,RcNT)
  colnames(rc) <- c("AA","Rc","RcNt")
  try(rc <- data.frame(rc,stringsAsFactors = F),silent = T)
  rc$Rc <- as.numeric(rc$Rc)
  rc$RcNt <- as.numeric(rc$RcNt)
  
  if(!exists("rc")){
    "RcConstantsTable is Missing!"
  }
  SeqL <- unlist(strsplit(Seq,""))
  rcm <- rc[match(SeqL,rc$AA),]
  
  N <- nchar(Seq)
  
  K_L = 1
  if(N< 10){K_L = 1-0.027*(10-N)}
  if(N> 20){K_L = 1-0.014*(N-20)}
  
  
  
  if(N< 38){
    H = K_L*(sum(rcm$Rc)+0.42*rcm$RcNt[1]+0.22*rcm$RcNt[2]+0.05*rcm$RcNt[3])
  }
  if(N>= 38){
    H = K_L*(sum(rcm$Rc)+0.42*rcm$RcNt[1]+0.22*rcm$RcNt[2]+0.05*rcm$RcNt[3])-0.3*(K_L*(sum(rcm$Rc)+0.42*rcm$RcNt[1]+0.22*rcm$RcNt[2]+0.05*rcm$RcNt[3])-38)
    
  }
  
  return(H)
}
ModelRT <- function(a,b,sel = rep(T,length(a)),windows = 10,AutoWindowQuantile = 70,alpha = 0.1,QuantileLimit = c(0.05,0.95),...){
  lo <- windows
  if(length(a) <50*15){
    lo <- length(a)/2
  }

  windows <- seq(min(a),max(a),length.out = lo)
  Wi <- cbind(windows[-length(windows)],windows[-1])
  
  Outlier <- rep(0,length(a))
  BiquanTable <- c()
  BiquanTableAuto <- c()
  
  for(i in 1:dim(Wi)[1]){
    TempI <- Wi[i,]
    bi <- b[sel <-(a >= TempI[1]&a<=TempI[2])]
    ai <- a[(b >= TempI[1]&b<=TempI[2])]
    bi <- bi
    
    Biquan <- quantile(bi,probs = QuantileLimit,na.rm = T)
    AutoWindowQuantile2 <- (100-AutoWindowQuantile)/2/100
    BiquanAuto <- quantile(ai,probs = c(0+AutoWindowQuantile2,1-AutoWindowQuantile2),na.rm = T)
    
    selb <- bi <Biquan[1]|bi> Biquan[2]
    BiquanTable <- rbind(BiquanTable,c(Biquan,TempI))
    BiquanTableAuto <- rbind(BiquanTableAuto,c(BiquanAuto,TempI))
    
    Outlier[sel][selb] <- 1
  }
  colnames(BiquanTable) <- c(paste("Quantile",QuantileLimit,sep = "_"),"UpperWindow","LowerWindow")
  dco <- densCols(a,b)
  dco[Outlier==1] <- "red"
  hui <- loess(a[Outlier == 0]~b[Outlier == 0],...)
  return(list(Model = hui,BiquanTable = BiquanTable,BiquanTableAuto = BiquanTableAuto,outlier = dco,Windows = Wi))
}
FeatureRetentionPlot<- function(initTable,SpecFilFil,UpperLimit,RTrange,SelectedSequence = NULL,TypeName = "Features",SILAC = "none",...){
  RTrange[is.infinite(RTrange)] <- 0
  
  leFu <- apply(initTable,2,function(x){length(x[!is.na(x)])})
  par(mai = c(0.8,1,0.5,1),cex = 1.5,mgp = c(1.5,0.5,0))
  
  TypeNameFirstPlot <- TypeName
  if(SILAC != "none"){
    TypeNameFirstPlot <- paste(TypeNameFirstPlot,"(SILAC-Pairs)")
  }
  plot(1,type = "n",xlim = c(0,1),ylim = c(0,1),xlab = "elution time in min",axes = F,ylab = TypeNameFirstPlot,...)
  ElutionDensity <- matrix(leFu,dim(initTable)[1],dim(initTable)[2],byrow = T)
  image(t(ElutionDensity),col = colorRampPalette(c("white","orange"))(300),axes = F,add = T)
  
  if(dim(initTable)[1] == 1){
    image_input <- (initTable)
  }else{
    image_input <- t(initTable)
    
  }
  cbPalette <- rev(c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","tomato3","#999999")[-1])
  
  cols <- colorRampPalette(cbPalette)(max((SpecFilFil$ProteinColor)))
  cols <- cols#
  seqlength <- dim(initTable)[1]
  if(seqlength == 0|seqlength == 1){
    seqlength <- 1
    yset <- 0.65
    image(image_input,axes = F,col = cols,add = T)
    
  }else{
    yset <- (1:seqlength)/seqlength
    image(image_input,axes = F,col = cols,add = T,y = yset)
    
  }
  

  
  pr <- pretty(0:1 )

  Pos <- apply((initTable),1,function(x){min(grep(1,as.numeric(!is.na(x))),na.rm = T)})

  
  le <- dim(initTable)[1]
  if(le == 0){
    le <- 1
  }
  yp <- seq(0,1,length.out =  le )
  if(!any(as.logical(SpecFilFil$PredictionRange))){
    text((Pos/dim(initTable)[2])[!SpecFilFil$PredictionRange],yp[!SpecFilFil$PredictionRange],"*",pos = 4,col = "white",cex = 2)
  }
  
  if(!all(RTrange == 0)){
    axis(1,at = pr,labels =seq(round(RTrange[1],1),round(RTrange[2],1),length.out = length(pr)) , tck = -0.01)
  }

  ypr <- seq(0,1,length.out = seqlength)
  abline(h = (1:seqlength)/seqlength, col = "lightgrey",lty = "dotted")
  if(length(SelectedSequence) > 0){
    sel <- SelectedSequence == SpecFilFil$Modified_sequence
    Pos <- ((1:seqlength)/seqlength)[(sel)]#/dim(SpecFilFil)[1]
    if(any(sel)){
      abline(h = Pos,col = "red")
    }
  }
  
  axis(2,at = yset,label = 1:seqlength, tck = -0.01)
  fufi <- unique(cbind(SpecFilFil$Proteins,SpecFilFil$ProteinColor))
  if(length(fufi[,1]) <= 20){
    # legend("topleft",legend = fufi[,1],col = cols[unlist(fufi[,2])],pch = 20,bg = "#FFFFFF80",border = "transparent",box.col = "#FFFFFF80",cex = 0.4)
  }
  if(!all(RTrange == 0)){
    par(new = T)
    if(SILAC!= "none"){
      leFu <- leFu*2
    }
    plot(seq(0,1,length.out = length(leFu)),leFu,xlim = c(0,1),type = "n",col = "black",lwd = 3,xlab = "",ylab = "",axes = F)
    
    points(seq(0,1,length.out = length(leFu)),leFu,xlim = c(0,1),type = "s",col = "black",lwd = 3)
    
    axis(4, tck = -0.01)
    mtext(paste("# of ",TypeName ," monitored in parallel",sep = ""),4,line = 1.5,cex = 1.5)
    abline(h = UpperLimit,lty = "dashed") 
  }
}
ExtractInfoFromDB <- function(dbpath,Proteins,IDType,Fragmentation,MassAnalyzer,CurrentWD = getwd(),Charge,MissCleave = F,species = "mouse",RescaledHydrophobicities = T,RescaledHPType = "HP_all",ExcludeKPRPfragments = F,mzrange = c(300,1700),SILAC = "none"){ 
  setwd(CurrentWD)
  if(!exists("db")){
    db <- dbConnect(SQLite(),dbname = dbpath)
  }
  
  if(species == "human"){
    if(IDType == "Uniprot accession"){
      SeqID <- dbGetQuery(db,Wrapper("*","ProteinID","Unique_identifier",Proteins))
    }else{
      SeqID <- dbGetQuery(db,Wrapper("*","GeneID","Gene_name",toupper(Proteins)))
      
    }
  }
  if(species == "mouse"){
    if(IDType == "Uniprot accession"){
      SeqID <- dbGetQuery(db,Wrapper("*","ProteinID_Mouse","Unique_identifier",Proteins))
    }else{
      SeqID <- dbGetQuery(db,Wrapper("*","GeneID_Mouse","Gene_name",toupper(Proteins)))
      
    }
  }

  
  Seq  <<- dbGetQuery(db,Wrapper("*","Sequences","SeqID",unique(SeqID$SeqID)))
  if(MissCleave == "Include"){
    MissCleave <- 1
  }else{
    MissCleave <- 0
  }
  # Additional Filters
  Seq <- Seq[Seq$MissingCleavages <= MissCleave,]
  if(ExcludeKPRPfragments){
    Seq <- Seq[Seq$RPKP_Fragment_TrypticPeptid ==0,]
  }
  if(species == "mouse"){
    Seq$Unique_identifiers <- Seq$Unique_identifiers_Mouse
  }
  
  if(SILAC != "none"){
    Seq <- Seq[grepl("K$|R$",Seq$Sequence),]
  }
  
  Spec <- dbGetQuery(db,Wrapper("*","Spectra","Sequence",unique(Seq$Sequence)))
  Spec <- Spec[Spec$mz>=min(mzrange)&Spec$mz <= max(mzrange),]
  
  
  Spec$Proteins <- Seq$Unique_identifiers[SEL <- match(Spec$Sequence,Seq$Sequence)]
  Spec$GeneSymbol <- Seq$Gene_names[SEL]
  Spec$DwellTimes <- Seq$DwellTimes[SEL]
  Spec$RT_sd[is.na(Spec$RT_sd)] <- 0

  
  if(RescaledHydrophobicities){
    Spec$Hydrophobicity <- Spec$Rescaled_Hydrophobicity_Lo 
  }else{
    Spec$Hydrophobicity <- Seq$Hydrophobicity[SEL]
  }
  
  SpecFil <- FilterFun(Spec,Fragmentation,MassAnalyzer,Charge)
  
  SpecFil$Hydrophobicity_Rescaled <- NULL
  SpecFil$Hydrophobicity_Rescaled_Moving_Average <- NULL
  SpecFil$CollisionEnergy[is.na(SpecFil$CollisionEnergy)] <- 0
  SpecFil$CollisionEnergy[SpecFil$CollisionEnergy == 0] <- 20

  #---------- Check if All Proteins Are Listed
  Missing <- Proteins[is.na(match(Proteins,SpecFil$Proteins))]
  if(length(Missing) > 0){
    # cat("Warning, missing",Missing,"\nRelax filtering to include more spectra.")
  }
  dbDisconnect(db)
  return(SpecFil)
}
SpecFilFilSRM <- function(SpecFilFil,Transitions,convertToSRM = F,    totali = 3,totaliabove = 2,unmodifiedOnly = T,singlechargedonly = T,session = NULL,mzrange = c(300,1700),SILAC = "none"){
  
  if(convertToSRM){
    Transitions <- Transitions[!is.na(Transitions$Matches),]
    if(SILAC != "none"){
      # To get only Transitions with a native c terminus (to have the label!)
      Transitions <- Transitions[grep("^[xyz]",Transitions$Matches),]
    }
    if(unmodifiedOnly){
      Transitions <- Transitions[grep("-",Transitions$Matches,fixed = T,invert = T),]
      
    }
    if(singlechargedonly){
      Transitions <- Transitions[grep("+",Transitions$Matches,fixed = T,invert = T),]
    }
    Transitions <- Transitions[Transitions$Masses > min(mzrange) & Transitions$Masses<max(mzrange),]
    
    SpecIDTransitionCount <- table(unlist(Transitions$SpecID))
    ToLowTransitionCount <- SpecIDTransitionCount > totali
    Transitions <- Transitions[!is.na(match(Transitions$SpecID,names(SpecIDTransitionCount)[ToLowTransitionCount])),]

    ML    <- match(Transitions$SpecID,SpecFilFil$SpecID)
    MLd   <-  SpecFilFil[ML,]

    
    if(length(session) > 0){
      try(TansitionExtraction$close(),silent = T)
      TansitionExtraction <- Progress$new(session,min = 1,max = length(SpecFilFil$SpecID))
      on.exit(TansitionExtraction$close())
    }
    
    AllOut <- c()
    it <- 0
    for(i in SpecFilFil$SpecID){
      it = it+1
      if(length(session) > 0){
        try(TansitionExtraction$set(message = 'Assembling Transition List',value = it),silent = T)
      }
      SpecFilFilB.i <- SpecFilFil[i == SpecFilFil$SpecID,]
      tempi <- MLd[MLd$SpecID == i,]
      tempTransitions <- Transitions[Transitions$SpecID == i,]
      tempiAbove <- tempTransitions[tempTransitions$Masses > unique(unlist(SpecFilFilB.i$mz)),]
      
      if(dim(tempiAbove)[1]>=totaliabove ){
        SelectAbove <- tempiAbove[order(tempiAbove$Intensities,decreasing = T),][1:totaliabove,]
      }else{
        SelectAbove <- NULL
      }
      RestI <- totali-totaliabove
      if(length(RestI) > 0){
        tempibelow <- tempTransitions[tempTransitions$Masses < unique(unlist(SpecFilFilB.i$mz)),]
        SelectBelow <- tempibelow[order(tempibelow$Intensities,decreasing = T),][1:RestI,]
      }else{
        SelectBelow <- NULL
      }
      
      SelectAll <- rbind(SelectAbove,SelectBelow)
      SelectAll <- cbind(SpecFilFilB.i[rep(1,dim(SelectAll)[1]),],SelectAll)
      AllOut <- rbind(AllOut,SelectAll)
      
      
    }
    
    SpecFilFil <- AllOut
    SpecFilFil <- SpecFilFil[!is.na(SpecFilFil$Matches),]
    Below <- (tabname <- table(unlist(SpecFilFil$SpecID))) <totali
    if(any(Below)){
      SpecFilFil <- SpecFilFil[is.na(match(SpecFilFil$SpecID,names(tabname)[Below])),]
    }
    if(dim(SpecFilFil)[1] == 0){
      return(NULL)
    }
  }
  return(SpecFilFil)
}
OptimizeTable <- function(MD,SpecFilFil,RT_Template,Windowsize = NULL,UpperLimit,InputProteins = NULL,ProteouniqueSeparate = F,ProteouniqueSeparateType = "accession",OnlyIsoFormSpecific = F, scheduled = T,session = NULL,UseMovingAverage = F,SILAC = "none",Rate_model = "PEP"){
  
  try(TansitionExtraction$set(message = 'Predicting Retention Times',value = length(SpecFilFil$SpecID)/2),silent = T)
  
  if(length(MD) > 0 & scheduled){
    
    if(UseMovingAverage){
      

      RTCorrected <- sapply(SpecFilFil$Hydrophobicity,function(x){
        tx <- MD$Model$x-x
        median(MD$Model$fitted[abs(tx) == min(abs(tx))])
      })
      RTpred <- RTCorrected
      RTpred[RTpred <0] <- 0
      
    }else{
      
      md <<- MD$Model
      HyBa <<- SpecFilFil$Hydrophobicity
      SpecFilFil$Hydrophobicity[SpecFilFil$Hydrophobicity > max(md$x)] <- max(md$x)
      SpecFilFil$Hydrophobicity[SpecFilFil$Hydrophobicity < min(md$x)] <- min(md$x)
      

      RTpred <- predict(md,unlist(SpecFilFil$Hydrophobicity))
      RTpred[RTpred <0] <- 0
      RTpred[RTpred>max(md$y)] <- max(md$y)
    }

    MatchedRetention_Time <-  RT_Template$`Retention time`[match(unlist(SpecFilFil$Sequence),RT_Template$Sequence)]
    if(any(!is.na(MatchedRetention_Time))){
      RTpred[!is.na(MatchedRetention_Time)] <- MatchedRetention_Time[!is.na(MatchedRetention_Time)]
    }

    if(length(Windowsize) == 0){
      
      Windowsize <- max(Diffis<- apply(MD$BiquanTableAuto[,1:2],1,diff),na.rm = T)/2
      if(is.na(Windowsize)){
        Windowsize <- 5
      }
    }else{
      Windowsize <- Windowsize/2
    }
    
    try(TansitionExtraction$set(message = 'Predicting Retention Times',value = length(SpecFilFil$SpecID)/2),silent = T)
    
    RTpredPlus <- RTpred+Windowsize
    RTpredMinus <- RTpred-Windowsize
    SpecFilFil$PredictionRange <- as.numeric(!is.na(RTpred))
    OutOfRange <- NULL
    if(any(is.na(RTpred))){
      if( any(sel1 <- (unlist(SpecFilFil$Hydrophobicity) < min(md$x,na.rm = T)))){
        try(RTpredMinus[is.na(RTpred&sel1)] <- 0,silent = T)
        try(RTpredPlus[is.na(RTpred)&sel1] <- min(md$x,na.rm = T)+Windowsize,silent = T)
      }
      
      if(any(sel2 <- (unlist(SpecFilFil$Hydrophobicity) > max(md$x,na.rm = T)))){
        RTpredMinus[is.na(RTpred)&sel2] <- max(md$x,na.rm = T)-Windowsize
        RTpredPlus[is.na(RTpred)&sel2] <- max(md$x,na.rm = T)+Windowsize*2

      }
      if(any(c(sel1|sel2))){
        OutOfRange <-  range(unlist(SpecFilFil$Hydrophobicity)[sel1|sel2],unlist(SpecFilFil$Hydrophobicity)[sel1|sel2])
        OutOfRange <- unique(OutOfRange)
      }
      

    }
    
    
    try(TansitionExtraction$set(message = 'Predicting Retention Times',value = length(SpecFilFil$SpecID)),silent = T)
    
    
    SpecFilFil$Start <- RTpredMinus
    SpecFilFil$End <- RTpredPlus
    RTrange   <- range(RT_Template$`Retention time`)
    Scheduled <- T
    
  }else{
    OutOfRange <- NULL
    
    SpecFilFil$Start <- 0
    
    SpecFilFil$End <- 0
    RTrange <- c(0,0)
    Scheduled <- F
    SpecFilFil$PredictionRange <- 1
    
    
  }
  
  # Acquisition Time Calculator
  #--------------------------
  
  initTable <- c()
  Gradient <- seq(round(RTrange[1],1),round(RTrange[2],1),by = 0.1)
  template <- rep(NA,length(Gradient))
  SpecFilFil$ProteinColor <- as.numeric(as.factor(sapply(SpecFilFil$Proteins,unique)))
  SpecFilFil <- SpecFilFil[order(SpecFilFil$Start),]
  
  for(i in 1:dim(SpecFilFil)[1]){
    tempx <- SpecFilFil[i,]
    # tempx$Start[tempx$Start == ""] <- 0
    # tempx$End[tempx$End == ""] <- 1
    M <- match(c(round(unlist(tempx$Start),1),round(unlist(tempx$End),1)),round(Gradient,1))
    if(is.na(M[1])){M[1] <- 1}
    if(is.na(M[2])){M[2] <- length(Gradient)}
    tempxi <- template
    tempxi[M[1]:M[2]] <- tempx$ProteinColor
    initTable <- rbind(initTable,tempxi)
    
  }
  
  
  SpecFilFil$ID <- 1:dim(SpecFilFil)[1]
  
  # Apply Additional Filters on Protein Level
  SpecFilFil$ProteinsTemp <- SpecFilFil$Proteins
  
  Seli <- unique(unlist(strsplit(InputProteins," |;|,|\t")))
  Seli <- Seli[Seli!= ""]
  
  if(any(!is.na(match(tolower(Seli),tolower(SpecFilFil$GeneSymbol))))& length(intersect(tolower(Seli),unlist(strsplit(as.character(SpecFilFil$Proteins),";")))) == 0){
    Seli <- unique(SpecFilFil$Proteins)
  }
  
  
  Pinit <- unlist(SpecFilFil$ProteinsTemp)
  
  if(!ProteouniqueSeparate){
   
      SpecFilFil$ProteinsTemp <- sapply(SpecFilFil$Proteins,function(x){unique(unlist(strsplit(x,"-|;")))[1]})
  }
  if(OnlyIsoFormSpecific){
    SpecFilFil <- SpecFilFil[selProteo <- !grepl(";",SpecFilFil$Proteins),]
    initTable <- initTable[selProteo,]
    if(is.vector(initTable)){
      initTable <- t(as.matrix(initTable))
    }
  }
  
  
  
  SpecFilFilB <- SpecFilFil
  initTableB <- initTable
  ContinueOptimization <- T
  it <- 0
  leFu <- apply(initTable,2,function(x){length(x[!is.na(x)])})
  # stop()
  if(length(session) > 0){
    try(TansitionExtraction$close(),silent = T)
    ma <- max(leFu) -UpperLimit
    ma[ma< 0] <- 0
    optimization <- Progress$new(session,min = 1,max = ma)
    on.exit(optimization$close())
    
  }
  
  if(SILAC!= "none"){
    UpperLimit <- UpperLimit/2
  }
  
  while(any(leFu >UpperLimit)&ContinueOptimization ){
    leFu <- apply(initTableB,2,function(x){length(x[!is.na(x)])})
    
    it <- it+1
    if(length(session) >0){
      optimization$set(value = it,message = "Filtering Table")
    }
    cat("\roptimization round",it)
    
    ToMany <- Gradient[leFu >  UpperLimit & leFu == max(leFu,na.rm = T)]
    ToMany <- unique(round(ToMany,1))[1]
    
    if(length(ToMany) == 0){
      ContinueOptimization <- F
    }else{
    SpecFilFilSel <- SpecFilFilB[sel <- round(SpecFilFilB$Start,1) <= min(ToMany) & round(SpecFilFilB$End,1) >= min(ToMany), ]
    
    P <- table(unlist(SpecFilFilB$ProteinsTemp))
    if(any(P) > 2){
      P <- P[P > 2]
    }
    P <- names(P)[P == max(P)]
    sel <- !is.na(match(SpecFilFilSel$ProteinsTemp,P))
    
    if(!any(sel)&length(sel) > 0){
      P <- table(unlist(SpecFilFilSel$ProteinsTemp))
      if(length(P) != 0){
        P <- names(P)[P == max(P)]
      }else{
        P <- NULL
      }
      
      sel <- !is.na(match(SpecFilFilSel$ProteinsTemp,P))
    }
    
    SpecFilFilSel <- SpecFilFilSel[sel,]
    if(Rate_model == "Score"){
      sc <- unlist(SpecFilFilSel$Score)
      if(length(sc) ==0){
        Remove <- SpecFilFilSel[NULL,]
      }else{
        Remove <- SpecFilFilSel[SpecFilFilSel$Score == min(sc,na.rm = T),]
      }
    }

    if(Rate_model == "PEP"){
      sc <- unlist(SpecFilFilSel$PEP)
      if(length(sc) ==0){
        Remove <- SpecFilFilSel[NULL,]
      }else{
        Remove <- SpecFilFilSel[SpecFilFilSel$PEP == max(sc,na.rm = T),]
      }
    }
    
    
    
    if(length(Remove$ID) == 0){
      ContinueOptimization <- F
      
      
    }else{
      
      rem <- is.na(match(SpecFilFilB$ID,Remove$ID))
      initTableB <- as.matrix(initTableB[rem,])
      SpecFilFilB <- SpecFilFilB[rem,]
      
    }#
    }
  }
  if(is.vector(initTableB)){
    initTableB <- t(as.matrix(initTableB))
  }
  OutList <- list(SpecFilFilB = SpecFilFilB,initTableB = initTableB,SpecFilFil = SpecFilFil,initTable = initTable,OutOfRange = OutOfRange)

  return(OutList)
}
summaryTable <- function(SpecFilFilB,TypeName = "Feature"){
  Proteins <- unique(unlist(SpecFilFilB$Proteins))
  ProteinSummary <- sapply(Proteins,function(x){
    tempx <- SpecFilFilB[SpecFilFilB$Proteins == x,]
    GS <- paste(unique(tempx$GeneSymbol),collapse = ";")
    UniSeq <- length(unique(unlist(tempx$Sequence)))
    TotalSeq <- length(tempx$Sequence)
    AverageScore <- round(mean(unlist(tempx$Score),na.rm = T))
    Seq <- paste(unique(sort(unlist(as.character(tempx$Modified_sequence)))),collapse = ";\n")
    return(c(GS,TotalSeq,AverageScore,Seq))
    
  })
  ProteinSummary <- t(ProteinSummary)
  ProteinSummary <- cbind(rownames(ProteinSummary),ProteinSummary)
  colnames(ProteinSummary) <- c("Protein_Accession","Gene_Symbol",paste(TypeName,"Count",sep = "_"),"Average_Score","Peptides")
  rownames(ProteinSummary) <- NULL
  try(ProteinSummary <- data.frame(ProteinSummary,check.names = F,check.rows =F),silent = T)
  ProteinSummary <- ProteinSummary[order(ProteinSummary$"Gene_Symbol"),]
  return((ProteinSummary))
}

ExtractMSMS <- function(SpecFilFil,tab){
  
  TransiMap <- cbind(SpecFilFil[match(tab$SpecID,SpecFilFil$SpecID),],tab)
  return(TransiMap)
}
PlotFun <- function(x,type = "both",SpecID = 877977){
  
  msms <- x[c("Masses","Intensities","Matches")]#[sel,]
  msms <- msms[!is.na(msms$Intensities),]
  rawmsms <- x[c("Masses","Intensities_RAW","Matches")]#[sel,]
  rawmsms <- rawmsms[!is.na(rawmsms$Intensities_RAW),]
  colnames(rawmsms) <- colnames(msms)
  
  
  if(type == "Both"){
    
    msms$Intensities <- msms$Intensities*-1
    msms <- rbind(rawmsms,msms)
    
  }
  if(type == "Raw File Spectrum"){
    msms <- rawmsms
  }
  return(msms)
}

# ui function --------

triggerType <- "hover"
ui <- fluidPage(
  
  fluidRow(
    column(9,
           titlePanel("Picky - Your Magic PRM and SRM Designer Tool",tags$head(tags$link(rel = "icon", type = "image/png", href = "./favicon.png"),tags$link(rel = "shortcut icon", type = "image/png", href = "./favicon.png"),tags$title("Picky"))) # Title
           
           ),
    column(3,
           conditionalPanel("output.showexportbutton == true",
                            downloadButton("ExportButton", "Export inclusion list",width = '100%',style = "background-color: #0096ff;margin-top: 25px")
           )
    )
  ),
  fluidRow(
    
    column(9,
           checkboxInput("triggerType","Display Tooltips",value = T,width = '100%')
           
           ),
    column(3,
           conditionalPanel("output.showexportbutton == true",
                            uiOutput("ExportTypes",width = '100%')
           )
    )
  ),
  
  fluidRow(

    
  ),
  
  tabsetPanel(id = "tabs",
    #-------- Settings
    tabPanel("Settings",
             wellPanel(style = "background-color: #ffffff;",
                       fluidRow(
                         column(12,
                                # Query
                                wellPanel(
                                  

                                  fluidRow(
                                    column(9,mainPanel(h3("1. Database Query"), width = "100%"))
                                    
                                  ),
                                  wellPanel(
                                    
                                  
                                  fluidRow(
                                    column(8,
                                           radioButtons(inputId = "IDType","Type of ID",c("Gene symbol","Uniprot accession"),
                                                        selected = "Gene symbol",inline = T,width = '100%')
                                    ),
                                    
                                    conditionalPanel(condition = includeMouse,
                                                     column(4, radioButtons(inputId = "Species","Species",c("Human","Mouse"),inline = T,
                                                                            selected = "Human",width = '100%'))
                                    )
                                  ),
                                  fluidRow(
                                    
                                    column(12,
                                           textInput(inputId = "Proteins","List of proteins",value = "TP53 MYC BCL2 BRCA1 BRCA2 PIK3CA GRID2 NFKB1 APOE",width = '100%') # Protein Selecter
                                          
                                    )
                                  )
                                  
                                  
                                  )
                                  ,
                                  fluidRow(
                                    
                                    column(4,
                                           radioButtons(inputId = "SRM","Method type",c("SRM","PRM"),selected = "PRM",inline = T,width = '100%')  ,  
                                           radioButtons(inputId = "MissCleave","Miss-cleaved peptides",choices = c("Include","Exclude"),selected = c("Exclude"),inline = T)
                                           
                                    ),
                                    
                                    column(4,
                                           checkboxGroupInput(inputId = "MassAnalyzer","Detector types",choices = c(c("ITMS","FTMS")),selected = c("FTMS"),inline = T),
                                           checkboxGroupInput(inputId = "Charge","Allowed charge state",choices = c(c(".*","1","2","3","4","5")),selected = c(2:4),inline = T)
                                           
                                           
                                    ),
                                    column(4,
                                           checkboxGroupInput(inputId = "Fragmentation","Fragmentation types",choices = c("CID","HCD","ETHCD","ETCID"),selected = c("HCD","CID"),inline = T),
                                           
                                           sliderInput("mzrange","m/z range",0,2000,value = c(300,1700),step = 10)
                                           
                                            
                                    )
                                    
                                  )
                              
                                  
                                )
                                  
                                  
                                
                         ),
                         # RT estimation/RT Template
                         #--------------
                         column(12,
                                
                                wellPanel(
                                  mainPanel(h3("2. Retention Time Estimation (optional)"), width = "100%"),
                                  
                                  fluidRow(
                                    column(2,   radioButtons(inputId = "PRMType","Acquisition type",c("static","scheduled"),selected = "scheduled",inline = T,width = '100%')    ),
                                    
                                    column(6,
                                           fileInput(inputId = "RetInput","Retention time calibration file",accept = c(
                                             ".txt"),width = "100%"))
                                    ,
                                    column(2,
                                           checkboxInput("ExampleFile","Example file",width = "100%",value = T)
                                           ,style = "margin-top: 15px"
                                           
                                    ),
                                    column(2,
                                           downloadButton("ExampleFileDL", "Template",width = '100%',style = "margin-top: 25px")
                                           
                                           )
                                    
                                    
                                    
                                  ),
                                  fluidRow(textOutput("HydrophobicityRangeCheck"), style="color:#3C8CC5;margin-bottom: 25px",width = "100%"),
                                 
                                    fluidRow(
                                      column(4,
                                             checkboxInput("RetWinAuto",label = "Auto adjust retention time window",value = T),
                                             conditionalPanel(condition = "input.RetWinAuto",
                                                              fluidRow(
                                                                column(2),
                                                                column(6,sliderInput("QuantileAuto", "Quantile limits", 80, min = 10, max = 100, step = 1,
                                                                                      width = NULL)),
                                                                column(4,uiOutput("QuantileAutoResult"),style = "margin-top: 25px")
                                                              
                                                              
                                                              )
                                                              ),
                                             conditionalPanel(condition = "input.RetWinAuto==false",
                                                              sliderInput(inputId = "RetWin",# id for identifiying the input later on
                                                                          label = "Retention time window", # Title or Explanation
                                                                          value = 15,min = 1,max = 100,step = 1 # input specific arguments
                                                              )
                                             ),
                                             checkboxInput("AdvancedSettings",label = "Advanced curve adjustment") ,
                                             conditionalPanel(
                                               condition = "input.AdvancedSettings",
                                               wellPanel(
                                                 sliderInput(inputId = "loessWin",# id for identifiying the input later on
                                                             label = "# Windows for subsetting", # Title or Explanation
                                                             value = 30,min = 1,max = 100 # input specific arguments
                                                 ),
                                                 sliderInput(inputId = "loessSpan",# id for identifiying the input later on
                                                             label = "Span of loess", # Title or Explanation
                                                             value = 30,min = 1,max = 100 # input specific arguments
                                                 ),
                                                 sliderInput(inputId = "loessQuan",# id for identifiying the input later on
                                                             label = "Quantile limit", # Title or Explanation
                                                             value = 0.1,min = 0,max = 0.2,step = 0.01 # input specific arguments
                                                 )
                                               )
                                               
                                             )
                                             
                                      ),
                                      column(8,
                                             # mainPanel(h3("Hydrophobicity"),width = "100%"),
                                             plotOutput("ret", width = "100%",height = "400px")
                                      )
                                      
                                    )
                                    
                                    
                                    
                                  #)
                                )
                                
                                
                         )
                       ),
                       # Additional Settings
                       #--------------
                       wellPanel(
                         mainPanel(h3("3. Additional Settings"),width = "100%"),
                         
                         fluidRow(
                           column(6,
                                  mainPanel(h4("Modifications"),width = "100%"),
                                  
                                  selectInput("ModType",NULL,c("all","only modified peptides","only unmodified peptides"),selected = "only unmodified peptides"),
                                  
                                  checkboxInput("SeparateModifications","Consider modified and unmodified peptides separately")
                                  ,mainPanel(h4("Proteins"),width = "100%"),
                                  checkboxInput("ProteouniqueSeparate","Consider protein isoforms separately",value = T),
                                  checkboxInput("OnlyIsoFormSpecific","Allow only proteotypic peptides",value = F),
                                  mainPanel(h4("Stable isotope label"),width = "100%"),
                                  selectInput("SILAC",NULL,c("none","lysine and arginine"),selected = "none")),
                           
                           column(6,
                                  # mainPanel(h4("Resolution"),width = "100%"),
                                  fluidRow(column(9,uiOutput("UpperLimit")),column(3,actionButton("show", "Advanced"))),
                                  
                                  conditionalPanel(
                                    condition = "input.SRM == 'SRM'",
                                    # mainPanel(h4("Transitions"),width = "100%"),
                                    sliderInput("SRMtotalcount",label = "# of transitions",min = 3,max = 20,step = 1,value = 3),
                                    uiOutput("SRMabovecountInitiate"),
                                    fluidRow(mainPanel(h4("Dwell Time"),width = "100%")),
                                    fluidRow(
                                      column(8,checkboxInput("AbuBinDwell","Abundance weighted dwell times",value = T)),
                                      column(4,
                                        conditionalPanel(
                                          condition = "input.AbuBinDwell==false",
                                          fluidRow(
                                            column(8,numericInput("FixedDwellTime",label = NULL,value = 50)),
                                            column(3,mainPanel("ms",style = "margin-top: 5px",width = "100%"))
                                            
                                          )
                                            
                                        )
                                        
                                      )
                                    )
                                    
                                  )
                           )
                         )
                       )
             )
             ),
    # Table Panel
    #--------------
    tabPanel("Table",
             wellPanel(style = "background-color: #ffffff;",
             tableOutput("ExportTable")
             )
             ),
    # Spectra Panel
    #--------------
    tabPanel("Spectra",
             wellPanel(style = "background-color: #ffffff;",
             wellPanel(
             fluidRow(
               
               column(4,
                      selectInput("spectraPlotType",label = "MSMS Type",choices = c("MaxQuant Spectrum","Raw File Spectrum","Both"),selected = "Both")
                      ),
               column(4,
                      uiOutput("MSMS_pro")
                      ),
               column(4,
                      uiOutput("MSMS_seq")
                      )
             ))
             
             ,
             plotOutput("FragmentPlot",height = "100px"),
             
             plotOutput("MSMSplot")
             )
             
             
             )
    
  ),style='width: 100%; height: 100%',
  
  
  # Output
  
  
  wellPanel(style = "background-color: #ffffff;",
            # mainPanel(h3("Acquisition plots"), width = "100%"),
            
            fluidRow(
              
              column(6,
                     plotOutput("PRM_Schedule"),
                     tableOutput("ScheduledSummary")
              ),
              column(6,
                     plotOutput("PRM_ScheduleFil"),
                     tableOutput("ScheduledSummaryFil")
                    
                     
              )
            )
  )
  ,
  hr(),
  
  helpText("This application was developed by the" ,
           a("Selbach",     href="https://www.mdc-berlin.de/1151285/en/research/research_teams/intrazellul_re_signalwege_und_massenspektrometrie"),
           "lab. Please contact us (henrik.zauber(at)mdc-berlin.de) if you have feature requests or problems."
           ,"If you are using Picky, please cite us: TO BE ANNOUNCED.",
  "The Picky spectrum library is based on ProteomeTools (Zolg et al (2017). Nature Methods. http://doi.org/10.1038/nmeth.4153). At its current state, ",AllGenes,"proteins covered by",AllSeq,"sequences and",AllSpec,"spectra are implemented in the Picky database.","This application uses the",a("shiny,",href = "http://shiny.rstudio.com")," ",a("RSQLite,",href = "https://cran.r-project.org/web/packages/RSQLite/index.html")," ",a("data.table",href = "https://cran.r-project.org/web/packages/data.table/index.html"),"and",a("shinyBS",href ="https://cran.r-project.org/web/packages/shinyBS/index.html"),"R-packages.",
  "Copyright 2017",
  h5("Legal Disclosure"),
  "This homepage is presented by the MDC, MAX-DELBRÜCK-CENTRUM FÜR MOLEKULARE MEDIZIN in der Helmholtz-Gemeinschaft, a corporation under Public Law.",
  h5("Responsible for the content"),
  "Henrik Zauber",br(),
  "Matthias Selbach",
  # "...",br(),
  h6("Address"),
  "MAX-DELBRÜCK-CENTRUM FÜR MOLEKULARE MEDIZIN in der Helmholtz-Gemeinschaft",br(),
  "Robert-Rössle-Straße 10",br(),
  "13125 Berlin-Buch",br(),
  "Tel: +49 / (0)30 / 94 06 - 0",br(),
  "Fax: +49 / (0)30 / 94 94 161",br(),
  "http://www.mdc-berlin.de",br(),
  "VAT identification number: DE81 12 61 930",br(),
  "The MAX-DELBRÜCK-CENTRUM FÜR MOLEKULARE MEDIZIN in der Helmholtz-Gemeinschaft is a corporation under Public Law of the state of Berlin. It was founded as a foundation under Public Law in 1992.",
  h5("Disclaimer"),
  h6("Accountability for content"),
  "The content and code of this page bas been created with the utmost care.
  However, we cannot guarantee the functionality, contents' accuracy, completeness or topicality. 
  According to statutory provisions, we are furthermore responsible for our own content on these web pages. 
  In this context, please note that we are accordingly not obliged to monitor merely the transmitted or saved information of third parties, or investigate circumstances pointing to illegal activity. 
  Our obligations to remove or block the use of information under generally applicable laws remain unaffected by this as per $$ 8 to 10 of the Telemedia Act (TMG)",
  h6("Accountability for links"),
  "Responsibility for the content of external links (to web pages of third parties) lies solely with the operators of the linked pages.",
  "No violations were evident to us at the time of linking. Should any legal infringement become known to us, we will remove the respective link immediately.",
  h6("Copyright"),"Our web pages and their contents are subject to German copyright law. Unless expressly permitted by law ($ 44a et seq. of the copyright law), every form of utilizing, reproducing or processing works subject to copyright protection on our web pages requires the prior consent of the respective owner of the rights. Individual reproductions of a work are allowed only for private use, so must not serve either directly or indirectly for earnings. Unauthorized utilization of copyrighted works is punishable ($ 106 of the copyright law).",
  "Quelle:",a( "twigg.de",href="https://www.twigg.de/haftungsausschlussimpressumenglisch.htm"),
  h6("Data protection"),
  "A visit to our website can result in the storage on our server of information about the access (date, time, page accessed). This does not represent any analysis of personal data (e.g., name, address or e-mail address). If personal data are collected, this only occurs – to the extent possible – with the prior consent of the user of the website. Any forwarding of the data to third parties without the express consent of the user shall not take place.
  We would like to expressly point out that the transmission of data via the Internet (e.g., by e-mail) can offer security vulnerabilities. It is therefore impossible to safeguard the data completely against access by third parties. We cannot assume any liability for damages arising as a result of such security vulnerabilities.
  The use by third parties of all published contact details for the purpose of advertising is expressly excluded. We reserve the right to take legal steps in the case of the unsolicited sending of advertising information; e.g., by means of spam mail.",
  "source:",a("Mustervorlage.net",href="http://www.mustervorlage.net/disclaimer-muster#Englisch"),
  h6("Access-Data/Server-Logfiles"),
  "The provider (respectively webspace-provider) imposes data about every access to the offer. Access-Data are: Name of the called up website, file, date and time of the demand, amount of transferred data, notice about successful demand, type of browser plus version, operation-system of user, referrer URL (site that has been visited before), IP-address and inquiring provider.
  These log-data are collected only insofar as required for statistical analyses, technical purposes, security and improvement of the offer. Provider reserves his right to review the log-data subsequently, there is concrete evidence for unlawful use."


  ),
  
  #,
  # tableOutput("FilteredTable")
  # verbatimTextOutput("ret")
  
  ### Adding HoverInfo
  bsPopover(id = "", title = "",placement = "bottom", trigger = "hover")#,title = "Picky")
)

# Server --------

server <- function(input, output, session){
  

  
  
  options(shiny.maxRequestSize=30*1024^2) 
  startTime <- Sys.time()
  # Sets the value for SRMabovecountInitiate
  MinSet <- reactive({
    if(input$SRMtotalcount >=2){
      MinSet <- 2
    }else{
      MinSet <- 0
    }
    MinSet
  })
  # For SRM condition dependent slider
  output$SRMabovecountInitiate <- renderUI({
    
    
    sliderInput("SRMabovecount",label = "# of transitions > precursor m/z",min = 0,max = input$SRMtotalcount,value  = MinSet(),step = 1)
    
  })
  # Creating Protein List
  SpecFilFilPre <- reactive({
    invalidateLater(1000000,session)
    if(input$Proteins != ""){
      
 
      excludeString <- c(":","'","\"","&","\\$","\\|","\\\\","\\(","\\)")
      excludeString <- paste(excludeString,collapse = "|")
      ProteinsCor <- input$Proteins
      ProteinsCor <- gsub(excludeString,"",ProteinsCor)
      
      Seli <- unique(unlist(strsplit(ProteinsCor," |;|,|\t")))

      Seli <- Seli[Seli!= ""]
      if(length(Seli) > 100){
        Seli <- Seli[1:100]
        warning1 <- Progress$new(session)
        on.exit(warning1$close())
        
        warning1$set(message = 'Warning!',
                     detail = 'You exceeded the limit of 100 IDs. Only the first 100 IDs will be queried.')
      }
      if(length(Seli) == 0){
        return(NULL)
      }
      
      progress <- Progress$new(session)
      on.exit(progress$close())
      
      progress$set(message = 'Searching Database',
                   detail = 'This may take a while...')
      dbpath <- dbpath
      Seli <- Seli
      
      speciesControl <- (input$Species)
      if(length(speciesControl) == 0){
        speciesControl <- "human"
      }else{
        speciesControl <- tolower(speciesControl)
      }
      SpecFil <- NULL
      try(SpecFil <- ExtractInfoFromDB(dbpath,Seli,IDType = input$IDType,Fragmentation = input$Fragmentation,MassAnalyzer= input$MassAnalyzer,Charge = input$Charge,MissCleave = input$MissCleave,species = speciesControl,RescaledHydrophobicities = T,mzrange = input$mzrange,RescaledHPType = "HP_rescaled_RT",SILAC = input$SILAC))      
      if(dim(SpecFil)[1] == 0){
        return(NULL)
      } 
      Mods <- grepl("(",SpecFil$`Modified_sequence`,fixed = T) # Filter for Modified Sequences
      if(input$ModType == "all"){
      }
      if(input$ModType == "only modified peptides"){
        SpecFil <- SpecFil[Mods,]
      }
      if(input$ModType == "only unmodified peptides"){
        SpecFil <- SpecFil[!Mods,]
        
      }
      if(dim(SpecFil)[1] == 0){
        return(NULL)
      }
      
      # SpecFiIl <<- SpecFil
      if(input$SeparateModifications){
        InitSeq <- SpecFil$`Modified_sequence`
      }else{
        InitSeq <- SpecFil$`Sequence`
        
      }
      progress$set(message = 'Prefiltering Table',
                   detail = 'This may take a while...')
      SpecFilFil<- sapply(unique(InitSeq),function(x){
        tempx <- SpecFil[InitSeq == x,]
        tempx <- tempx[tempx$Score == max(tempx$Score),]
        return(tempx[1,])
      })
      SpecFilFil <- t(SpecFilFil)
      rownames(SpecFilFil) <- NULL
      try(SpecFilFil <- data.frame(SpecFilFil,stringsAsFactors = F),silent = T)
      return(SpecFilFil)
    }else{return(NULL)}
    
  })
  # SRM only Extract Transitions, called by SpecFilFilSRM
  Trans <- reactive({

    if(is.null(SpecFilFilPre())){
      return(NULL)
    }
    db <- dbConnect(SQLite(),dbname = dbpath)
    # SpecFilFil <- FilteredTable()$SpecFilFilB
    if(length(session) > 0){
      try(TansitionExtraction$close(),silent = T)
      TansitionExtraction <- Progress$new(session,min = 1,max = length(SpecFilFilPre()$SpecID))
      TansitionExtraction$set(message = 'Searching Transitions',
                              detail = 'This may take a while...')    
      on.exit(TansitionExtraction$close())
      
    }
    hui <- SpecFilFilPre()#$SpecID
    TransiPos <- strsplit(as.character(hui$SpecIDPosition)," ")

    rows <- unique(unlist(lapply(TransiPos,function(x){x[1]:x[2]})))

    Transitions <- dbGetQuery(db,Wrapper("*","Transitions","rowid",rows))
    
    Transitions <- Transitions[!is.na(Transitions$Intensities_RAW),]
    dbDisconnect(db)
    rm(db)

    return(Transitions)
  })
  
  # For Testing Purposesonly
  output$path  <- renderText({paste(dbpath,getwd(),sep = ";")})
  
  
  RetTab  <- reactive({
  
    
    
    inFile <- input$RetInput
    
   
    if (is.null(inFile)&!input$ExampleFile)
      return(NULL)

    if(input$ExampleFile&is.null(inFile)){
      RT_Template <- sapply(unlist(strsplit(evidence(),";")),function(x){unlist(strsplit(x," ",fixed = T))})
      RT_Template <- t(RT_Template)
      colnames(RT_Template) <- c("Sequence","Retention time","PH")
      RT_Template <- data.table(RT_Template)
      RT_Template$"Retention time" <- as.numeric(RT_Template$"Retention time")
    }else{
      Er <- class(try(RT_Template <- fread(inFile$datapath,sep = "\t"),silent = T))
      if(Er == "try-error"){
        RT_Template <- matrix()
      }
      if(dim(RT_Template)[2] == 2){
        NACHECK  <- RT_Template[!is.na(as.numeric(unlist(RT_Template[,2]))),]
        NACHECK2 <- tolower(unique(unlist(strsplit(as.character(unlist(NACHECK[,1])),""))))
        NACHECK2 <- NACHECK2[NACHECK2 != " "]
        NACHECK2 <- setdiff(NACHECK2,letters)
        if(length(NACHECK2) != 0 | dim(NACHECK)[1] < 50){
          RT_Template <- matrix()
          
        }
      }
      if(dim(RT_Template)[2] <2){
        progress <- Progress$new(session)
        on.exit(progress$close())
        
        progress$set(message = 'Warning',
                     detail = "Incorrect Format. Returned to Template File.")
        
        RT_Template <- sapply(unlist(strsplit(evidence(),";")),function(x){unlist(strsplit(x," ",fixed = T))})
        RT_Template <- t(RT_Template)
        colnames(RT_Template) <- c("Sequence","Retention time","PH")
        RT_Template <- data.table(RT_Template)
        RT_Template$"Retention time" <- as.numeric(RT_Template$"Retention time")
        
      }
    }
    if(!all(any(colnames(RT_Template) == "Sequence")&any(colnames(RT_Template) == "Retention time"))){
      
      progress <- Progress$new(session)
      on.exit(progress$close())
      
      progress$set(message = 'Warning',
                   detail = "Table contains no \'Sequence\' or \'Retention time\' column. First and second column will be interpreted as \'Sequence\' and \'Retention time\' column respectively.")
      colnames(RT_Template)[1:2] <- c("Sequence","Retention time")
    }
    RT_Template <- RT_Template[order(RT_Template$`Retention time`)]
    if(dim(RT_Template)[1] > 8000){
      RT_Template <- RT_Template[sample(1:dim(RT_Template)[1],8000),]
      
    }

    cat("\rCalculating Hydrophobicity of template sequences\n")
    uniSeq <- unique(RT_Template$Sequence)
    
    progress <- Progress$new(session, min=1, max=length(uniSeq))
    on.exit(progress$close())
    progress$set(message = 'Calculating Hydrophobicity',
                 detail = 'This may take a while...')
    progress$set(value = 0)
    
    
    if(!is.null(inFile)){
    
    it <<- 0 
    if(Use_parallel){
      cl <- makeCluster(getOption("cl.cores", 4))
      progress$set(value = length(uniSeq)/2)
      
      PH <- parSapply(cl,uniSeq,function(x){PeptideHydrophobicity(x)})
      progress$set(value = length(uniSeq))
      stopCluster(cl)
    }else{
      PH <- sapply(uniSeq,function(x){it <<- it+1;      progress$set(value = it); cat("\rCalculated",it,"from",length(uniSeq));PeptideHydrophobicity(x)})
      
    }
    
    PH <- PH[match(RT_Template$Sequence,uniSeq)]
    RT_Template$PH <- PH
    }
    return(RT_Template)
  })
  fileUploaded <- reactive({
    return(!is.null(RetTab()))
  })
  output$fileUploaded <- reactive({
    fileUploaded()
  })
  # Loess Model
  md <- reactive({
    if (is.null(input$RetInput)&is.null(RetTab()))
      return(NULL)
    ah <- RetTab()$`Retention time`
    bh <- as.numeric(RetTab()$PH)
    MD <<- ModelRT(ah,bh,windows = input$loessWin,span = input$loessSpan,QuantileLimit = c(0+input$loessQuan,1-input$loessQuan),AutoWindowQuantile=input$QuantileAuto)
    
    

 

    
    return(MD)
  }
  )
  # Retention Hydrophobicity Plot Output
  output$ret <- renderPlot({
    
    if (!is.null(input$RetInput)|!is.null(RetTab())){
      a <- RetTab()$PH
      b <- RetTab()$"Retention time"
      mdf <- md()$Model
      col = densCols(a,b)
      col <- md()$outlier
      par(mai = c(0.8,0.8,0.1,0.1),bg = "#f5f5f5")
      plot(b,a,col = col,frame = F,ylab = "Hydrophobicity",xlab = "Retention Time [min]",pch = 20,type = "n",bg = "grey")
      grid(col = "black",lwd = 1,lty = "dotted")
      points(b,a,col = col,ylab = "Hydrophobicity",xlab = "Retention Time [min]",pch = 20)

      points(mdf$fitted[order(mdf$x)],mdf$x[order(mdf$x)],type = "l",col =2)
    }else{return(NULL)}
    
  })
  
  # SRM Prefiltering, returns the initial Table if SRM == F
  SpecFilFil <- reactive({
    if(is.null(SpecFilFilPre())){
      return(NULL)
    }
    return(SpecFilFilSRM(SpecFilFil = SpecFilFilPre(),Transitions = Trans(),input$SRM == "SRM",input$SRMtotalcount,input$SRMabovecount,session = session,mzrange = input$mzrange))
    
  })
  
  # Applies FIlters to the Table
  FilteredTable <- reactive({
    
    if(input$RetWinAuto){
      win <- NULL
    }else{
      win <- input$RetWin
    }
    
    if(is.null(SpecFilFilPre())){
      return(NULL)
    }
    OptimizeTable(md(),SpecFilFil = SpecFilFil(),RT_Template = RetTab(),UpperLimit = input$UpperLimit,InputProteins = input$Proteins,ProteouniqueSeparate = input$ProteouniqueSeparate
                  ,OnlyIsoFormSpecific = input$OnlyIsoFormSpecific,Windowsize = win,scheduled = input$PRMType == "scheduled",UseMovingAverage = F,SILAC = input$SILAC,session = session)

  })
  

  output$HydrophobicityRangeCheck <- renderText({
    if(is.null(FilteredTable())){
      return(NULL)
    }
    
    if(!is.null(FilteredTable()$OutOfRange)){
      paste("Hydrophobicity is out of range (",paste(FilteredTable()$OutOfRange,collapse = " and "),") for one or more peptides. Maybe try a different retention time template.")
      
    }else{return(NULL)}
   
  })
  # Image Plot of aquistion lists
 
  output$ScheduledSummary <- renderTable({
    SpecFilFil <- FilteredTable()$SpecFilFil
    Entries <- table(unlist(SpecFilFil$Proteins))
    Pis <-unique(unlist(strsplit(input$Proteins," |;|,|\t")))
    Pis <- Pis[""!=Pis]
    
    Seli <- length(Pis)
    if(input$IDType == "Uniprot accession"){
      Matched <- match(tolower(Pis),tolower(unlist(strsplit(as.character(SpecFilFil$Proteins),";"))))
    }else{
      Matched <- match(tolower(Pis),tolower(unlist(strsplit(as.character(SpecFilFil$GeneSymbol),";"))))
    }
    
    Missing <- Pis[is.na(Matched)]
    if(length(Missing) == 0){
      Missing <- ""
    }
    
    EntriesSum <- c(length(Pis),length(Matched[!is.na(Matched)]),paste(Missing,collapse = "; "))
    
    names(EntriesSum)[1:3] <- c("Queried Proteins","Proteins before Filtering","Missing Proteins")
    EntriesSum <- t(as.matrix(EntriesSum))

    return(EntriesSum)
  })
  output$ScheduledSummaryFil <- renderTable({
    SpecFilFil  <-  FilteredTable()$SpecFilFilB
    SpecFilFil$Proteins <- tolower(SpecFilFil$Proteins)
    Entries <- table(unlist(SpecFilFil$Proteins))
    Pis <-unique(unlist(strsplit(input$Proteins," |;|,|\t")))
    Pis <- Pis[Pis!= ""]
    Seli <- length(Pis)
    Pis <- tolower(Pis)
    SpecFilFil$GeneSymbol <- tolower(SpecFilFil$GeneSymbol)
    if(input$IDType == "Uniprot accession"){
      Matched <- match(tolower(Pis),tolower(unlist(strsplit(as.character(SpecFilFil$Proteins),";"))))
    }else{
      Matched <- match(tolower(Pis),tolower(unlist(strsplit(as.character(SpecFilFil$GeneSymbol),";"))))
    }
    
    Missing <- Pis[is.na(Matched)]
    if(length(Missing) == 0){
      Missing <- ""
    }

    EntriesSum <- c(length(Pis),length(Matched[!is.na(Matched)]),paste(Missing,collapse = "; "))
    
    names(EntriesSum)[1:3] <- c("Queried Proteins","Proteins after Filtering","Missing/Excluded Proteins")
    EntriesSum <- t(as.matrix(EntriesSum))
    
    return(EntriesSum)
  })
  

  # Summary Of Type
  output$ExportTable <- renderTable({
    if(is.null(SpecFilFilPre())){
      return(NULL)
    }
    
    if(input$SRM == "SRM"){
      TypeName <- "Transition"
    }else{
      TypeName <- "Feature"
    }
    hui <-summaryTable(FilteredTable()$SpecFilFilB,TypeName = TypeName)
    rownames(hui) <- NULL
    try(hui <- data.frame(hui),silent = T)
    return(hui)
  })
  # Template File
  evidence <- reactive({
    "ENGVEPEVVLYLETPADAATLR 30.04 50.339880;YLNFNQLSQYTEK 24.577 28.210000;HDLAVEPPAPTVLQK 18.475 30.410000;LAVPTLIH 23.049 26.785990;IATDPFVGNLTFFR 30.318 48.960000;KPLVTTAGPR 8.1571 14.862000;AIQELAEDIEK 15.38 22.360000;FADSESLFR 20.098 24.500140;EDGTIDFDDGSK 14.684 17.970000;NGEDPVAIRK 8.4106 12.150000;GPFIDLHLLK 26.199 45.930000;LHGYLPSR 10.949 15.561700;NPDMPADELEK 14.572 15.930000;TFEVTEDESK 12.229 14.010000;EQEVYMGEIPLMTDNGTFVINGTER 28.23 55.437300;NLDITDTR 13.242 15.013020;LPIVVYTPDNVDVK 24.655 35.400000;MSEIAIVLDEAER 27.923 35.340000;SGITFSQELK 18.428 26.400000;AIPNYNVMGLAK 22.432 27.210000;AAGAELVGMEDLADQIK 27.595 35.290000;PLTLDQLQQQNGK 17.93 21.250000;KGDVIIGANQQAVK 11.258 22.782000;YETVTAEELLSPMADK 26.877 33.930000;NFDNMREDEGLADR 11.625 20.310000;VALQGNMDPSMLYAPPAR 24.174 34.110000;LAALVDPGGDAEK 16.322 17.915000;LIAGESDGLPGITIDR 24.476 34.385000;TFVHPTSNK 6.5559 11.272205;ADGINPEELLGNSSAAAPR 24.627 25.200000;DALLEVTR 17.731 25.163600;KYTIDLTER 14.682 20.065206;IAEDLGLVTAK 19.462 27.960000;GNVIDPLDMVDGISLPELLEK 33.847 66.185250;MVEEDPAHPR 7.6034 5.730000;GVTAIITMTESGR 22.529 26.840000;ITFNAPTVPVVNNVDVK 25.307 37.290000;DNMAFSPEER 14.935 18.105000;VVADDQAPAEQSLR 12.359 13.705000;LEPYFTEGR 17.544 19.907580;YSYEASLMALHDR 18.014 29.590000;VTAEDVLKPGTADILVPR 24.605 39.015000;DAEALYGLLK 27.138 34.900000;SDILNAVGITLDETTTR 29.338 41.880000;HQKPVPALNQPGGIVEK 12.4 21.630000;VGDNDNLSALAAILAGADK 31.212 37.040000;IGYPVIIK 20.834 27.613740;TGNTDSLALQNIGLETDSR 23.251 33.650000;TMLYAINGGVDEK 19.645 29.040000;LEVVVNER 13.752 18.527410;IAEAAWQVNESTENIGAR 19.663 26.760000;VGFFNPIASEK 23.85 29.040000;AVGVHQSAYK 6.5711 8.210000;EHEDTLAGIEATGVTQR 17.131 24.170000;EFGYAQVEVVNDSALVR 25.461 41.750000;ALAPATTK 11.203 7.960590;TAAQDGDMIFFGADNK 24.218 32.395000;NIDVLIADTAGR 23.065 29.790000;ALDAIIASVTESLK 31.591 37.640000;VVYDVFDSNEVLEGK 26.183 34.630000;NVNPQIVENYR 16.358 17.940000;WVAMGEEYK 17.086 20.272455;AAMSGMLSPELK 20.546 26.965000;AANGVFDDANVQNR 14.434 14.490000;HPLVMGNWK 15.381 27.642930;EIETRPGSIVR 11.306 19.630000;IAAGQPLSIK 15.311 19.935000;AGSGALTLGQPNSPGVPADFAK 22.428 29.082240;EAGAVDATFVTK 15.832 23.260000;VQLLGSGSILR 23.076 34.060000;VEVETPEENT 23.073 8.755000;MIPGFEDGIK 19.443 25.230000;LLSEQLGEGEIELR 23.552 36.490000;IVADAMGALR 18.036 23.855000;MTLQQQIIK 17.824 25.628820;LGLDTLGIETVER 26.726 36.270000;VMGFIGGTSDRPAPISDK 18.752 27.930000;IGVMFGNPETTTGGNALK 22.528 31.065000;FNQIGSLTETLAAIK 28.892 41.410000;VLLEAADK 12.71 17.652360;HVEPEQALR 8.6586 12.571160;LEGVQTQDGYDEMVK 17.445 21.310000;MGLVEPLK 18.847 24.350040;FEVGEGIEK 15.645 18.024825;LPIYLDYSATTPVDPR 26.53 36.400000;VIDQGTGETMGDMTK 14.037 16.080000;MYPEEPVYLPPR 21.536 25.130000;YSSTISDPDTNVK 12.153 11.990000;IVQVIGAVVDVEFPQDAVPR 29.914 49.780000;EEVVGAEMMR 16.375 24.405000;VKPNQQVTIIDSEGK 12.729 17.302000;GTLADILK 22.459 26.724500;GEVDDIDHLGNR 13.939 19.265000;AGMVAGVIVNR 17.478 26.695000;ASYTLTEPQLTAPLK 23.075 32.720000;VLAEQALAQPTTDELMTLVNK 29.332 47.559710;TESFAQLFEESLK 33.895 40.890000;GLIIPDLR 22.11 31.426120;PQGYQQAVTVK 11.351 12.050000;DLSDVTLGQFAGK 24.73 31.550000;DEDGNYLVDVILDEAANK 32.328 41.770000;PVNYHFAPR 12.436 17.144260;SHLMEELK 7.8232 22.259380;TLLGAAAR 12.311 18.418620;LYHEEEVTVYDPQDVEFK 21.398 34.360000;QVWAGAIDTVGDK 19.728 24.810000;MLHDMTGADSSVSK 10.663 15.810000;AYGLFSEESELAQTLR 28.76 39.560000;GLNIFNSK 19.375 21.635020;LIGSDDQSDIALLQIQNPSK 26.146 39.110000;SFVVIIPAR 24.054 35.499905;LPQPVAEQAHK 9.0988 8.750000;ITNTEGGSFNSHYGVK 11.945 16.790000;GWLLNEPNYR 22.682 30.570000;LNPQIDMGQR 14.705 15.820000;GNTGENLLALLEGR 30.299 36.650000;SMQDPIADMLTR 26.838 27.940000;TLWVPALK 24.719 32.750520;MEFPEPVISIAVEPK 29.107 40.680000;DLGSMDVNEVIEK 22.539 30.050000;PVTGEITYGLER 18.936 25.620000;TIDLGLER 18.609 24.681140;ANTGVTLEPINSQNAPK 16.771 18.420000;QVVVGLPDVR 21.022 25.135000;VLDLIAHISK 21.18 33.460000;GYYPSYVLNEWER 26.61 34.890000;AAILNTAR 10.449 16.025240;ELGMDPTIHDIR 17.76 29.510000;LDTTGLIDRK 14.246 21.750000;AFIEENALK 17.392 23.770390;LAVEAGLLAR 20.853 29.615000;MNFSHGSPEDHK 7.3995 6.240000;DWLQAIESVIR 34.269 42.150000;AALAGAVNTLMR 21.746 28.790000;LDEFETVGNTIR 20.459 28.950000;AQFTDAAIK 13.353 17.163720;LHQSPTGNDIR 7.0491 8.850000;AEGDGEGFIVIDDLVDTGGTAVAIR 31.665 52.879800;SVTGMVAR 6.3097 14.322440;TDYQAVVSGIAEGYK 24.063 26.430000;GESAPTTAFGAAVVGGDNGR 19.458 20.290000;LREEIAEQHR 6.9846 12.330000;AATALQLEQVLR 24.523 32.790000;LVTDEAELAGMPESALAAAK 26.254 36.160000;GALDCSGVK 27.12 10.625160;GPAAVNVTAI 21.248 22.205000;IENQEGLNNFDEILEASDGIMVAR 31.097 53.175520;ELAAATSSADEGASVAYK 14.486 18.985000;LHVDPENFR 14.791 17.343725;LDGVIPGWTEGLK 26.682 37.550000;SLVLGISGGQDSTLAGK 22.024 33.945000;MQLIGGTWYGGEMK 25.37 36.060000;MLQEAVDALLDNGR 28.174 31.660000;RPEMTYEK 7.285 11.153340;TIVVYGTR 13.577 19.927490;VGNFMDDSAITAK 19.645 24.940000;TCVADESHAGCEK 3.8829 2.705000;ATLEDLGQAK 14.244 17.020000;ESVLFGDSTLALR 26.688 44.965000;TMLFDAPLQMK 26.803 39.840000;EVAVDDGILGGLK 24.322 34.455000;GYDHLELNGK 12.479 18.390000;HLDALVADEDLSR 17.716 30.650000;VDDLAVDLVER 23.516 32.320000;NGTFVTTYLSPR 22.465 29.750000;ALQEASGFIR 17.281 23.640000;DASESDVEASLDR 15.16 14.700000;INPIEDMPYDENGTPVDIVLNPLGVPSR 31.645 52.871520;DIAVLEDAGVPYQLLESSR 30.781 49.245000;ISNEESISAMFEH 23.842 26.990000;VEGGQHLNVNVLR 15.434 25.980000;LDTVVVPTNR 15.428 21.050000;DNTPMFVK 13.021 22.164780;LPIDLSQLK 24.429 29.190000;TVTYDFER 15.29 18.579440;INGQAFDMNK 13.824 16.890000;VAFQAVIK 18.281 23.375660;SALDSQQGEPWQTIR 20.74 24.720000;LTSENPIDLVR 22.49 26.970000;AITDYGVTGR 13.758 13.260000;NLDNDDIEYK 13.502 16.770000;YNEEHGITPQGLNK 11.948 14.490000;GIVDIPVYPLK 27.392 37.565000;QVIDASHAEGK 7.1914 7.510000;QGTHQAMDDIR 8.2186 10.870000;AQGIEVLLDDR 23.923 29.140000;VLEVSNLR 15.817 20.774160;EIGVETGGSNVQFAVNPK 20.084 28.030000;TLQLLAAGAR 19.528 28.670000;TLTISDNGVGMTR 18.302 24.770000;TLLGTALR 18.551 26.365020;DDDTADILTAASR 20.274 23.310000;HETISEDELR 10.722 17.970000;HAVTEASPMVK 8.6034 15.075000;LADAEFFFNTDR 26.614 35.940000;ETIDEGLWDTLSAR 28.298 39.840000;TPGNIAAIRPMK 14.362 22.930000;FGNMSGQMR 7.7721 14.118230;EMLIADGIDPNELLNSLAAVK 33.071 61.999680;FGAIGIGSR 20.251 20.807605;FGVPTVTDIIK 27.361 32.885000;GAMEESGAVLIK 17.275 27.395000;VFDMGVDVR 20.383 24.033100;SAVNAIETR 10.232 13.227935;AIAGDSVSGLGTGGMSTK 16.618 19.435000;VEAVNMDER 10.636 12.410615;SDIQMIFQDPLASLNPR 31.24 50.280000;WNAIMTVLR 26.807 36.823185;LDAVVFTGGIGENAAMVR 26.706 44.525000;PMDFSGIIPALQTK 29.363 40.220000;QNVQGVDVR 9.5117 10.260285;AALQNTILK 16.214 23.536870;VQFGDYQANGMMAVAK 22.319 28.860000;IQISHISR 10.018 15.665760;DVFMGVDELQVGMR 27.011 43.920000;YAGGAEENAPLHK 8.3877 8.660000;AWGLEPESR 17.326 16.871820;QGTHQAMDDIR 8.2188 10.870000;HDSGDPVEWGEK 12.88 13.310000;GVLAVMPDMDK 22.655 28.840000;HIINNTQDSVVNVDK 12.714 22.220000;IVELEAPQLPR 22.466 27.380000;IIAVLEPR 19.498 24.506130;ANDIDVPAALIDSEIDVLR 31.805 52.620000;FGQGEAAPVVAPAPAPAPEVQTK 19.179 23.384780;SISEITQDSLVQGLGK 24.17 33.590000;LVGLTGIDDAAR 20.619 25.760000;HLNATIINEGDINTR 16.457 28.750000;IDAAIAATR 12.014 16.682085;EEPQSIEVHPR 9.8575 14.180000;GSSMASDAFFPFR 26.154 37.450000;EENGLLTYIGTR 24.303 32.830000;ARGEDIEPLR 10.855 16.380000;GEGDIDNAPWPVADEK 20.879 23.890000;VVGDVVFEDAAK 19.859 25.030000;VSGAHLIK 7.7682 16.971240;VLTGGVDANALHRPK 12.121 19.360000;EQEVLAMGADK 14.563 21.310000;RDDEVIVLTGK 15.599 29.090000;AGLIGFSK 18.225 23.479720;EAPLAIELDHDK 17.894 28.010000;DYFGAHTYK 12.825 17.387510;EGDWAVAEMR 19.191 25.190000;LDYQGIYTPDQVER 20.4 24.650000;NQIADLVGADPR 19.54 21.420000;GEILGGMAAVEQPEK 15.201 28.340000;DAPDYQYLAAR 17.948 20.650000;LAASAEFIER 17.38 24.615000;ETLVEYGFR 20.577 30.055970;QQLGLSPQDSLGTR 17.8 21.690000;TLYPDAETIK 16.284 21.370000;SYFIPPPQMK 20.943 26.790000;LPTPEEYQTYVAQVDK 23.627 22.750000;AFLEAENPQR 14.445 15.880000;GQFAAVPLNILGDK 28.173 40.570000;VLEASGLR 7.567 17.084760;ELVADDATNTVYISGIGK 23.099 36.785000;VEMNVQYR 8.1519 14.525830;APLVMVVDHQR 15.64 27.200000;NDAADLGAAGGMGGMGGMGGMM 24.912 34.510860;PLTEVEQK 8.6831 11.683100;NHFASEYIYNAYK 17.787 28.730000;AFVEYLNK 18.171 23.039830;SQAIEGLVK 15.095 20.963285;LLDNAAADLAAISGQK 26 29.190000;YNPDALMTDLPK 24.055 26.540000;FDVIVMDPPK 23.789 30.965000;AGTEVAAAHAGWR 11.202 16.720000;AFLPGSLVDVRPVR 24.241 38.480000;IVDLLTER 20.595 25.617680;FLFDEYVR 26.641 30.773380;EEFLADNPGIDAEDANVQQFNAQK 23.72 38.637920;LWIDSPAR 18.953 21.133640;VLTEAAVAGK 11.451 15.560000;LAELGPQGLITTLK 27.135 40.240000;TVEVGPLANMLVK 27.568 38.640000;ANQFLGMEPLPTFIANDVIK 31.855 64.120000;LQNADLVVWVGPEMEAFMQK 32.29 59.490000;FLGFEQTFK 25.225 32.235490;RLDEVYALYADPDADFDK 25.16 41.430000;SFFENDLR 19.951 25.551460;NGEYTYIINTTSGR 19.511 23.250000;FHTLSGGKPQVEGAEDYTDSDD 15.028 18.166680;SPHYIVMNDK 8.0977 20.880000;ELGLGTAR 7.8887 17.510460;AGEITAMAHTR 9.994 14.720000;AIWEQELTDMR 21.074 30.710000;FGPDAFTVQATR 20.8 23.360000;AWHSSSETIAK 8.2557 14.490000;ALEMIDMHGGDLFSEE 28.462 43.340000;ADVQGSVEAISDSLLK 28.255 33.175000;HPEIGDVR 8.5854 11.834460;ITGTPFEALVR 25.123 31.090000;SMVMTEAQK 9.7249 14.414995;GEQQAMQVATR 9.6698 11.590000;ESFLNPGFR 21.307 29.472170;NPDIVAGVAALK 21.752 29.530000;LIPAADISEQISTAGK 22.431 26.660000;DVTTGDTLCDPDAPIILER 26.044 41.220000;PGNVVLTPTILR 23.571 38.430000;TIANLLTA 26.256 27.495490;GIDTVLAELR 27.386 31.890000;IVGQADPVAWVSLQDIQGK 28.728 41.680000;DNTPMFVK 13.228 22.164780;GVNVLADAVK 18.48 22.840000;NSETLENFSEK 14.223 18.150000;NYFDLVAR 22.929 26.478540;LGAVPGGTER 9.5374 9.245000;RTDAQNTAAYIGNGDR 10.415 11.810000;IETLVPDFR 23.706 29.900290;KAGNVAADGVIK 8.7591 15.352000;GLIEEMASAYEDPK 26.794 25.220000;EQFPQEQAER 10.047 10.310000;VTITIAADSIETAVK 25.23 33.990000;ELLPLIPAEK 26.107 37.010000;ETDIVTLPIQHDSAAELAQPLDVK 26.936 53.232160;EVDAIMNR 8.1501 18.144280;TALAIDAIINQR 25.523 35.020000;IGVLGLNGAGK 19.052 24.165000;HNGAPAFSPDGSK 8.8823 7.230000;ELLDQIAELLR 32.98 45.410000;PAPTPQAPAQNTTPVTK 11.425 6.650000;ADEQILDIGDASAQELAEILK 32.71 51.272000;FAEVVFER 20.04 26.090680;IYPGQVLR 16.339 19.184880;EMIISEVR 17.461 25.381180;HLNSGGELR 2.9675 12.405750;GLENLSGDLYEK 20.242 27.070000;TQPLYDTQVSVSHTTR 14.242 21.420000;ANPEQLEEQREETR 9.3342 6.070000;LIAETAPDANNLLR 21.178 30.385000;EGGIDAYPVDPK 15.785 18.190000;MEPFFPSGLK 25.177 32.630000;PNIPATGEFK 15.588 18.880000;AALELAEQR 13.939 17.504270;LDVLDGLK 22.836 26.417050;GLIEEMASAYEDPK 17.6 25.220000;AFFANPVLTGAVDK 25.939 36.580000;AQAGGVLTR 8.6369 11.787895;DDYNPETEQYTLTISQR 21.042 26.910000;ANDAAGDGTTTATVLAQAIITEGLK 31.518 38.148600;LPQVASPLTFATEEEIR 28.264 40.050000;DLDEIITIAGQELNEK 31.805 41.950000;SIVVPAPEAQNDPR 15.992 15.965000;FTDMIDGQTITR 16.477 28.810000;DVVIGETITVGELANK 26.059 39.895000;EVDALGGLMAK 21.611 29.580000;KFEELVQTR 12.912 21.544166;YGVVEFDK 17.024 19.643690;EAVFQYALYEPGLLEPEGVSK 29.8 55.447710;AFSIGGGQADALMLEK 25.25 37.580000;LATLLSDASR 17.83 23.540000;LAPVVVDVPDDVLVLR 30.223 53.290000;GQGIVLNEPSVVAIR 23.448 38.970000;VAAPESLAVLFNPAR 29.182 37.285000;AEIEGEIGDSHMGLAAR 14.149 28.510000;VEEINDDNHAEGLIGYVR 20.62 32.980000;PLGIGVLTTAEK 22.948 30.650000;DSLGETAFNMLLDR 31.098 45.530000;VDGVFTADPAK 15.509 18.820000;AGDNAPMAYIELVDR 27.614 31.920000;QNQLEQLAEQYDELK 27.231 27.970000;PLLPEVLFDNGQAR 29.246 39.950000;LEFNNDNRK 7.7682 9.934330;ADIDYNTSEAHTTYGVIGVK 18.642 30.550000;EHPVLLNR 10.084 23.290520;HAVTEASPMVK 6.9421 15.075000;GDIDANAFQHK 11.651 16.980000;AVGDSLEAQQYGIAFPK 24.35 31.710000;AQTFTLVAK 17.006 23.877420;TFTAKPETVK 8.4622 14.710000;TVHSLTQALAK 12.968 22.790000;SLNFLDFEQPIAELEAK 34.803 55.070000;SVDPNTASPYASYLQYGHIAGIDEILEGK 29.996 47.667960;IVIAGEITEADK 18.366 23.930000;TGWLDTVAVR 23.479 32.500000;SLIGPDGEQYK 14.971 16.020000;TGFYMSLIGTPDEQR 27.067 36.850000;VPYIAQVMNDAPAVASTDYMK 24.117 39.065320;VNAGHGLTYHNVK 7.9563 14.815000;QNLAQVER 9.0518 11.512820;KLVAASQAALGL 20.619 32.177000;VPYIAQVMNDAPAVASTDYMK 27.145 39.065320;DASGTINVDIDHK 14.04 19.100000;EIGNGFSELNDAEDQAQR 19.826 22.430000;VAVVLGESEVANGTAVVK 21.173 33.785000;AGYAEDEVVAVSK 15.503 17.620000;FTNTSGFANK 10.483 15.010000;SEQNNTEMTFQIQR 21.52 21.990000;LLESAGIAYTVNQR 19.784 28.090000;LIGNPETR 9.6081 9.374860;AFALAGLR 21.286 26.161630;ALTEANGDIELAIENMR 27.856 37.140000;AVDLISNVAGDR 19.742 23.010000;IPYVESFPTGTPQSPYGK 23.361 24.670000;LDAVVFTGGIGENAAMVR 24.777 44.525000;QEDANFSNNAMAEAFK 15.64 24.710000;EGLIDPNHSVQIGIR 20.84 35.790000;LLNDTDMAIIDK 17.556 32.890000;QIQEINTGILIANGADMK 25.255 41.310000;ATALILDK 15.744 26.104870;GEMPSDFDAK 9.4811 16.065000;DLGSLIFIDMR 31.875 50.350000;VAVTELAHIVDETLAANNLDR 31.079 46.623010;FVVEGDLR 16.826 21.734350;VVMTADAVK 11.671 16.156665;AVMDDFAAFVEK 23.98 34.185000;IAAHAADLAK 7.9627 14.535000;ITLNMGVGEAIADK 23.314 31.090000;ATVILAHTIK 13.88 30.795000;DVNDLPELLK 25.789 32.720000;HNVPAGSESHFK 6.5181 10.805000;TTLTAAITTVLAK 28.482 37.850000;AVASVLMSK 15.29 21.002205;AVAAALEQMPR 17.627 19.685000;SYAIDPITLPSAR 24.057 30.765000;GASEAFISTLR 20.584 28.720000;DPEFQNIMK 19.847 24.723930;AQPDWSIALLR 27.807 36.690000;MLDLTPAEFR 24.895 31.260000;MGHAGAIIAGGK 7.2823 17.190000;FDGTVEVK 12.115 15.977940;IDDDLTLLSETLEEVLR 36.155 57.670000;KLEHAVPMAK 7.6455 17.402000;ELGSVNDFIVDGNR 23.048 32.410000;AWHSSSETIAK 8.2539 14.490000;IAQFNVVSEAHNEGTIVSVSDGVIR 25.184 45.067800;DAIPTQSVLTITSNVVYGK 29.501 43.650000;SDLFNVNAGIVK 23.628 36.430000;VSTGSLIMVFEVAGEAGAAAPAAK 30.191 45.632960;ISLDEAMASGHVK 17.15 21.490000;EGDFLLLQK 23.72 38.521070;NEEIAQAAAAVK 14.602 17.390000;THGQPATPSTIGK 7.4606 7.330000;LLEIPFNPVYR 28.751 39.590000;DVDVLYFR 25.127 32.939720;EFLDWMPK 27.931 33.724900;DSVGAVVMGPYADLAEGMK 26.646 41.005000;AISTIAESK 9.5982 12.804680;DAHFIGLR 16.138 27.386700;FTGNVIGDIVDAQYK 27.676 34.110000;DGVGLLPTVLDVVENPK 32.149 48.305000;DDDPSFDELVALAVETGR 32.575 43.110000;NHIATLQER 7.8644 16.618840;SLNLVSEQLLAANGLK 28.048 47.270000;MSTGLALDSEGK 15.524 19.440000;ELVSGLTQSATGK 16.098 20.285000;VTVPVLWDK 24.592 31.442495;THGQPATPSTIGK 7.4674 7.330000;GGLGNLMK 16.154 20.764700;IATSYPHLLK 16.078 25.560000;FGAPSYNR 11.441 8.878210;FINAESAK 14.183 11.493900;TIGFPTANVPLR 23.55 32.690000;FHNVTNGITPR 11.452 17.190000;LVIESGDSAQSR 11.326 12.410000;EGDPNLGVIAETLTEHGTR 26.548 32.290000;QQAQVEQVLK 13.386 15.465000;EIQLSPTLIGGSK 22.809 31.630000;DNVVIYEGELESLR 26.347 41.505000;LYDQMLEPK 13.426 21.318430;TKPHVNVGTIGHVDHGK 7.8771 15.612000;QMLMSWVER 24.488 32.070080;SNILLIGPTGSGK 21.401 32.000000;STMDHYAASNPLNK 13.387 16.825000;NNGSEVQSLDPHK 9.2908 8.750000;VTAIPSDFITVASTESCPFAIMANEEK 31.222 52.960930;YLGLPSEEAFK 23.476 28.110000;NNQHDVAIVR 8.317 16.050000;YNAVIDQMLEEGTAYK 27.665 33.165000;LNDSNLFR 16.798 21.256620;QAGTPIDVVWPK 24.308 27.140000;IFLDASSEER 16.908 20.550000;DEDGNYLVDVILDEAANK 32.329 41.770000;MFIVGKPTIMGER 20.883 35.250000;MEQSVANLVDMR 22.226 26.580000;DNISLDLGNNAEAVILR 28.652 49.080000;GFSVNFER 19.108 21.389060;LNAPVDEQGR 9.6937 8.045000;EHLSQEVLGK 11.805 21.770000;TQVVVLGAGPAGYSAAFR 24.904 39.245000;EMLQNSPMALR 18.89 28.880000;VMSLLEPTKK 15.368 22.730000;ISGAGIQESHVHDVTITK 12.789 25.490000;ELVHNVALR 12.833 26.159105;QTWCTIERK 18.71 16.171260;DIIAILGMDELSEEDK 31.402 47.520000;SLADIGEALK 21.839 25.145000;AYEDAETVTGVINGK 18.724 18.460000;VHGGILGR 8.7512 15.627920;NAIDWLSR 22.491 25.892020;LAEDDVNLPAQLK 20.594 26.840000;LADIASGQPLPVAGK 19.377 25.840000;QTVYAFLGDGEMDEPESK 24.813 31.145000;RAVGLPGDK 8.2134 12.507915;MNQNLLVTK 14.757 23.390920;ALSELEQIVTR 24.295 28.440000;GVEPIYETMPGWSESTFGVK 28.916 45.940000;AGLNEINLPELQAGSSIMPAK 28.255 44.981320;ASSGLNEDEIQK 10.788 12.520000;ELADQLVPYAK 21.053 27.585000;NGVMIVNTSR 14.227 22.125000;INVTAYQGK 12.03 10.863545;SNQMTGLFSTIDEK 20.286 31.350000;IVDAQGNDVLIPGTDMPAQYFLPGK 28.425 51.782400;TQDLLYDAR 19.054 23.322810;TAHQLVLSK 8.9074 21.863310;QINLPTEVSEAIYNR 26.961 30.310000;DAPYFGFEIPNAPGK 28.313 33.850000;AAAAPVTGPLADDPIQETITFDDFAK 30.14 47.874740;IAGINIPDHK 15.433 16.960000;EISNPENLMLSEELR 25.663 39.330000;DDVIVSPQTVQVK 18.803 28.985000;QQLPDDATLR 13.955 16.690000;GPIEFSNQELDDLIIR 29.842 51.580000;SLTPGQEFVK 16.3 21.570000;YHVNFMGGDLGK 18.845 26.645000;LYFCLFSKR 8.6527 33.675530;LYAIQPEETLTLDVK 26.537 40.185000;AGVDVLGISTDKPEK 18.275 22.195000;EMNIADYDAELWQAMEQEK 29.349 43.080000;WKEGEATLAPSLDLVGK 23.193 39.782000;FQVFGADAMR 22.09 26.805000;APVVVPAGVDVK 18.222 23.675000;THAPVDFDTAVASTITSHDAGYINK 23.293 36.739650;MKPNIHPEYR 8.7236 11.802000;NVLQPIGWDAFGLPAEGAAVK 30.451 55.156840;VVADIAGVPAQINIAEVR 27.138 40.205000;RVEIDPSLLEDDK 20.079 29.800000;AVVVTSGTTSEVLLNK 19.618 33.485000;ATIVNAGSSDFLEGEQVEYSR 23.24 35.959420;YAEIADHLGLSAPGDR 19.639 26.760000;LHEEAMALPSEEEFAER 20.317 32.150000;FMASEIIR 19.713 26.085950;TGSLIMIFEVEGAAPAAAPAK 29.232 51.518500;QGGTLQLFR 20.666 26.241810;DIDGEVTTLEK 18.225 22.470000;TIIGFGSPNK 17.93 22.140000;EEYPQSAAIDLR 18.67 24.630000;AHPADVLPQEMDK 14.354 17.750000;LFIDNFDK 23.074 27.604280;DLLAPWVPDAPSR 28.869 35.650000;LDFFDDEIDSLR 29.542 42.350000;AGLGLGIK 15.119 22.533720;TNPLLTPFELPPFSK 31.374 50.900000;FGGNLTAAEEK 12.916 15.610000;LTNIDNTWNQR 16.044 21.070000;LLEDPVEVSANPSTR 19.07 20.790000;DGIASAILPGR 22.12 28.280000;LIPITSQNR 13.978 16.988580;MFTGSIVAIVTPMDEK 29.008 43.000000;AGPLAGYPVVDMGIR 22.52 35.970000;MIQEQTMLNVADNSGAR 20.058 24.880000;QADGTMEPLPK 13.033 13.340000;GEVLDIVAEPR 21.352 29.665000;LEEGLNPNAILVPQQGVTR 24.642 39.610000;TGNAVILR 11.94 22.940500;ALLMHGTEGEVYANPQR 15.823 26.640000;QALTYAVNK 12.096 15.801520;LIVAVNSDASTK 14.554 18.985000;GEPLPLLSQVAER 26.293 34.940000;VTTVVIDDDVAVK 19.207 30.640000;TPETHPLTWR 13.676 22.930000;SIGTLSAFEQNALEGMLDTLK 33.443 57.375340;NVALEEQAVEAVLAK 27.856 34.215000;MQDYTLEADEGR 14.815 14.860000;IQPSGTFYDYEAK 18.216 20.360000;MAVAESGMGIVEDK 16.19 23.685000;DGNLNLTGQAPYK 16.963 22.230000;VIVTDMDGTFLNDAK 24.572 36.055000;AFTDFFAR 25.365 28.171880;AGQTFTFTTDK 16.537 19.720000;IHPGEAMDK 7.1982 8.095360;AMEIVYDEADLRR 18.875 29.310000;EGSGMVEFSR 14.896 20.890000;MVGGVTPGK 8.9094 8.299690;AGLVDILGASGAENVQGEVQQK 26.477 33.359040;QFDYLVNSMR 22.673 29.530000;YVSPDEFDEMK 19.259 19.580000;EDLAVLAK 16.427 26.365020;EPAQSLAESAAK 10.869 12.345000;VTIDFTGSVDGEEFEGGK 24.333 31.290000;LMEEEGILAGISSGAAVAAALK 31.164 48.949920;FVMPGAALNEDEFR 25.846 34.875000;LATQQSHIPAK 7.3938 10.540000;KGDIVGPIR 10.906 19.442486;DSWLALLDEAGMK 31.691 47.280000;LESLVEDLVNR 27.158 33.010000;IFSFTALTVVGDGNGR 29.161 39.450000;IVAQAQAEIEAERK 12.886 15.355000;DAGVMAASWAQYQAQDALIK 28.228 45.000000;SNPLGEDFDYR 19.838 22.500000;FICIADASK 16.829 19.849200;NAGDTASIPTIEAILNEEK 29.414 33.720000;APVALLDFAASAR 27.951 36.275000;DQLLENLQEGMEVK 26.424 36.850000;MNSLIAQYPLVK 24.334 37.640000;TPASFEPSIDYVVTK 23.706 33.805000;ILVIAADER 17.289 24.894205;SVGTQVDDGTLEVR 15.903 20.740000;VAEPGEVLAAGGR 15.83 15.910000;QEAAPAAAPAPAAGVK 14.582 10.185000;LVHTLNGSGLAVGR 14.821 27.910000;LTALGDELR 17.576 23.882285;YLSLLPYTDR 25.635 31.810000;EGVNSTESGLQFR 16.278 22.465000;ADAFAVIVK 20.983 29.360275;GVINILSAEK 22.037 27.890000;GITINTSHVEYDTPTR 15.952 22.190000;IYGDTIPGHQADK 11.065 11.530000;TYDLEVWIPAQNTYR 28.294 41.790000;HLAEQVQPLR 11.333 20.425000;DGLPWSVK 21.977 24.718980;TTDVTGTIELPEGVEMVMPGDNIK 27.916 48.946400;GTEITLHLR 16.004 27.681850;VENFSLWSLK 9.2548 40.480000;VPGYLEEEGANK 13.77 13.520000;FIAATNVNDTVPR 17.033 22.525000;YGSGGLVQGKK 7.2199 9.590000;PAVLEVITPTQVR 23.579 33.875000;SIGVGQYQHDVSQTQLAR 16.156 23.690000;VDLAQLVK 22.106 24.520320;QGITSYSLVPFPGH 27.195 33.820000;IAEQEGIAEDGYR 13.397 14.560000;GEILGGMAAVEQPEK 20.356 28.340000;ILELAGFLDSYIPEPER 32.347 52.810000;LANAVGIGAVK 15.396 21.640000;SNNTEPTWFPDSQNLAFTSDQAGRPQVYK 25.424 38.150100;DAGIEASQIGYVNAHGTSTPAGDK 16.698 20.956800;GQIGAYEAALMNTK 21.911 26.720000;MVVTLIHPIAMDDGLR 24.571 50.705000;LEIEAIAVR 21.111 28.275380;SMQDPIADMLTR 30.739 27.940000;IEAQLNDVIADLDAVR 29.914 41.605000;LSYDTEASIAK 14.505 17.170000;AVDALAGQSAK 9.2772 11.410000;MQVSVETTQGLGR 16.358 18.135000;AMAVIDPVK 13.23 21.099505;AENMLWHK 13.252 22.666160;EGSSLLGSDAGELAGAGK 19.845 26.190000;PTSLPFEK 16.934 20.083580;VAGQSVIFEGTEGLPAQISR 25.393 39.310000;LAAGANALDTSR 12.143 13.915000;VLNNGDLGENK 9.7039 12.360000;ALESFQGTGR 12.365 14.240000;APVVVPAGVDVK 18.219 23.675000;VTIPGSDNEYYK 16.622 15.090000;VFQTHSPVVDSISVK 17.885 28.700000;LMTEFNYNSVMQVPR 22.174 36.160000;DYYAIMGVKPTDDLK 22.745 34.970000;TTDVTGTIELPEGVEMVMPGDNIK 24.884 48.946400;FDSTDTQLR 11.829 15.752870;ITVPVDATEEQVR 17.791 19.965000;ANGLVSDVFEAR 22.996 27.720000;GYVVTNNHVVDNATVIK 16.263 31.965000;TMGKPPVDDDIPAEEQVR 14.175 20.640000;MNIAGASFANIEGVK 24.387 31.990000;AGGNYLSSLLVGSEAR 27.147 32.520000;GQYGHVVIDMYPLEPGSNPK 23.15 35.270000;FGGESVLAGSIIVR 25.525 40.510000;TVVADGVGQGYK 12.19 14.815000;AQQLTLTAPEMTALVGGMR 29.076 43.940000;AVIPVAGLGTR 20.247 24.860000;VILVGNLGQDPEVR 22.539 32.580000;DMGGLPDALFVIDADHEHIAIK 28.14 61.061040;IEVPEIGEEVIEIK 28.445 38.905000;GQHLLAMK 6.5032 22.060720;LINYLVEEFKK 23.481 36.810000;YQGLVAQIELLK 27.993 41.610000;DEDGAVVTFTGK 17.355 23.170000;LVDLPGYGYAEVPEEMK 26.202 35.460000;LLTTCNIPVPSDVR 22.223 28.690000;QYGEAFEK 11.285 10.983060;SNQPVVFVDSTGR 17.37 22.950000;NTLTPTQMR 11.692 15.713950;SQTVHFQGNPVTVANSIPQAGSK 18.976 25.645660;VAAEMATGAIER 14.101 17.885000;VIDGADGAR 5.8587 7.083440;DGDGFAWIER 24.852 31.930000;ASDFVLAMGQGR 21.554 28.120000;HPDTLLAAR 10.776 20.734630;NTSPEIAEAIFEVAGYDEK 30.413 37.450000;SSPLDVSVYEPK 21.009 22.600000;RAEITPANADTVTR 10.532 16.280000;EKPAVEVR 6.5918 11.306592;DLGMVGAEELR 17.801 29.250000;ATGIDFDVR 18.741 22.009260;DLTGDPWGGR 18.641 18.250000;NYLGIDGIPEFGR 27.8 36.390000;DIHGAPVGDTLTLAR 19.321 32.020000;GFEGGQMPLYR 14.785 26.110000;VVITATQMLDSMIK 29.786 41.580000;YALVGDVGGTNAR 16.623 17.860000;LALMAMEQAR 20.632 26.440000;IVTGVTASQALLDEAVR 26.364 37.080000;SPIAGMPVLEVWK 29.393 44.980000;GVALVDGK 11.129 16.096190;LGELLEALK 25.285 33.150110;DLHDDAEWMAK 17.466 25.300000;SGAFVPGIRPGEQTAK 14.807 21.825000;EIGPATLK 11.165 17.056380;SLDVMGSPIR 19.823 25.270000;YTAGSGDPLEHESVK 12.117 12.665000;QFPNIDNAYMELGTNR 26.453 31.580000;ATVLGHIQR 10.449 19.455135;YLLLPNNPNMVIVDTK 28.39 48.110000;TVYSTENPDLLVLEFR 29.483 51.240000;ALVIVESPAK 17.146 24.215000;VLTQEMVK 12.648 17.936160;LVPGYEAPVMLAYSAR 27.232 38.010000;AAFNQMVQGHK 7.3994 14.490000;QDDAMVAFQNFIK 29.146 38.750000;LAANWPLEQDELLTR 27.706 43.515000;IVEENTYGIVK 16.095 22.980000;LLSENGYDPVYGAR 18.927 22.090000;NQIFLINER 21.28 33.296060;VMVTSHLGR 8.678 17.226965;SIFILPPSK 23.317 32.877670;AIISDVNASDEDR 13.876 15.010000;IVNPIMGVK 19.544 23.819040;LAFEPGVIPELLK 30.375 46.640000;LGESLWAPR 16.34 25.171510;LSSTDPNLQNIPVR 18.682 24.270000;IDMASPFIVMR 24.532 40.645000;DPALSELYLVEGDSAGGSAK 26.776 36.985000;SAFDEFSTPAAR 19.259 22.020000;EFGYAQVEVVNDSALVR 25.45 41.750000;VANLAEAQPDREIEK 12.984 17.410000;EMNIADYDAELWQAMEQEK 30.407 43.080000;VLYEMDGVPEELAR 24.165 32.860000;VQSMPEINDADK 13.985 12.560000;MGSTGIGNGIAIPHGK 16.435 21.740000;EDEIVIDR 14.393 23.432420;DGSMIGLNK 14.613 22.116290;DYYAIMGVK 22.117 26.922910;LLADDIVPSR 19.787 24.665000;EHQATVCK 17.197 4.890820;VGTPAITR 6.9283 11.389840;TMNTPHGDAITVFDLR 21.935 36.940000;VDSSGFQTEPVAADGK 14.011 15.420000;TPEGQAPEEIIMDQHEEIEAVEPEASAEQVDPR 24.663 33.317140;FSPDMTPEDPIVMESIK 26.962 38.760000;GNYSMGVR 5.4824 12.061500;VYFAADEQTLLK 22.965 34.080000;ANPDDTFEAQLFYGDLK 28.613 40.670000;GHEDYQIVGNQR 10.209 12.830000;TTPSIIAYTQDGETLVGQPAK 24.961 35.101600;FATHGGYLLQGK 14.72 25.080000;FSLEGGDALIPMLK 30.071 47.210000;VVVYDMPGTTR 17.229 19.305000;YYLATGGGDISQAEVLLK 26.888 43.030000;GPAIAQAFDAEGKPSK 15.426 17.605000;AYPQEAAEFTR 15.738 14.210000;NNSLSQEVQNAQHQR 12.2 7.450000;VINELTEK 11.369 15.590080;AAFNQMVQGHK 11.568 14.490000;ELVELFFEEIR 31.938 52.985000;ETILSLPGSEWEK 26.134 38.440000;YTDTPAGAALVAAGPK 18.244 19.890000;DEAPQAPWEGTTIPAYDGWSDDGK 27.703 33.176880;DYQEIDDVLFK 28.026 36.370000;MIAVTTTSGTGSEVTPFAVVTDDATGQK 27.704 33.615240;FIATNPDTHGR 8.8383 10.725000;ANDIVVALLQEK 28.158 35.620000;SELEAFEVALENVRPTVEVK 29.726 49.890000;DEWQAVAPSWR 23.24 30.420000;QTTFNDMIK 17.368 22.155210;TVVPYTSEIYGR 19.895 25.015000;ENSMLPEAEQQR 13.56 15.590000;AGENVGVLLR 18.553 26.520000;PSDLIPELQGR 21.678 26.830000;DAQSALTVSETTFGR 21.212 26.200000;DDVIIGVNR 16.698 27.132105;SAEHEVSLQSAK 8.1922 11.820000;TEFLFMDR 25.182 36.222340;ADQIPVLK 16.698 20.906600;LNDMVNWGR 19.103 23.420110;EFYEKPTTER 9.2772 13.550000;AATYEQIK 9.2607 10.301940;MFEPMELTNDAVIK 27.425 39.400000;LINDAYDSEYFATK 20.836 28.310000;PVTGGYYFAPSLDK 22.674 30.420000;IAELAGFSVPENTK 22.688 27.560000;VGAATEVEMK 8.0368 13.715000;QSVDQPVQTGYK 10.655 9.045000;PGNVVLTPTILR 23.571 38.430000;ELSFLNSGVSIR 25.562 38.510000;LPQLGIEFSGPGAK 24.926 30.150000;GDIEQIPSMYSALK 24.735 35.080000;QYTTVVADTGDIAAMK 20.836 27.110000;GPHVPNMR 8.8995 9.157280;HPTFVEGDIR 13.15 23.310000;DALTGELFR 23.864 30.357600;LLEQEMVNFLFEGK 31.875 50.290000;RQDIESNLQYDAGDK 13.446 17.430000;ETEAIGHVVR 9.5446 20.490000;TEFYADLNR 17.399 24.899070;ILIGEVTVVDESEPFAHEK 25.206 42.160000;AMSAQDLLNSR 18.919 19.510000;NELFTATR 14.516 21.180940;VIGIGSPR 13.152 13.508880;MGPEPIIK 15.79 18.153740;MGMNIINDDITGR 20.642 29.215000;GMNTAVGDEGGYAPNLGSNAEALAVIAEAVK 30.995 40.980240;QDYYEILGVSK 23.087 28.250000;TLAEGQNVEFEIQDGQK 19.409 26.645000;ELTPAAVTGTLTTPVGR 23.854 31.310000;GTGNMELHLSR 13.37 22.450000;QEAAPAAAPAPAAGVK 11.199 10.185000;QAIVAEVSEVAK 19.807 21.890000;KNSNFYLK 7.9646 19.281372;EQFSDGVGYSWIDTLK 28.696 44.510000;ILVDENMPYAR 19.746 25.185000;EFASEPVWFGDSDR 25.326 34.525000;TSLDDWLR 24.918 28.616500;SSQTPWTLDPVAVGALDYR 29.314 44.550000;VEGWENAEAAK 12.79 13.880000;NAMPLLNR 17.365 22.510070;MTLTEIVAK 16.316 26.115320;YDSTHGRFDGTVEVK 11.606 18.170000;QNTFFVTNSGVQNR 18.154 25.170000;SAEQLAQAWK 17.277 20.620000;GFEELDTSK 13.957 16.842630;YIQDQHPEVK 8.0907 9.130000;PLSNLQNGQMVK 14.35 22.350000;WTQPGNIVTNGAYTLK 22.464 32.670000;RQLAYPINK 10.933 20.267590;QGNALGWATAGGSGFR 21.992 26.970000;PVDIVNAALK 19.728 27.120000;LLMAIGEPVMVDTAQEPR 27.686 41.965000;IGQVVEGYR 12.821 15.752870;ANLQPVLITGMEK 23.26 34.920000;DTVEIQGEVDK 13.193 18.705000;VEFFEGNFTEYEEYK 27.029 39.780000;AGDNAPMAYIELVDR 23.72 31.920000;FPGVDPLLGPEMR 27.504 35.690000;VTVQSLDVVR 18.695 25.015000;QYDINEAIALLK 28.993 37.310000;IATDPFVGNLTFFR 30.316 48.960000;LAEEIIYGPEHVSTGASNDIK 20.003 33.464840;DDLPEGVYNEQFK 20.416 28.710000;QNLATFCQTWDDENVHK 19.942 29.570000;LIGAPPGYVGFDQGGLLTDAVIK 30.17 56.435780;AEIDEEQLAAAPVIIR 24.82 41.010000;DYGVGSDVYSVTSFTELAR 29.352 41.570000;VDQLSNDVNAMR 15.79 19.520000;ISDEQLDQAIANIAK 27.538 28.690000;VVTFRPGQK 9.1867 12.483590;EHFAIYKD 12.958 22.202620;STVLAPTTVVTR 17.453 28.025000;GALALVER 16.788 23.574320;LDAPLIVVATQGGK 25.355 32.825000;FIWVQQPGS 25.058 25.979100;SLQEAVMEEIIK 28.348 34.570000;REGTDLFLK 15.678 29.822450;PVVGYIAGVTAPK 20.641 26.295000;TLAEGQNVEFEIQDGQK 19.398 26.645000;MTISVEGK 13.065 13.612940;QQEIAFTDK 13.61 16.433970;NQGDHLLHSTR 6.7476 13.970000;EVTTTPLAADDWK 20.185 27.680000;EISSHDSSTNGLINR 11.173 18.630000;NINLVIPR 21.383 26.951540;SNSDFVPSTTR 12.655 14.550000;VLVVEDNALLR 22.965 37.035000;IGNPEYFTDIEFR 26.676 36.490000;AEITPANADTVTR 12.767 14.710000;IGEEYELK 14.247 16.356340;GYGMGDAAEGK 9.513 7.090000;QPLEPMDITR 18.562 22.750000;TLELMEEVAK 22.285 29.770000;RPEMTYEK 7.2924 11.153340;YQLTALEAR 18.14 20.929230;LGTDGLQLYSSGK 19.361 23.170000;NLPIETNIMDLEAAK 28.181 39.520000;AEAGDVANAILDGTDAVMLSGESAK 30.412 37.697550;TPAQIVIR 14.688 22.614130;IEIPVLPLR 28.388 36.857240;VEDVLFAMVK 22.435 38.980000;SINAEVTNSAIELK 19.614 28.090000;LVDIVEPTEK 17.272 21.660000;IIAIDNSPAMIER 21.323 32.205000;ALSGGVGAEELK 14.063 18.040000;NELNEAAETLANFLK 32.215 39.790000;ETPADAEVISHQLMLR 22.724 40.040000;RNDVNPEITDR 8.8903 12.910000;NGLLWDNSLNVDGIK 28.083 46.850000;IDAILVDR 18.148 26.908970;ELSAEGFNFIGTGVSGGEEGALK 27.937 44.748180;GKPFAPLLEK 17.296 29.512000;ITTVQAAIDYINGHQA 28.175 29.690000;SENLYSAAR 10.662 14.585270;ANVSQVMHIIGDVAGR 26.723 31.595000;ELDAETAQAIISR 21.334 27.110000;GMPIATPVFDGAK 22.818 28.990000;ELINTTGDYAHIAE 20.591 29.460000;STLFNALTK 21.629 29.335950;HAVIALTSIYGVGK 21.988 37.675000;SPAGAIQFEAVDAPEIIPDPFDPSK 31.212 47.899650;AAATQPAEALR 10.501 10.965000;EALEWGTTGAGLR 20.947 30.760000;PLADVQEQVK 12.844 16.025000;MDHFMDWLAK 25.433 39.170000;RIEALAEDFSDK 17.265 27.050000;RYAGVGDIIK 15.201 23.725000;VVISGLQK 13.525 18.522680;TDVLIDEVR 18.872 29.097565;NTDIAEELELPPVK 25.649 34.250000;DQPWPEGGK 12.55 10.216500;EITLEAAR 13.866 19.042980;LAGITVPDR 16.927 17.358320;SATPAQAQAVHK 5.0708 5.020000;YEVISTLSK 16.976 23.259565;VVSLPSTDIFDAQDEEYR 26.891 35.330000;VYGDDAPQAWQK 14.694 13.480000;HHITADGYYR 8.221 14.460000;VTKPEAGHFAK 6.7811 10.720000;GNNVVVLGTQWGDEGK 22.032 30.150000;LSGFGNFDLR 24.326 32.170000;DAASFAPLHNPAHLIGIEEALK 26.795 54.407700;MPLVSLNTK 17.499 24.149860;FNVPVSDADIEK 19.826 22.885000;IDRPEEYADIATK 14.827 17.820000;DPATLDYVVPFK 27.172 37.385000;MAPPQISAEVLK 21.562 24.760000;HAVEFVASNAR 10.421 19.375000;EEQQVAESIALTDDTLVPFLAGETVR 32.365 59.750680;TFAHVDPVK 11.104 18.180505;MLQHTPIR 8.8123 15.760360;NHEAGGIYLFTDEK 20.158 30.330000;GNYSMGVR 10.047 12.061500;LDTWFSEMQPVEETQDGE 29.812 36.850000;QIVLNLYEK 21.422 30.148405;PVILAADGSEPDLSQQALTEK 22.896 37.734220;LLANQEEGTQIR 12.773 17.765000;FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK 27.396 53.315080;TDITELEAFRK 19.5 30.080000;GLNNPDLDAAVGEDLAQQLR 27.805 37.570000;VYVTFLNDKDEDAVK 20.13 30.555000;LGSQIYAMAGMQTR 21.164 27.470000;FESEVYILSK 21.515 32.850000;EGADVAISYLPVEEEDAQDVK 25.306 37.039090;LQTWIGVPGVDQSR 24.778 30.090000;FDGTVEVK 8.2128 15.977940;ENNFALPAVNCVGTDSINAVLETAAK 30.679 48.264040;HGAVSVMQYR 11.669 18.105000;GMPLYEHIAELNGTPGK 19.518 32.990000;LGVALATAESVVDAIER 31.215 39.745000;YHMEDVHR 7.3231 9.218770;WKDEYLDQGNQTQP 16.494 17.282000;VEEFPSEPPFDGVISR 26.055 35.280000;QDVPSFRPGDTVEVK 17.194 21.625000;ENVTSIIGNGVVLSPAALMK 30.273 55.065000;NPEAMAASLK 9.6805 17.430000;IMNVLGEPVDMK 16.607 30.980000;DLPDLVYMTEAEK 28.196 34.500000;YMLVDGQGNFGSIDGDSAAAMR 25.095 37.402560;GLMLNVTDPASIESVLEK 29.496 48.145000;LPSQPLPIIGSGK 22.581 27.750000;AAVIEAMTK 14.227 18.355645;LFPEITIK 24.456 28.834080;LEQAAYEMTALR 19.796 27.410000;FDLAFANDPDYDR 24.471 30.490000;VIDTTAAGDSFSAGYLAVR 24.433 37.180000;HQGTFDVAR 8.7591 13.573350;GAGPGVIK 10.988 11.181720;QIGIYSPNGQQYTPQDR 17.905 15.310000;NLMAETYPR 15.577 17.655085;TVLIPFENKR 18.893 30.140000;EGVHLHDEDPAR 7.8814 13.565000;HLESVVTNK 7.7333 14.449050;TIIMPVDVFEMELSDK 31.64 55.340000;DGFSLYDR 19.196 23.299980;MGLQNYLQAQIREEG 23.032 29.540000;VSVPTMDAAEAFK 22.575 24.615000;NVDLALAGITITDER 27.691 40.340000;AESFTTTNR 8.2595 10.080280;ANDAAGDGTTTATVLAQAIITEGLK 31.517 38.148600;NDHLMLAR 11.223 24.387880;GLEVGATGFDPK 18.763 21.870000;TMNLGTVSEER 10.309 18.140000;WEHNQDAMAVEK 7.8092 17.660000;IFTTFDSVAQDAAEK 22.594 28.350000;GPTDFVENYAK 19.604 19.630000;IGAGPWVVK 19.689 23.707145;GLMLNVTDPASIESVLEK 30.439 48.145000;NQPELSEDTIKK 10.657 13.920000;LLDGGWQTWSDAGLPVER 28.729 45.190000;HIVGAIANEGDISSR 13.779 23.745000;DFSNVGPMVR 18.703 24.290000;FVLLTSGATVADYNDAPADAQQSEVLK 26.989 44.919600;YYLGNADEIAAK 17.69 21.530000;VNIEIDPQTQAVVDTVER 25.27 31.690000;SSDYDAIIK 15.671 20.189750;DHLDDPVIGELR 20.037 32.010000;VLDLASPIGR 21.721 26.260000;VNALLADK 12.808 19.785590;GVEPIYETMPGWSESTFGVK 25.876 45.940000;TVAVEHAEPVYLR 15.491 29.015000;MESLASLYK 20.277 25.083940;LMTEFNYNSVMQVPR 22.204 36.160000;MITGIQITK 27.116 23.624440;DAIAAAIDVLNEER 29.655 34.150000;GDQDMILLLSK 26.498 41.630000;IHAEVPLSEMFGYATQLR 29.009 47.345000;ITEQEAQEMVDHLVMK 27.782 33.390000;QIEAGAPADLFISADQK 25.004 32.910000;DGHLIVNGK 9.2951 20.802740;GFGFITPEDGSK 21.896 24.910000;IVDLLTER 20.904 25.617680;VADGATVVSTSTR 10.543 11.310000;GGVVIEPMLTDQWYVR 29.286 52.825000;QDGPTALILSR 20.771 27.250000;ESWQDLPQNK 15.017 18.440000;IVIGELLK 24.325 31.435580;AVVTSMAETAR 12.49 14.785000;MMSPMSAEATVVR 20.133 23.430000;SSTAVNLALALAAEGAK 28.634 36.850000;IGVVSADGASTLDALEAK 23.954 32.365000;ILTLLGPNGAGK 21.578 27.710000;ALEGDAEWEAK 14.824 17.440000;VEGIDVVSLPFYK 29.107 43.080000;LAATIAQLPDQIGAK 22.849 30.515000;VETTDGVVQLSGTVDSQAQSDR 19.106 20.781360;AVSEATAEVDVISAAALSEQQLAK 29.204 38.147040;HAEQENMTLTELK 14.935 23.000000;LTGLEPGELFVHR 23.796 38.070000;HDGEAEDVAVALNEQYQPR 19.449 22.610000;GLALLDEELAK 25.282 36.745000;TDTLLEISVLPLDSYAK 31.136 57.430000;RFQDEEVQR 7.6632 12.133310;AVYEAIGFVAKP 23.425 30.710000;SQSTSLIER 10.899 15.928010;TKPHVNVGTIGHVDHGK 7.8761 15.612000;LPGVNDTR 10.365 8.182900;KPITDLGVK 11.618 20.250076;QLTAQAPVDPIVLGK 23.168 33.590000;KQVEEAGDKLPADDK 7.4915 9.277000;AADMTGADIEAMTR 10.77 20.290000;SLDDFLIK 25.562 33.081620;IGVMFGNPETTTGGNALK 22.538 31.065000;QYLIAPSILSADFAR 29.804 48.110000;GEDIEPLR 13.548 18.532140;DWPFFSTR 28.615 31.312600;GFMEVETPMMQVIPGGASAR 28.779 43.285000;LNHLVDDK 8.9496 15.438720;PWNSTWFANTK 23.295 30.250000;VNDEGIIEDAR 15.121 17.440000;ALGVGEVK 11.748 13.565640;LNNAPADSWR 12.894 14.970000;AIDDHTLEVTLSEPVPYFYK 28.671 51.160000;AVTAEVEAALGNR 20.239 18.210000;HESGVVTDPQTVLPTTTLR 20.33 34.370000;DVAGYAAGLELFDR 29.182 39.495000;VTLPEFER 19.915 21.512040;YGSDKPDLR 7.6205 7.774270;YETLHADR 7.3765 11.569580;ATPAVGFAMGLER 23.671 31.470000;ATAAEVFGLPLETVTSEQR 29.222 39.995000;LFPDTDPAFK 20.904 24.180000;TQTQLQQQHLENQINNNSQR 12.134 16.470000;PLIDQGGLSR 15.476 21.800000;VWVVEGSK 14.676 17.912510;DAFIDEMQR 19.119 24.714200;QVIGQVAADLR 26.931 23.910000;LLAEHNLDASAIK 15.289 28.065000;MLNETPALAPDGQPYR 19.012 22.760000;NVMGVPAVFVNGK 24.661 31.815000;EDLLASGR 12.46 19.837620;GFLIGGTSGR 18.016 23.310000;ANPWQQFAETHNK 18.336 16.970000;EHLIQFGGAEPVEGK 17.674 31.070000;DAGVLANYETPK 17.301 20.000000;SSAIFINAGR 17.148 27.025000;ATEFSPFELTK 24.608 30.520000;GITDILVVDNLK 28.24 41.990000;ASLEANVR 9.2954 12.222320;AQYVLAEQVTR 16.919 21.940000;YHFEQSSTTTQPAR 9.4461 9.970000;FNSSLSEDGQR 10.44 11.310000;NVALEEQAVEAVLAK 27.856 34.215000;ITEQEAQEMVDHLVMK 26.585 33.390000;RGIEEAYR 8.2348 13.016960;ISFLIDADGK 23.789 30.490000;QVTIAQLEDVKPLLMK 25.279 47.760000;GMLTTDVESYDK 18.159 21.040000;AIAEALGWEDK 20.718 25.635000;VLVVDDLLATGGTIEATVK 30.587 50.735000;EGDALLQGGSLTGNGSVEK 19.099 28.890000;EGVAGWAGENLPLVR 26.722 40.565000;QILEGVEYK 16.546 20.929230;MVVTLIHPIAMDDGLR 26.577 50.705000;LMPSEQWK 15.298 17.226660;NVLQPIGWDAFGLPAEGAAVK 30.454 55.156840;LNELLEFPTPFTYK 30.554 49.570000;ETDLVTIGNDAWATGNPVFK 28.932 49.890000;GFIDVEQVR 19.815 24.772580;LASTEWVDIVNEENEVIAQASR 30.641 43.390080;ANGLSEAMLDR 17.56 22.620000;EGTRPAVVIPTNEELVIAQDASR 23.595 40.609620;GPQLPAPNMLMMDR 26.352 36.230000;EAFATIAVAADK 19.665 28.160000;GDLGLVIACLPYA 32.406 48.530000;LLADGMESFNK 19.329 25.165000;DRFNVPVSDADIEK 18.098 28.390000;RVVEPLITLAK 19.499 37.675000;VNYQGIGSSGGVK 12.644 11.840000;YGIPQISTGDMLR 21.365 30.940000;DQNEDIDREDFAQR 12.858 16.350000;PLSLEEEK 9.3585 15.750900;VATEFSETAPATLK 17.574 23.810000;FTVLISPHVNK 18.62 31.785000;KVYAGVLK 6.1604 19.177312;VMLSELFR 27.027 32.570780;SQNGAAMSFGR 8.1644 13.670000;MQAASGQLQQSHLLK 12.329 24.735000;AIAQVGTISANSDETVGK 16.19 19.635000;VAEETPHLIHK 9.1819 17.610000;RFPLHEMR 10.168 22.912120;FFVASDVHPQTLDVVR 22.855 41.945000;SDILNAVGITLDETTTR 29.333 41.880000;HAVTEASPMVK 9.0334 15.075000;VGLLPNLINR 26.379 36.140000;IPAQPQLMLDMEQAR 25.56 34.145000;EFAQLAFK 21.037 29.065850;MNIIEANVATPDAR 20.121 24.790000;SNNALQTIINAR 22.309 25.250000;EMLEHMASTLAQGER 20.233 28.880000;IMLDNEDITHVPAENR 18.648 27.580000;FFETVNPLK 22.283 28.187810;NMITGAAQMDGAILVVAATDGPMPQTR 30.548 46.804780;EMYTFEDR 16.008 19.374080;LFGSIGTR 16.698 18.853780;DSLFWGEQTIER 25.732 39.930000;YGQAIGHIGK 9.7633 14.690000;IPGFEQPWEEDFGKPER 24.833 32.370000;AFDLPHIK 18.996 22.779680;FIFLEEGDIAEITR 29.691 49.250000;VTITIAADSIETAVK 25.221 33.990000;AEAVSIMTDAVR 19.882 25.235000;IVGIDAVVLSTQHSEEIDQK 23.56 38.180000;IAVMWSEK 15.859 23.020910;VGLFGGAGVGK 19.355 22.840000;ETGEIHYGR 7.8661 12.347370;ILFYTGVNHK 17.245 26.910000;TLEFLGHK 15.813 25.040620;YNMELTLEEAVR 26.045 32.265000;ISHGQVDLSELGPNADELLSK 24.623 39.676640;DVNQLTPR 11.335 13.357520;TQFTGDTASLLVENGR 21.407 32.170000;SGDPLLSPSAVATLR 24.459 34.750000;HSADEITVR 8.165 13.043065;IDDIDLNLEDFVQR 29.647 45.170000;AGQTSMIAR 8.6547 12.473860;KLTPEQAEQIK 9.9583 14.002000;IVGDDEADFK 13.61 16.080000;VNGIAPGAILTDALK 27.969 37.640000;LLLENNTK 8.2418 17.680740;VTDGSGAAVQEGQGDLWVK 19.215 27.740000;IHSEEDERPIGR 7.0905 8.670000;VALVQPHEPGATTVPAR 14.403 20.910000;DGMNLVAK 13.719 20.154530;LRDELPGVK 12.614 17.737790;ASDTLLAGGTMNNLGGEDSDTIVENGSIYR 26.35 43.705200;MYDLGFIK 24.863 30.726080;PANQLLTQPVK 15.291 22.500000;VDLVFAPSVK 23.297 32.120000;ELPGLIADGR 19.828 26.160000;YIAETFLEDAR 24.054 29.605000;EFLDANLA 25.02 28.900300;EYLPASYHEGSK 11.534 14.730000;ALIAAPLR 18.125 24.491940;SLSTEATAK 6.5 8.630510;GGTTVIAFVEDPDGYK 25.808 32.550000;SANMTLEQLESLNAEQK 23.232 30.920000;ELQLVYNK 16.881 23.848660;TEQQLTAMK 6.6189 16.531270;GVELAPGESVPMVGVVEK 25.12 37.440000;ANISVHEPR 8.4575 9.603510;QVFGGQVVGQALYAAK 23.761 33.060000;GAELELLR 20.365 27.925920;TYHYYSLPLAAK 19.22 31.040000;NDYGMPEQVQFWSVADNVR 28.723 43.030000;RPFSMASTPDEK 11.595 16.490000;HTEALGELTR 11.955 20.430000;VLHTLAGK 7.3277 16.280660;NQNLLVTK 12.225 20.783620;LMVEQILTDLQK 28.388 38.935000;GITINTSHVEYDTPTR 15.946 22.190000;YRAEELAEER 10.324 13.025000;DIADAVTAAGVEVAK 23.867 25.945000;MNPETNSIANR 8.0116 8.890000;LESGDLPLEEALNEFER 31.079 43.210000;HGYAFNELDLGK 19.7 30.430000;HGLITGATGTGK 8.3783 15.530000;PLLGTLEYGLPHK 23.551 37.450000;TVTHMQDEAANFPDPVDR 18.096 24.040000;VPEGIGETAIVQIR 22.61 31.520000;PVMIIGHQK 11.545 22.374135;ELIPGSELWPR 26.943 36.560000;SEDQVELVEK 13.248 19.590000;KVEIPGVATTASPSSEVGR 17.154 19.872000;DTIGDIIILPR 28.888 45.080000;QIIANTVDFGASDAPLSDEK 23.47 35.160000;AAGAELVGMEDLADQIK 21.347 35.290000;QLYQHMAK 7.4843 12.761540;IGSLGMDVYENER 16.151 25.190000;LTDDDMTIIEGK 14.59 26.470000;GFTSEITVTSNGK 16.557 19.410000;LVALEDGVSLLNR 26.681 39.435000;ELLSQYDFPGDDTPIVR 26.933 43.010000;ELTQPDDVRK 7.832 11.310000;ASIPEGVWTK 18.666 22.270000;VPTPNVSVVDLTVR 24.622 30.820000;TFTDAAEVIGEAWESR 29.237 35.410000;TSLGSLLDHTGAFGESEK 23.881 34.850000;VVVETPVGLNER 18.17 23.005000;ELGASDEADLQR 12.681 17.110000;FVEDPHTVVK 11.898 18.600000;ELAAFSQFASDLDDATRK 27.852 37.485000;LDGPVTGNGK 7.9527 7.350000;AEIEGEIGDSHMGLAAR 17.397 28.510000;EVQEISPNLR 14.667 20.780000;MLGFGYDWSR 26.731 32.460000;ALQAIAGPFSQVR 23.161 28.740000;YMPESMDIVHYVDK 22.999 30.930000;AALANLFSELPSK 28.183 35.790000;APLAAEAMGIIAPR 27.004 34.200000;EIPFLYASSAATYGGR 27.004 36.180000;THLTEDVINAAEK 15.477 23.030000;SVQTVTGQPDVDQVVLDEAIK 25.795 37.408840;FVNILMVDGK 25.695 35.900000;IMIDLDGTENK 20.004 22.630000;VIGTATSENGAQAISDYLGANGK 24.03 23.643440;AANGVFDDANVQNR 14.434 14.490000;LLADDIVPSR 19.814 24.665000;AVLDELNK 10.851 18.078060;LMVEQILTDLQK 28.404 38.935000;HADNTLTFGPR 13.595 18.500000;NMNVPGEDQYR 9.8345 10.140000;SVGGEVIEQPR 12.543 14.940000;QLLEYDDVANDQR 18.069 21.190000;AVVESIQR 16.025 13.324410;FSEEAASWMQEQR 19.653 22.710000;LIDDAVAWAK 21.353 28.410000;GTAMNPVDHPHGGGEGR 7.1567 5.425000;GELLVPER 15.704 24.018940;SPHYIVMNDK 8.0796 20.880000;EYASFTQEQVDK 14.901 17.905000;LVGPLIDVMGSAGEDLK 29.443 42.860000;GEVNPPLR 11.146 14.156890;YFTEAGVGFVVPDDSR 25.554 33.950000;AFDQIDNAPEEK 13.873 13.980000;SGLYEDGVR 12.076 16.492350;ILFYTGVNHK 17.252 26.910000;GTVSTESGVLNQQPYGFNTR 20.18 28.625000;VFLGTDSAPHAR 12.645 18.200000;YSASGARIPR 7.5949 8.865000;EAQAPIVAITGSNGK 18.05 20.860000;ALYTVVTEGK 16.984 19.640000;NEQDGGDLVYFQGHISPGVYAR 23.249 39.064680;YHVSNYQPSPMVR 13.633 17.245000;GETPFEINSR 16.748 19.190000;FLDADIQNNPQK 15.669 17.930000;FDMSEYMER 21.643 21.663845;MVSNASALGR 10.414 13.730000;TVINFDNAIIAAGSR 26.03 38.090000;MTLTEIVAK 20.835 26.115320;FVDEAFEQK 15.466 19.265400;VVADAIAK 9.2866 14.005530;EDLLISGPGK 17.593 27.270000;NVGENALAISR 15.147 20.340000;EIVFTSGATESDNLAIK 22.556 39.105000;LFDDAGLILVDFK 31.509 54.830000;AIQQQIENPLAQQILSGELVPGK 29.956 46.424680;GQAHWEGDIK 10.099 16.045000;VMLQAYDEGR 15.163 17.630000;TGFHNGEPVK 7.0361 12.750000;QLIQVNPDILMR 24.771 40.240000;EANNLGIPVFAIVDTNSDPDGVDFVIPGNDDAIR 32.82 62.599440;FNLMLETK 22.209 28.862460;WDLGDIIGAR 29.526 33.950000;QDLDQLQAGAR 13.193 16.850000;ANDIDVPAALIDSEIDVLR 31.805 52.620000;AFTGVGGTPLFIEK 25.594 39.280000;SQPEVNDAITK 10.869 12.020000;QTDLVEAMAK 17.206 21.070000;NIPVGSTVHNVEMK 14.256 23.040000;VTIEGWNGPVEIDQIK 27.14 38.890000;LEHNIIELQAK 15.768 29.460000;GHVMNALPEDAK 8.1111 19.105000;GEAIGVIAAQSIGEPGTQLTMR 26.621 43.803180;VLSGPQAQPAGDK 8.7096 5.760000;IPAQPQLMLDMEQAR 25.568 34.145000;MQLNSTEISELIK 23.923 34.060000;AVDGVTLR 12.841 15.807660;EAYELVAPILTK 26.321 39.360000;KGDEIAAVVLQVDAER 23.631 32.882000;AVVEELAR 13.37 16.824610;NESHGIIVK 8.3389 19.061070;QDPLAYLER 23.346 23.449300;IIVVTSGK 11.786 17.221930;NINMTSYASSK 13.552 13.290000;VAFSGGEFQLNADR 21.687 27.710000;VVGQLGQVLGPR 19.861 25.330000;SGFQYHGR 7.8717 10.642500;KAGNVAADGVIK 8.7584 15.352000;ALLNSMVIGVTEGFTK 26.537 46.740000;QAYEQGFSQR 10.542 9.540000;YFGTSDMEYGK 17.024 17.150000;TMNLGTVSEER 14.832 18.140000;SGNANTDYNAAIALVQDK 22.745 24.950000;VIIETGELK 16.679 23.867690;ESDLALMTNAGTEIGVASTK 25.255 39.490000;ARLPATDGQVK 8.3423 11.980000;VQNAAGDIVSLR 17.893 22.160000;PGMDSLAPEDGSHRPAAEPTPPGAQPTAPGSLK 16.788 17.591090;AVEAMPGVMR 16.282 19.710000;RWEQLIEGK 12.825 26.105590;DVFLGLDK 23.518 31.331520;GAFVSEVLPGSGSAK 20.51 26.320000;EGFEVTYLAPQR 22.714 31.090000;EHGTQHSDAILTR 7.4517 16.670000;ILVAVEADK 16.782 21.002205;PAIVGEPEMPSTIADVAQEK 24.837 31.950000;QYAPMSVAQQSLVLFAAER 30.393 46.985000;GANVTVPFK 16.608 19.382160;ALGLATSIR 9.3646 23.585520;TFTDAAEVIGEAWESR 29.231 35.410000;LVSPEGLAPEWD 27.844 29.460000;ILEITPER 16.726 19.591660;ANPQPELLK 13.871 16.609110;GFYAEHDGK 7.3881 9.934330;EMTPVEQK 13.538 10.008680;KFEELVQTR 12.912 21.544166;ILVAVGNISGIASVDDQVK 27.129 42.385000;AGNTIGQLFR 21.283 25.220000;RQQEIAFSDK 10.088 17.030000;AVAAVNGPIAQALIGK 25.342 32.885000;ATEVNYHDSGATIR 10.029 14.620000;LSVFKPIAQPR 17.297 26.845000;ILDVLIPPAK 24.777 34.910000;IALVTGASR 13.699 17.377780;GGDTVTLNETDLTQIPK 22.697 31.250000;IGNTEMPSR 9.1414 8.358070;TSLGQSIAK 10.236 16.200450;AVEDLVNTQIR 18.805 23.810000;NYTPYEGDESFLAGATEATTTLWDK 30.698 44.909700;ALQSSINEDK 8.7868 10.140000;LVADSITSQLER 20.637 24.435000;NVWLGDISR 21.606 29.180270;FNSLTPEQQR 12.061 14.010000;LGGRPEYR 6.9672 6.593620;IDGVFGDTAVVTEWLK 31.308 50.270000;GVVVAIDK 13.719 20.353190;FEDNALPWDSIDAATYVK 29.05 45.150000;DRLVPIIADEAR 19.936 34.690000;VEGKDDVTGEELTTR 11.344 14.280000;VFTYYTPK 16.599 17.974000;LFDLGQVPK 22.977 25.132590;INPGNIGNEER 11.38 8.240000;TLGADALEPK 15.264 18.470000;GITDILVVDNLK 28.239 41.990000;GALIDSQAAIEALK 24.314 35.520000;DLPLIASNFR 26.776 37.300000;EAGVQEADFLANVDK 23.075 29.460000;VGDIVVSDEAR 13.878 19.440000;APVAGALLIWDK 28.007 42.875000;LADAIAEPLLDK 22.891 33.240000;GDVAALNVDALTENQK 21.784 28.705000;DLHVSVEDK 11.545 16.151800;VNDGPFADFNGVVEEVDYEK 29.146 37.640000;SLHLNGPIVDK 15.368 26.520000;IHLTPAER 9.1318 13.121020;LVTDELVIALVK 30.299 46.260000;ITAYADELLNDLDK 28.306 36.265000;QAFLDFFHSK 25.675 36.240000;VINGVVELPK 19.377 25.380000;YSLLGAMR 21.092 24.775740;HNDLENVGYTAR 12.071 16.630000;GEGMVLTGPK 9.9506 20.290000;IKLPADSTR 10.235 11.775246;MAGMAIDGK 7.5779 14.604730;TTPSIIAYTQDGETLVGQPAK 24.97 35.101600;EAGVQEADFLANVDK 23.075 29.460000;TAEDYLGEPVTEAVITVPAYFNDAQR 30.724 48.657920;SLNFLDFEQPIAELEAK 34.806 55.070000;AVYADTEGFVSEMDTR 22.965 26.710000;EAMEPEFK 10.116 16.777310;VIGITNEEAISTAR 18.143 24.580000;MLHLTVEQAR 15.008 23.610000;SPAFDSIMAETLK 26.895 35.405000;AVGDSLEAQQYGIAFPK 24.357 31.710000;NYITESGK 8.6035 8.551840;LMGMTPTGNGR 13.907 12.260000;FPYNTPTSK 12.608 9.331070;YVESQVYQGVVENLASEQAAR 25.92 29.757480;FAEGLETVGDNFLR 27.272 38.680000;DLHVSFKVK 11.791 25.006100;SDWMEMEK 17.85 22.022880;HEIGSAYPGDEVR 11.278 15.520000;PLDANQMAALVELLK 32.05 47.450000;DGDHYATLNTAQK 9.4786 14.130000;TLQEQYPNDTDLLGR 20.403 27.770000;GDGAILTATQNGYGK 15.377 20.630000;VLAEQALAQPTTDELMTLVNK 29.333 47.559710;NGEFIEITEK 18.872 27.250000;GDVVYLDLR 22.281 32.892265;ATLESIAYQTR 17.771 21.420000;DMTAAELK 12.53 16.857720;AILGSMER 14.415 17.746960;IVQSPDVIPADSEAGR 16.755 18.580000;DLSPLWEMYK 28.898 39.050000;TQVVVLGAGPAGYSAAFR 24.89 39.245000;GQAIIATADGR 11.778 17.345000;STCTGVEMFR 12.452 21.700000;QANDLDYASLPDSVVEQVR 26.405 32.940000;VPQDYVTQSGPLR 16.506 17.720000;ALANFMFDSDEAMVR 24.922 42.615000;HVDPAAAIQQGK 10.571 10.520000;REEESAAAAEVEER 9.2408 10.050000;FQDDILAGR 18.349 22.310890;TQLNEAEANVTVAK 14.119 19.470000;ETVDFVDNYDGTEK 18.901 23.565000;GMNTEGVLPGPLR 21.56 26.740000;DVLLDINK 20.27 28.966520;AGETLGLVGESGSGK 15.487 18.220000;AADMTGADIEAMTR 13.859 20.290000;STPFAAQVAAER 16.287 19.700000;GYSLTIPIAQEDGAPVSTILDETYK 30.563 51.233700;LFGKPEIDGSR 13.674 17.330000;AGVLAEVR 13.793 18.158470;EFVETALR 16.796 24.714250;IDLTIVGPEAPLVK 27.747 42.470000;QILPEANSQIVGFR 23.609 35.110000;FGNMSGQMR 8.2737 14.118230;VVDLNNLTFDR 23.636 32.330000;VVEEANNVDIR 12.424 16.830000;RSDVIEIR 11.661 21.578260;GDEIAAVVLQVDAER 27.259 35.730000;IMLDGEPYAVEASEFVK 27.197 40.980000;SEVAVPGIDASTFDGIIQK 28.199 44.065000;LAASIAFK 17.343 22.718190;APVVVPAGVDVK 18.814 23.675000;VAALNGLNR 13.222 17.985905;FFINPTGR 20.15 20.547120;AMMDNLQTETVINR 20.131 28.985000;GIANSILIK 19.861 30.226245;VGDAVQADVDEARR 11.534 10.740000;ITEQEAQEMVDHLVMK 22.58 33.390000;VGSLGWANK 14.829 19.207020;DDYLPLIR 25.387 33.214060;VHVHVEEGSPK 6.7584 7.295000;QAIVAEVSEVAK 19.13 21.890000;VMLQAYDEGR 13.194 17.630000;DWADRPGEENK 9.1588 8.525000;DIGEPSVLNR 16.699 20.770000;QGGALVTSTAATVTGINR 20.04 26.870000;LMTEFNYNSVMQVPR 24.508 36.160000;GMPIATPVFDGAK 18.638 28.990000;IIAATHQNLEQR 9.2579 17.605000;VEEHIADALEK 13.509 20.480000;AYTPAWAEQITGVSR 24.377 26.860000;LTIVPAQTSAEDVLK 23.312 33.020000;ELNEQLEENLGYK 19.987 27.910000;DYYAIMGVKPTDDLK 22.745 34.970000;LGPVYSVR 15.024 17.330720;EHTTEHLRAELK 7.767 19.070000;GIAVADTAR 9.8118 13.393345;DVSIMPFK 24.508 27.925920;EIGNGFSELNDAEDQAER 20.383 23.330000;LSYTGEVK 11.362 11.891220;GLGEMNPEQLWETTMDPESR 30.308 36.870000;GAVSVLDNLSPIK 26.361 33.995000;GHSSAQYSGEIK 7.258 8.130000;KYEQEIDVR 10.008 14.032606;LAHEAQLDVAPLGK 16.714 26.990000;KVEVMNTDAEGR 8.8472 12.272000;ISVVGMYR 17.944 21.441090;AAVVVPDNVLFEGGK 25.261 35.265000;VIYQAGFTDR 15.975 21.280000;AIEQTAITR 10.445 14.556080;DSGSDMVLVGLLR 30.071 44.030000;VGSFDGGWGASYMAR 22.924 29.940000;GETDELAALVAQQR 23.967 26.390000;WLLPDPIETLK 30.488 42.490000;VLVTTLTK 15.291 22.642510;VVNVGDVVEVMVLDIDEER 29.718 50.630000;GTGNMELHLSR 9.6174 22.450000;NNEVYLIEVNPR 21.858 30.850000;YREPVLVSGTDGVGTK 15.269 23.250000;GTLEDPNLFIR 23.695 37.550000;TVAAMDVLAPGIGEIIGGSQR 30.508 46.652590;DVLETLGTDK 19.117 23.320000;IINEPTAAALAYGLDK 24.926 34.430000;FQAAMLAADDDR 18.002 22.605000;LPAAQSFAALSGMEEGGK 24.942 28.825000;YMANAMGPEGVR 12.46 15.555000;ESTHEYLK 6.8954 13.518340;MAIPEQPLEILR 28.058 37.760000;EQVLAAISLVR 26.041 39.085000;MYAVFQSGGK 16.932 17.655000;GDVLNYDEVMER 17.171 28.805000;MIQEQTMLNVADNSGAR 18.509 24.880000;AAFQPVFLEVVDESYR 30.995 46.590000;SAYALGASLGR 17.303 22.520000;GYAGDTATTSEIK 11.019 11.565000;IFSSAAYNAGPGR 14.447 14.950000;GIEGSSLDVPENIVHSGK 19.905 27.290000;ILGGDLEPTLGNVSLDPNER 27.277 39.310000;IMNVLGEPVDMK 23.516 30.980000;TLLVTTGSEAVENAVK 22.677 31.870000;LGEELDAAK 12.014 14.760410;TPEAALNFMR 21.743 28.930000;DDANIHAIQR 9.0974 18.385000;LKPDPNTLCDEFK 18.225 21.132000;TSDIHETIIK 11.762 25.150000;VFDEFKPLVEEPQNLIK 27.848 49.200000;MGHAGAIIAGGK 8.9764 17.190000;WNGVTVTPK 14.657 17.192910;LDNKPGYVNTGFIDAEVLK 24.174 42.950000;QLVNGTVDGR 10.608 12.465000;MGFNNLGVDNLVENVK 27.603 39.940000;GGLVPGALLAR 23.398 34.350000;DNPPQDLIDLNPNQSVPTLVDR 28.411 42.748560;EEQQDNAIFSAK 14.572 18.830000;VLLPAFPDIR 28.599 37.760000;KHESGVVTDPQTVLPTTTLR 18.116 32.162000;TSAESILTTGPVVPVIVVK 29.31 52.525000;IHPGEAMDK 6.8547 8.095360;YHVNQYTGDESR 8.6731 5.845000;YAEIADHLGLSAPGDR 19.631 26.760000;DDLQAVMAMVR 28.807 34.910000;SVPLAEVQPGMLLR 27.557 42.990000;LNEMGIVLEDGPQGTTWR 27.01 41.470000;GEMPQTIGGGIGQSR 15.93 18.165000;SHDALTAVTSLSVDK 19.631 28.430000;AGGGSATLSMGQAAAR 8.3175 13.120000;SFLESLGSDQAK 20.547 24.010000;LNAQLAQAEEK 11.771 13.945000;RAVASVLMSK 12.459 24.955000;ENTLQQAVGLPDQK 18.716 22.690000;GAPDVFEQFNTAVQK 24.768 28.870000;NGIIHTTIGK 11.568 22.600000;VDAYAGDPILTLMER 30.042 40.995000;LAEPAPTGEQLQNILR 24.38 30.640000;VANLGSLGDQVNVK 18.98 24.310000;LDQALAEMFPDYSR 27.632 35.550000;TTPTPQPGSDEINR 10.637 8.000000;VGAITAANR 8.0703 11.204095;ALNAAGFR 13.115 15.836040;MSADLLIVPFIDK 31.365 53.315000;TPAAEITPR 10.325 12.459265;HYGALQGLNK 10.977 18.170000;VEAVMAGMDITPER 21.651 28.655000;DGVGLLPTVLDVVENPK 32.159 48.305000;ADNPFDLLLPAAMAK 32.009 46.200000;HTSDTPFDVSK 12.061 13.730000;DTPINWTTK 17.472 21.678440;YMANAMGPEGVR 14.588 15.555000;VIGQNEAVDAVSNAIR 22.442 23.580000;LQTYIDQVEGK 17.481 18.990000;MNFSHGSPEDHK 7.4059 6.240000;KEAAPAAAPAAAAAK 7.2284 7.997000;NILNELQK 19.271 21.086340;EITPVNIEEELK 24.467 30.330000;NQASNDLPN 9.4819 6.854785;TGMTFDAQPGR 14.486 16.525000;ELPYMNFPK 24.054 27.594280;ETTFNELMNQQA 23.812 27.590000;QPVDVTPVWMMR 26.846 36.725000;ANEAYLQGQLGNPK 15.572 17.920000;AVQLGGVALGTTQVINSK 23.306 34.810000;RAEITPANADTVTR 10.523 16.280000;DLVNQDGGENWAK 16.863 19.925000;NLLIGNVIK 23.479 34.512310;AITGIFFGSDTGNTENIAK 26.359 38.560000;FGLLPDPAR 23.338 25.210430;INALETVTIASK 20.673 27.565000;AAMEDVLK 15.416 17.373290;KEAAPAAAPAAAAAK 7.2264 7.997000;VLSSIADK 10.324 14.908960;LEPIVSVLQSDPEQFEQLK 29.691 49.360000;TLLTQVAPPGVTAHVVDVAK 23.955 42.370000;NTQGVILIR 17.61 29.919750;AEGQQLVNQAMGGILQDSINEMGAK 29.984 40.324800;QIEDAGLR 10.428 13.442660;IAPDLSEEELIQVADSLVR 32.768 49.110000;TVLEADPMPVIK 22.629 33.440000;IGTDTTYAPFSSK 17.808 18.190000;LRPSPVDEAK 7.8278 9.480000;LAAHAAIK 5.2497 13.068990;ANPDMSAMVEGIELTLK 31.002 41.770000;ESGALFVDR 16.244 25.969370;NDLQGATLAIVPGDPDR 23.154 31.930000;PGSLLINASR 16.418 26.630000;MQQLQNIIETAFER 30.122 36.660000;WSVPLTVR 22.741 26.719770;TVYDNVAIPLIIAGASGDDIR 30.912 53.579240;MKPNIHPEYR 8.7264 11.802000;SAEHEVSLQSAK 8.1944 11.820000;LMEVEQVLESAR 25.732 28.160000;PSSEVSMIHAR 10.843 18.230000;TQTGALIMIFDSADGAADAAPAQAEEK 29.691 40.472740;WTADVEAGK 12.372 13.665785;WTDNAEPTEDSSDYHVTTSQHAR 12.203 12.808460;GFQASTEQQNNPPAK 8.6829 4.810000;EQTQAAVSLLEK 18.372 25.110000;TFGMEGLFR 20.469 33.578230;VYLASAAPEIR 18.41 24.080000;FLNDFPGAETIR 23.884 31.830000;VLVVDGGGSVR 14.006 19.735000;QASITQELTEIVSGAAAV 34.299 36.440000;VLHANISR 7.8001 14.294060;APEVVEEQR 8.6643 9.146200;TNVLFFTK 23.249 34.268850;AATFAASLGLK 22.074 28.990000;GVNLPGVSIALPALAEK 29.767 45.440000;TVGFKPAGGVR 10.381 17.140000;ATAVAIHGK 1.7882 12.157635;EGDAVQLVGFGTFK 26.971 40.290000;PAGQVIAQYYEFLR 28.517 40.900000;ALEADLAR 9.8554 16.592840;TPPAAVLLK 17.976 27.127240;IVQVPFAELVK 28.192 37.180000;VFEGNRPTNSILLR 17.444 33.500000;ISDQWFPQ 25.502 23.829740;GPIPLPTR 17.308 18.333480;EALMGVMGDK 18.39 24.960000;AVDGEMITVTVEGK 20.46 25.210000;KTNDTLAVTGEAFSR 16.378 24.582000;IVQAMEEIK 15.485 20.510840;QPGLDFQSAK 14.928 16.750000;IDNEEVLIMSESDILAIVEA 35.73 65.270000;PLDIDLPQLIVK 29.98 47.350000;ADFLSDPNK 14.363 16.735600;EGNDFYHEMTDSNVIDK 20.504 29.590000;AVITGDVTQIDLPR 23.298 31.860000;ASAPGQISVNDLR 17.403 18.795000;MILVTGHR 11.944 20.036280;YSTEMMNVISAGLDK 27.159 33.790000;VATYDLQPEMSSAELTEK 18.383 29.310000;PHGAPVPENFR 11.945 15.010000;PGMSVEAIQGIIASMK 27.542 39.905000;SALDSQQGEPWQTIR 20.724 24.720000;IATVTNFPHGNDDIDIALAETR 24.573 42.048720;VHIINLEK 13.974 23.810820;TIVSDGKPQTDNDTGMISYK 15.264 22.265000;VAVTGAGQSPALDVTVHAIGK 20.508 30.748410;EAEAYTNEVQPR 11.622 10.760000;DSILEAIDAGIK 29.031 36.280000;LTIPVDMATEFMTQVAK 32.01 46.720000;EHEGEIITGVVK 14.887 26.370000;GGVNVAAGTGISNYINTIPVEEQPEYPGNLELER 29.105 49.064100;LEFLAADALR 25.629 37.310000;LDMLNEELSDK 20.591 27.625000;LMLLAPIIK 28.88 43.843380;RSDVIEIR 11.672 21.578260;WNMLHPLETPR 22.105 32.145000;ELHEETVR 6.3097 12.827760;MQQQLGDQYSELAANEGYMK 21.783 28.560000;TLGLYPDEVVLR 25.995 41.070000;GGNFFQPTILVDVPANAK 28.616 47.350000;MLEDQNLISAHGK 14.717 22.760000;DFNPSGIILSGGPESTTEENSPR 25.882 29.401020;VLPAVAMLEER 27.288 31.410000;DALPTEEEQFAAYK 20.651 26.800000;VWTGGGDEAALAR 15.787 20.960000;MVYAPTQEHGK 8.479 7.530000;DGEDPGYTLYDLSER 24.071 27.630000;ISDNTVMTTSR 11.155 13.990000;NTADGAAVNAMR 6.7753 12.425000;FQAMAAEGVK 13.431 17.205000;NIIQHVENNR 10.813 14.640000;IINIPSAEAAR 17.201 19.630000;MPIITFIDTPGAYPGVGAEER 29.804 46.509620;GVELAPGESVPMVGVVEK 21.889 37.440000;TLTPDVPVLSEEDPPGWEVR 28.35 43.170000;GLEDFYYSVGK 24.119 28.570000;IPVGHVEAGLR 13.237 21.845000;NDAAVLLVDHQAGLLSLVR 29.057 62.005000;ASLSAFDYLIR 28.633 40.320000;GADVGITMIAR 20.01 28.220000;YFPDATILALTTNEK 28.885 38.900000;NNDTVHDFTK 8.7346 13.350000;GYEIHISDEALK 18.348 27.590000;FLDQVAAK 13.993 17.907780;AAQVPVVVAVNK 17.471 22.690000;YMANAMGPEGVR 9.2429 15.555000;EMLIADGIDPNELLNSLAAVK 33.062 61.999680;EMLIADGIDPNELLNSLAAVK 32.187 61.999680;DDVMGESLYNPMLPGIVADLK 32.185 61.413010;GEVVASTFDEPASR 15.661 21.765000;MESIEQQLTELR 25.515 29.380000;TSHITVVVSDR 11.447 22.900000;GGDFAALAK 15.005 21.357350;AAVLPANLIQAQR 22.117 30.465000;KVVMTGPSK 6.3367 11.040631;QTHQTPVIMLTAR 17.189 28.320000;METDLQAK 8.8939 12.089880;SVANSQDIDLVVDANAVAYNSSDVK 26.08 36.655950;NHFASEYIYNAYK 17.789 28.730000;IFVDEGPSMK 13.526 20.625000;FLTGYDLR 20.186 25.286580;KQTEQALLK 8.9785 16.737546;LILPLAIGK 27.744 36.789130;ATLPAIGR 15.139 17.046920;VGAATEVEMK 7.6099 13.715000;HILTSNIEK 14.847 20.501110;LIEVPVEYIAGK 26.208 33.410000;ADLNVPVK 13.817 16.933400;KDEDGNYLVDVILDEAANK 30.038 38.362000;PVLLEPIMK 19.821 35.923160;INPAGAPTYVPGEYK 17.607 17.040000;AHPQLAEEFTR 14.133 19.250000;LPLPTQYQTAIK 21.331 26.550000;MELLLLSNSTLPGK 28.547 47.580000;SSPHLELLR 16.116 27.438600;MNTEATHDQNEALTTGAR 9.545 11.140000;VALYGIDYLMK 28.405 40.810000;NETNELFIPPGPR 21.734 28.790000;LGADGNALFR 18.57 23.745000;VIESLIHSGEPLGLEAGSK 22.209 39.580000;EDDVGNIVHK 10.356 17.470000;IVGNNVHPGTAK 6.8794 8.080000;IAYNVEAAR 12.57 13.096580;GITVNVVAPGFIETDMTR 28.596 46.090000;ELDAQPTGFLDSR 21.117 27.610000;ADKPLVGTGMER 10.945 17.880000;ALMVHVGGDNMSDQPK 15.114 20.415000;GYRPQFYFR 18.15 25.920720;EADAALGR 6.7878 11.692560;VVGDDFAK 11.354 13.461580;ELGEMDDAPDENEALK 17.54 22.510000;VDYSTFLQEVNNDQVR 26.989 31.720000;ADIAVHSMK 8.6243 16.784250;QLPLNFYQIQTK 25.537 35.540000;FDFSILDK 27.205 32.154540;LVDAINQLR 18.196 24.286080;TPPGAGDEIQLTDAIDMLIEK 31.423 49.378880;VTYDPVSK 10.805 9.497840;ALAESIGITVEK 19.935 27.115000;LEDAQVQLENNR 13.409 17.210000;REEESAAAAEVEER 9.2406 10.050000;IEVPEIGEEVIEIK 28.443 38.905000;AIAEATEAGLLGK 19.142 25.035000;LTTDAEALAR 13.612 18.170000;AEGQPDAAMEDGDALIFMNFR 30.502 47.781560;DMALIGQALIR 28.051 42.295000;APMILALANPEPEILPPLAK 30.113 62.075000;FKDDVNEVR 8.9606 13.254206;YMMDNNYQYSK 14.386 14.355000;NIVNDPSVVFDDIVTNEQIQK 29.761 46.110290;SNTFVAELK 17.672 25.152050;ALGEYLEK 15.541 17.822640;DAAPDATFR 11.892 14.181475;MTIGEIIAEPLR 28.352 38.190000;ITWEGSQNQDADVSSDGK 14.497 12.540000;TWFELAPK 24.165 29.770620;FEFFIGGR 26.437 33.157300;ALVEQEPSDNDLAEEELLSQGATQR 23.758 32.377950;TLVDGIRK 10.44 17.638170;ISAAAHNELTK 7.7269 13.265000;SGEQTAVAQDSVAAHLR 14.351 19.350000;INIIDTPGHVDFTIEVER 26.458 46.740000;PILPAEWLTSATK 29.05 38.170000;MAVAESGMGIVEDK 11.843 23.685000;SGTGSVDYAK 7.0903 7.850000;QINADNVHNLK 9.38 14.510000;VSELQIPVNAGSVDMLDQLGQVSPGHR 28.805 45.767480;LDNLHVAMVGDLK 20.619 37.550000;YYDELPTEGNEHGQAFR 16.651 21.030000;AAVEAETLK 9.9662 13.879845;VANLTAFTPDFK 25.206 32.210000;DGHLIVNGK 8.3702 20.802740;LFGVTTLDIIR 29.499 44.530000;DEVIDHLGTIAK 19.532 32.545000;GQYTAGFAQGK 11.742 12.270000;MKPFIFGAR 19.672 30.164946;LIAAAPTAVAPEESGFYAR 22.748 34.485000;DLDLEASAAAHPVR 16.62 25.450000;LAPQELEQK 12.077 12.250070;NTNITGVIVNK 16.042 24.650000;SALLVLEDGTQFHGR 22.065 40.220000;TDKPQPVNALLK 15.061 23.710000;EANNLGIPVFAIVDTNSDPDGVDFVIPGNDDAIR 32.824 62.599440;GPSIMPGGQK 11.573 11.530000;GEYQYLNPNDHVNK 13.344 17.190000;TLASVQTLDESR 15.958 21.545000;VEDALHATR 7.4798 12.824140;QEVHGNGLSSYPHPK 9.0924 10.685000;GMVLTGGGALLR 22.351 36.515000;HVDSLITIPNDK 17.823 26.720000;LSDEVTDSPMVDK 15.376 18.770000;ILLINPTDSDAVGNAVK 24.581 35.410000;QYVAFASQR 13.939 16.234505;QLADGTAKPEVQDETLVTYAEK 19.004 27.084780;VIAAGANVVR 12.165 18.255000;EVENRPAVSLK 10.204 17.480000;DITADVLK 18.873 22.297220;VMQAQGSQLTNK 8.9255 10.630000;FSGNYGNMTEVSYQVAK 19.916 25.310000;TSPSIDVAFQAVDQDALK 27.415 38.100000;ASDESELESLLDATAYEALLEDE 35.206 50.409960;VADVINNDPLVGDK 19.529 25.110000;AAVEEGVVAGGGVALIR 22.527 35.065000;MAPPQISAEVLK 21.553 24.760000;ASAQLETIK 11.315 16.341535;DGIPLSVIEAAQQSPVYK 28.468 40.180000;REDVVVATK 7.1331 17.368050;DTFWWADK 26.173 33.422180;DQLLENLQEGMEVK 26.422 36.850000;RGIEEAYR 8.2298 13.016960;DLEHPIEVPVGK 18.584 25.950000;IKEELFYPSNEEK 16.475 24.902000;GMVIGTTGFDEAGK 16.258 25.915000;RLGPDPDQFGG 16.575 16.230000;MFDNLTDR 17.147 18.352400;LPSPQVVGAESEEEDASHAA 17.807 14.550000;AGDNAPMAYIELVDR 27.6 31.920000;PGNTLPMR 11.079 14.975180;LMVEQILTDLQK 30.202 38.935000;TTDVTGTIELPEGVEMVMPGDNIK 24.879 48.946400;KLVPYYTVK 14.572 23.426921;LAAEGHFAEAR 9.8118 15.615000;DFGYAVFGK 23.852 28.304570;ELAAFSQFASDLDDATR 29.452 39.385000;DAIPTQSVLTITSNVVYGK 29.501 43.650000;DGIPVLLETEAR 25.964 35.780000;QIADDYQQALR 16.117 17.185000;GGHSGGEIHVGLGNANK 9.4911 14.100000;ELAASGVMLTSDENVVR 21.953 37.385000;IVGLQTEAPLKR 15.444 24.680000;LIMGLADGEVLVDGR 27.856 42.985000;DQAVLIEPFDTVTVQGFYR 30.839 58.725000;QGWMVLNGPK 21.789 27.120000;AREVPAAIQK 8.1111 13.380000;TGVSNTLENEFK 18.763 22.825000;NVEASFELNDASK 17.665 21.540000;GSAAHGYFQPYK 12.645 16.725000;SPNAAEEHLK 6.8296 9.430000;ESVAELAYYR 19.717 25.865000;YIASLADELR 21.874 27.505000;ELIPGGVNSPVR 18.726 23.860000;LIDQATAEIVETAK 21.095 26.210000;ASVEIDRK 5.9711 9.360670;GATAYQITEAR 12.912 14.520000;GDEIAAVVLQVDAER 27.264 35.730000;GHNVTVIDPVEK 13.217 21.430000;KGNTSLYDHNNNTSDYSK 7.8006 7.482000;IELSSAQQTDVNLPYITADATGPK 26.349 36.183520;IGENINIR 14.312 17.680740;QQIEEATSDYDREK 9.9931 7.040000;NSTAMLTTFNEVNMK 20.6 35.850000;AEAGDVANAILDGTDAVMLSGESAK 28.034 37.697550;VMLLFTNPTDVER 25.001 39.530000;INPAGAPTYVPGEYK 17.611 17.040000;SPFVTSGIR 16.107 22.700090;SLYEADLVDEAK 19.973 26.070000;MYAVFQSGGK 27.508 17.655000;YGDMINHR 8.5607 12.572340;VGINELLR 21.826 26.100140;AMVSNTATVLR 15.015 22.385000;IGYSGDSSSDISLKPLNYEQK 20.173 28.682740;ERVYLAEEGGR 10.571 19.425000;NDDYLILDAVNNQVYVNPTNEVIDK 29.539 54.990900;VMSLLEPTK 19.352 23.964990;APLVEELYR 21.288 26.952100;ISGAGIQESHVHDVTITK 12.785 25.490000;IVGLQTEAPLK 18.613 25.980000;AFYPTDAK 12.548 11.333080;GEMPSDFDAK 13.702 16.065000;VAIAGAGGR 7.9856 8.523480;IFVDEGPSMK 16.862 20.625000;LQTMVSHFTIDPSR 20.379 31.390000;VADITLDPK 15.822 18.788630;GAEAMQAAGLK 7.6972 16.120000;TSAESILTTGPVVPVIVVK 29.308 52.525000;DRGEDALIIYDDLSK 25.272 39.690000;LLMEEADFQK 16.87 26.965000;TVNGLEPVQQGEITVDGIVVNDK 26.54 41.902920;GATVLPHGTGR 8.3384 13.120000;ALVTGGDSGIGR 12.947 15.115000;YEVPVIIEAFPETLAGEK 31.082 50.505000;PANPPSMPNDPSK 10.855 2.800000;KGDLAAVAVK 10.042 20.282000;MTSVGSQDTTGPMTR 12.132 11.340000;LETIEGSK 8.64 12.118260;VGEHTLEK 7.292 9.497840;LASLTSLI 30.469 31.161240;MQDLSLEAR 15.615 20.102180;GVIVQYGGQTPLK 18.392 26.690000;DADNLLQHR 11.744 17.708600;ATGVPLAK 9.8556 13.925120;ASTPLGVGGFGAAR 19.545 23.120000;LVLADEPTGNLDAR 18.989 26.160000;GQNEDQNVGIK 8.0473 7.570000;TGVLDQVLDALK 31.423 38.225000;QTAFSQYDR 12.266 12.400885;MIAPILDEIADEYQGK 29.526 37.555000;FGVSAAAAVAVAAGPVEAAEEK 27.652 30.409020;GGIYLYPSTASHPDGK 16.419 22.400000;NFLRCANWEEK 18.99 25.610000;AVIESENSAERDQLLENLQEGMEVK 21.836 39.022800;GAVLASEQALR 14.664 24.195000;SAGGIVLTGSAAAK 14.397 22.420000;QLVAQLEK 13.951 18.413890;DLTPEVSFR 20.552 24.470950;IDLPAADPER 16.519 15.870000;AVSTLPADVQAK 13.772 16.910000;GSIPENTVIAGVPAK 18.498 26.300000;VYVMSEAK 7.0903 13.863630;GLMLNVTDPASIESVLEK 30.427 48.145000;ISYISTGGGAFLEFVEGK 30.508 48.090000;VMVMVDDK 14.1 18.830130;GVTYSAPLR 13.772 17.650220;AAIERDFGSVDNFK 18.059 27.540000;ELIVASSYSK 15.341 24.060000;IIGEQLGVK 15.858 21.337890;AHPQLAEEFTR 14.133 19.250000;EVPAAIQK 9.1619 13.556180;VVLGLVTPDEGVIK 25.86 40.030000;GVSADQISIVSYGK 20.122 25.240000;PAHLAQSIR 7.9637 16.589650;GDVLEMNIR 14.395 27.929965;KFTGEVSLTGQPFVMEPSK 22.888 40.142000;DSAIVPVYYYVNAR 27.061 38.205000;HNEEDPIFVR 14.765 22.930000;DTFLMIDK 23.028 34.368180;AYEAYDFHEVVQR 18.011 25.160000;TEFDVILK 24.274 32.911340;VGWSADYAEALK 22.263 26.490000;QADGVVIVTPEYNYSVPGGLK 25.801 39.578040;KLEHAVPMAK 7.6493 17.402000;EMLLDAMENPEK 20.09 29.480000;ANAVVMATGGAGR 12.294 14.595000;LADPEELEFMGIR 27.851 38.540000;IFWDPATDTLTIK 29.008 42.100000;LGQSATTLSGGEAQR 10.639 11.570000;VIIAAHGNSLR 11.121 23.030000;TQMSAAHTPEQITR 10.137 14.045000;LDHQAPEEMR 8.07 10.800000;LVSEALAER 13.549 18.448080;NPSISTHLLGSNASSVIR 19.733 35.730000;FTAEGVQEIDYK 17.883 22.985000;TTDVTGTIELPEGVEMVMPGDNIK 27.492 48.946400;TFEAIQYYLLK 28.917 45.410000;AGALIAPIQIVEER 27.592 39.695000;GGAGFSTGLK 12.779 18.125000;ALAEAGCSAVK 9.1928 10.715000;LYLNSFNQTR 18.695 23.810000;PAVEAVENGDIQFVPK 22.288 31.575000;NLVELVQK 18.412 23.692570;SVALEEGAILAR 20.785 32.815000;GPGETQLETDRR 7.8671 8.030000;KRVSQASDSYYYR 8.1776 12.917000;VGNLAFLDVTGR 26.204 35.140000;TQTGALIMIFDSADGAADAAPAQAEEK 28.55 40.472740;GIILAGGSGTR 14.234 21.840000;SILLTALAR 24.449 35.991270;LGTSTVSPIELENAVR 24.183 31.970000;PMLSPQAELELLETDER 29.088 42.320000;LVIPPELAYGK 24.191 29.610000;AQNLQLDAEDKK 9.4423 12.840000;ILADIAVFDK 25.634 35.185000;HQKPVPALNQPGGIVEK 12.419 21.630000;LGTLAAGK 8.811 14.256220;YSPELDSHGQYSLPASGK 16.151 18.840000;NNLTSAYIR 15.374 20.870850;MYALTQGR 13.086 13.674430;TGTLTEGKPQVVAVK 12.526 23.750000;GGDTVTLNETDLTQIPK 22.709 31.250000;MFEPMELTNDAVIK 25.411 39.400000;KFVQAYQSDEVYEAANK 15.071 22.917000;IVDTLTASALNEPSAMAAALEL 33.14 54.412560;MTVQDYLLK 23.181 28.815395;MSYVEGLLSSNQK 22.624 25.340000;ELPPEERPAAGAVINEAK 15.186 22.660000;AQIHVEDTER 8.0723 9.790000;GGTPYGATTIAGGDGSR 9.0741 10.850000;AAGYELGK 9.6564 10.207340;ASFGGQIITVK 19.472 28.220000;TVDDFINEVIEPNK 28.364 33.940000;PLTEVEQK 8.9209 11.683100;VVANQIIHMLESN 28.596 34.005000;DDDPSFDELVALAVETGR 32.568 43.110000;GLEEDAEGLR 15.39 16.870000;EGQSLPVGVGQPTLK 19.009 26.890000;NAVWAAIK 18.615 23.834470;GLVEPQQER 9.6021 10.357585;TELFLVEGDSAGGSAK 21.969 33.190000;ESIEEAVSEVVNALK 34.038 33.540000;VPDMSAYR 13.493 11.560120;GDGTTGEDITSNVR 11.797 11.930000;SGETEDATIADLAVGTAAGQIK 26.857 31.930200;GPIEFSNQELDDLIIR 29.825 51.580000;VNEDEMYPGEAGIDIYNLTK 28.064 38.640000;TQLINLFEVADGK 29.295 40.770000;VFDVDSQR 12.261 12.392600;MGMNIINDDITGR 23.067 29.215000"
  })
  # Final Extraction of the Transitions
  Transitions <- reactive({
    if(!exists("db")){
      db <- dbConnect(SQLite(),dbname = dbpath)
      
    }
    SpecFilFil <- FilteredTable()$SpecFilFilB
    
    
    TransiPos <- strsplit(as.character(SpecFilFil$SpecIDPosition)," ")
    
    rows <- unique(unlist(lapply(TransiPos,function(x){x[1]:x[2]})))
    
    Transitions <- dbGetQuery(db,Wrapper("*","Transitions","rowid",rows))
    # Transitions <- dbGetQuery(db,Wrapper("*","Transitions","SpecID",SpecFilFil$SpecID))
    dbDisconnect(db)
    rm(db)
    return(Transitions)
  })
  # Export Button
  ExportEvi <- reactive({
    print("Start")
    RT_Templatehumpe <- sapply(unlist(strsplit(as.character(evidence()),";")),function(x){unlist(strsplit(as.character(x)," ",fixed = T))})
    RT_Templatehumpe <- t(RT_Templatehumpe)
    RT_Templatehumpe <- RT_Templatehumpe[,1:2]
    colnames(RT_Templatehumpe) <- c("Sequence","Retention time")
    print("Write")
    return(RT_Templatehumpe)
  })
  output$ExampleFileDL <- downloadHandler(
    filename = function() {
      paste("ExampleFile_RetentionTime_Template.txt", sep="")
    },
    content = function(file) {
      EE <- ExportEvi()
      write.table(EE,file,quote = F,row.names = F,sep = "\t")
    },
    contentType = "text/txt",  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
    
    
  )
  output$ExportButton <- downloadHandler(
    filename = function() {
      paste(paste("Picky_Export",input$SRM,gsub(" ","_",input$ExportTypes,fixed = T),sep = "_"), ".zip", sep="")
    },
    content = function(file) {
      wd <- getwd()
      try(ExportWrapper(FilteredTable()$SpecFilFilB,input,Transitions(),file,session))
      setwd(wd)
    }
    ,
    contentType = "application/zip",  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

  )
  
  
  
  ###### 
  # Prepare MSMS

  TransiTemp <- reactive({
    
    ExtractMSMS(FilteredTable()$SpecFilFilB,Transitions())
  })
  
  
  output$MSMS_pro <- renderUI({
    if(is.null(SpecFilFilPre())){
      return(NULL)
    }
    
    selectInput("MSMS_pro","Protein",sort(unique(unlist(paste(FilteredTable()$SpecFilFilB$GeneSymbol,FilteredTable()$SpecFilFilB$Proteins,sep = " # ")))))
  })
  output$QuantileAutoResult <- renderUI({
    testfufi <- md()$BiquanTableAuto
    if(length(testfufi) > 0){
      ht <- paste(round(max(Diffis<- apply(md()$BiquanTableAuto[,1:2],1,diff),na.rm = T)),"min")
      
    }else{ht <- ""}
    helpText(ht,label = "HUMP")
  })
  output$MSMS_seq <- renderUI({
    if(is.null(input$MSMS_pro)){
      return(NULL)
    }
    SpecFilFil <- FilteredTable()$SpecFilFilB
    inp <- input$MSMS_pro
    inp <- sapply(strsplit(inp," # "),function(x){x[2]})
    Sel <- SpecFilFil$Proteins == inp
    Seq <- SpecFilFil$Modified_sequence[Sel]
    Seq <- sort(unique(unlist(Seq)))
    Seq <- Seq[order(nchar(Seq))]
    
    
    selectInput("MSMS_seq","Sequence",Seq)
  })
  
  
  msmsSelection <- reactive({
    if(is.null(SpecFilFilPre())){
      return(NULL)
    }
    
    SpecFilFil <- FilteredTable()$SpecFilFilB
    if(is.null(FilteredTable())){
      return(NULL)
    }

    se <- input$MSMS_seq
    SpecID <- SpecFilFil$SpecIDPosition[SpecFilFil$Modified_sequence==se]
    TransiPos <- strsplit(as.character(SpecID)," ")
    rows <- unique(unlist(lapply(TransiPos,function(x){x[1]:x[2]})))
    
    db <- dbConnect(SQLite(),dbname = dbpath)
    Transitions <- dbGetQuery(db,Wrapper("*","Transitions","rowid",rows))
    dbDisconnect(db)
    
    msms <- PlotFun(Transitions,type = input$spectraPlotType)
    return(msms)
  })
  output$MSMSplot <- renderPlot({
    if(is.null(msmsSelection())|is.null(input$MSMS_pro)){
      return(NULL)
    }
    if(is.null(SpecFilFilPre())){
      return(NULL)
    }
    msms <- msmsSelection()
    if(dim(msms)[1] == 0){
      return(NULL)
    }
    SE <-  input$MSMS_seq

    plot(msms[,1],msms[,2],type = "h",xlab = "m/z",ylab = "Relative Intensity",frame = F) 
    abline(h = 0,lty = "dotted")
    pos <- rep(1,dim(msms)[1])
    pos[msms[,2] >=0] <- 3
    
    col <- rep("dodgerblue3",dim(msms)[1])
    col[grep("a",msms[,3])] <- "red"
    col[grep("b",msms[,3])] <- "red"
    col[grep("c",msms[,3])] <- "red"
    
    text(msms[,1],msms[,2],msms[,3],srt = 90,pos = pos,cex = 0.8,col = col,xpd = NA,offset = 1)
    if( input$spectraPlotType == "Both"){
      legend("bottomright",legend = "MaxQuant m/z deconvoluted Spectrum",bty = "n")
      legend("topright",legend = "Raw File Spectrum",bty = "n")
      
    }
  })
  output$FragmentPlot <- renderPlot({
    if(is.null(SpecFilFilPre())|is.null(input$MSMS_seq)|is.null(msmsSelection())|is.null(input$MSMS_pro)){
      return(NULL)
    }
    xl <- c(1,SeSeq(input$MSMS_seq))
    xl[is.na(xl)] <- 2
    par(mai = c(0,0,0,0))
    
    plot(1,type = "n",ylim = c(-1,4),xlim = xl,frame = F,axes = F,xlab = "",ylab = "")
    
    try(FragmentPlot(msms = msmsSelection(),se = input$MSMS_seq),silent = T)
    
  })
  triggerTypos <- observe({
    
    
    if(fileUploaded()&!input$triggerType){
    }
    if(!input$triggerType){
      for(i in names(input)){
        removeTooltip(session,i)
      }
    }else{
      
      addTooltip(session, id = "Proteins", title = "Protein search query field: Paste here uniprot accessions or gene symbols separated by Space, Comma or Semicolon. You can also have look at http://www.uniprot.org to find potential proteins of interest.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "Charge", title = "Set the charge states of precursors. These will be included in the initial precursor list before filtering.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "ExampleFile", title = "Example Retention Time Template. Required is a tab delimited text file, with >Sequence< and >Retention time< columns.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "Fragmentation", title = "Fragmentation Types to be included based on the ProteomeTools data. HCD is the default for beam-type fragmentation.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "IDType", title = "Defines Protein ID Type used in the Protein SearchQuery Field. Choose between uniprot accession and gene symbol. You can also have look at http://www.uniprot.org to find potential proteins of interest.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "loessQuan", title = "Sets the outer quantile borders. All points in a window outside these borders will be considered as outliers.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "loessSpan", title = "Controls the degree of smoothing.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "loessWin", title = "Defines the number of windows the data will be split into. This has a direct effect on the accuracy of the outlier detection procedure.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "MassAnalyzer", title = "Detector Types to be included based on the ProteomeTools data. FTMS yields high resolution fragmentation spectra.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "MissCleave", title = "Missing cleavages (i.e. when trypsin does not cleave after a lysine or arginine) to be included in the initial precursor set. It is recommended to select complete tryptic peptides. Nevertheless for specific biological questions or experimental set-ups the inclusion of miss-cleaved peptides can be meaningful.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "ModType", title = "Includes available modified peptidoforms. In the current version of ProteomeTools only Oxidation (M) was considered as variable modification and is therefore available in Picky.",placement = "right", trigger = "hover")
      addTooltip(session, id = "OnlyIsoFormSpecific", title = "Include only protein isoform specific peptides. This will reduce the intial set of protein representing peptides. Please note that some proteins/isoforms are only represented by non-proteounique peptides and will therefore be excluded if this option is enabled.",placement = "right", trigger = "hover")
      addTooltip(session, id = "PRMType", title = "Choose between a scheduled or a static method. A static method triggers all selected precursor masses during the whole data-acquisition. 
In contrast to this, a scheduled method triggers the selected precursor masses only within predefined retention time windows. 
Therefore the scheduled method can monitor more precursors than the static method while obtaining similar cycle times (recommended). 
                 The scheduled method requires a template table with a >Sequence< and >Retention time< column from any complex sample analyzed on your current LC-MS system. 
                 Picky uses these data to correlate experimentally observed RTs (from the uploaded list) with their calculated hydrophobicities (Krokhin et al., 2004, doi: 10.1074/mcp.M400031-MCP200). Finally, this fit is used to predict RTs of peptides to be targeted via their adjusted hydrophobicities.",placement = "right", trigger = "hover")
      addTooltip(session, id = "ProteouniqueSeparate", title = "Enables seperate filtering of Protein Isomers. If enabled Picky will filter for the best scoring peptides of each protein but will keep all protein isomers in the list. If disabled, Picky will group all isomeric proteins as one ProteinGroup and filter this entire list.",placement = "right", trigger = "hover")
      addTooltip(session, id = "RetInput", title = "Tab delimited txt with >Sequence< and >Retention time< columns.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "SeparateModifications", title = "Enables separate filtering of modified peptides (in the current version of ProteomeTools, this includes only Oxidation (M)).",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "spectraPlotType", title = "Display the corresponding MaxQuant m/z deconvoluted and/or the raw file spectrum.",placement = "right", trigger = "hover")
      addTooltip(session, id = "SRM", title = "Choose between single reaction monitoring (SRM) and parallel reaction monitoring (PRM). In case of PRM Picky will export an inclusion list for the best peptides of the selected proteins. For SRM Picky will export inclusion lists with a number of transitions per peptide (to be specified below).",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "SRMtotalcount", title = "Total number of selected transitions per peptide.",placement = "left", trigger = "hover")
      # y(# addTooltip(session, id = "tabs", title = "This is an input.",placement = "right", trigger = "hover")
      addTooltip(session, id = "UpperLimit", title = "The maximal number of peptides (in case of PRM) or transitions (in case of SRM) to be monitored in parallel during the HPLC run. This value should be chosen by considering both the injection times (i.e. dwell times) and the duty cycle of the mass spectrometer: Longer injection times increase sensitivity but result in longer cycle times. For example, an injection time of 100 ms allows for the analysis of ~10 features per second. With a maximum cycle time of 2 sec this would enable the analysis of up to 20 coeluting features.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "ExportButton", title = "Exports all resulttables as one zip file. Provided is the inclusion acquisition list and spectral libraries in txt format (according to msms.txt import into Skyline).",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "RetWin", title = "Retention time window width in min. This value should correspond to the deviation in the retention time prediction plot. The bigger this value, the less precursors can be added to the acquisition list.",placement = "bottom", trigger = "hover")
      addTooltip(session, id = "RetWinAuto", title = "Sets the retention time window automatically based on the set quantile deviation observed in the provided retention time template.",placement = "bottom", trigger = "hover")
      addTooltip(session,id = "QuantileAuto",title = "Quantile deviation in % for automated setting of the retention time window. This value defines the quantile deviation width of retention times from uploaded peptides. The maximal observed retention time window is used.",trigger = "hover")
      addTooltip(session,id = "RescaledHydrophobicities", title = "Enables the usage of rescaled hydrophobicities based on observed retention times in the ProteomeTools dataset. In most cases, this will result in an up to two-fold more accuarate retention time prediction.",trigger = "hover")
      addTooltip(session,id = "mzrange", title = "Sets the range of mass to charge ratios of peptides to be included in the query. This range must match the specification of your machine and method.",trigger = "hover")
      addTooltip(session,id = "SRMabovecount", title = "Defines the number of transitions that should be higher than the corresponding precursor mass to charge ratio.",trigger = "hover",placement = "left")
      addTooltip(session,id = "Species", title = "Selected species for the query. All sequences provided by Picky are based on Proteome Tools data and therefore represent the human proteome. The mouse dataset is a subset of this data and contains approximately 80,000 sequences that are conserved between mouse and human.",trigger = "hover")
      addTooltip(session,id = "SILAC", title = "Enables the inclusion of stable isotopically labeled peptides into the precursor/transition list. If set to \"lysine and arginine\" precursor and corresponding fragment masses will be shifted for each lysine by 8.01419 and arginine by 10.0082 Da and added to the list of selected sequences. This will increase the complexity of the selected peptide set by two and less peptides can be included in the final acquisition list.",trigger = "hover")
      addTooltip(session,id = "ExportTable","List of selected peptides to be included in the final targeted acquisition list. Not all peptides faithfully represent the abundance of the corresponding protein. To minimize this problem, it is generally recommended to target at least two peptides per protein. Please adjust the number of proteins to be targeted and other parameters accordingly.",placement = "center",trigger = "hover")
      addTooltip(session,id = "AbuBinDwell",title = "For SRM methods Picky can preset dwell-times based on their average abundance as are provided by ProteomicsDB. To this end, the abundance range (based on iBAQ) was split into three equal windows of low, average and high abundant proteins which were assigned to the dwell times 100, 50 and 10 ms respectively. Proteins not identified in ProteomicsDB are considered to be low abundant and therefore assigned to the 100 ms dwell time fraction.")
      addTooltip(session,id = "FixedDwellTime",title = "Sets a fixed dwell time for all sequences.")
    }
  })
  # 
  
  # Adding Tooltips

  output$PRM_Schedule <- renderPlot({
    if(is.null(SpecFilFilPre())){
      return(NULL)
    }

    HUI <- input$tabs == "Spectra"
    
    reti <- RetTab()$`Retention time`
    if(is.null(reti)){
      RTrange <- c(0,1)
    }else{
      RTrange <- range(reti)
      
    }

    par(mai = c(0.8,0.8,0.1,0.1))
    if(HUI){
      SeSe <- input$MSMS_seq
    }else{
      SeSe <- NULL
    }
    if(input$SRM == "SRM"){
      TypeName <- "Transitions"
    }else{
      TypeName <- "Features"
    }
    FeatureRetentionPlot(initTable = FilteredTable()$initTable,SpecFilFil = FilteredTable()$SpecFilFil,UpperLimit = input$UpperLimit,RTrange = RTrange,main = "Unfiltered",SelectedSequence = input$MSMS_seq,SILAC = input$SILAC,TypeName = TypeName)
    # }
  })
  # Image plot of filtered Acquisition list
  output$PRM_ScheduleFil <- renderPlot({
    if(is.null(SpecFilFilPre())){
      return(NULL)
    }

    HUI <- input$tabs == "Spectra"

    if(is.null(RetTab()$`Retention time`)){
      RTrange <- c(0,1)
    }else{
      RTrange <- range(RetTab()$`Retention time`)
    }

    par(mai = c(0.8,0.8,0.1,0.1))
    if(HUI){
      SeSe <- input$MSMS_seq
    }else{
      SeSe <- NULL
    }
    if(input$triggerType){
      
      addTooltip(session, id = "PRM_ScheduleFil", title = "These plots provide a graphical representation of the peptides to be targeted within certain retention time windows. Every peptide is shown as a horizontal bar that is colored based on the corresponding protein. The x-axis shows the retention time. The left y-axis shows the total number of features (peptides for PRM, transitions for SRM) during the entire HPLC gradient. The right y-axis (scheduled acquisition only) shows the number of coeluting features that is indicated by the black line in the plot. The horizontal dashed line represents the user defined maximum of co-eluting features. The left plot shows the unfiltered list (all possible peptides/transitions matching to the specified proteins). The right plot shows the optimized set of peptides/transitions selected by Picky that does not exceed the user-defined maximum of coeluting features. Not all peptides faithfully represent the abundance of the corresponding protein. To minimize this problem, it is generally recommended to target at least two peptides per protein. Please adjust the number of proteins to be targeted and other parameters accordingly. Details about targeted peptides and proteins can be found by clicking the Table tab.",placement = "left", trigger = "hover")
      addTooltip(session, id = "PRM_Schedule", title = "These plots provide a graphical representation of the peptides to be targeted within certain retention time windows. Every peptide is shown as a horizontal bar that is colored based on the corresponding protein. The x-axis shows the retention time. The left y-axis shows the total number of features (peptides for PRM, transitions for SRM) during the entire HPLC gradient. The right y-axis (scheduled acquisition only) shows the number of coeluting features that is indicated by the black line in the plot. The horizontal dashed line represents the user defined maximum of co-eluting features. The left plot shows the unfiltered list (all possible peptides/transitions matching to the specified proteins). The right plot shows the optimized set of peptides/transitions selected by Picky that does not exceed the user-defined maximum of coeluting features. Not all peptides faithfully represent the abundance of the corresponding protein. To minimize this problem, it is generally recommended to target at least two peptides per protein. Please adjust the number of proteins to be targeted and other parameters accordingly. Details about targeted peptides and proteins can be found by clicking the Table tab.",placement = "top", trigger = "hover")
      
    }else{
      removeTooltip(session,"PRM_ScheduleFil")
      removeTooltip(session,"PRM_Schedule")
    }
    if(input$SRM == "SRM"){
      TypeName <- "Transitions"
    }else{
      TypeName <- "Features"
    }
    FeatureRetentionPlot(initTable = FilteredTable()$initTableB,SpecFilFil = FilteredTable()$SpecFilFilB,UpperLimit = input$UpperLimit,RTrange = RTrange,main = "Filtered",SelectedSequence = SeSe,SILAC = input$SILAC,TypeName = TypeName)
    # }
  })
  output$showexportbutton <- reactive({
   as.character(as.numeric(!is.null(SpecFilFilPre())&input$Proteins != ""))
  })
  
  output$UpperLimit <- renderUI({
    if(input$SRM == "SRM"){
      TypeName <- "Transitions"
      MaxVal <- 300*input$SRMtotalcount
      MinVal <- input$SRMtotalcount
      stepfun <- input$SRMtotalcount
      if(input$SILAC != "none"){
        stepfun <- stepfun*2
        MinVal = MinVal * 2
      }
    }else{
      TypeName <- "Features"
      MaxVal <- 300
      stepfun <- 1
      MinVal <- 1
      if(input$SILAC != "none"){
        stepfun <- 2
        MinVal <- 2
      }
    }
    if(length(input$PeakWidth) != 0){
      PeakWith  <- input$PeakWidth# om s
      MSSpeed   <- input$ScanSpeed
      RequiredPointsPerPeak <- input$PpP
    }else{
      PeakWith <- 30
      MSSpeed <- 60
      RequiredPointsPerPeak <- 20
    }
    
    
    CycleTime <- PeakWith/RequiredPointsPerPeak
    CycleTime <- CycleTime*1000/MSSpeed
    
    sliderInput(inputId = "UpperLimit",# id for identifiying the input later on
                                               label = paste("Maximal number of ",tolower(TypeName), " monitored in parallel",sep = ""), # Title or Explanation
                                               value = CycleTime,min = MinVal,max = MaxVal,step = stepfun # input specific arguments
  )})
  

  output$ExportTypes <- renderUI({
    
    if(input$SRM == "SRM"){
      machines = TypeList$SRM
      selected = "Thermo TSQ"
    }else{
      machines = TypeList$PRM
      selected = "Thermo Q Exactive"
    }
    
    selectInput("ExportTypes",NULL,machines,selected = selected)
    
    
  })
  
  
  observeEvent(input$show, {
    if(length(input$PeakWidth) != 0){
      PeakWith  <- input$PeakWidth# om s
      MSSpeed   <- input$ScanSpeed
      RequiredPointsPerPeak <- input$PpP
    }else{
      PeakWith <- 30
      MSSpeed <- 60
      RequiredPointsPerPeak <- 20
    }
    
    showModal(modalDialog(
      numericInput("PeakWidth",label = "Peak width [s]",PeakWith,min = 0),
      numericInput("ScanSpeed",label = "Scan speed [ms]",MSSpeed,min =0),
      numericInput("PpP",label = "Points per peak",value = RequiredPointsPerPeak,min = 0),
      
      footer = modalButton("Close")
    ))
  })
  
  outputOptions(output, "showexportbutton", suspendWhenHidden = FALSE)
  session$allowReconnect(FALSE)
  
}


shinyApp(ui = ui,server = server)
