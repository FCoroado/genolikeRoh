#Functions to calculate Runs Of Homozigosity.
#Input: Beagle File
#more information on how to get this data at http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods
#output: Length of each ROH detected and total sum of ROHs
#---
#Get a matrix with the chromossome, position and a column with the prob of
#being homozygous P(hom) for each individuals
#---
#Input: beagle file, nIndividuals
#Output: Matrix with Chr number, ChrPosition and P(hom) for Ind 1 and 2 

get_position <- function(beagle,nind){
  chr_pos <- strsplit(beagle[,1], split = '_')
  chr_mx <- matrix(unlist(chr_pos), nrow=2)
  chr_mx <- apply(chr_mx, 1, function(x) as.numeric(x))
  
  chr_mx_2<-cbind(chr_mx,1-as.numeric(beagle[,5])) #Ind 1
  vector_Phet <- seq(from = 5, to = 9999, by = 3)
  var <- chr_mx_2
  for (i in 1:nind){
    
    
    var <- cbind(var,1-as.numeric(beagle[,vector_Phet[i]])) #add remaining Ind
  }  
  chr <- var
} 

#-------------------------------------------------------------------
#For an Individual Chr get a Matrix with the Position and and P(hom) 
#-------------------------------------------------------------------
#Input: Chr Number, Individual (1, 2, etc...) and get_position (beagle file)
#Output: Matrix with Position and PHom
get_chr_ind_mx <- function(chr_number,ind,get_pos_beagle){
  chr <- get_pos_beagle #mx with position on chr and P(hom)
  chr_n <- which(chr[,1]==chr_number) #Chr Number to get 
  
  mySeq <- c(3:99)
  prob_ind <- chr[chr_n,mySeq[ind]]
  positions <- chr[chr_n,2] 
  data <- positions
  data <- c(data, prob_ind)
  chr_n_mx <- matrix(data = c(positions,prob_ind), ncol = 2, nrow = length(positions))
  chr_n_mx
}

#---
#Get a true or false matrix where True = for sites with Phom > Y (Sites
#with P > X but are between 2 sites with Phom > Y are included)
#Example: X = 0.5 and Y = 0.9
#---
#Input: Matrix with Positions for one Chr and Phom
#Output: T/F matrix 
        # 1st Col = True(Phom>50 2nd Col) = True(PHom>90) 
        # 3rd Col = Phom condition of being True stated above 
get_tf_matrix<-function(m, PhomMin, PhomMax){
  v_90 <- m[,2]>=PhomMax
  v_50 <- m[,2]>PhomMin & m[,2]< PhomMax
  v_tf <- v_90 | v_50

  matrix_tf<-matrix(data=c(c(v_50),c(v_90),c(v_tf)),nrow=length(v_90),ncol=3)
  colnames(matrix_tf) <- c('v_50','v_90','v_tf')
  x<-which(matrix_tf[,3] == T) #values of v_tf that are true
  y<-which(matrix_tf[,1] == T) #values of v_50 that are true
  #the ones that dont lay between v_90==T 
  #are to convert from T to F on column 3
  r<-which(matrix_tf[,1] == T & matrix_tf[,2] == F & matrix_tf[,3] == T)
  #See if a v_50 == TRUE lays between two sites of v_90==TRUE 
  x<-c()
  for (i in 1:length(r)-1){
    x<-c(x,matrix_tf[r[i],1] == matrix_tf[r[i]-1,2] & matrix_tf[3,1] == matrix_tf[r[i]+1,2])
  }
  
  #update the matrix col 3
  row_to_change<-c()
  for (i in 1:length(x)){
    if(length(matrix_tf[r[x],3]) > 0){
      if(matrix_tf[r[x],3] != x[i]) row_to_change <- c(row_to_change,r[i])  
    }
  
  }
  matrix_tf[row_to_change,3]<-F
  
  #previous method doesnt consider the last position, so lets correct it
  ifelse (matrix_tf[nrow(matrix_tf),3]==matrix_tf[nrow(matrix_tf),2], matrix_tf[nrow(matrix_tf),3]<-matrix_tf[nrow(matrix_tf),3], matrix_tf[nrow(matrix_tf),3]<-FALSE)
  
  matrix_tf
}

#Get a Matrix with position and 3rd col of the TF Matrix
#Input: TF matrix and output from get_chr_ind_mx
#Output: Matrix with Col 1 = Positions Col2 = Phom(T/F) 
get_tf_positions_mx<-function(matrix_tf,chr_n_mx){
  positions<-chr_n_mx[,1]
  tf<-matrix_tf[,3]
  position_tf <- matrix(data = c(positions,tf), ncol = 2, nrow = length(positions))
  position_tf
}

#Calculate ROH size
#Input Matrix with Position and TF
#Output: Vector with ROH sizes

get_roh <- function(position_tf){
  out <- rle(position_tf[,2] == 1)
  out$lengths[out$values]
  start<-c(0, cumsum(out$lengths))[which(out$values)] + 1
  end<-cumsum(out$lengths)[which(out$values)]
  index_roh<-sort(c(start,end),decreasing=F)
  position_roh<-position_tf[index_roh]
  keep_singles <- function(v){
    v[!(v %in% v[duplicated(v)])] 
  }
  single_positions<-keep_singles(position_roh)
  odd_sites <- single_positions[seq(2, length(single_positions), 2)]
  even_sites <- single_positions[seq(1, length(single_positions), 2)]
  roh <- sort(odd_sites-even_sites,decreasing=F)
  roh_MB <- roh/1000000 #get the ROHs in Mb
  roh_MB
}

calcROH <- function(beagle, PhomMin, PhomMax){
  beagle_1 <- get_position(beagle, nind=1)
  beagle_2 <- get_chr_ind_mx(20,ind = 1, get_pos_beagle = beagle_1)  
  beagle_3 <- get_tf_matrix(beagle_2, PhomMin, PhomMax)
  beagle_4 <- get_tf_positions_mx(matrix_tf = beagle_3, chr_n_mx = beagle_2)
  roh <- get_roh(beagle_4)
  roh
}

calcSNroh <- function(roh, minSize){
  roh_min <- roh[which(roh > minSize)]
  sroh <- sum(roh_min)
  nroh <- length(roh_min)
  out <- c(sroh, nroh)
  out
} 
