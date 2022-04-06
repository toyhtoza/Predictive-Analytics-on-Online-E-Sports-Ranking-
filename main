#Messay Method

erg3020 <- read.csv("erg3020_data.csv")
names(erg3020)
Diff_KDA <- erg3020$WKDA - erg3020$LKDA
erg3020_data <- cbind(erg3020,Diff_KDA)
B <- as.matrix(erg3020[,1:16])
Bt <- t(B)
BtB <- Bt %*% B
new_BtB <- rbind(BtB[-16,],rep(1,16))
v <- as.matrix(erg3020[,17])
Btv <- Bt %*% v
new_Btv <- rbind(as.matrix(Btv[-16,1]),c(0))
#r * BtB = Btv
r <- solve(new_BtB) %*% new_Btv
rank(-r[,1]) 


B <- as.matrix(erg3020[,1:16])
Bt <- t(B)
BtB <- Bt %*% B
new_BtB <- rbind(BtB[-16,],rep(1,16))
v <- as.matrix(erg3020_data[,20])
Btv <- Bt %*% v
new_Btv <- rbind(as.matrix(Btv[-16,1]),c(0))
#r * BtB = Btv
r <- solve(new_BtB) %*% new_Btv
rank(-r[,1]) 


Messay Method (n=100)
all_matrix <- matrix(c(16:1),nrow=1) 
colnames(all_matrix) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")

k <- 1 

while (k <= 100){
  game_initial <- matrix(rep(0,1920),nrow = 120,ncol = 16)
  colnames(game_initial) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
  n <- 1
  for (i in 1:16){
    for (j in 1:16){
      if (i < j){
        win_prob <- j/(i+j)
        lose_prob <- i/(i+j)
        result <- sample(c(1,-1),1,prob = c(win_prob,lose_prob))
        game_initial[n,j] <- result
        game_initial[n,i] <- 0-result
        n <- n+1
      }
    }
  }
  v_random <- matrix(sample(c(1,2),120,replace = T),ncol = 1)
  B <- as.matrix(game_initial[,1:16])
  Bt <- t(B)
  BtB <- Bt %*% B
  new_BtB <- rbind(BtB[-16,],rep(1,16))
  v <- v_random
  Btv <- Bt %*% v
  new_Btv <- rbind(as.matrix(Btv[-16,1]),c(0))
  #r * BtB = Btv
  r <- solve(new_BtB) %*% new_Btv
  all_matrix <- rbind(all_matrix,rank(-r[,1]))
  k <- k+1
}
#all_matrix <- all_matrix[-1,]

average_matrix <- matrix(c(16:1),nrow=1)
colnames(average_matrix) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
for (i in 1:16){
  average_matrix[1,i] <- (sum(all_matrix[-1,i])/100)-all_matrix[1,i]
}
sum(abs(average_matrix))


#Bradley Terry 
data1<-read.csv("BT_data1.csv", header=T)
row.names(data1) = data1[,1]
data1 = data1[,-1]

library(CVXR)
sita = Variable(16)
function1 = 0

for (i in 1:16){
  for (j in 1:16){
    part1 =  data1[i,j]*sita[i]
    sel = matrix(rep(0,16))
    sel[i] = 1
    sel[j] = 1
    part2 = log_sum_exp(sel*sita)
    part3 = data1[i,j]*part2
    part= part3-part1
    function1 = function1+part
  }
}

constraint1 = list(sum(exp(sita))<=1)
problem1 <- Problem(Minimize(function1),constraints = constraint1)
result1 <- solve(problem1)
rank1 = result1$getValue(sita) 
row.names(rank1) = row.names(data1)  #??
#rank1[order(rank1[,1],decreasing = T),]
rank(-rank1[,1]) 


all_matrix <- matrix(c(16:1),nrow=1)
colnames(all_matrix) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
k <- 1
while (k <= 2){
  game_initial <- matrix(rep(0,256),nrow=16)
  colnames(game_initial) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
  rownames(game_initial) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
  for (i in 1:16){
    for (j in 1:16){
      if (i<j){
        win_prob <- j/(i+j)
        lose_prob <- i/(i+j)
        result <- sample(c(1,-1),1,prob = c(win_prob,lose_prob))
        if (result == 1){
          game_initial[i,j] <- sample(c(1,2),1)
        }
        else{
          game_initial[j,i] <- sample(c(1,2),1)
        }
      }
    }
  }
  data1 <- game_initial
  library(CVXR)
  sita = Variable(16)
  function1 = 0
  
  for (i in 1:16){
    for (j in 1:16){
      part1 =  data1[i,j]*sita[i]
      sel = matrix(rep(0,16))
      sel[i] = 1
      sel[j] = 1
      part2 = log_sum_exp(sel*sita)
      part3 = data1[i,j]*part2
      part= part3-part1
      function1 = function1+part
    }
  }
  constraint1 = list(sum(exp(sita))<=1)
  problem1 <- Problem(Minimize(function1),constraints = constraint1)
  result1 <- solve(problem1)
  rank1 = result1$getValue(sita) 
  row.names(rank1) = row.names(data1)  #??
  #rank1[order(rank1[,1],decreasing = T),]
  #rank(-rank1[,1])
  all_matrix <- rbind(all_matrix,rank(-rank1[,1]))
  k <- k+1
}
#all_matrix <- all_matrix[-1,]
average_matrix <- matrix(c(16:1),nrow=1)
colnames(average_matrix) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
for (i in 1:16){
  average_matrix[1,i] <- (sum(all_matrix[-1,i])/2)-all_matrix[1,i]
}
sum(abs(average_matrix))


#elo algorithm
erg3020 <- read.csv("erg3020_data.csv")
names(erg3020)
elo_data <- erg3020[,1:16]
initial_matrix <- matrix(rep(1500,16),nrow=1)
colnames(initial_matrix) <- names(erg3020)[1:16]
for (i in 1:120){
  for (j in 1:16){
    if (elo_data[i,j]==1){
      win_location <- j
    }
    if (elo_data[i,j]==-1){
      lose_location <- j
    }
  }
  E_win <- 1/(1+10^((initial_matrix[lose_location]-initial_matrix[win_location])/400))
  E_lose <- 1/(1+10^((initial_matrix[win_location]-initial_matrix[lose_location])/400))
  point_win <- initial_matrix[win_location] + 32*(1-E_win)
  point_lose <- initial_matrix[lose_location] + 32*(0-E_lose)
  initial_matrix[win_location] <- point_win
  initial_matrix[lose_location] <- point_lose
}
rank(-initial_matrix)
rbind(rank(initial_matrix),rank(initial_matrix))

#ģ?????? Elo (n=100)
all_matrix <- matrix(c(16:1),nrow=1)
colnames(all_matrix) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
k <- 1
while (k <= 100){
  game_initial <- matrix(rep(0,1920),nrow = 120,ncol = 16)
  colnames(game_initial) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
  n <- 1
  for (i in 1:16){
    for (j in 1:16){
      if (i < j){
        win_prob <- j/(i+j)
        lose_prob <- i/(i+j)
        result <- sample(c(1,-1),1,prob = c(win_prob,lose_prob))
        game_initial[n,j] <- result
        game_initial[n,i] <- 0-result
        n <- n+1
      }
    }
  }
  initial_matrix <- matrix(rep(1500,16),nrow=1)
  colnames(initial_matrix) <- names(game_initial)[1:16]
  for (i in 1:120){
    for (j in 1:16){
      if (game_initial[i,j]==1){
        win_location <- j
      }
      if (game_initial[i,j]==-1){
        lose_location <- j
      }
    }
    E_win <- 1/(1+10^((initial_matrix[lose_location]-initial_matrix[win_location])/400))
    E_lose <- 1/(1+10^((initial_matrix[win_location]-initial_matrix[lose_location])/400))
    point_win <- initial_matrix[win_location] + 32*(1-E_win)
    point_lose <- initial_matrix[lose_location] + 32*(0-E_lose)
    initial_matrix[win_location] <- point_win
    initial_matrix[lose_location] <- point_lose
  }
  all_matrix <- rbind(all_matrix,rank(-initial_matrix))
  k <- k+1
}
#all_matrix <- all_matrix[-1,]
average_matrix <- matrix(c(16:1),nrow=1)
colnames(average_matrix) <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16")
for (i in 1:16){
  average_matrix[1,i] <- (sum(all_matrix[-1,i])/100)-all_matrix[1,i]
}
sum(abs(average_matrix))




