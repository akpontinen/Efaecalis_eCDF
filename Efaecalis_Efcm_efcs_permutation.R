
## Empirical cumulative distribution functions (ECDFs) for 
## relative accessory gene frequencies (0.01-0.99) 
## per each major host type.
## ECDFs are further compared using permutation tests. 
## Test statistic is the maximum difference of the two ECDFs compared.
## NB when permutations are run without seed number, p-values may vary slightly from one run to another

################ Codes for E. faecalis
### FUNCTIONS

as.numeric.character.fun <- function(x) {  
  as.numeric(as.character(x)) }

# as.numeric.factor.fun <- function(x) {
#   as.numeric(levels(x))[x] }


### DATA

# import Panaroo output: gene presence/absencee
#presence_absence <- read.delim("gene_presence_absence.Rtab", header=FALSE, skip = 1)
presence_absence <- read.delim("gene_presence_absence.Rtab", header=FALSE, stringsAsFactors = FALSE)

# import metadata 
meta.host <- read.delim("meta_host.txt", header = FALSE)
colnames(meta.host) <- c("strain.id", "host")

# merge meta info with gene frequencies:
# transpose the presence/absence data such that strains are on rows and genes in columns
presence_absence.t <- t(presence_absence)

# strain.id and gene names as column names
colnames(presence_absence.t) <- c("strain.id", presence_absence.t[1,-1])

#merge host metadata with presence absence
merged.data.host <- merge(meta.host, presence_absence.t[-1,], by.x = "strain.id", by.y = "strain.id")



### RELATIVE FREQUENCIES

## for whole collection
freq.data <- sapply(presence_absence[-1,-1], as.numeric.character.fun)
rel.freq <- rowMeans(freq.data)

# histogram
hist(rel.freq, density = FALSE, breaks = 100, 
     main = 'a. Pan-genome gene frequency distribution \n for the whole E. faecalis collection',
     xlab = "Relative gene frequency")

# only accessory (defined by rel.freq < 0.99), with 1%-99% cutoff
hist(rel.freq[rel.freq > 0.01 & rel.freq < 0.99], density = FALSE, col = "blue", border = "black", breaks = 100, 
     main = "b. Accessory gene frequency distribution \n for the whole E. faecalis collection \n (gene frequencies 0.01-0.99)",
     xlab = "Relative gene frequency")


## for different hosts

# HOSPITALISED PATIENT
patient.data <- subset(merged.data.host, host == "Hospitalised patient")
patient.freq.data <- sapply(patient.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.patient <- colMeans(patient.freq.data)


# WILD ANIMAL (BIRD)
bird.data <- subset(merged.data.host, host == "Wild animal (Bird)")
bird.freq.data <- sapply(bird.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.bird <- colMeans(bird.freq.data)


# NON-HOSPITALISED PERSON
healthy.data <- subset(merged.data.host, host == "Non-hospitalised person")
healthy.freq.data <- sapply(healthy.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.healthy <- colMeans(healthy.freq.data)


# ENVIRONMENT
env.data <- subset(merged.data.host, host == "Environment")
env.freq.data <- sapply(env.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.env <- colMeans(env.freq.data)


# FARM ANIMAL
farm.data <- subset(merged.data.host, host == "Farm animal")
farm.freq.data <- sapply(farm.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.farm <- colMeans(farm.freq.data)


# OLD ISOLATES (HUMAN)
old.data <- subset(merged.data.host, host == "Human (Old isolate)")
old.freq.data <- sapply(old.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.old <- colMeans(old.freq.data)


# FOOD PRODUCT
food.data <- subset(merged.data.host, host == "Food product")
food.freq.data <- sapply(food.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.food <- colMeans(food.freq.data)


### ECDFs

# empirical cumulative distribution function for accessory genome for different host samples
Fn.patient <- ecdf(rel.freq.patient[rel.freq.patient > 0.01 & rel.freq.patient < 0.99])
Fn.bird <- ecdf(rel.freq.bird[rel.freq.bird > 0.01 & rel.freq.bird < 0.99])
Fn.healthy <- ecdf(rel.freq.healthy[rel.freq.healthy > 0.01 & rel.freq.healthy < 0.99])
Fn.env <- ecdf(rel.freq.env[rel.freq.env > 0.01 & rel.freq.env < 0.99])
Fn.farm <- ecdf(rel.freq.farm[rel.freq.farm > 0.01 & rel.freq.farm < 0.99])
Fn.old <- ecdf(rel.freq.old[rel.freq.old > 0.01 & rel.freq.old < 0.99])
Fn.food <- ecdf(rel.freq.food[rel.freq.food > 0.01 & rel.freq.food < 0.99])

# plot
plot(Fn.env, col = "yellowgreen", pch=20, main = "c. Empirical cumulative distribution function (CDF) \n for major E. faecalis hosts \n (gene frequencies 0.01-0.99)", xlab = "Empirical CDF", ylab = "ecdf(x)")
lines(Fn.farm, col = "hotpink", pch=20)
lines(Fn.patient, col = "red2", pch=20)
lines(Fn.healthy, col = "darkturquoise", pch=20)
lines(Fn.bird, col = "mediumblue", pch=20)
legend(0.5, 0.5, legend=c("Environment","Farm animal","Hospitalised patient", "Non-hospitalised person", "Wild bird"), col=c("yellowgreen", "hotpink", "red2", "darkturquoise", "mediumblue"), pch=20, cex=1.0)


#### PERMUTATION TESTS

# hospitalised vs bird

# NB:
# evaluate ecdf for both sets at same points
# we have chosen the points to always be the unique relative frequencies of the smaller group (here bird)
new.x <- sort(unique(rel.freq.bird[rel.freq.bird > 0.01 & rel.freq.bird < 0.99]))
bird.points <- Fn.bird(new.x)
patient.points <- Fn.patient(new.x)

# observed value of the test statistic
# maximum difference between the actual observations, NB absolute values
# (doesn't matter which one is higher at any given point)
T.obs <- max(abs(bird.points - patient.points))

# concatenate the point values for both sets
all <- c(bird.points, patient.points)

#label the bird and patient values as such
group.label <- rep(c("bird","patient"), each = length(new.x))

# permutation loop
n.loop <- 10000
T.perm <- c()
for(i in 1:n.loop){
	perm <- group.label[sample(length(all))]
  bird.perm <- all[which(perm == "bird")]
  patient.perm <- all[which(perm == "patient")]
  bird.perm.sort <- sort(bird.perm)
  patient.perm.sort <- sort(patient.perm)
  T.perm[i] <- max(abs(bird.perm.sort - patient.perm.sort))	
}
# p-value
p <- sum(T.perm >= T.obs)/n.loop

# histogram of the test statistics from the permutations vs the observed test statistic
hist(T.perm, main = "Histogram of permutation test \n Hospitalised patient vs. Bird \n p=0.6928", breaks=100)
abline(v = T.obs, col="red", lwd=2, lty=2)


# Hospitalised vs environment

new.x <- sort(unique(rel.freq.env[rel.freq.env > 0.01 & rel.freq.env < 0.99]))
env.points <- Fn.env(new.x)
patient.points <- Fn.patient(new.x)

# observed value of the test statistic
T.obs <- max(abs(env.points - patient.points))

# concatenate the point values for both sets
all <- c(env.points, patient.points)

#label the environment and patient values as such
group.label <- rep(c("env","patient"), length(new.x))

# permutation loop
n.loop <- 10000
T.perm <- c()
for(i in 1:n.loop){
	perm <- group.label[sample(length(all))]
  env.perm <- all[which(perm == "env")]
  patient.perm <- all[which(perm == "patient")]
  env.perm.sort <- sort(env.perm)
  patient.perm.sort <- sort(patient.perm)
  T.perm[i] <- max(abs(env.perm.sort - patient.perm.sort))	
}

# p-value
p <- sum(T.perm >= T.obs)/n.loop

# histogram of the test statistics from the permutations vs the observed test statistic
hist(T.perm, main = "Histogram of permutation test \n Hospitalised patient vs. Environment \n p=0.2564", breaks=100)
abline(v = T.obs, col="red", lwd=2, lty=2)


# Hospitalised vs healthy

new.x <- sort(unique(rel.freq.healthy[rel.freq.healthy > 0.01 & rel.freq.healthy < 0.99]))
healthy.points <- Fn.healthy(new.x)
patient.points <- Fn.patient(new.x)

# observed value of the test statistic
T.obs <- max(abs(healthy.points - patient.points))

# concatenate the point values for both sets
all <- c(healthy.points, patient.points)

#label healthy and patient values as such
group.label <- rep(c("healthy","patient"), each = length(new.x))

#permutation loop
n.loop <- 10000
T.perm <- c()
for(i in 1:n.loop){
	perm <- group.label[sample(length(all))]
  healthy.perm <- all[which(perm == "healthy")]
  patient.perm <- all[which(perm == "patient")]
  healthy.perm.sort <- sort(healthy.perm)
  patient.perm.sort <- sort(patient.perm)
  T.perm[i] <- max(abs(healthy.perm.sort - patient.perm.sort))	
}

# p-value
p <- sum(T.perm >= T.obs)/n.loop

# histogram of the test statistics from the permutations vs the observed test statistic
hist(T.perm, main = "Histogram of permutation test \n Hospitalised patient vs. Non-Hospitalised \n p=0.2008", breaks=100)
abline(v = T.obs, col="red", lwd=2, lty=2)


# Hospitalised & farm animal

new.x <- sort(unique(rel.freq.farm[rel.freq.farm > 0.01 & rel.freq.farm < 0.99]))
farm.points <- Fn.farm(new.x)
patient.points <- Fn.patient(new.x)

# observed value of the test statistic
T.obs <- max(abs(farm.points - patient.points))

# concatenate the point values for both sets
all <- c(farm.points, patient.points)

#label the farm and patient values as such
group.label <- rep(c("farm","patient"), each = length(new.x))

# permutation loop
n.loop <- 10000
T.perm <- c()
for(i in 1:n.loop){
	perm <- group.label[sample(length(all))]
farm.perm <- all[which(perm == "farm")]
patient.perm <- all[which(perm == "patient")]
farm.perm.sort <- sort(farm.perm)
patient.perm.sort <- sort(patient.perm)
T.perm[i] <- max(abs(farm.perm.sort - patient.perm.sort))	
}

# p-value
p <- sum(T.perm >= T.obs)/n.loop

# histogram of the test statistics from the permutations vs the observed test statistic
hist(T.perm, main = "Histogram of permutation test \n Hospitalised patient vs. Farm animal \n p=0.4787", breaks=100)


# Hospitalised vs food product

new.x <- sort(unique(rel.freq.food[rel.freq.food > 0.01 & rel.freq.food < 0.99]))
food.points <- Fn.food(new.x)
patient.points <- Fn.patient(new.x)

# observed value of the test statistic
T.obs <- max(abs(food.points - patient.points))

# concatenate the point values for both sets
all <- c(food.points, patient.points)

# label the food and patient values as such
group.label <- rep(c("food","patient"), each = length(new.x))

# permutation loop
n.loop <- 10000
T.perm <- c()
for(i in 1:n.loop){
	perm <- group.label[sample(length(all))]
  food.perm <- all[which(perm == "food")]
  patient.perm <- all[which(perm == "patient")]
  food.perm.sort <- sort(food.perm)
  patient.perm.sort <- sort(patient.perm)
  T.perm[i] <- max(abs(food.perm.sort - patient.perm.sort))	
}

# p-value
p <- sum(T.perm >= T.obs)/n.loop

# histogram of the test statistics from the permutations vs the observed test statistic
hist(T.perm, main = "Histogram of permutation test \n Hospitalised patient vs. Food product \n p=0.7748", breaks=100)
abline(v = T.obs, col="red", lwd=2, lty=2)




################ Codes for E. faecium

### FUNCTIONS

as.numeric.character.fun <- function(x) {  
  as.numeric(as.character(x)) }

# as.numeric.factor.fun <- function(x) {
#   as.numeric(levels(x))[x] }


### DATA

# import Panaroo output: gene presence/absencee
#presence_absence <- read.delim("gene_presence_absence.Rtab", header=FALSE, skip = 1)
presence_absence <- read.delim("gene_presence_absence_Efcm.Rtab", header=FALSE, stringsAsFactors = FALSE)

# import metadata 
meta.host <- read.delim("meta_host_Efcm.txt", header = FALSE)
colnames(meta.host) <- c("strain.id", "host")

# merge meta info with gene frequencies:
# transpose the presence/absence data such that strains are on rows and genes in columns
presence_absence.t <- t(presence_absence)

# strain.id and gene names as column names
colnames(presence_absence.t) <- c("strain.id", presence_absence.t[1,-1])

#merge host metadata with presence absence
merged.data.host <- merge(meta.host, presence_absence.t[-1,], by.x = "strain.id", by.y = "strain.id")



### RELATIVE FREQUENCIES

## for whole collection
freq.data <- sapply(presence_absence[-1,-1], as.numeric.character.fun)
rel.freq <- rowMeans(freq.data)

# histogram
hist(rel.freq, density = FALSE, breaks = 100, 
     main = 'a. Pan-genome gene frequency distribution \n for the E. faecium collection',
     xlab = "Relative gene frequency")

# only accessory (defined by rel.freq < 0.99), with 1%-99% cutoff
hist(rel.freq[rel.freq > 0.01 & rel.freq < 0.99], density = FALSE, col = "blue", border = "black", breaks = 100, 
     main = "b. Accessory gene frequency distribution \n for the E. faecium collection \n (gene frequencies 0.01-0.99)",
     xlab = "Relative gene frequency")


## for different hosts

# HOSPITALISED PATIENT
patient.data <- subset(merged.data.host, host == "Hospitalised patient")
patient.freq.data <- sapply(patient.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.patient <- colMeans(patient.freq.data)


# NON-HOSPITALISED PERSON
healthy.data <- subset(merged.data.host, host == "Non-hospitalised person")
healthy.freq.data <- sapply(healthy.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.healthy <- colMeans(healthy.freq.data)


# OTHERS
other.data <- subset(merged.data.host, host == "Other")
other.freq.data <- sapply(other.data[,-(1:2)], as.numeric.character.fun)
# note, now genes are in columns as we transposed the data for merging
rel.freq.other <- colMeans(other.freq.data)


### ECDFs

# empirical cumulative distribution function for accessory genome for different host samples
Fn.patient <- ecdf(rel.freq.patient[rel.freq.patient > 0.01 & rel.freq.patient < 0.99])
Fn.healthy <- ecdf(rel.freq.healthy[rel.freq.healthy > 0.01 & rel.freq.healthy < 0.99])
Fn.other <- ecdf(rel.freq.other[rel.freq.other > 0.01 & rel.freq.other < 0.99])

# plot

plot(Fn.patient, col = "red2", main = "c. Empirical cumulative distribution function (CDF) \n for E. faecium hosts \n (gene frequencies 0.01-0.99)", xlab = "Empirical CDF", ylab = "ecdf(x)")
lines(Fn.other, col = "gray")
lines(Fn.healthy, col = "darkturquoise")
legend(0.4, 0.4, legend=c("Hospitalised patient", "Non-hospitalised person", "Others"), col=c("red","darkturquoise", "gray"), pch=20, cex=1.0)


#### PERMUTATION TESTS

# Hospitalised vs healthy

new.x <- sort(unique(rel.freq.healthy[rel.freq.healthy > 0.01 & rel.freq.healthy < 0.99]))
healthy.points <- Fn.healthy(new.x)
patient.points <- Fn.patient(new.x)

# observed value of the test statistic
T.obs <- max(abs(healthy.points - patient.points))

# concatenate the point values for both sets
all <- c(healthy.points, patient.points)

#label healthy and patient values as such
group.label <- rep(c("healthy","patient"), each = length(new.x))

#permutation loop
n.loop <- 10000
T.perm <- c()
for(i in 1:n.loop){
	perm <- group.label[sample(length(all))]
  healthy.perm <- all[which(perm == "healthy")]
  patient.perm <- all[which(perm == "patient")]
  healthy.perm.sort <- sort(healthy.perm)
  patient.perm.sort <- sort(patient.perm)
  T.perm[i] <- max(abs(healthy.perm.sort - patient.perm.sort))	
}

# p-value
p <- sum(T.perm >= T.obs)/n.loop

# histogram of the test statistics from the permutations vs the observed test statistic
hist(T.perm, main = "Histogram of permutation test \n Hospitalised patient vs. Non-Hospitalised \n p=0.0187", breaks=100)
abline(v = T.obs, col="red", lwd=2, lty=2)


# Hospitalised vs other product

new.x <- sort(unique(rel.freq.other[rel.freq.other > 0.01 & rel.freq.other < 0.99]))
other.points <- Fn.other(new.x)
patient.points <- Fn.patient(new.x)

# observed value of the test statistic
T.obs <- max(abs(other.points - patient.points))

# concatenate the point values for both sets
all <- c(other.points, patient.points)

# label the other and patient values as such
group.label <- rep(c("other","patient"), each = length(new.x))

# permutation loop
n.loop <- 10000
T.perm <- c()
for(i in 1:n.loop){
	perm <- group.label[sample(length(all))]
  other.perm <- all[which(perm == "other")]
  patient.perm <- all[which(perm == "patient")]
  other.perm.sort <- sort(other.perm)
  patient.perm.sort <- sort(patient.perm)
  T.perm[i] <- max(abs(other.perm.sort - patient.perm.sort))	
}

# p-value
p <- sum(T.perm >= T.obs)/n.loop

# histogram of the test statistics from the permutations vs the observed test statistic
hist(T.perm, main = "Histogram of permutation test \n Hospitalised patient vs. other product \n p=0.0047", breaks=100)
abline(v = T.obs, col="red", lwd=2, lty=2)