set.seed(2)
y <- rpois(5,4)
n <- length(y)

lambda <- 1000
lambda.new <- 3.5
j <- 1
while(abs(lambda.new - lambda) > 0.001){
    lambda <- lambda.new
    lambda.new <- lambda + (lambda)^2/(n*mean(y)) *(n *(mean(y)/lambda - 1))
    print(lambda.new)
    print(j)
    j <- j + 1
}

lambda.new
mean(y)

