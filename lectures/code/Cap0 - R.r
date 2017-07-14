

# Presentación


1 + 1 		# Aritmetica
2 + 3 *4 		# Orden
3^2			# Exponencial
exp(0.5)		# Funciones pre establecidas
sqrt(9)
pi			# Constantes pre establecidas
2*pi


x <- 2
x
y <- 3
z <- 4
x + y*z
X + Y*Z

estado.civil <- "casado"
estado.civil

a <- 2 > 3
a



x <- c(1,2,3,4,5)
x
y <- rep(4,5)
y
z <- seq(4,5,0.1)
z
x^2
x+y
sqrt(x)

z <- seq(4,5,0.2)
z
z[1]
z[2:4]
z[3] <- 20
z
z < 4.5
z[-1]

dat <- read.csv(file.choose())

dat$t1[1:6]
dat$t2[1:6] + dat$t1[1:6]
dat[1,]
dat[1:4,1:3]


fun1 <- function(x){
 x*x + x
}

fun1(3)
fun1(c(1,2,3,4))

fun2 <- function(x){
 sum(x)
}

fun2(c(1,2,3,4))

fun3 <- function(a,b){
  aux <- sum(a)
  z   <- b + aux
  return(z)
}

x <- seq(0,10,1)
y <- 11
h <- c(1,2,3)

fun3(x,y)
fun3(x,h)


fun  <- function (x, a) (x - a)^2
xmin <- optimize(fun, c(0, 1), a = 1/3)
xmin


fr <- function(x) {   
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}

res <- optim(c(-1.2,1), fr)
res

