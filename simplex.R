#Simplex method for solving Linear Programms in Standard Equality Form.
#Input: A,b,c defining the Linear Programming problem in standrt equality form,
# x,B: x a basic feasible solution determined by the basis B.
# Output: If there exists an optimal solution: x', B' an optimal basic feasible solution.
# If the LP is unbounded: x,d such that x+ad for a in R+ is a certificate.
# Note: The smallest subscript rule is used and the method in this form can not handel degenerate solutions.

# number of decimal places. Increase round_coef for more accuracy.
round_coef = 10

simplex <- function(A,b,c,x,B){
  n = length(c)
  y = solve(t(A[,B]), c[B])
  s = c - t(y)%*%A
  # if all entries in s are smaller equal 0, we found an optimal solution.
  if (sum(round(s, round_coef) > rep(0, n)) == 0){return(list("optimal", x, y))
  } else{
    for (i in 1:length(s)){
      if (round(s[i], round_coef) > 0){k = i; break}
    }
  }
  d = solve(A[,B], A[,k])
  d_hat = rep(0, n)
  d_hat[B] = - d
  d_hat[k] = 1
  
  # if all entries in d are smaller equal 0, then the Linear Programming Problem is unbounded
  if (sum(round(d_hat, round_coef) < rep(0,n)) == 0){return(list("unbounded", x, d))}
  
  tem = x[B]/d
  idx = which(d > 0)
  min = min(tem[idx])
  min_idx = B[idx[which.min(tem[idx])]]
  
  x_new = x + min*d_hat
  B_new = sort(setdiff(union(B, k), min_idx))
  return(list("new bfs", x_new, B_new))
}

vec = c(1,1,1,1,2,1,1,2,2,1,2,1,2,2,1,1,2,1,4,-1,2,1,-2,4,2,1,0,0)
A = matrix(vec, nrow=4)
b = c(7,5,5,5)
c = c(1,2,2,2,6,4,-4)
x = c(1,1,1,1,0,0,0)
B = c(1,2,3,4)

repeat{
  output = simplex(A,b,c,x,B)
  print(output)
  conclusion = output[1]
  if (conclusion == "new bfs"){
    x = unlist(output[2])
    B = unlist(output[3])
  }
  if (conclusion == "unbounded"){
    print("The LP is unbounded. x and d provide a certifiate.")
    print("x = ")
    print(output[2])
    print("d = ")
    print(output[3])
    break
  }
  if (conclusion == "optimal"){
    print("The LP has an optimal solution x. y is the optimal solution of the dual and therefore provides a certificate.")
    print("x = ")
    print(output[2])
    print("y = ")
    print(output[3])
    break
  }
  
}

