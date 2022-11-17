## Project 4
## Baoyan Deng (s2402728); Yani Shi (s2308441); Qiming Xiong (s2442309)

## Github repository: 

## Work contribution:
## Baoyan Deng (33%): Write the approximate framework of the newt function
## Yani Shi (33%): Write the functions: hessian_matrix, is_positive_definite and 
##                 transfer, write part of comments
## Qiming Xiong (33%): Write all 4 warnings and test cases to debugging, write 
##                     part of comments

## General Description:
## In this project, we write the function called newt, which is used to implement
## Newton's method for minimization of functions. We start from a guess at the 
## parameter theta, and find the quadratic function f, then we generate f(theta),
## the gradient and the second derivative at the particular guess. Then we minimize 
## the quadratic to find an improved guess, and repeat the process, until we arrive
## at the convergence.
## We compute the second derivative as Hessian matrix. And we check whether the 
## Hessian matrix is positive definite or not. If it is not positive definite, we
## need to perturb it. Because positive definiteness ensures that the specific
## direction will decrease the function.
## Although the direction is specified, the distance moved in the direction in 
## each iteration is not. The corresponding function value is constantly 
## calculated to obtain a delta in that direction that makes the function value 
## smaller by constantly decrease the value of delta to its half.
## We separate 6 steps to finish newton optimizer. And it contains three based 
## functions, hessian_matrix(to generate the hessian matrix), is_positive_definite
## (to check whether it is positive definite or not), and transfer(to make matrix
## which is non-positive definite to positive definite).

## Also, the function issues some warnings:
## 1. If the objective is not finite or not define at the initial theta.
## 2. If the derivatives are not finite or not define at the initial theta.
## 3. If the step fails to reduce the objective despite trying max.half step
##    halvings.
## 4. If maxit is reached without convergence.
## 5. If the Hessian is not positive definite at convergence.
## 6. If the function value of new theta is not numeric variable, as non-numeric
##    variable cannot compare with numeric variable
## 7. If the function value of new theta is too large, for example: inf.

library(MASS)
newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,
                 max.half=20,eps=1e-6) {
  ## This function is used to execute the newton method to minimize the function.
  
  ## It includes 10 inputs:
  ## (theta) is a vector of initial values for the optimization parameters.
  
  ## (func) is the objective function to minimize. Its first argument is the vector
  ## of optimization parameters. Remaining arguments will be passed from newt 
  ## using ‘...’.
  
  ## (grad) is the gradient function. It has the same arguments as func but 
  ## returns the gradient vector of the objective w.r.t. the elements of parameter
  ## vector.
  
  ## (hess) is the Hessian matrix function. It has the same arguments as func but 
  ## returns the Hessian matrix of the objective w.r.t. the elements of parameter
  ## vector. If not supplied then newt should obtain an approximation to the
  ## Hessian by finite differencing of the gradient vector. 
  
  ## (...) any arguments of func, grad and hess after the first (the parameter 
  ## vector) are passed using this.
  
  ## (tol) the convergence tolerance.
  
  ## (fscale) a rough estimate of the magnitude of func near the optimum - used 
  ## in convergence testing.
  
  ## (maxit) the maximum number of Newton iterations to try before giving up.
  
  ## (max.half) the maximum number of times a step should be halved before 
  ## concluding that the step has failed to improve the objective.
  
  ## (eps) the finite difference intervals to use when a Hessian function is not 
  ## provided
  
  ## Return a list includes:
  ## (f) the value of the objective function at the minimum.
  
  ## (theta) the value of the parameters at the minimum.
  
  ## (iter) the number of iterations taken to reach the minimum.
  
  ## (g) the gradient vector at the minimum (so the user can judge closeness to 
  ## numerical zero).
  
  ## (Hi) the inverse of the Hessian matrix at the minimum 
  
  ## step 1: we have already get f(theta) and grad(theta) in the parameters of
  ## function.
  
  ## warning 1: If the objective is not finite or not define at the initial theta.
  if (all(is.finite(func(theta))==FALSE) | anyNA(func(theta)) == TRUE){
    stop('Objective is not finite at the initial theta.')
  }
  ## warning 2: If the derivatives are not finite or not define at the initial 
  ## theta.
  if (all(is.finite(grad(theta))==FALSE) | anyNA(grad(theta)) == TRUE) {
    stop('Derivatives is not finite at the initial theta.')
  }
  test_value = max(abs(grad(theta)))
  # store the length of theta, n will be used to transfer hessian matrix 
  # to positive definite
  ## Initial the Newton algorithm
  iteration = 0
  while (iteration <= maxit) {
    ## step 2: Test whether theta is minimum, and terminate if it is.
    ## tol*(abs(f) + fscale), which gives us a value against which to judge 
    ## whether a gradient is zero. 
    ## 'convergence' relates to whether the method has 'converged' to the optimum.
    if (test_value<tol*(abs(func(theta))+fscale)) {
      break
    }
    else{
      ## step 3: Check whether the function of calculating hessian matrix is 
      ## provided in input or not.
      ## If it is not provided, we use hessian_matrix() to generate one.
      if (is.null(hess) == TRUE) {
        hessian = hessian_matrix(theta,grad)
      }
      else {
        hessian = hess(theta)
      }
      ## Check if the matrix is positive definite and transform it if it is not
      ## by calling transfer().
      if (is_positive_definite(hessian) == FALSE) {
        hessian = transfer(hessian)
      }
      
      ##step 4: Get the direction in each iteration
      delta = -chol2inv(chol(hessian))%*%grad(theta) ## Get the direction
      f_theta = func(theta)   ## Get the value of function 
      theta = theta + delta    ## Update the value of theta
      f_theta_del = func(theta)   ## Get the value of function of the new theta
      
      ##step 5: Find the reasonable value to move in the specified direction
      half_time = 0     ## Initial the half time of direction
      flag = TRUE
      while (f_theta_del > f_theta) {
        if (half_time >= max.half) {
          ## if the half time exceed the max.half, stop the halving process.
          flag = FALSE
          break
        }
        delta = delta/2   ## Halve the delta.
        theta = theta-delta    ## Update theta.
        f_theta_del = func(theta)     ## Calculate the new point value
        
        ## New warning: Check whether the new point value is numeric or not and 
        ## finite or not
        if (!is.numeric(f_theta_del)) {
          stop('The function value of new theta is not numeric, cannot compare 
               non-numeric variable with numeric variable.')
        }
        if(all(is.finite(f_theta_del)==FALSE)) {
          stop('The function value of new theta is too large.')
        }
        half_time = half_time + 1
      }
      ## warning 3: If the step fails to reduce the objective despite trying 
      ## max.half step halvings.
      if (flag == FALSE) {
        stop('Fails to reduce the objective despite trying max.half step halvings.')
      }
      ## step 6: Update the new value of theta, Hessian and test value
      ## If the function of calculating hessian matrix is not given, generate it
      ## by hessian_matrix()
      if (is.null(hess) == TRUE) {
        hessian = hessian_matrix(theta,grad)
      }
      ## Hessian matrix is given in the input
      else {
        hessian = hess(theta)
      }
      ## Transfer the hessian matrix to a positive definite matrix if it's not
      if (is_positive_definite(hessian) == FALSE) {
        hessian = transfer(hessian)
      }
      test_value = max(abs(grad(theta)))    ## Get the new test value
      iteration = iteration + 1   ## Update the iteration
    }
  }
  ## warning 4: If maxit is reached without convergence.
  if (test_value>=tol*(abs(func(theta))+fscale)) {
    if (iteration >= maxit) {
      stop('maxit is reached without convergence.')
    }
  }
  else{
    ## Get the return value in the function
    result = list(f=func(theta),theta=theta,iter=iteration,g=grad(theta))
    ## warning 5: If the Hessian is not positive definite at convergence.
    if (is_positive_definite(hessian_matrix(theta,grad)) == FALSE) {
      warning('The heissian is not positive definite at convergence.')
    }
    else {
      ## Compute the inverse hessian matrix
      inverse_hessian = ginv(hessian_matrix(theta,grad))
      result$HI = inverse_hessian
    }
    return(result)
  }
}

hessian_matrix <- function(theta,grad,...,eps=1e-7){
  ## Create the hessian matrix. The inputs are get same from newton function.
  ## (theta) is a vector of initial values for the optimization parameters.
  ## (grad) is the gradient function.
  ## (eps) the finite difference intervals to use when a Hessian function is not 
  ## provided
  ## Return the hessian matrix
  
  n <- length(theta)
  Hfd <- matrix(0,n,n)
  for (i in 1:n) {
    th1 <- theta
    th1[i] <- th1[i]+eps
    grad1 <- grad(th1)
    Hfd[i,] <- (grad1 - grad(theta))/eps
  }
  return(Hfd)
}

is_positive_definite <- function(x) {
  ## Check if the matrix is positive definite by attempting Choleski decomposition
  ## of the matrix, if it's possible then the matrix is a positive definite matrix,
  ## if chol(x) returns error the it's not.
  ## There is one input:
  ## x is the given matrix
  ## If the matrix is positive definite, return True. If it isn't, return False.
  if(inherits(try(chol(x), silent=TRUE), "try-error")) {
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

transfer <- function(x,...) {
  ## Transfer the non-positive definite matrix to positive definite.
  ## There is one input:
  ## x is the given matrix
  ## Return the matrix of x being transferred into a positive definite matrix
  n <- length(theta)
  x_norm <- norm(x)
  i = -6
  while (is_positive_definite(x) == FALSE){
    x <- x + diag(n)*x_norm*10^i
    i = i+1
  }
  return(x)
}

