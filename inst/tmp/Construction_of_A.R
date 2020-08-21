# Relationship between L and A matrices.


K = 4

L = cbind(rep(-1,K-1),diag(1,K-1))

A.0 = kronecker(L,L)

# Covariance matrix Simga_K_curl is symmetric. There are duplicate entries corresponding to the covariances. We need to remove the corresponding rows to reduce the dimension of the system of equations. So the row dimension corresponds to the number of quantities available from the estimator of the singular matrix. 

Index.matrix0 = matrix(seq(1:(K-1)^2),c(K-1,K-1))
Index.vec0 = 0
for (i in 1:K-1){
	for (j in 1:K-1){
		if (i < j){Index.vec0 = c(Index.vec0,Index.matrix0[i,j])}
	}
}

Index.vec0 = Index.vec0[-1]

A.1 = A.0[-Index.vec0,]

# Now we use the symmetry in Sigma_K to reduce the column dimension. We need to add the corresponding columns. For example, [1,2] and [2,1] are identical so the corresponding columns should be added and one of them need to removed.

Index.matrix1 = matrix(seq(1:K^2),c(K,K))
Index.mat1 = c(0,0)
for (i in 1:K){
	for (j in 1:K){
		if (i < j){Index.mat1 = rbind(Index.mat1,c(Index.matrix1[i,j],Index.matrix1[j,i]))}
	}
}

Index.mat1 = Index.mat1[-1,]

A.2 = A.1
A.2[,Index.mat1[,2]] = A.1[,Index.mat1[,2]] + A.1[,Index.mat1[,1]]
A.2 = A.2[,-Index.mat1[,1]]

# Now remove the columns that correspond to the structural 0s in the model covariance matrix. The indices for the model covariance matrix should go from 1 to K(K+1)/2. This should equal the number of columns in A.2

check = ncol(A.2) - K*(K+1)/2
if (check == 0){cat("Good to go")}


# Indices for the variances/covariances that are non-zero. These are the indices in the vectorized form of the lower triangle including the diagonal form of Sigma_K.

Index.nonzero = which(Sigma_K[lower.tri(Sigma_K,diag=T)] != "NA")   # This changes as per the model for covariance.  


A = A.2[,Index.nonzero]

# Number of columns of A should be smaller than or equal to the number of rows. Furthermore, svd(A) should give positive singular values (eigenvalues). This assures that t(A)%*%A is invertible, of the same rank as the number of unknowns. The solution is unique. 
# Sigma_K_curl = L %*% Sigma_K_star %*% t(L)  # This will be a (K-1) by (K-1) symmetric matrix. We only need the upper diagonal of this matrix (including the diagonals) in the appropriate order.

# Sigma_hat = solve(t(A) %*% A) %*% t(A) %*% vec.upper(Sigma_K_curl)
