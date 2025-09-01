#!/Users/randhal/Documents/Interpreters/xptopics/bin/python
import numpy as np
import time 


# Solution to part 2.1
def matmul(A, B):

	na = len(A)
	ma = len(A[0]) # Matrices must be evenly distributed.
	
	nb = len (B)
	mb = len(B[0])

	# Check that matrices are compatible
	if ma == nb:

		# Create a new matrix to store the product 
		C = [[0 for i in range(mb)] for j in range(na)]
		
		for i in range(na):
			for j in range(mb):
				for k in range(ma):
					C[i][j] = C[i][j] + A[i][k]*B[k][j]

		return C	

	else:

		print("Matrices are not compatible.")


# Solution to part 2.2 
def npmatmul(A, B):
	na = len(A)
	ma = len(A[0]) # Matrices must be evenly distributed.
	
	nb = len (B)
	mb = len(B[0])

	# Check that matrices are compatible
	if ma == nb:
	# Use numpy for compare results. 
		return np.dot(np.array(A),np.array(B))
	else:
		print("Matrices are not compatible.")


if __name__ == "__main__":
	# Define the matrices
	A = [[1, 3, 7], [5, 1, 2], [0, 1, 2]]
	B = [[2, 0, 1], [4, 1, 3], [3, 5, 2]]
	C = [[1, 3, 7], [5, 1, 2]]
	D = [[2], [4], [3]]
	
	# Print AxB
	# This case will be used for the speed test
	
	print("AxB")	
	
	time_start = time.perf_counter()

	# NOTE: np.array is only used to stylize the output.	
	print(np.array(matmul(A, B)))
	print(npmatmul(A, B))
	
	for i in range(1000):
		matmul(A, B)
	time_end  = time.perf_counter()
	print(f"For 1000 iterations the custom algorithm takes: {time_end - time_start}")

	time_start = time.perf_counter()
	for i in range(1000):
		npmatmul(A, B)
	time_end  = time.perf_counter()
	print(f"For 1000 iterations the np.dot algorithm takes: {time_end - time_start}")

	# Print AxC
	print("\nAxC")	
	print(np.array(matmul(A, C)))
	print(npmatmul(A, C))

	# Print CxD
	print("\nCxD")	
	print(np.array(matmul(C, D)))
	print(npmatmul(C, D))

	# Print BxD
	print("\nBxD")	
	print(np.array(matmul(B, D)))
	print(npmatmul(B, D))

	print("""\nUsing NumPy's np.dot function is faster than the custom Python implementation for matrix multiplication. \n
		This is because NumPy uses highly optimized, pre-compiled code, often written in C or Fortran, to perform the \n 
		calculations. In contrast, the custom code uses standard Python loops, which are interpreted line by line. """)
