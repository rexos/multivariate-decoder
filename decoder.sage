## 
# This is a very naive implementation of the MLD to MQ reduction.
# There are a number of improvements that can be implemented in 
# order to make it more efficient, both time- and memory-wise.
#
##

reset()
n = 14
k = 3
l = ceil(log(n)) + 1
m = n-k
P = PolynomialRing(GF(2),n+n*l, var_array=['x'])
P.inject_variables()

C =  codes.random_linear_code(GF(2),n,k)
t = floor((C.minimum_distance()-1)/2)
H = C.parity_check_matrix()

def to_binary(val):
	bits = val.bits()
	return bits + [0 for i in range(l-len(bits))]
	
def gen_error(weight):
	p = Permutations(n).random_element()
	etmp = [1 for i in range(weight)]+[0 for i in range(n-weight)]
	return [etmp[p(i)-1] for i in range(1,n+1)]
print(t)
e = gen_error(t)
t = to_binary(t)
s = H*vector(e)


def parity_check(H,s):
	S = []
	j=0
	for row in H.rows():
		S.append(sum([P(var('x'+str(i))) for i in range(n) if row[i]==1])+s[j])
		j +=1
	return S

def weight_computation():
	S = []
	A = [0 for i in range(n)]
	for i in range(n):
		B = []
		for j in range(l):
			eq = var('x'+str(n+ i*l + j)) + A[j] + prod(A[:j])*var('x'+str(i))
			S.append(eq)
			B.append(var('x'+str(n+ i*l + j)))
		A = copy(B)
	return S
	
def weight_constraint():
	f = []
	for i in range(l):
		f.append((var('x'+str(n + (n-1)*l + i)) + t[i]) * prod([(var('x'+str(n+ (n-1)*l + h)) + t[h] + 1) for h in range(i+1,l)]))
		
	return P(sum([f[i]*(t[i]+1) for i in range(l)]))
	
	
def field_equations():
	return [P(var('x'+str(i))^2 + var('x'+str(i))) for i in range(n+n*l)]



S = parity_check(H,s)
S += weight_computation()
S.append(weight_constraint())
S += field_equations()
