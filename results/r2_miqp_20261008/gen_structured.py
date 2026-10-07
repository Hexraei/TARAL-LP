# Structured R2 families, expressible in the harness tuple format (ni,nc,int_hi,cont_hi,rows,q,c).
# All diagonal convex Q (prototype scope). Seeds 26119000+idx (26119 family). Feasible by construction.
import random

def portfolio(idx):
    # min sum d_j x_j^2 - mu_j x_j + pen sum y : sum x = 1, 0 <= x_j <= y_j, sum y <= k, y binary
    rng = random.Random(26119000+idx)
    n = rng.randint(4,6); k = rng.randint(2, n-1)
    d  = [rng.uniform(0.2,2.0) for _ in range(n)]
    mu = [rng.uniform(0.0,1.5) for _ in range(n)]
    pen = rng.uniform(0.01,0.1)
    ni, nc_ = n, n
    int_hi, cont_hi = [1]*ni, [1.0]*nc_
    c = [-mu[j] for j in range(n)] + [-mu[j]+pen for j in range(n)]  # x then y... objective uses only x quad; y linear pen
    # variables: y_0..y_{n-1} (int) then x_0..x_{n-1} (cont)
    q = [0.0]*ni + d
    c = [pen]*ni + [-mu[j] for j in range(n)]
    rows = []
    # sum x = 1
    rows.append(([0.0]*ni + [1.0]*nc_, 1.0, 1.0))
    # sum y <= k
    rows.append(([1.0]*ni + [0.0]*nc_, -1e6, float(k)))
    # x_j - y_j <= 0  (x_j <= y_j)
    for j in range(n):
        a=[0.0]*(ni+nc_); a[j]=-1.0; a[ni+j]=1.0
        rows.append((a, -1e6, 0.0))
    return ni, nc_, int_hi, cont_hi, rows, q, c

def feature_select(idx):
    # min sum d_j (x_j - t_j)^2 + lam sum y_j  s.t. -M y_j <= x_j <= M y_j, y binary
    # expand: sum d_j x_j^2 - 2 d_j t_j x_j + const ; const dropped (shift invariant for comparison? NO:
    # both sides must drop the same const; objective shift is constant per case so argmin identical,
    # but reported objectives differ by the const. Keep const=0 by absorbing: we cannot add constants
    # in this format, so both proto and ref see the shifted objective identically -> comparable.)
    rng = random.Random(26119000+1000+idx)
    n = rng.randint(4,7); M = rng.choice([2.0,5.0]); lam = rng.uniform(0.05,0.5)
    d = [rng.uniform(0.3,2.0) for _ in range(n)]
    t = [rng.uniform(-1.5,1.5) for _ in range(n)]
    ni, nc_ = n, n
    int_hi, cont_hi = [1]*ni, [M]*nc_
    q = [0.0]*ni + d
    c = [lam]*ni + [-2*d[j]*t[j] for j in range(n)]
    rows = []
    # x_j - M y_j <= 0 ; -x_j - M y_j <= 0  -> |x_j| <= M y_j ; x_j >= -M handled via cont lower? cont bounds are 0..M
    # format has cont lower bound 0; allow negative x by shifting? Instead restrict t_j >= 0 region:
    for j in range(n):
        a=[0.0]*(ni+nc_); a[j]=-M; a[ni+j]=1.0
        rows.append((a, -1e6, 0.0))
    return ni, nc_, int_hi, cont_hi, rows, q, c

def unit_commit(idx):
    # min sum a_i p_i^2 + b_i p_i + f_i u_i : sum p = D, Pmin_i u_i <= p_i <= Pmax_i u_i, u binary
    rng = random.Random(26119000+2000+idx)
    n = rng.randint(3,6)
    a=[rng.uniform(0.1,1.0) for _ in range(n)]; b=[rng.uniform(0.5,2.0) for _ in range(n)]
    f=[rng.uniform(5.0,20.0) for _ in range(n)]
    pmin=[rng.uniform(0.5,2.0) for _ in range(n)]; pmax=[pmin[i]+rng.uniform(1.0,6.0) for i in range(n)]
    D = rng.uniform(0.4*sum(pmin), 0.8*sum(pmax))
    ni, nc_ = n, n
    int_hi, cont_hi = [1]*ni, pmax[:]
    q = [0.0]*ni + a
    c = f[:] + b[:]
    rows = []
    rows.append(([0.0]*ni + [1.0]*nc_, D, D))             # sum p = D
    for i in range(n):                                     # p_i - Pmax_i u_i <= 0
        r=[0.0]*(ni+nc_); r[i]=-pmax[i]; r[ni+i]=1.0
        rows.append((r, -1e6, 0.0))
    for i in range(n):                                     # p_i - Pmin_i u_i >= 0
        r=[0.0]*(ni+nc_); r[i]=-pmin[i]; r[ni+i]=1.0
        rows.append((r, 0.0, 1e6))
    return ni, nc_, int_hi, cont_hi, rows, q, c

FAMS=[portfolio, feature_select, unit_commit]
def gen_structured(idx):
    return FAMS[idx % 3](idx // 3)
