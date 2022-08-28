from math import lgamma
from math import log
from math import exp

def compute_upper_bound(n, p_positive, p_at_least_1):
    return log(1 - p_at_least_1) - lgamma((1-p_positive)*n + 1) + lgamma(n + 1)
    
    
def inequality_left_side(m, n, p_positive, p_at_least_1):
    return lgamma(n - m + 1) - lgamma((1-p_positive)*n - m + 1)
    
    
def dgamma_approx(m, n, p_positive, p_at_least_1):
    return -log(n - m + 1) + log((1-p_positive)*n - m + 1) + 1/(n - m + 1) - 1/((1-p_positive)*n - m + 1)
    
    
def step_size(upper_bound, m, n, p_positive, p_at_least_1):
    change_needed = upper_bound - inequality_left_side(m, n, p_positive, p_at_least_1)
    dy_over_dx = dgamma_approx(m, n, p_positive, p_at_least_1)
    return change_needed / dy_over_dx 
  
  
def compute_sample_size(n, p_positive, p_at_least_1):
    
    asymptotic_limit = int(log(1-p_at_least_1) / log(1-p_positive)) + 1
    m = min([asymptotic_limit, int((1-p_positive)*n)])
    
    if n <= int(1/p_positive) + 1 or m >= n or n <= asymptotic_limit:
        return {"n" : n, "m" : n, "p_detection" : 1}
    
    else:
    
        upper_bound = compute_upper_bound(n, p_positive, p_at_least_1)    
        delta_m = 0
        while inequality_left_side(m, n, p_positive, p_at_least_1) < upper_bound and m <= n:
            if (1-p_positive)*n - m + 1 <= 0:
                break
            delta_m = 0.5*(upper_bound - inequality_left_side(m, n, p_positive, p_at_least_1)) / dgamma_approx(m, n, p_positive, p_at_least_1)
            if delta_m == 0:
                break
            else:
                m += delta_m

        
        m = min([int(m) + 1, n])
        while inequality_left_side(m, n, p_positive, p_at_least_1) > upper_bound and m < asymptotic_limit and m < n:
            m += 1
        
        p = 1 - exp(lgamma((1-p_positive)*n + 1) + lgamma(n - m + 1) - lgamma((1-p_positive)*n - m + 1) - lgamma(n + 1))
        p = p if p >= p_at_least_1 and p < 1 else 1 if m == n else p_at_least_1
        
        return {"n" : n, "m" : m, "p_detection" : p}
