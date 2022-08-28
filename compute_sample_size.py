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
    
    m = int(log(1-p_at_least_1) / log(1-p_positive))
    asymptotic_limit = m + 1
    
    if n <= 1/p_positive or m >= n:
        return {"n" : n, "m" : n, "p_detection" : 1}
    
    else:
    
        upper_bound = compute_upper_bound(n, p_positive, p_at_least_1)
        if inequality_left_side(m, n, p_positive, p_at_least_1) > upper_bound:
            return {"n" : n, "m" : asymptotic_limit, "p_detection" : p_at_least_1}
        
        delta_m = 0
        while inequality_left_side(m, n, p_positive, p_at_least_1) < upper_bound:
            delta_m = 0.5 * step_size(upper_bound, m, n, p_positive, p_at_least_1)
            if delta_m == 0 or -delta_m > 0.5*m:
                break
            else:
                m += delta_m
        
        m = int(m) + 1
        while inequality_left_side(m, n, p_positive, p_at_least_1) > upper_bound and m < asymptotic_limit:
            m += 1
        
        p = 1 - exp(lgamma((1-p_positive)*n + 1) + lgamma(n - m + 1) - lgamma((1-p_positive)*n - m + 1) - lgamma(n + 1))
        p = p if p !=0 and p != 1 else p_at_least_1
        

        return {"n" : n, "m" : m, "p_detection" : p}
