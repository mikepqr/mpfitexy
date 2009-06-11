;-------------------------------------------------------------------------------
function plineschi2, p, x1 = x1, y1 = y1, x2 = x2, y2 = y2, $
    e_x1 = e_x1, e_y1 = e_y1, e_x2 = e_x2, e_y2 = e_y2, e_int = e_int
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Given two sets of data, returns the total unreduced chi-squared found from
    ; data-model residuals weighted by both error bars and optional intrinsic 
    ; scatter when they are fitted to two lines with a common slope. 
    ;---------------------------------------------------------------------------
    ; INPUTS
    ; x1, y1, x2, y2: independent and dependent variables for two samples
    ;     e_x1, etc.: corresponding error bars
    ;          e_int: intrinsic scatter. 
    ;              p: [slope, intercept, delta], common slope, intercept of 
    ;                 line for sample 1, and difference between intercepts of
    ;                 two samples.
    ;---------------------------------------------------------------------------
    ; OUTPUT
    ; chi-squared of data from models with these data and choice of paramters
    ;---------------------------------------------------------------------------
    slope = p[0]
    intercept = p[1]
    delta = p[2]
    f1 = slope * x1 + intercept
    f2 = slope * x2 + intercept + delta

    if n_elements(e_int) eq 0 then e_int = 0.0
    
    resid1 = (y1 - f1)/sqrt((e_y1^2 + slope^2 * e_x1^2 + e_int^2))
    resid2 = (y2 - f2)/sqrt((e_y2^2 + slope^2 * e_x2^2 + e_int^2))

    return, total([resid1, resid2]^2)
end
;-------------------------------------------------------------------------------
function twolinefit, x1, y1, x2, y2, e_x1, e_y1, e_x2, e_y2, e_int = e_int, $
    guess = guess, p0range = p0range, p1range = p1range, p2range = p2range, $
    chired = chired
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Explore chi-squared in parameter space in extremly crude brute force grid 
    ; fashion to find optimum paramters for Eqn. (3) of Bureau et al. 1996 (i.e.
    ; minimizes the same quantity as MPFITEXY2
    ;---------------------------------------------------------------------------
    ; INPUTS
    ; x1, y1, x2, y2: independent and dependent variables for two samples
    ;     e_x1, etc.: corresponding error bars
    ;          guess: starting points for common slope, intercept of sample 1, 
    ;                 difference between intercepts of the two samples. No
    ;                 default, must be set
    ;  p1range, etc.: distance to explore either side each parameter
    ;          e_int: intrinsic scatter in data. Should be adjusted to ensure
    ;                 sqrt(minchi2/dof) ~= 1.0 for meaningful errors.
    ;---------------------------------------------------------------------------
    ; OUTPUTS
    ;         chired: redued chi of final model (~= 1.0 for meaningful errors)
    ;         return: best parameters of model: slope, intercept, delta
    ;---------------------------------------------------------------------------
    n = 100
    p0 = guess[0] + range(0, p0range*2, n) - p0range
    p1 = guess[1] + range(0, p1range*2, n) - p1range
    p2 = guess[2] + range(0, p2range*2, n) - p2range
    meshgrid3d, p0, p1, p2, pp0, pp1, pp2
    chi2grid = dblarr(n , n , n)
    dof = n_elements(x1) + n_elements(x2) - 3

    for i = 0, n - 1 do begin
        for j = 0, n - 1 do begin
            for k = 0, n - 1 do begin
                chi2grid[i, j, k] = plineschi2([p0[i], p1[j], p2[k]], $
                    x1 = x1, y1 = y1, x2 = x2, y2 = y2, e_x1 = e_x1, $
                    e_y1 = e_y1, e_x2 = e_x2, e_y2 = e_y2, e_int = e_int)
            endfor
        endfor
    endfor

    chired = sqrt(min(chi2grid, minhere)/dof)
    return, [pp0[minhere], pp1[minhere], pp2[minhere]] 
end
