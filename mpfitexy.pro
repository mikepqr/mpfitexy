;-------------------------------------------------------------------------------
function lineresid, p, x = x, y = y, e_x = e_x, e_y = e_y, e_int = e_int
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Utility function called by mpfitexy. Given a set of data, returns the
    ; residuals weighted by both error bars and optional intrinsic scatter 
    ; when fitted with a straight line
    ;---------------------------------------------------------------------------
    ; INPUTS
    ;           x, y: independent and dependent variables
    ;       e_x, e_y: corresponding error bars
    ;          e_int: intrinsic scatter
    ;              p: [slope, intercept]
    ;---------------------------------------------------------------------------
    ; OUTPUT
    ; Residual of data from models with these data and choice of paramters
    ;---------------------------------------------------------------------------
    slope = p[0]
    intercept = p[1]
    f = slope * x + intercept
    if n_elements(e_int) eq 0 then e_int = 0.0
    resid = (y - f)/sqrt((e_y^2 + slope^2*e_x^2 + e_int^2))
    return, resid
end
;-------------------------------------------------------------------------------
function mpfitexy, x, y, e_x, e_y, fixslope = fixslope, errors = perror, $
    guess = guess, minchi2 = minchi2, dof = dof, quiet = quiet, e_int = e_int, $
    fixint = fixint, x0 = x0, reduce = reduce, inv = inv, silent = silent
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Uses MPFIT to determine straight line fit to data with errors in both
    ; variables using the weighting defined in Numerical Recipes, i.e. emulates
    ; the ASTROIDL library function FITEXY.
    ;---------------------------------------------------------------------------
    ; INPUTS
    ;           x, y: independent and dependent variables
    ;       e_x, e_y: corresponding error bars
    ;             x0: fit the relationship y = a(x - x0) + b
    ;          guess: starting points for slope, intercept Slope can be fixed 
    ;                 (see below). Default = [0.0, 0.0, 0.0]
    ;      /fixslope: fix the slope to guess[0] 
    ;         /quiet: Suppress MPFIT's text output
    ;        /silent: Do not print MPFIT status code (see mpfit.pro for docs)
    ;          e_int: intrinsic scatter in data. Should be adjusted to ensure
    ;                 sqrt(minchi2/dof) ~= 1.0
    ;        /fixint: fix the intercept to guess[0]
    ;        /reduce: adjust intrinsic scatter e_int to ensure chired ~= 1.0
    ;           /inv: fits the inverse relation x = a' y + b' (still returns
    ;                 slope of forward relation, i.e. slope = 1/a' and
    ;                 intercept = -b'/a' 
    ;---------------------------------------------------------------------------
    ; OUTPUTS
    ;         errors: 1 sigma fitting errors in paramters slope and intercept.
    ;                 Dubious if sqrt(minchi2/dof) != 1.0
    ;        minchi2: chi-squared of final model
    ;            dof: degrees of freedom
    ;         return: best parameters of model: slope, intercept
    ;---------------------------------------------------------------------------

    ;---------------------------------------------------------------------------
    ; DEFAULTS
    ;---------------------------------------------------------------------------
    if n_elements(e_int) eq 0 then e_int = 0.d
    if n_elements(guess) eq 0 then guess_ = [1.d, 1.d] else $
        guess_ = double(guess)
    if n_elements(x0) eq 0 then x0 = 0.d
    if keyword_set(reduce) then e_int = 0.1d

    ;---------------------------------------------------------------------------
    ; RESCALE X-COORDS TO X0
    ;---------------------------------------------------------------------------
    x_ = double(x) - x0
    y_ = double(y)
    e_x_ = double(e_x)
    e_y_ = double(e_y)

    ;---------------------------------------------------------------------------
    ; SWAP X AND Y AND INVERT GUESS IF FITTING INVERSE FUNCTION
    ;---------------------------------------------------------------------------
    if keyword_set(inv) then begin
        xtemp = x_
        e_xtemp = e_x_
        x_ = y_
        e_x_ = e_y_
        y_ = xtemp
        e_y_ = e_xtemp
        forwardslopeguess = guess_[0]
        forwardinterceptguess = guess_[1]
        guess_[0] = 1 / forwardslopeguess
        guess_[1] = - forwardinterceptguess/forwardslopeguess
    endif

    ;---------------------------------------------------------------------------
    ; FIX SLOPE/LABEL PARAMETERS
    ;---------------------------------------------------------------------------
    pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], parname:''},2)
    if keyword_set(fixslope) then pi(0).fixed = 1
    if keyword_set(fixint) then pi(1).fixed = 1
    pi(0).parname = '    Slope'
    pi(1).parname = 'Intercept'

    ;---------------------------------------------------------------------------
    ; CHECK ALL VALUES FINITE
    ;---------------------------------------------------------------------------
    ok = where(finite(x_) and finite(y_))

    ;---------------------------------------------------------------------------
    ; CALL MPFIT ONCE
    ;---------------------------------------------------------------------------
    result = mpfit('lineresid', guess_, functargs = {x:x_[ok], y:y_[ok], $
        e_x:e_x_[ok], e_y:e_y_[ok], e_int:e_int}, parinfo = pi, $
        status = status, errmsg = errmsg, bestnorm = minchi2, dof = dof, $
        perror = perror, quiet = quiet, covar = covar)
    chired = sqrt(minchi2/dof)
    if ~keyword_set(quiet) then print, chired, e_int
    if ~keyword_set(silent) and status ne 1 then print, status
    ;---------------------------------------------------------------------------
    ; CALL MPFIT UNTIL REDUCED CHI2 ~= 1.0 IF REQUIRED
    ;---------------------------------------------------------------------------
    if keyword_set(reduce) then begin
        while abs(chired - 1.) gt 0.01 do begin
            e_int = e_int * chired^(4./3)
            result = mpfit('lineresid', guess_, functargs = {x:x_[ok], y:y_[ok], $
                e_x:e_x_[ok], e_y:e_y_[ok], e_int:e_int}, parinfo = pi, $
                status = status, errmsg = errmsg, bestnorm = minchi2, dof = dof, $
                perror = perror, quiet = quiet, covar = covar)
            chired = sqrt(minchi2/dof)
            if ~keyword_set(quiet) then print, chired, e_int
            if ~keyword_set(silent) and status ne 1 then print, status
        endwhile
    endif

    ;---------------------------------------------------------------------------
    ; FLIP BEST-FITTING PARAMS AND THEIR ERRORS IF FITTING INVERSE FUNCTION
    ; See PhD notes VIII.3 for derivation of error propagation
    ;---------------------------------------------------------------------------
    if keyword_set(inv) then begin
        forwardslope = result[0]
        forwardintercept = result[1]
        forwardslopeerror = perror[0]
        forwardintercepterror = perror[1]
        result[0] = 1 / forwardslope
        result[1] = - forwardintercept/forwardslope
        perror[0] = forwardslopeerror / forwardslope^2
        perror[1] = sqrt(forwardintercepterror^2/forwardslope^2 + $
            forwardslopeerror^2*forwardintercept^2/forwardslope^4 - $
            2 * covar[0,1] * forwardintercept/forwardslope^3)
        e_int = abs(e_int * result[0])
    endif
    return, result
end

;-------------------------------------------------------------------------------
pro testmpfitexy, shuffle = shuffle
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Code to test mpfitexy2 and mpfitxy
    ; Result should be a plot of two lines which overlap. Blue line is fit from
    ; mpfitexy, dashed red line from fitexy. Outputs fit parameters, errors 
    ; in these and sqrt(chi-squared/DOF) for both methods, which should all 
    ; match
    ;---------------------------------------------------------------------------
    
    ;---------------------------------------------------------------------------
    ; SETUP
    ;---------------------------------------------------------------------------
    x = [2.43184, 2.34696, 2.21818, 2.23748, 2.16919, 2.12780, 2.15438, $
        2.24805, 2.29982, 2.14897, 2.17739, 2.22395, 2.12592, 2.44776, $
        2.19388, 2.37343, 2.39492, 2.23249, 2.17179, 2.40253, 2.38844, $
        2.37377, 2.22880, 2.21212, 2.26511, 2.32707, 2.33750, 2.15525]
    y = [-25.3406, -24.8878, -22.6971, -22.8706, -22.7128, -22.4252, -22.9278, $
        -23.8840, -24.7333, -23.0504, -23.5204, -24.0678, -23.0894, -25.5364, $
        -23.5926, -25.3422, -25.0512, -24.2071, -24.1777, -25.1833, -24.8495, $
        -24.5180, -24.4620, -23.2211, -25.2567, -24.8611, -25.2406, -23.2333]
    e_x = [0.00571056, 0.0249614, 0.00348239, 0.0192876, 0.0228633, $
        0.00863316, 0.0357660, 0.00380802, 0.00545385, 0.00564655, $
        0.00300506, 0.0340787, 0.0335453, 0.00403909, 0.0522591, $
        0.0109578, 0.00863588, 0.0209927, 0.0510806, 0.0656906, 0.0100137, $
        0.00491871, 0.0204541, 0.00346392, 0.0107160, 0.00327621, $
        0.00947749, 0.0233549]
    e_y = [0.0434049, 0.0429708, 0.208040, 0.159271, 0.0503973, 0.103848, $
        0.0504721, 0.0643180, 0.0422348, 0.0273715, 0.0239415, 0.0999366, $
        0.0457358, 0.0461065, 0.0563985, 0.0450141, 0.0479903, 0.0541895, $
        0.0601011, 0.0521960, 0.0512739, 0.0498970, 0.0448522, 0.0457561, $
        0.0708393, 0.0523630, 0.0460620, 0.0519679]
    
    if keyword_set(shuffle) then begin
        i = randperm(n_elements(x))
        x = x[i]
        y = y[i]
        e_x = e_x[i]
        e_y = e_y[i]
    endif

    ;---------------------------------------------------------------------------
    ; FITTING
    ;---------------------------------------------------------------------------
    ; Fit relationship y = a * (x - x0) + b
    x0 = 2.5

    ;---------------------------------------------------------------------------
    ; Fit using MPFITEXY
    print, "### MPFITEXY, reduce for intrinsic scatter"
    a = mpfitexy(x, y, e_x, e_y, x0 = x0, errors = mpfiterrors, $
        minchi2 = minchi2, dof = dof, e_int = e_int, guess = [-9., -26.], $
        /reduce, /quiet)
    print, "Fit = ", a
    print, "Errors = ", mpfiterrors
    print, "Reduced chi = ", sqrt(minchi2/dof), " for scatter ", e_int
    print, ""

    ;---------------------------------------------------------------------------
    ; Fit using FITEXY
    ; Use intrinsic scatter found in MPFITEXY, /reduce
    print, "### FITEXY, forward , adopt intrinsic scatter used in mpfitexy"
    fitexy_e_int = e_int
    fitexy, x - x0, y, c, d, x_sig = e_x, y_sig = sqrt(e_y^2 + fitexy_e_int^2), $
        fitexyerrors, minchi2exy, qexy
    ; Reverse because fitexy returns parameters backwards
    print, "Fit = ", [d, c]
    print, "Errors = ", reverse(fitexyerrors)
    print, "Reduced chi = ", sqrt(minchi2exy/dof), " for scatter ", fitexy_e_int
    print, ""

    ;---------------------------------------------------------------------------
    ; Fit using MPFITEXY, /inv
    print, "### MPFITEXY, /inv, reduce for intrinsic scatter"
    a_inv = mpfitexy(x, y, e_x, e_y, x0 = x0, errors = mpfiterrors_inv, $
        minchi2 = minchi2_inv, dof = dof, e_int = e_int_inv, $
        guess = [-9., -26.], /reduce, /quiet, /inv)
    print, "Fit = ", a_inv
    print, "Errors = ", mpfiterrors_inv
    print, "Reduced chi = ", sqrt(minchi2_inv/dof), " for scatter ", e_int_inv
    print, ""

    ;---------------------------------------------------------------------------
    ; Fit using MPFITEXY
    print, "### MPFITEXY, with x and y reversed, reduce for intrinsic scatter"
    a = mpfitexy(y, x - x0, e_y, e_x, x0 = 0., errors = mpfiterrors, $
        minchi2 = minchi2, dof = dof, e_int = e_int, guess = [-1/9., 26./9.], $
        /reduce, /quiet)
    print, "Fit = ", [1/a[0], -a[1]/a[0]]
    print, "Errors (IGNORE) = ", mpfiterrors
    print, "Reduced chi = ", sqrt(minchi2/dof), " for scatter ", e_int
    print, ""

    ;---------------------------------------------------------------------------
    ; Fit using FITEXY
    ; Use no intrinsic scatter
    print, "### FITEXY, x and y swapped, no intrinsic scatter"
    fitexy_e_int_inv = 0
    fitexy, y, x - x0, c_inv, d_inv, x_sig = e_y, $
        y_sig = sqrt(e_x^2 + fitexy_e_int_inv^2), $
        fitexyerrors_inv, minchi2exy_inv, qexy_inv
    print, "Fit = ", [1/d_inv, -c_inv/d_inv]
    ; Convert errors in inverse to forward. FITEXY does not return full 
    ; covariance matrix, so these derived errors are severely biased.
    e_slope = fitexyerrors_inv[1]/d_inv^2
    e_inter2 = fitexyerrors_inv[1]^2 * c_inv^2/d_inv^4 + $
        fitexyerrors_inv[0]^2/d_inv^2
    print, "Errors (BIASED) = ", [e_slope, sqrt(e_inter2)]
    print, "Reduced chi = ", sqrt(minchi2exy_inv/dof), $
        " for scatter ", fitexy_e_int_inv

    ;---------------------------------------------------------------------------
    ; PLOTTING
    ;---------------------------------------------------------------------------
    ; Plot data
    ;---------------------------------------------------------------------------
    xrange = [2.0, 2.5]
    yrange = [-22, -26]
    plot, xrange, yrange, yrange = yrange, /nodata
    oploterror, x, y, e_x, e_y, psym = 1
    ;---------------------------------------------------------------------------
    ; Plot fit from mpfitexy (blue line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, a[0]*(!x.crange - x0) + a[1], $
        color = fsc_color("Blue"), thick = 4
    ;---------------------------------------------------------------------------
    ; Plot fit from mpfitexy (cyan line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, a_inv[0]*(!x.crange - x0) + a_inv[1], $
        color = fsc_color("Cyan"), thick = 4
    ;---------------------------------------------------------------------------
    ; Plot fit from fitexy (dashed red line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, c + d*(!x.crange - x0), color = fsc_color("Red"), $
        thick = 2
    ;---------------------------------------------------------------------------
    ; Plot inverse fit from fitexy (dashed red line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, -c_inv/d_inv + (1/d_inv)*(!x.crange - x0), $
        color = fsc_color("Brown"), thick = 2
end
