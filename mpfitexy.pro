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
    fixint = fixint, x0 = x0, reduce = reduce
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
    ;          e_int: intrinsic scatter in data. Should be adjusted to ensure
    ;                 sqrt(minchi2/dof) ~= 1.0
    ;        /fixint: fix the intercept to guess[0]
    ;        /reduce: adjust intrinsic scatter e_int to ensure chired ~= 1.0
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
    if n_elements(e_int) eq 0 then e_int = 0.0
    if n_elements(guess) eq 0 then guess = [0., 0.] else guess = float(guess)
    if n_elements(x0) eq 0 then x0 = 0
    if keyword_set(reduce) then e_int = mean(e_x)

    ;---------------------------------------------------------------------------
    ; RESCALE X-COORDS TO X0
    ;---------------------------------------------------------------------------
    x_ = x - x0

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
    ok = where(finite(x_) and finite(y))

    ;---------------------------------------------------------------------------
    ; CALL MPFIT ONCE
    ;---------------------------------------------------------------------------
    result = mpfit('lineresid', guess, functargs = {x:x_[ok], y:y[ok], $
        e_x:e_x[ok], e_y:e_y[ok], e_int:e_int}, parinfo = pi, $
        status = status, errmsg = errmsg, bestnorm = minchi2, dof = dof, $
        perror = perror, quiet = quiet)
    chired = sqrt(minchi2/dof)
    if ~keyword_set(quiet) then print, chired, e_int
    ;---------------------------------------------------------------------------
    ; CALL MPFIT UNTIL REDUCED CHI2 ~= 1.0 IF REQUIRED
    ;---------------------------------------------------------------------------
    if keyword_set(reduce) then begin
        while abs(chired - 1.) gt 0.01 do begin
            e_int = e_int * chired^(4./3)
            result = mpfit('lineresid', guess, functargs = {x:x_[ok], y:y[ok], $
                e_x:e_x[ok], e_y:e_y[ok], e_int:e_int}, parinfo = pi, $
                status = status, errmsg = errmsg, bestnorm = minchi2, dof = dof, $
                perror = perror, quiet = quiet)
            chired = sqrt(minchi2/dof)
            if ~keyword_set(quiet) then print, chired, e_int
        endwhile
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
    a = mpfitexy(x, y, e_x, e_y, x0 = x0, errors = mpfiterrors, $
        minchi2 = minchi2, dof = dof, e_int = e_int, guess = [-9., -26.], $
        /reduce, /quiet)
    print, "Fit = ", a
    print, "Errors = ", mpfiterrors
    print, "Reduced chi = ", sqrt(minchi2/dof), " for scatter ", e_int

    ;---------------------------------------------------------------------------
    ; Intrinsic scatter for chi-squared ~= 1.0
    e_int = 0.45
    fitexy, x - x0, y, c, d, x_sig = e_x, y_sig = sqrt(e_y^2 + e_int^2), $
        fitexyerrors, minchi2exy, qexy
    ; Reverse because fitexy returns parameters backwards
    print, [d, c], reverse(fitexyerrors), sqrt(minchi2exy/dof)

    ;---------------------------------------------------------------------------
    ; PLOTTING
    ;---------------------------------------------------------------------------
    ; Plot data
    ;---------------------------------------------------------------------------
    ploterror, x, y, e_x, e_y, psym = 1
    ;---------------------------------------------------------------------------
    ; Plot fit from mpfitexy (thick blue line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, a[0]*(!x.crange - x0) + a[1], $
        color = fsc_color("Blue"), thick = 4
    ;---------------------------------------------------------------------------
    ; Plot fit from fitexy (dashed red line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, c + d*(!x.crange - x0), color = fsc_color("Red"), $
        linestyle = 2, thick = 2
end
;-------------------------------------------------------------------------------
pro testmpfitexyorder, ps = ps
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Code to test whether shuffling arrays makes a difference
    ;---------------------------------------------------------------------------
    
    ;---------------------------------------------------------------------------
    ; SETUP
    ;---------------------------------------------------------------------------
    ; Set up two samples of data.
    s0 = ['NGC_0128', 'ESO_151G04', 'NGC_1381', 'NGC_1596', 'NGC_2310', $
        'ESO_311G12', 'NGC_3203', 'NGC_4469', 'NGC_4710', 'NGC_6771', $
        'ESO_597G36', 'NGC_1032', 'NGC_3957', 'NGC_5084']
    sp = ['NGC_3390', 'PGC_44931', 'ESO_443G42', 'NGC_5746', 'IC_4767', $
        'NGC_6722', 'ESO_185G53', 'IC_4937', 'IC_5096', 'ESO_240G11', $
        'NGC_1886', 'NGC_4703', 'NGC_7123', 'IC_5176']

;   ; Old galaxy order (FIX: this gives different answers for log mass fits,
;   ; and possibly others. This may be a bug in mpfitexy)
;   s0 = ['NGC_0128', 'NGC_1381', 'NGC_1596', 'NGC_2310', 'ESO_311G12', $
;           'NGC_3203', 'NGC_4469', 'NGC_4710', 'NGC_6771', 'ESO_597G36', $
;           'NGC_1032', 'NGC_3957', 'NGC_5084', 'ESO_151G04']
;   sp = ['NGC_3390', 'ESO_240G11', 'ESO_443G42', 'IC_4767', 'IC_4937', $
;           'IC_5096', 'IC_5176', 'NGC_1886', 'ESO_185G53', 'NGC_4703', $
;           'NGC_5746', 'NGC_6722', 'NGC_7123', 'PGC_44931']

    ; Fix distance and apparent magnitude errors
    fixderror = 0.0
    fixapperror = 0.0
    logmasssp = logmass(sp, error = e_logmasssp, fixapperror = fixapperror, $
        fixderror = fixderror) * 1.d
    logmasss0 = logmass(s0, error = e_logmasss0, fixapperror = fixapperror, $
        fixderror = fixderror) * 1.d
    ; Evaluate model circular velocity and error
    vcsp = modelmeanvc(sp, error = e_vcsp) * 1.d
    logvcsp = alog10(vcsp) * 1.d
    e_logvcsp = e_vcsp/(vcsp * alog(10))
    vcs0 = modelmeanvc(s0, error = e_vcs0) * 1.d
    logvcs0 = alog10(vcs0) * 1.d
    e_logvcs0 = e_vcs0/(vcs0 * alog(10))

    x0 = 2.4
    vmassguess = [3.5, 11.0] * 1.d
    print, mpfitexy([logvcsp, logvcs0], [logmasssp, logmasss0], $
        [e_logvcsp, e_logvcs0], [e_logmasssp, e_logmasss0], /quiet)

end
