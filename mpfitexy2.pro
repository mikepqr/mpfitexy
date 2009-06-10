;-------------------------------------------------------------------------------
function parallellineresid, p, x1 = x1, y1 = y1, x2 = x2, y2 = y2, $
    e_x1 = e_x1, e_y1 = e_y1, e_x2 = e_x2, e_y2 = e_y2, e_int = e_int
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Utility function called by mpfitexy2. Given two sets of data, returns the
    ; residuals weighted by both error bars and optional intrinsic scatter 
    ; when they are fitted to two lines with a common slope. 
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
    ; Residual of data from models with these data and choice of paramters
    ;---------------------------------------------------------------------------
    slope = p[0]
    intercept = p[1]
    delta = p[2]
    f1 = slope * x1 + intercept
    f2 = slope * x2 + intercept + delta

    if n_elements(e_int) eq 0 then e_int = 0.0
    
    resid1 = (y1 - f1)/sqrt((e_y1^2 + slope^2 * e_x1^2 + e_int^2))
    resid2 = (y2 - f2)/sqrt((e_y2^2 + slope^2 * e_x2^2 + e_int^2))

    return, [resid1, resid2]
end
;-------------------------------------------------------------------------------
function mpfitexy2, x1, y1, x2, y2, e_x1, e_y1, e_x2, e_y2, guess = guess, $
    fixslope = fixslope, errors = perror, dof = dof, minchi2 = minchi2, $
    quiet = quiet, e_int = e_int, chired = chired, x0 = x0, reduce = reduce, $
    latex = latex
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Uses MPFIT to determine a common slope and two intercepts for two sets of
    ; data. Uses the method of Bureau (1996), i.e. minimizes their Equation (3).
    ;---------------------------------------------------------------------------
    ; INPUTS
    ; x1, y1, x2, y2: independent and dependent variables for two samples
    ;     e_x1, etc.: corresponding error bars
    ;          guess: starting points for common slope, intercept of sample 1, 
    ;                 difference between intercepts of the two samples. Slope
    ;                 can be fixed (see below). Default = [0.0, 0.0, 0.0]
    ;      /fixslope: fix the slope to guess[0] 
    ;         /quiet: Suppress MPFIT's text output
    ;         /latex: LaTeX output of fit params
    ;          e_int: intrinsic scatter in data. Should be adjusted to ensure
    ;                 sqrt(minchi2/dof) ~= 1.0
    ;---------------------------------------------------------------------------
    ; OUTPUTS
    ;         errors: 1 sigma fitting errors in paramters slope, intercept and 
    ;                 delta. Dubious if sqrt(minchi2/dof) != 1.0
    ;        minchi2: chi-squared of final model
    ;            dof: degrees of freedom
    ;         return: best parameters of model: slope, intercept, delta
    ;---------------------------------------------------------------------------

    ;---------------------------------------------------------------------------
    ; DEFAULTS
    ;---------------------------------------------------------------------------
    if n_elements(e_int) eq 0 then e_int = 0.0
    if n_elements(guess) ne 3 then guess = [1.0, 1.0, 1.0] else guess = float(guess)
    if n_elements(x0) eq 0 then x0 = 0
    if keyword_set(reduce) then e_int = 0.1

    ;---------------------------------------------------------------------------
    ; RESCALE X-COORDS TO X0
    ;---------------------------------------------------------------------------
    x1_ = x1 - x0
    x2_ = x2 - x0

    ;---------------------------------------------------------------------------
    ; FIX SLOPE/LABEL PARAMETERS
    ;---------------------------------------------------------------------------
    pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], parname:''},3)
    if keyword_set(fixslope) then pi(0).fixed = 1
    pi(0).parname = '    Slope'
    pi(1).parname = 'Intercept'
    pi(2).parname = '    Delta'

    ;---------------------------------------------------------------------------
    ; CHECK ALL VALUES FINITE
    ;---------------------------------------------------------------------------
    ok1 = where(finite(x1_) and finite(y1))
    ok2 = where(finite(x2_) and finite(y2))

    ;---------------------------------------------------------------------------
    ; CALL MPFIT ONCE
    ;---------------------------------------------------------------------------
    result = mpfit('parallellineresid', guess, functargs = {x1:x1_[ok1], $
        y1:y1[ok1], x2:x2_[ok2], y2:y2[ok2], e_x1:e_x1[ok1], e_y1:e_y1[ok1], $
        e_x2:e_x2[ok2], e_y2:e_y2[ok2], e_int:e_int}, parinfo = pi, $
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
            result = mpfit('parallellineresid', guess, functargs = {x1:x1_[ok1], $
                y1:y1[ok1], x2:x2_[ok2], y2:y2[ok2], e_x1:e_x1[ok1], e_y1:e_y1[ok1], $
                e_x2:e_x2[ok2], e_y2:e_y2[ok2], e_int:e_int}, parinfo = pi, $
                status = status, errmsg = errmsg, bestnorm = minchi2, dof = dof, $
                perror = perror, quiet = quiet)
            chired = sqrt(minchi2/dof)
            if ~keyword_set(quiet) then print, chired, e_int
        endwhile
    endif

    ;---------------------------------------------------------------------------
    ; EVALUATE TOTAL SCATTER (see Bedregal et al. 2006, eqn. 18 or Verheijen
    ; et al. 2001, section 7)
    ;---------------------------------------------------------------------------
    w1 = 1/(e_y1[ok1]^2 + result[0]^2 * e_x1[ok1]^2 + e_int^2)
    w2 = 1/(e_y2[ok2]^2 + result[0]^2 * e_x2[ok2]^2 + e_int^2)
    numerator = total(w1 * (y1[ok1] - (result[0] * x1_[ok1] + result[1]))^2) + $
        total(w2 * (y2[ok2] - (result[0] * x2_[ok2] + result[1] + result[2]))^2)
    denominator = total(w1) + total(w2)
    scatter = sqrt(numerator/denominator)

    if keyword_set(latex) then begin
        sep = "&"
        math = "$"
        pm = "\pm"
        newline = "\\ "
        print, math, result[0], pm, perror[0], math, sep, $
            math, result[1], pm, perror[1], math, sep, $
            math, result[2], pm, perror[2], math, sep, $
            chired, sep, $
            e_int, sep, $
            scatter, newline, $
            format = '(A1,F5.2,A3,F5.2,A1,A2,A2,F8.2,A3,F5.2,A1,A2,A2,F5.2,' + $
             'A3,F5.2,A1,A2,F5.2,A2,F5.2,A2,F5.2,A4)'
    endif
    return, result
end

;-------------------------------------------------------------------------------
pro testmpfitexy2, ps = ps
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Code to test mpfitexy2 and mpfitxy
    ;---------------------------------------------------------------------------
    
    ;---------------------------------------------------------------------------
    ; SETUP
    ;---------------------------------------------------------------------------
    ; Set up two samples of data. Populate sample 1 with smaller values
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

    sample1 = where(y ge -8.3 * x - 5.4)
    sample2 = where(y lt -8.3 * x - 5.4)

    x1 = x[sample1]
    y1 = y[sample1]
    x2 = x[sample2]
    y2 = y[sample2]

    x0 = 2.5

    e_x1 = e_x[sample1]
    e_y1 = e_y[sample1]
    e_x2 = e_x[sample2]
    e_y2 = e_y[sample2]

    ;---------------------------------------------------------------------------
    ; FITTING
    ;---------------------------------------------------------------------------
    ; Intrinsic scatter ~= 0.25 mag found to give sqrt(chi2/DOF) ~= 1.0
    e_int = 0.25
    ;---------------------------------------------------------------------------
    ; Fit the whole sample at once, totally free, using FITEXY
    fitexy, x, y, a, b, x_sig = e_x, y_sig = e_y, fitexyerrors, minchi2exy
    print, "# FIT ALL (FITEXY, no intrinsic error)"
    print, "Fit = ", b, a
    print, "Errors = ", reverse(fitexyerrors)
    print, "Reduced chi2 = ", sqrt(minchi2exy/(n_elements(x1) + n_elements(x2) - 2))
    ;---------------------------------------------------------------------------
    ; Fit the whole sample at once, totally free, using MPFITEXY
    fitall = mpfitexy(x, y, e_x, e_y, e_int = 0., errors = errorsall, $
        /quiet, dof = dofall, minchi2 = minchi2all)
    print, "# FIT ALL (MPFITEXY, no intrinsic error)"
    print, "Fit = ", fitall
    print, "Errors = ", errorsall
    print, "Reduced chi2 = ", sqrt(minchi2all/dofall)
    ;---------------------------------------------------------------------------
    ; Fit the two samples individually using mpfitexy2
    fit = mpfitexy2(x1, y1, x2, y2, e_x1, e_y1, e_x2, e_y2, errors = errors, $
        dof = dof, minchi2 = minchi2, x0 = x0, /reduce, e_int = e_int_reduced, $
        /quiet, /latex)
    print, "MPFITEXY2"
    print, "Fit = ", fit
    print, "Errors = ", errors
    print, "Reduced chi2 = ", sqrt(minchi2/dof), $
        ", e_int for this reduced chi2 = ", e_int_reduced
    ;---------------------------------------------------------------------------
    ; Fit the two samples seperately, constrain the slope to be that found with
    ; mpfitexy2
    fit1 = mpfitexy(x1, y1, e_x1, e_y1, guess = [fit[0], 0.], /fixslope, $
        errors = errors1, dof = dof1, minchi2 = minchi21, e_int = e_int, /quiet)
    fit2 = mpfitexy(x2, y2, e_x2, e_y2, guess = [fit[0], 0.], /fixslope, $
        errors = errors2, dof = dof2, minchi2 = minchi22, e_int = e_int, /quiet)
    print, "MPFITEXY of two samples with fixed slope " + strtrim(fit[0], 2)
    print, "Fits = ", fit1, fit2
    print, "Errors = ", errors1, errors2
    print, "Reduced chi2s = ", sqrt(minchi21/dof1), sqrt(minchi22/dof2)
    ;---------------------------------------------------------------------------

    ;---------------------------------------------------------------------------
    ; PLOTTING
    ;---------------------------------------------------------------------------
    if keyword_set(ps) then begin
        openpsdev, file = '~/Desktop/mpfitexy2.ps', /mnras, xsize = 13, $
            ysize = 10
        !x.margin = [5, 1]
        !y.margin = [3, 1]
    endif

    color1 = fsc_color("Red")
    color2 = fsc_color("Blue")
    xrange = [2, 2.6]
    yrange = [-22, -26]

    ;---------------------------------------------------------------------------
    ; Plot data
    plot, xrange, yrange, yrange = yrange, /nodata, /xstyle, /ystyle
    oploterror, x1, y1, e_x1, e_y1, color = color1, errcolor = color1, psym = 1
    oploterror, x2, y2, e_x2, e_y2, color = color2, errcolor = color2, psym = 1
    ;---------------------------------------------------------------------------
    ; Plot single fit for whole sample from mpfitexy (solid white line)
    oplot, xrange, fitall[0] * xrange + fitall[1]
    ;---------------------------------------------------------------------------
    ; Plot two lines from mpfitexy2 (colored lines)
    oplot, xrange, fit[0] * (xrange - x0) + fit[1], color = color1, thick = 2
    oplot, xrange, fit[0] * (xrange - x0) + fit[1] + fit[2], color = color2, thick = 2
    ;---------------------------------------------------------------------------
    ; Plot two lines from running mpfitexy on each sample seperately 
    ; (dashed lines)
    oplot, xrange, fit1[0] * xrange + fit1[1], linestyle = 2
    oplot, xrange, fit2[0] * xrange + fit2[1], linestyle = 2
    if keyword_set(ps) then closepsdev
end

;-------------------------------------------------------------------------------
pro testmpfitexy2order, ps = ps, shufs0 = shufs0, shufsp = shufsp, old = old, $
    rev = rev, srt = srt
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

    if keyword_set(rev) then begin
        s0 = reverse(s0)
        sp = reverse(sp)
    endif

    if keyword_set(srt) then begin
        i = sort(s0)
        s0 = s0[i]
        j = sort(sp)
        sp = sp[j]
    endif

    if keyword_set(old) then begin
        ; Old galaxy order (FIX: this gives different answers for log mass fits,
        ; and possibly others. This may be a bug in mpfitexy)
        s0 = ['NGC_0128', 'NGC_1381', 'NGC_1596', 'NGC_2310', 'ESO_311G12', $
                'NGC_3203', 'NGC_4469', 'NGC_4710', 'NGC_6771', 'ESO_597G36', $
                'NGC_1032', 'NGC_3957', 'NGC_5084', 'ESO_151G04']
        sp = ['NGC_3390', 'ESO_240G11', 'ESO_443G42', 'IC_4767', 'IC_4937', $
                'IC_5096', 'IC_5176', 'NGC_1886', 'ESO_185G53', 'NGC_4703', $
                'NGC_5746', 'NGC_6722', 'NGC_7123', 'PGC_44931']
    endif

    if keyword_set(shufs0) then begin
        i = randperm(n_elements(s0))
        s0 = s0[i]
    endif
    if keyword_set(shufsp) then begin
        j = randperm(n_elements(sp))
        sp = sp[j]
    endif


    ; Fix distance and apparent magnitude errors
    fixderror = 0.0
    fixapperror = 0.0
    logmasssp = logmass(sp, error = e_logmasssp, fixapperror = fixapperror, $
        fixderror = fixderror)
    logmasss0 = logmass(s0, error = e_logmasss0, fixapperror = fixapperror, $
        fixderror = fixderror)
    ; Evaluate model circular velocity and error
    vcsp = modelmeanvc(sp, error = e_vcsp)
    logvcsp = alog10(vcsp)
    e_logvcsp = e_vcsp/(vcsp * alog(10))
    vcs0 = modelmeanvc(s0, error = e_vcs0)
    logvcs0 = alog10(vcs0)
    e_logvcs0 = e_vcs0/(vcs0 * alog(10))

    x0 = 2.4
    vmassguess = [3.0, 11.0, 0.2] * 1.d
    logvc_logmass = mpfitexy2(logvcsp * 1.d, logmasssp * 1.d, logvcs0 * 1.d, $
        logmasss0 * 1.d, e_logvcsp * 1.d, e_logmasssp * 1.d, e_logvcs0 * 1.d, $
        e_logmasss0 * 1.d, /quiet, x0 = x0, guess = vmassguess, /latex)
end
