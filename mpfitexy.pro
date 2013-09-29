;+
; NAME:
;   MPFITEXY
;
; AUTHOR:
;   Michael Williams <williams@astro.columbia.edu>, Columbia University
;
; PURPOSE
;   Uses MPFIT to determine straight line fit to data with errors in both
;   variables using the weighting defined in Numerical Recipes, i.e. 
;   emulates the ASTROIDL library function FITEXY with some bells and 
;   whistles.
;
;   MPFITEXY is merely a wrapper to MPFIT. In the absence of intrinsic scatter
;   it minimizes the same figure of merit at the ASTROIDL FITEXY routine. It
;   adds support for intrinsic scatter (see Tremaine et al. 2002) and provides a
;   number of other convenient features. See http://purl.org/mike/mpfitexy for
;   an overview, links to background reading, and important usage warnings.
;   
;   The features of the code are discussed, and relevant literature cited, in
;   Section 4 of Williams, Bureau & Cappellari, 2010, MNRAS, 409, 1330. In the
;   interests of reproducibility, you should either cite that publication (which
;   I would prefer) or refer to the permanent URL at which MPFITEXY can always
;   be found: http://purl.org/mike/mpfitexy. MPFITEXY is dependent on the MPFIT
;   package, which you must separately acknowledge by citing Markwardt C. B.,
;   2009, in Astronomical Data Analysis Software and Systems XVIII, Bohlender,
;   D., Dowler P., Durand D., eds., Astronomical Society of the Pacific
;   Conference Series. [http://adsabs.harvard.edu/abs/2009ASPC..411..251M]
;
;   The /reduce keyword, which adjusts the intrinsic scatter to ensure the
;   reduced chi^2 is ~= 1.0 uses the simple iterative procedure described in
;   Section 3.2.1 of Bedregal et al. 2006, MNRAS, 373, 1125.
; 
; CALLING SEQUENCE:
;  result = mpfitexy(x, y, e_x, e_y, x0 = x0, guess = guess, 
;                    e_int_guess = e_int_guess, /fixslope, /fixint, /quiet, 
;                    /silent, /reduce, /inv, dof = dof, errors = errors, 
;                    minchi2 = minchi2)
;
; INPUTS:
;           x, y: independent and dependent variables [required]
;                 Non-finite values are silently ignored.
;       e_x, e_y: corresponding error bars [required]
;                 Must be non-zero.
;
;             x0: fit the relationship y = a(x - x0) + b. Default: x0 = 0.
;          guess: starting point guesses for slope, intercept.
;                 fixed (see below). Default = [1., 1.]
;    e_int_guess: intrinsic scatter in data in units of y1, y2. For fit and 
;                 errors to be meaningful, this quantity should be adjusted 
;                 to ensure sqrt(minchi2/dof) ~= 1.0. This can be done either by
;                 manually adjusting e_int_guess, or by using the /reduce 
;                 keyword, which will adjust the intrinsic scatter.
;
;      /fixslope: fix the slope to guess[0]. Cannot be used with fixint.
;        /fixint: fix the intercept to guess[1]. Cannot be used with fixslope.
;        /reduce: adjust intrinsic scatter e_int to ensure chired ~= 1.0
;           /inv: fits the inverse relation x = a' y + b' (still returns
;                 slope of forward relation, i.e. slope = 1/a' and
;                 intercept = -b'/a'. Output variables (slope, intercept, 
;                 errors, etc.) are returned for the forward relation.
;                 See Section 4 of Williams et al. (2010) for details.
;         /quiet: Suppress MPFIT's text output
;        /silent: Do not print MPFIT status code (see mpfit.pro for docs)
;         /latex: Print result and errors in LaTeX table format
;
; OUTPUTS:
;         result: two-element array containing best parameters of model: 
;                 [slope, intercept]
;         errors: 1-sigma fitting errors in slope and intercept. Not 
;                 meaningful if sqrt(minchi2/dof) != 1.0 (which implies
;                 observational uncertainties and/or intrinsic scatter not
;                 well-estimated) or if used with /fixint or /fixslope.
;        minchi2: unreduced chi-squared of final model
;            dof: degrees of freedom
;   e_int_reduce: Intrinsic scatter used to ensure sqrt(minchi2/dof) ~= 1.
;        scatter: "Total" scatter, i.e. RMS distance of points from model.
;                 The e_int_reduce/scatter gives a useful idea of what 
;                 fraction of the scatter is observational and what is 
;                 intrinsic.
;
; KNOWN ISSUES:
; Vulnerable to rounding errors if input data are not prenormalized.
; 
;-
; MODIFICATION HISTORY:
; Pre-2009.08 - Initial private releases
;  2009.08.05 - Correctly propagate covariance term into intercept error after 
;               inverse fit (thanks to Tim Davis, University of Oxford)
;  2010.05.15 - Initial public release (v1.0, hg revision 24)
;  2011.06.20 - Update contact details, references, acknowledgment instructions
;               (Release tagged v1.0.1)
;  2013.09.29 - Update contact details
;               (Release tagged v1.0.2)
;- 
; Copyright (C) 2009-2013, Michael Williams <williams@astro.columbia.edu>
; This software is provided as is without any warranty whatsoever. Permission 
; to use, copy, modify, and distribute modified or unmodified copies is 
; granted, provided this copyright notice and disclaimer are included unchanged.
; All other rights reserved.
;-

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
    ; Residual of data from models with these data and choice of parameters
    ;---------------------------------------------------------------------------
    slope = p[0]
    intercept = p[1]
    f = slope * x + intercept
    if n_elements(e_int) eq 0 then e_int = 0.0
    resid = (y - f)/sqrt((e_y^2 + slope^2*e_x^2 + e_int^2))
    return, resid
end
;-------------------------------------------------------------------------------
function mpfitexy, x, y, e_x, e_y, x0 = x0, guess = guess, $
    fixslope = fixslope, fixint = fixint, e_int_guess = e_int_guess, $
    reduce = reduce, inv = inv, quiet = quiet, silent = silent, latex = latex, $
    errors = perror, minchi2 = minchi2, dof = dof, $
    e_int_reduce = e_int_reduce, scatter = scatter
    ;---------------------------------------------------------------------------
    ; DEFAULTS
    ;---------------------------------------------------------------------------
    if n_elements(e_int_guess) eq 0 then e_int = 0.d else e_int = e_int_guess
    if n_elements(guess) eq 0 then guess_ = [1.d, 1.d] else $
        guess_ = double(guess)
    if n_elements(x0) eq 0 then x0 = 0.d
    if keyword_set(fixint) and keyword_set(fixslope) then $
        message, "MPFITEXY cannot be used with both fixint and fixslope"

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
    if ~keyword_set(silent) and status ne 1 then print, "MPFIT error: ", status

    if chired lt 1.0 and keyword_set(reduce) then begin
        print, "chi-squared already less than 1.0. Not attempting to adjust e_int."
        print, chired
        reduce1 = 0
    endif else if abs(chired -1.) ge 0.01 and keyword_set(reduce) then begin
        reduce1 = 1 
    endif else reduce1 = 0

    ;---------------------------------------------------------------------------
    ; CALL MPFIT UNTIL REDUCED CHI2 ~= 1.0 IF REQUIRED
    ;---------------------------------------------------------------------------
    if keyword_set(reduce1) then begin
        e_int = mean(e_y_)/10.
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
    ; EVALUATE TOTAL SCATTER (see Williams et al. 2010, equations 6, 7 and 10.
    ; See also Bedregal et al. 2006, eqn. 18, Verheijen et al. 2001, section 7)
    ;---------------------------------------------------------------------------
    w = 1/(e_y_[ok]^2 + result[0]^2 * e_x_[ok]^2 + e_int^2)
    numerator = total(w * (y_[ok] - (result[0] * x_[ok] + result[1]))^2)
    denominator = total(w)
    scatter = sqrt(numerator/denominator)

    ;---------------------------------------------------------------------------
    ; FLIP BEST-FITTING PARAMS AND PROPAGATE ERRORS IF FITTING INVERSE FUNCTION
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
        scatter = abs(scatter * result[0])
    endif

    if keyword_set(latex) then begin
        sep = "&"
        math = "$"
        pm = "\pm"
        newline = "\\ "
        print, math, result[0], math, sep, math, perror[0], math, sep, $
            math, result[1], math, sep, math, perror[1], math, sep, $
            chired, sep, $
            e_int, sep, $
            scatter, newline, $
            format = '(A1,F6.2,A1,A1,A1,F5.2,A1,A2,A2,F8.2,A1,A1,A1,F5.2,' + $
                'A1,A2,F5.2,A2,F5.2,A2,F5.2,A4)'
    endif

    e_int_reduce = e_int
    return, result
end
