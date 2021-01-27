;+
; NAME:
;   MPFITEXY2
;
; AUTHOR:
;   Michael Lee Williams <mike@mike.place>
;
; PURPOSE
;   Uses MPFIT to determine a common slope and two intercepts for two sets of
;   data. Uses the method of Williams (2010), i.e. minimizes that paper's
;   Equation 3 (or Equation 8 if called with /inv).
;                     y1 = a x1 + b 
;                     y2 = a x2 + b + d
;   MPFITEXY2 returns the array [a, b, d]
;
;   MPFITEXY2 is an implementation of the minimization described in Section 
;   4 of Williams M.J., Bureau M., Cappellari M., 2010, MNRAS, 409, 1330. You 
;   should read that short section before using the code. It adds support for
;   intrinsic scatter (see Tremaine et al. 2002) and provides a number of other
;   convenient features. See http://purl.org/mike/mpfitexy for an overview,
;   links to background reading, and important usage warnings.
;
;   The features of the code are discussed, and relevant literature cited, in
;   Section 4 of Williams, Bureau & Cappellari, 2010, MNRAS, 409, 1330. In the
;   interests of reproducibility, you should either cite that publication (which
;   I would prefer) or refer to the permanent URL at which MPFITEXY2 can always
;   be found: http://purl.org/mike/mpfitexy. MPFITEXY2 is dependent on the MPFIT
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
;  result = mpfitexy2(x1, y1, x2, y2, e_x1, e_y1, e_x2, e_y2, x0 = x0, 
;       guess = guess, e_int_guess = e_int_guess, /fixslope, /fixint, /reduce, 
;       /inv, /quiet, /silent, /latex, errors = perror, minchi2 = minchi2, $
;       dof = dof, e_int_reduce = e_int_reduce, scatter = scatter
;
; INPUTS:
; x1, y1, x2, y2: independent and dependent variables [required]. Note order.
;                 Non-finite values are silently ignored.
;    e_x1, etc. : corresponding error bars [required]
;                 Must be non-zero.
;
;             x0: fit the relationships 
;                   y = a(x - x0) + b, 
;                   y = a(x - x0) + b + d.
;                 Default: x0 = 0.
;          guess: starting point guesses for [slope, intercept, offset].
;                 fixed (see below). Default = [1., 1., 1.]
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
;         result: three-element array containing best parameters of model: 
;                 [slope, intercept, offset]
;         errors: 1-sigma fitting errors in slope, intercept and offset.
;                 Not meaningful if sqrt(minchi2/dof) != 1.0 (which implies
;                 observation uncertainties and/or intrinsic scatter not
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
;  2016.06.28 - Fix minor bug with NaNs in errors
;               https://github.com/williamsmj/mpfitexy/issues/1
;               (Release tagged v1.0.3)
;  2021.01.26 - Update name and email
;               (Release tagged v1.0.4)
;- 
; Copyright (C) 2009-2021, Michael Lee Williams <mike@mike.place>
; This software is provided as is without any warranty whatsoever. Permission 
; to use, copy, modify, and distribute modified or unmodified copies is 
; granted, provided this copyright notice and disclaimer are included unchanged.
; All other rights reserved.
;-

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
    ; Residual of data from models with these data and choice of parameters
    ;---------------------------------------------------------------------------
    slope = p[0]
    intercept = p[1]
    delta = p[2]
    f1 = slope * x1 + intercept
    f2 = slope * x2 + intercept + delta

    if n_elements(e_int) eq 0 then e_int = 0.d
    
    resid1 = (y1 - f1)/sqrt((e_y1^2 + slope^2 * e_x1^2 + e_int^2))
    resid2 = (y2 - f2)/sqrt((e_y2^2 + slope^2 * e_x2^2 + e_int^2))

    return, [resid1, resid2]
end
;-------------------------------------------------------------------------------
function mpfitexy2, x1, y1, x2, y2, e_x1, e_y1, e_x2, e_y2, x0 = x0, $
    guess = guess, fixslope = fixslope, fixint = fixint, $
    e_int_guess = e_int_guess, reduce = reduce, inv = inv, quiet = quiet, $
    silent = silent, latex = latex, errors = perror, minchi2 = minchi2, $
    dof = dof, e_int_reduce = e_int_reduce, scatter = scatter
    ;---------------------------------------------------------------------------
    ; DEFAULTS
    ;---------------------------------------------------------------------------
    if n_elements(e_int_guess) eq 0 then e_int = 0.d else e_int = e_int_guess
    if n_elements(guess) ne 3 then guess_ = [1.d, 1.d, 1.d] else $
        guess_ = double(guess)
    if n_elements(x0) eq 0 then x0 = 0.d
    if keyword_set(fixint) and keyword_set(fixslope) then $
        message, "MPFITEXY cannot be used with both fixint and fixslope"

    ;---------------------------------------------------------------------------
    ; RESCALE X-COORDS TO X0 AND CONVERT EVERYTHING TO DOUBLE PRECISION
    ;---------------------------------------------------------------------------
    x1_ = double(x1) - x0
    x2_ = double(x2) - x0
    e_x1_ = double(e_x1)
    e_x2_ = double(e_x2)
    y1_ = double(y1)
    y2_ = double(y2)
    e_y1_ = double(e_y1)
    e_y2_ = double(e_y2)

    ;---------------------------------------------------------------------------
    ; SWAP X AND Y AND INVERT GUESS IF FITTING INVERSE FUNCTION
    ;---------------------------------------------------------------------------
    if keyword_set(inv) then begin
        x1temp = x1_
        x2temp = x2_
        e_x1temp = e_x1_
        e_x2temp = e_x2_

        x1_ = y1_
        x2_ = y2_
        e_x1_ = e_y1_
        e_x2_ = e_y2_
        
        y1_ = x1temp
        y2_ = x2temp
        e_y1_ = e_x1temp
        e_y2_ = e_x2temp

        forwardslopeguess = guess_[0]
        forwardinterceptguess = guess_[1]
        forwarddeltaguess = guess_[2]
        guess_[0] = 1 / forwardslopeguess
        guess_[1] = - forwardinterceptguess/forwardslopeguess
        guess_[2] = - forwarddeltaguess/forwardslopeguess
    endif

    ;---------------------------------------------------------------------------
    ; FIX SLOPE/LABEL PARAMETERS
    ;---------------------------------------------------------------------------
    pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], parname:''},3)
    if keyword_set(fixslope) then pi(0).fixed = 1
    if keyword_set(fixint) then pi(1).fixed = 1
    pi(0).parname = '    Slope'
    pi(1).parname = 'Intercept'
    pi(2).parname = '    Delta'

    ;---------------------------------------------------------------------------
    ; CHECK ALL VALUES FINITE
    ;---------------------------------------------------------------------------
    ok1 = where(finite(x1_) and finite(y1_))
    ok2 = where(finite(x2_) and finite(y2_))

    ;---------------------------------------------------------------------------
    ; CALL MPFIT ONCE
    ;---------------------------------------------------------------------------
    result = mpfit('parallellineresid', guess_, functargs = {x1:x1_[ok1], $
        y1:y1_[ok1], x2:x2_[ok2], y2:y2_[ok2], e_x1:e_x1_[ok1], $
        e_y1:e_y1_[ok1], e_x2:e_x2_[ok2], e_y2:e_y2_[ok2], e_int:e_int}, $
        parinfo = pi, status = status, errmsg = errmsg, bestnorm = minchi2, $
        dof = dof, perror = perror, quiet = quiet, covar = covar)
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
        e_int = mean([e_y1_[ok1], e_y2_[ok2]])/10.
        while abs(chired - 1.) gt 0.01 do begin
            e_int = e_int * chired^(4./3)
            result = mpfit('parallellineresid', guess_, functargs = {x1:x1_[ok1], $
                y1:y1_[ok1], x2:x2_[ok2], y2:y2_[ok2], e_x1:e_x1_[ok1], $
                e_y1:e_y1_[ok1], e_x2:e_x2_[ok2], e_y2:e_y2_[ok2], e_int:e_int}, $
                parinfo = pi, status = status, errmsg = errmsg, bestnorm = minchi2, $
                dof = dof, perror = perror, quiet = quiet, covar = covar)
            chired = sqrt(minchi2/dof)
            if ~keyword_set(quiet) then print, chired, e_int
            if ~keyword_set(silent) and status ne 1 then print, status
        endwhile
    endif

    ;---------------------------------------------------------------------------
    ; EVALUATE TOTAL SCATTER (see Bedregal et al. 2006, eqn. 18 or Verheijen
    ; et al. 2001, section 7)
    ;---------------------------------------------------------------------------
    w1 = 1/(e_y1_[ok1]^2 + result[0]^2 * e_x1_[ok1]^2 + e_int^2)
    w2 = 1/(e_y2_[ok2]^2 + result[0]^2 * e_x2_[ok2]^2 + e_int^2)
    numerator = total(w1 * (y1_[ok1] - (result[0] * x1_[ok1] + result[1]))^2) + $
        total(w2 * (y2_[ok2] - (result[0] * x2_[ok2] + result[1] + result[2]))^2)
    denominator = total(w1) + total(w2)
    scatter = sqrt(numerator/denominator)

    ;---------------------------------------------------------------------------
    ; FLIP BEST-FITTING PARAMS AND THEIR ERRORS IF FITTING INVERSE FUNCTION
    ; See PhD notes VIII.3 for derivation of error propagation
    ;---------------------------------------------------------------------------
    if keyword_set(inv) then begin
        forwardslope = result[0]
        forwardintercept = result[1]
        forwarddelta = result[2]
        forwardslopeerror = perror[0]
        forwardintercepterror = perror[1]
        forwarddeltaerror = perror[2]
        result[0] = 1 / forwardslope
        result[1] = - forwardintercept/forwardslope
        result[2] = - forwarddelta/forwardslope
        perror[0] = forwardslopeerror / forwardslope^2
        perror[1] = sqrt(forwardintercepterror^2/forwardslope^2 + $
            forwardslopeerror^2*forwardintercept^2/forwardslope^4 - $
            2 * covar[0,1] * forwardintercept/forwardslope^3)
        perror[2] = sqrt(forwarddeltaerror^2/forwardslope^2 + $
            forwardslopeerror^2*forwarddelta^2/forwardslope^4 - $
            2 * covar[0,2] * forwarddelta/forwardslope^3)
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
            math, result[2], math, sep, math, perror[2], math, sep, $
            chired, sep, $
            e_int, sep, $
            scatter, newline, $
            format = '(A1,F6.2,A1,A1,A1,F5.2,A1,A2,A2,F8.2,A1,A1,A1,F5.2,' + $
                'A1,A2,A2,F5.2,A1,A1,A1,F5.2,A1,A2,F5.2,A2,F5.2,A2,F5.2,A4)'
    endif

    e_int_reduce = e_int
    return, result
end
