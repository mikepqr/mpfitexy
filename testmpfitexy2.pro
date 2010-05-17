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
    fitall = mpfitexy(x, y, e_x, e_y, e_int_guess = 0., errors = errorsall, $
        /quiet, dof = dofall, minchi2 = minchi2all)
    print, "# FIT ALL (MPFITEXY, no intrinsic error)"
    print, "Fit = ", fitall
    print, "Errors = ", errorsall
    print, "Reduced chi2 = ", sqrt(minchi2all/dofall)
    ;---------------------------------------------------------------------------
    ; Fit the two samples individually using mpfitexy2
    fit = mpfitexy2(x1, y1, x2, y2, e_x1, e_y1, e_x2, e_y2, errors = errors, $
        dof = dof, minchi2 = minchi2, x0 = x0, /reduce, e_int_reduce = e_int_reduced, $
        /quiet)
    print, "MPFITEXY2"
    print, "Fit = ", fit
    print, "Errors = ", errors
    print, "Reduced chi2 = ", sqrt(minchi2/dof), $
        ", e_int = ", e_int_reduced
    ;---------------------------------------------------------------------------
    ; Fit the two samples seperately, constrain the slope to be that found with
    ; mpfitexy2
    fit1 = mpfitexy(x1, y1, e_x1, e_y1, guess = [fit[0], 0.], /fixslope, $
        errors = errors1, dof = dof1, minchi2 = minchi21, e_int_guess = e_int_guess, /quiet)
    fit2 = mpfitexy(x2, y2, e_x2, e_y2, guess = [fit[0], 0.], /fixslope, $
        errors = errors2, dof = dof2, minchi2 = minchi22, e_int_guess = e_int_guess, /quiet)
    print, "MPFITEXY of two samples with fixed slope " + strtrim(fit[0], 2)
    print, "Fits = ", fit1, fit2
    print, "Errors = ", errors1, errors2
    print, "Reduced chi2s = ", sqrt(minchi21/dof1), sqrt(minchi22/dof2)
    ;---------------------------------------------------------------------------

    ;---------------------------------------------------------------------------
    ; PLOTTING
    ;---------------------------------------------------------------------------
    if keyword_set(ps) then begin
        openpsdev, file = '~/Desktop/mpfitexy2.ps', /slide
        !x.margin = [8, 1]
        !y.margin = [4, 1]
    endif

    color1 = fsc_color("Red")
    color2 = fsc_color("Blue")
    xrange = [2, 2.6]
    yrange = [-22, -26]

    ;---------------------------------------------------------------------------
    ; Plot data
    plot, xrange, yrange, yrange = yrange, /nodata, /xstyle, /ystyle, $
        xtitle = textoidl('log (v / km s^{-1})'), ytitle = textoidl('M_K')
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
    ; Plot two lines from running mpfitexy on each sample separately 
    ; (dashed lines)
    oplot, xrange, fit1[0] * xrange + fit1[1], linestyle = 2
    oplot, xrange, fit2[0] * xrange + fit2[1], linestyle = 2
    if keyword_set(ps) then closepsdev
end

;-------------------------------------------------------------------------------
pro testmpfitexy2inv, ps = ps
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Code to test inverse fitting of mpfitexy2 and mpfitexy
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
    ; Fit the whole sample at once using MPFITEXY
    print, "# MPFITEXY"
    fitall = mpfitexy(x, y, e_x, e_y, e_int = e_int1, errors = errorsall, $
        /quiet, dof = dofall, minchi2 = minchi2all, guess = [-8.d, -26.d], $
        x0 = x0, /reduce)
    print, "Fit = ", fitall
    print, "Errors = ", errorsall
    print, "Reduced chi2 = ", sqrt(minchi2all/dofall), ", e_int = ", e_int1
    print, ""
    ;---------------------------------------------------------------------------
    ; Fit the whole sample at once, using MPFITEXY with inverse function
    print, "# MPFITEXY (inverse function)"
    fitallinv = mpfitexy(x, y, e_x, e_y, e_int = e_int2, errors = errorsallinv, $
        /quiet, dof = dofallinv, minchi2 = minchi2allinv, $
        guess = [-8.d, -26.d], /inv, x0 = x0, /reduce)
    print, "Fit = ", fitallinv
    print, "Errors = ", errorsallinv
    print, "Reduced chi2 = ", sqrt(minchi2allinv/dofallinv), ", e_int = ", e_int2
    print, ""
    ;---------------------------------------------------------------------------
    ; Fit the two samples individually using MPFITEXY2
    print, "MPFITEXY2"
    fit = mpfitexy2(x1, y1, x2, y2, e_x1, e_y1, e_x2, e_y2, errors = errors, $
        dof = dof, minchi2 = minchi2, x0 = x0, /reduce, e_int = e_int_reduced, $
        /quiet, /latex)
    print, "Fit = ", fit
    print, "Errors = ", errors
    print, "Reduced chi2 = ", sqrt(minchi2/dof), $
        ", e_int = ", e_int_reduced
    print, ""
    ;---------------------------------------------------------------------------
    ; Fit the two samples individually using MPFITEXY2 with inverse function
    print, "MPFITEXY2 (inverse)"
    fitinv = mpfitexy2(x1, y1, x2, y2, e_x1, e_y1, e_x2, e_y2, errors = errorsinv, $
        dof = dofinv, minchi2 = minchi2inv, x0 = x0, /reduce, $
        e_int = e_int_reducedinv, /inv, /quiet)
    print, "Fit = ", fitinv
    print, "Errors = ", errorsinv
    print, "Reduced chi2 = ", sqrt(minchi2inv/dofinv), $
        ", e_int = ", e_int_reducedinv
    ;---------------------------------------------------------------------------

    ;---------------------------------------------------------------------------
    ; PLOTTING
    ;---------------------------------------------------------------------------
    if keyword_set(ps) then begin
        openpsdev, file = '~/Desktop/mpfitexy2.ps'
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
    oplot, xrange, fitall[0] * (xrange - x0) + fitall[1]
    oplot, xrange, fitallinv[0] * (xrange -x0) + fitallinv[1], linestyle = 1
    ;---------------------------------------------------------------------------
    ; Plot two lines from mpfitexy2 (colored lines)
    oplot, xrange, fit[0] * (xrange - x0) + fit[1], color = color1
    oplot, xrange, fit[0] * (xrange - x0) + fit[1] + fit[2], color = color2
    oplot, xrange, fitinv[0] * (xrange - x0) + fitinv[1], color = color1, linestyle = 1
    oplot, xrange, fitinv[0] * (xrange - x0) + fitinv[1] + fitinv[2], color = color2, linestyle = 1
    if keyword_set(ps) then closepsdev
end
