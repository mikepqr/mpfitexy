;-------------------------------------------------------------------------------
pro testmpfitexy, shuffle = shuffle
    ;---------------------------------------------------------------------------
    ; PURPOSE
    ; Code to test and mpfitexy
    ; Result should be a plot of four lines. Blue is mpfitexy, red if fitexy 
    ; (with mpfitexy's intrinsic scatter). These should agree. 
    ; 
    ; Cyan is mpfitexy, brown is fitexy with the order of
    ; the arguments reversed (and mpfitexy's intrinsic scatter). These should
    ; agree.
    ; 
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
        minchi2 = minchi2, dof = dof, e_int_guess = e_int_guess, guess = [-9., -26.], $
        /reduce, /quiet, e_int_reduce = e_int_reduce)
    print, "Fit = ", a
    print, "Errors = ", mpfiterrors
    print, "Reduced chi = ", sqrt(minchi2/dof), " for scatter ", e_int_reduce
    print, ""

    ;---------------------------------------------------------------------------
    ; Fit using FITEXY
    ; Use intrinsic scatter found in MPFITEXY, /reduce
    print, "### FITEXY, forward , adopt intrinsic scatter used in mpfitexy"
    fitexy_e_int = e_int_reduce
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
        minchi2 = minchi2_inv, dof = dof, e_int_reduce = e_int_inv, $
        guess = [-9., -26.], /reduce, /quiet, /inv)
    print, "Fit = ", a_inv
    print, "Errors = ", mpfiterrors_inv
    print, "Reduced chi = ", sqrt(minchi2_inv/dof), " for scatter ", e_int_inv
    print, ""

    ;---------------------------------------------------------------------------
    ; Fit using MPFITEXY
    print, "### MPFITEXY, with x and y reversed, reduce for intrinsic scatter"
    a_swap = mpfitexy(y, x - x0, e_y, e_x, x0 = 0., errors = mpfiterrors, $
        minchi2 = minchi2, dof = dof, e_int_reduce = e_int_reduce, guess = [-1/9., 26./9.], $
        /reduce, /quiet)
    print, "Fit = ", [1/a_swap[0], -a_swap[1]/a_swap[0]]
    print, "Errors (IGNORE) = ", mpfiterrors
    print, "Reduced chi = ", sqrt(minchi2/dof), " for scatter ", e_int_reduce
    print, ""

    ;---------------------------------------------------------------------------
    ; Fit using FITEXY
    ; Use no intrinsic scatter
    print, "### FITEXY, x and y swapped, with intrinsic scatter"
    fitexy_e_int_inv = 0
    fitexy, y, x - x0, c_inv, d_inv, x_sig = e_y, $
        y_sig = sqrt(e_x^2 + e_int_inv^2), $
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
        color = fsc_color("Blue"), thick = 4, linestyle = 1
    ;---------------------------------------------------------------------------
    ; Plot fit from mpfitexy /inv (cyan line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, a_inv[0]*(!x.crange - x0) + a_inv[1], $
        color = fsc_color("Cyan"), thick = 4, linestyle = 1
    ;---------------------------------------------------------------------------
    ; Plot fit from fitexy (dashed red line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, c + d*(!x.crange - x0), color = fsc_color("Red"), $
        thick = 2, linestyle = 1
    ;---------------------------------------------------------------------------
    ; Plot inverse fit from fitexy (swapped) (dashed brown line)
    ;---------------------------------------------------------------------------
    oplot, !x.crange, -c_inv/d_inv + (1/d_inv)*(!x.crange - x0), $
        color = fsc_color("Brown"), thick = 2, linestyle = 1
end

;-------------------------------------------------------------------------------
pro testfixslopeint, shuffle = shuffle
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

    x0 = 2.5

    guess = [-8.5, -26.3]
    print, guess

    fit = mpfitexy(x, y, e_x, e_y, x0 = x0, /reduce, e_int_reduce = e_int_reduce, $
        /quiet, /latex, guess = guess)
    fit_int = mpfitexy(x, y, e_x, e_y, x0 = x0, /reduce, e_int_reduce = e_int_reduce, $
        /quiet, /latex, /fixint, guess = guess)
    fit_slope = mpfitexy(x, y, e_x, e_y, x0 = x0, /reduce, e_int_reduce = e_int_reduce, $
        /quiet, /latex, /fixslope, guess = guess)

    xrange = [2.0, 2.5]
    yrange = [-22, -26]
    plot, xrange, yrange, yrange = yrange, /nodata
    oploterror, x, y, e_x, e_y, psym = 1

    oplot, !x.crange, fit[0]*(!x.crange - x0) + fit[1], $
        color = fsc_color("Blue"), thick = 4, linestyle = 1
    oplot, !x.crange, fit_int[0]*(!x.crange - x0) + fit_int[1], $
        color = fsc_color("Red"), thick = 4, linestyle = 1
    oplot, !x.crange, fit_slope[0]*(!x.crange - x0) + fit_slope[1], $
        color = fsc_color("Green"), thick = 4, linestyle = 1
end
;-------------------------------------------------------------------------------
pro testmpfiterrors
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

    x0 = 2.5

    e_ = e_x * 0.

    guess = [-8.5, -26.3]
    e_int_guess = 0.1

    fit = mpfitexy(x, y, e_, e_, x0 = x0, /reduce, e_int_reduce = e_int_reduce, $
        /latex, guess = guess, e_int_guess = e_int_guess)

    xrange = [2.0, 2.5]
    yrange = [-22, -26]
    plot, xrange, yrange, yrange = yrange, /nodata
    oploterror, x, y, e_x, e_y, psym = 1

    oplot, !x.crange, fit[0]*(!x.crange - x0) + fit[1], $
        color = fsc_color("Blue"), thick = 4, linestyle = 1
end
