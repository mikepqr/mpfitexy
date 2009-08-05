;-------------------------------------------------------------------------------
; See http://www.dfanning.com/code_tips/randperm.html
function randperm, n, seed = seed
      if n_elements(seed) eq 0 then seed = systime(1)
      x = lindgen(n)
      return, x[sort(randomu(seed, n))]
end
;-------------------------------------------------------------------------------
pro testmpfit_c, shuf = shuf, dbl = dbl, mpresult = a, plt = plt, $
    compare = compare, seed=seed, linfitresult=b, $
               base_sub=base_sub
    ;---------------------------------------------------------------------------
    ; /shuffle: shuffle data
    ;     /dbl: cast all data in to double precision
    ;     /plt: plot data and fits
    ; /compare: also fit with FITEXY from astroidl and LINFIT from IDL
    ; mpresult: return result of MPFIT [slope, intercept]
    ;---------------------------------------------------------------------------
    ; Sample data
    ;---------------------------------------------------------------------------
    x = [2.43184, 2.34696, 2.21818, 2.23748, 2.16919, 2.12780, 2.15438, $
        2.24805, 2.29982, 2.14897, 2.17739, 2.22395, 2.12592, 2.44776, $
        2.19388, 2.37343, 2.39492, 2.23249, 2.17179, 2.40253, 2.38844, $
        2.37377, 2.22880, 2.21212, 2.26511, 2.32707, 2.33750, 2.15525]
    y = [-25.3406, -24.8878, -22.6971, -22.8706, -22.7128, -22.4252, -22.9278, $
        -23.8840, -24.7333, -23.0504, -23.5204, -24.0678, -23.0894, -25.5364, $
        -23.5926, -25.3422, -25.0512, -24.2071, -24.1777, -25.1833, -24.8495, $
        -24.5180, -24.4620, -23.2211, -25.2567, -24.8611, -25.2406, -23.2333]
    e_y = [0.0434049, 0.0429708, 0.208040, 0.159271, 0.0503973, 0.103848, $
        0.0504721, 0.0643180, 0.0422348, 0.0273715, 0.0239415, 0.0999366, $
        0.0457358, 0.0461065, 0.0563985, 0.0450141, 0.0479903, 0.0541895, $
        0.0601011, 0.0521960, 0.0512739, 0.0498970, 0.0448522, 0.0457561, $
        0.0708393, 0.0523630, 0.0460620, 0.0519679]

    guess = [-9., -26.]
    if keyword_set(dbl) then begin
        x = x * 1.d
        y = y * 1.d
        e_y = e_y * 1.d
        guess = guess * 1.d
    endif
    
    if keyword_set(shuf) then begin
        i = randperm(n_elements(x), seed=seed)
        x = x[i]
        y = y[i]
        e_y = e_y[i]
    endif

    baseline = [0,0]
    if keyword_set(base_sub) then begin
       baseline = [avg(x), avg(y)]
       x = x - baseline[0]
       y = y - baseline[1]
    endif

    ; Fit using MPFIT
    ;---------------------------------------------------------------------------
    line = 'p[0] * x + p[1]'
    a = mpfitexpr(line, x, y, e_y, guess, bestnorm = minchi2, dof = dof, /quiet, $
                  status=status, $
                  xtol=(1d-10 + 3e-6*(1-keyword_set(dbl))), $
                  ftol=(1d-10 + 3e-6*(1-keyword_set(dbl))))
    print, "MPFIT: " + "y = " + strtrim(a[0], 2) + " * x + " + strtrim(a[1], 2), $
           '    status=',status

    b = linfit(x, y, measure_errors = e_y)
    if keyword_set(compare) then begin
        ; Fit using LINFIT for comparison
        ;-----------------------------------------------------------------------
        b = linfit(x, y, measure_errors = e_y)
        print, "LINFIT: " + "y = " + strtrim(b[1], 2) + " * x + " + strtrim(b[0], 2)

        ; Fit using FITEXY with no x-error for comparison
        ;-----------------------------------------------------------------------
        fitexy, x, y, c, d, x_sig = x * 0., y_sig = e_y, fitexyerrors, minchi2exy
        ; Reverse because fitexy returns parameters backwards
        print, "FITEXY: " + "y = " + strtrim(d, 2) + " * x + " + strtrim(c, 2)
    endif

    if keyword_set(plt) then begin
        ploterror, x, y, e_y, psym = 1
        oplot, !x.crange, a[0]*(!x.crange) + a[1], thick = 4, color = 50
        oplot, !x.crange, c + d*(!x.crange), linestyle = 2, thick = 2
    endif
end
;-------------------------------------------------------------------------------
pro batchtestmpfit_c, dbl = dbl, n = n, base_sub=base_sub
    ;---------------------------------------------------------------------------
    ; /dbl: cast data to double precision
    ;    n: repeat n times (default = 10000)
    ;---------------------------------------------------------------------------
    if n_elements(n) eq 0 then n = 1000
    mpresults = dblarr(2, n) & linfitresults = mpresults
    for i = 0, n - 1 do begin
        testmpfit_c, /shuf, dbl = dbl, mpresult = a, seed=seed, linfitresult=b, $
                   base_sub=base_sub
        mpresults[*, i] = a
        linfitresults[*, i] = b
    endfor
    !p.multi = [0, 1, 2]
    plothist, mpresults[0, *], bin = 0.0001, xtitle = 'Slope'
    plothist, mpresults[1, *], bin = 0.0001, xtitle = 'Intercept'
    !p.multi = 0

    print, '==== MPFIT'
    print, 'Slope = ', avg(mpresults[0,*]), '+/-', stddev(mpresults[0,*])
    print, 'Int   = ', avg(mpresults[1,*]), '+/-', stddev(mpresults[1,*])

    print, '==== LINFIT'
    print, 'Slope = ', avg(linfitresults[1,*]), '+/-', stddev(linfitresults[1,*])
    print, 'Int   = ', avg(linfitresults[0,*]), '+/-', stddev(linfitresults[0,*])
end
