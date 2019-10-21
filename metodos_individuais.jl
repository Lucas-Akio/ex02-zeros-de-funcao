# Regula falsi
function regula_falsi(enl :: EquacaoNL{T};
                   atol = √eps(T),
                   rtol = √eps(T),
                   max_iter = 10_000,
                   max_time = 30.0,
                   max_eval = 10_000,
                  ) where T <: AbstractFloat

  start_time = time()
  x = 0
  f(x) = fun_val(enl, x)
  a = enl.x₀
  fa = f(a)
  ϵ = atol + rtol * abs(fa)
  if abs(fa) ≤ ϵ
    return a, fa, :resolvido, 0, time() - start_time
  end
  b = a + 1
  fb = f(b)
  if abs(fb) ≤ ϵ
    return b, fb, :resolvido, 0, time() - start_time
  end

  while fa * fb ≥ 0
    δ = b - a
    if δ > 1e5
      return a, fa, :falha, 0, time() - start_time
    end
    if abs(fa) < abs(fb)
      a -= δ
      fa = f(a)
    else
      b += δ
      fb = f(b)
    end
  end
    # agora, com o ponto b:
    x0 = a
    x1 = b
    y1 = f(x1)
    y0 = f(x0)
    Δt = 0.0
    iter = 0
    excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
    status = :desconhecido
    while !(excedido)
        x2 = x1 - y1 * (x1-x0)/(y1-y0)
        # x-tolerance.
        if min(abs(x2-x0),abs(x2-x1)) < atol
            status = :resolvido
            return x2, f(x2), status, iter, Δt
        end
        y = f(x2)
        # y-tolerance.
        if abs(y) < rtol
            status = :resolvido
            return x, f(x), status, iter, Δt
        end
        if sign(y0*y)== 1
            x0 = x2
            y0 = y
        else
            x1 = x2
            y1 = y
        end
        iter += 1
        Δt = time() - start_time
        excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
    end
    if excedido
            if Δt > max_time
              status = :max_time
            else
              status = :max_iter
            end
    end
    return x, f(x), status, iter, Δt
end

# Inverse Quadratic Interpolation
function iqi(enl :: EquacaoNL{T};
             atol = √eps(T),
             rtol = √eps(T),
             max_iter = 10_000,
             max_time = 30.0,
             max_eval = 10_000,
             ) where T <: AbstractFloat
  reset!(enl)
  x    = 0
  f(x) = fun_val(enl, x)
  x₀   = enl.x₀
  fx₀  = f(x₀)
  ϵ    = atol + rtol * abs(fx₀)
  Δt         = 0.0
  start_time = time()
  if abs(fx₀) ≤ ϵ
    return x₀, fx₀, :resolvido, 0, time() - start_time
  end
  x₁ = x₀ - 1.0
  fx₁ = f(x₁)
  if abs(fx₁) ≤ ϵ
    return x₁, fx₁, :resolvido, 0, time() - start_time
  end
  x₂ = x₀ - 2.0
  fx₂ = f(x₂)
  if abs(fx₂) ≤ ϵ
    return x₂, fx₂, :resolvido, 0, time() - start_time
  end
  iter = 0
  excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
  status = :desconhecido
  while !(excedido)
    x = x₀*fx₁*fx₂/((fx₀-fx₁)*(fx₀-fx₂)) +
        x₁*fx₀*fx₂/((fx₁-fx₀)*(fx₁-fx₂)) +
        x₂*fx₀*fx₁/((fx₂-fx₀)*(fx₂-fx₁))
    if x == Inf || x == -Inf || x == NaN || x <= eps()
      return x, NaN, :falha, iter, Δt
    end
    if min(abs(x-x₀),abs(x-x₁),abs(x-x₂)) < √eps(T)
      status = :resolvido
      return x, f(x), status, iter, Δt
    end
    fx = f(x)
    if abs(fx) < √eps(T)
      status = :resolvido
      return x, f(x), status, iter, Δt
    end
    x₀    = x₁
    fx₀   = fx₁
    x₁    = x₂
    fx₁   = fx₂
    x₂    = x
    fx₂   = fx
    iter += 1
    Δt = time() - start_time
    excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
  end
  if excedido
    if Δt > max_time
      status = :max_time
    else
      status = :max_iter
    end
  end
  return x, f(x), status, iter, Δt
end

# Brent
function brent(enl :: EquacaoNL{T};
	       atol = √eps(T),
               rtol = √eps(T),
               max_iter::Integer=10_000,
               max_time = 30.0,
               max_eval = 10_000,
               ) where T <: AbstractFloat

  reset!(enl)
  x  = 0
  fx = 0
  f(x) = fun_val(enl, x)
  x₀   = enl.x₀
  fx₀  = f(x₀)
  ϵ = atol + rtol * abs(fx₀)
  start_time = time()
  Δt = 0.0
  if abs(fx₀) ≤ ϵ
    return x₀, fx₀, :resolvido, 0, time() - start_time
  end

  x₁   = x₀ + 1
  fx₁  = f(x₁)
  if abs(fx₁) ≤ ϵ
    return x₁, fx₁, :resolvido, 0, time() - start_time
  end

  # Encontrando valores de função com sinais opostos.
  while fx₀ * fx₁ ≥ 0
    δ = x₁ - x₀
    if δ > 1e5
      return x₀, fx₀, :falha, 0, time() - start_time
    end
    if abs(fx₀) < abs(fx₁)
      x₀ -= δ
      fx₀ = f(x₀)
      if abs(fx₀) ≤ ϵ
        return x₀, fx₀, :resolvido, enl.contador, time() - start_time
      end
    else
      x₁ += δ
      fx₁ = f(x₁)
      if abs(fx₁) ≤ ϵ
        return x₁, fx₁, :resolvido, enl.contador, time() - start_time
      end
    end
  end

  resolvido = abs(fx₀) ≤ ϵ
  iter = 0
  excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
  EPS = eps(Float64)
  if abs(fx₀) < abs(fx₁)
     x₀,  x₁ = x₁,  x₀
    fx₀, fx₁ = x₁, fx₀
  end
  x₂   = x₀
  fx₂  = fx₀
  x₃   = x₂
  bisection = true
  while !(resolvido || excedido)
    iter += 1
    if iter ≠ 1
      fx₀   = f(x₀)
      fx₁   = f(x₁)
    end
    if abs(fx₂) ≤ ϵ
      return x₂, fx₂, :resolvido, iter, time() - start_time
    end
    if fx₀ ≠ fx₂ && fx₁ ≠ fx₂
      x = x₀*fx₁*fx₂/((fx₀-fx₁)*(fx₀-fx₂)) +
          x₁*fx₀*fx₂/((fx₁-fx₀)*(fx₁-fx₂)) +
          x₂*fx₀*fx₁/((fx₂-fx₀)*(fx₂-fx₁))
      #println("rodou iqi, iter = $iter, x = $x, f(x) = $(f(x))")
    else
      x = x₁ - fx₁ * (x₁-x₀)/(fx₁-fx₀)
      #println("rodou secante, iter = $iter, x = $x, f(x) = $(f(x))")
    end
    delta = abs(2EPS*abs(x₁))
    min1 = abs(x-x₁)
    min2 = abs(x₁-x₂)
    min3 = abs(x₂-x₃)
    if (x < (3x₀+x₁)/4 && x > x₁) ||
        (bisection && min1 >= min2/2) ||
        (!bisection && min1 >= min3/2) ||
        (bisection && min2 < delta) ||
        (!bisection && min3 < delta)
        x = (x₀+x₁)/2
        bisection = true
    else
        bisection = false
    end
    fx = f(x)
    resolvido = abs(fx) ≤ ϵ
    x3 = x₂
    x₂ = x₁
    if sign(fx₀) != sign(fx)
        x₁ = x
        y1 = fx
    else
        x₁ = x
        fx₀ = fx
    end
    if abs(fx₀) < abs(fx₁)
    # Swap lower and upper bounds.
        x₀, x₁   = x₁, x₀
        fx₀, fx₁ = fx₁, fx₀
    end
    Δt = time() - start_time
    excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
  end

  status = :desconhecido
  if resolvido
    status = :resolvido
  elseif excedido
    if Δt > max_time
      status = :max_time
    else
      status = :max_iter
    end
  end

  return x, fx, status, iter, Δt
end

# Brent M1
function brent_M1(enl :: EquacaoNL{T};
	       atol = √eps(T),
               rtol = √eps(T),
               max_iter::Integer=10_000,
               max_time = 30.0,
               max_eval = 10_000,
               ) where T <: AbstractFloat

  reset!(enl)
  x  = 0
  fx = 0
  f(x) = fun_val(enl, x)
  x₀   = enl.x₀
  fx₀  = f(x₀)
  ϵ = atol + rtol * abs(fx₀)
  start_time = time()
  Δt = 0.0
  if abs(fx₀) ≤ ϵ
    return x₀, fx₀, :resolvido, 0, time() - start_time
  end

  x₁   = x₀ + 1
  fx₁  = f(x₁)
  if abs(fx₁) ≤ ϵ
    return x₁, fx₁, :resolvido, 0, time() - start_time
  end

  # Encontrando valores de função com sinais opostos.
  while fx₀ * fx₁ ≥ 0
    δ = x₁ - x₀
    if δ > 1e5
      return x₀, fx₀, :falha, 0, time() - start_time
    end
    if abs(fx₀) < abs(fx₁)
      x₀ -= δ
      fx₀ = f(x₀)
      if abs(fx₀) ≤ ϵ
        return x₀, fx₀, :resolvido, enl.contador, time() - start_time
      end
    else
      x₁ += δ
      fx₁ = f(x₁)
      if abs(fx₁) ≤ ϵ
        return x₁, fx₁, :resolvido, enl.contador, time() - start_time
      end
    end
  end

  #println("Começando Brent")
  # Brent melhorado.
  resolvido = abs(fx₀) ≤ ϵ
  iter = 0
  excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval

  if abs(fx₀) < abs(fx₁)
     x₀,  x₁ = x₁,  x₀
    fx₀, fx₁ = x₁, fx₀
  end
  fx₀   = f(x₀)
  fx₁   = f(x₁)
  while !(resolvido || excedido)
    iter += 1
#    if iter ≠ 1
#    end
    x₂    = (x₀ + x₁)/2
    fx₂   = f(x₂)
    if abs(fx₂) ≤ ϵ
      return x₂, fx₂, :resolvido, iter, time() - start_time
    end
    if fx₀ ≠ fx₂ && fx₁ ≠ fx₂
      x = x₀*fx₁*fx₂/((fx₀-fx₁)*(fx₀-fx₂)) +
          x₁*fx₀*fx₂/((fx₁-fx₀)*(fx₁-fx₂)) +
          x₂*fx₀*fx₁/((fx₂-fx₀)*(fx₂-fx₁))
      #println("rodou iqi, iter = $iter, x = $x, f(x) = $(f(x))")
    else
      x = x₁ - fx₁ * (x₁-x₀)/(fx₁-fx₀)
      #println("rodou secante, iter = $iter, x = $x, f(x) = $(f(x))")
    end
    fx = f(x)
    resolvido = abs(fx) ≤ ϵ
    if x₂ > x 
       x, x₂ =  x₂, x
      fx,fx₂ = fx₂,fx
    end
    if (fx₂*fx) < 0 
      x₀ = x
      x₁ = x₂
    else
      if (fx*fx₁) < 0
        x₀ = x₂
      else
        x₁ = x
      end
    end
    Δt = time() - start_time
    excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
  end

  status = :desconhecido
  if resolvido
    status = :resolvido
  elseif excedido
    if Δt > max_time
      status = :max_time
    else
      status = :max_iter
    end
  end

  return x, fx, status, iter, Δt
end

# Brent modificado
function brent_modificado(enl :: EquacaoNL{T};
	       atol = √eps(T),
               rtol = √eps(T),
               max_iter::Integer=10_000,
               max_time = 30.0,
               max_eval = 10_000,
               ) where T <: AbstractFloat

  reset!(enl)
  x  = 0
  fx = 0
  f(x) = fun_val(enl, x)
  x₀   = enl.x₀
  fx₀  = f(x₀)
  ϵ = atol + rtol * abs(fx₀)
  start_time = time()
  Δt = 0.0
  if abs(fx₀) ≤ ϵ
    return x₀, fx₀, :resolvido, 0, time() - start_time
  end

  x₁   = x₀ + 1
  fx₁  = f(x₁)
  if abs(fx₁) ≤ ϵ
    return x₁, fx₁, :resolvido, 0, time() - start_time
  end

  # Encontrando valores de função com sinais opostos.
  while fx₀ * fx₁ ≥ 0
    δ = x₁ - x₀
    if δ > 1e5
      return x₀, fx₀, :falha, 0, time() - start_time
    end
    if abs(fx₀) < abs(fx₁)
      x₀ -= δ
      fx₀ = f(x₀)
      if abs(fx₀) ≤ ϵ
        return x₀, fx₀, :resolvido, enl.contador, time() - start_time
      end
    else
      x₁ += δ
      fx₁ = f(x₁)
      if abs(fx₁) ≤ ϵ
        return x₁, fx₁, :resolvido, enl.contador, time() - start_time
      end
    end
  end

  # Brent modificado.
  resolvido = abs(fx₀) ≤ ϵ
  iter = 0
  excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval

  if abs(fx₀) < abs(fx₁)
     x₀,  x₁ = x₁,  x₀
    fx₀, fx₁ = x₁, fx₀
  end
  while !(resolvido || excedido)
    iter += 1
    if iter ≠ 1
      fx₀   = f(x₀)
      fx₁   = f(x₁)
    end
    x₂    = (x₀ + x₁)/2
    fx₂   = f(x₂)
    if abs(fx₂) ≤ ϵ
      return x₂, fx₂, :resolvido, iter, time() - start_time
    end
    if fx₀ ≠ fx₂ && fx₁ ≠ fx₂
      x = x₀*fx₁*fx₂/((fx₀-fx₁)*(fx₀-fx₂)) +
          x₁*fx₀*fx₂/((fx₁-fx₀)*(fx₁-fx₂)) +
          x₂*fx₀*fx₁/((fx₂-fx₀)*(fx₂-fx₁))
    else
      x = x₁ - fx₁ * (x₁-x₀)/(fx₁-fx₀)
    end
    fx = f(x)
    resolvido = abs(fx) ≤ ϵ
    if x₂ > x 
       x, x₂ =  x₂, x
      fx,fx₂ = fx₂,fx
    end
    if (fx₂*fx) < 0 
      x₀ = x
      x₁ = x₂
    else
      if (fx*fx₁) < 0
        x₀ = x₂
      else
        x₁ = x
      end
    end
    Δt = time() - start_time
    excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
  end

  status = :desconhecido
  if resolvido
    status = :resolvido
  elseif excedido
    if Δt > max_time
      status = :max_time
    else
      status = :max_iter
    end
  end

  return x, fx, status, iter, Δt
end

# Brent nerfado
function brent_nerfado(enl :: EquacaoNL{T};
	       atol = √eps(T),
               rtol = √eps(T),
               max_iter::Integer=10_000,
               max_time = 30.0,
               max_eval = 10_000,
               ) where T <: AbstractFloat

  reset!(enl)
  x  = 0
  fx = 0
  f(x) = fun_val(enl, x)
  x₀   = enl.x₀
  fx₀  = f(x₀)
  ϵ = atol + rtol * abs(fx₀)
  start_time = time()
  Δt = 0.0
  if abs(fx₀) ≤ ϵ
    return x₀, fx₀, :resolvido, 0, time() - start_time
  end

  x₁   = x₀ + 3
  fx₁  = f(x₁)
  if abs(fx₁) ≤ ϵ
    return x₁, fx₁, :resolvido, 0, time() - start_time
  end

  # Encontrando valores de função com sinais opostos.
  while fx₀ * fx₁ ≥ 0
    δ = x₁ - x₀
    if δ > 1e5
      return x₀, fx₀, :falha, 0, time() - start_time
    end
    if abs(fx₀) < abs(fx₁)
      x₀ -= δ
      fx₀ = f(x₀)
      if abs(fx₀) ≤ ϵ
        return x₀, fx₀, :resolvido, enl.contador, time() - start_time
      end
    else
      x₁ += δ
      fx₁ = f(x₁)
      if abs(fx₁) ≤ ϵ
        return x₁, fx₁, :resolvido, enl.contador, time() - start_time
      end
    end
  end

  # Brent modificado.
  resolvido = abs(fx₀) ≤ ϵ
  iter = 0
  excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval

  if abs(fx₀) < abs(fx₁)
     x₀,  x₁ = x₁,  x₀
    fx₀, fx₁ = x₁, fx₀
  end
  while !(resolvido || excedido)
    iter += 1
    if iter ≠ 1
      fx₀   = f(x₀)
      fx₁   = f(x₁)
    end
    x₂    = (x₀ + x₁)/2
    fx₂   = f(x₂)
    if abs(fx₂) ≤ ϵ
      return x₂, fx₂, :resolvido, iter, time() - start_time
    end
    if fx₀ ≠ fx₂ && fx₁ ≠ fx₂
      x = x₀*fx₁*fx₂/((fx₀-fx₁)*(fx₀-fx₂)) +
          x₁*fx₀*fx₂/((fx₁-fx₀)*(fx₁-fx₂)) +
          x₂*fx₀*fx₁/((fx₂-fx₀)*(fx₂-fx₁))
    else
      x = x₁ - fx₁ * (x₁-x₀)/(fx₁-fx₀)
    end
    fx = f(x)
    resolvido = abs(fx) ≤ ϵ
    if x₂ > x 
       x, x₂ =  x₂, x
      fx,fx₂ = fx₂,fx
    end
    if (fx₂*fx) < 0 
      x₀ = x
      x₁ = x₂
    else
      if (fx*fx₁) < 0
        x₀ = x₂
      else
        x₁ = x
      end
    end
    Δt = time() - start_time
    excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
  end

  status = :desconhecido
  if resolvido
    status = :resolvido
  elseif excedido
    if Δt > max_time
      status = :max_time
    else
      status = :max_iter
    end
  end

  return x, fx, status, iter, Δt
end
