
# Computes the logarithmic mean: (aR-aL)/(LOG(aR)-LOG(aL)) = (aR-aL)/LOG(aR/aL)
# Problem: if aL~= aR, then 0/0, but should tend to --> 0.5*(aR+aL)
#
# introduce xi=aR/aL and f=(aR-aL)/(aR+aL) = (xi-1)/(xi+1)
# => xi=(1+f)/(1-f)
# => Log(xi) = log(1+f)-log(1-f), and for small f (f^2<1.0E-02) :
#
#    Log(xi) ~=     (f - 1/2 f^2 + 1/3 f^3 - 1/4 f^4 + 1/5 f^5 - 1/6 f^6 + 1/7 f^7)
#                  +(f + 1/2 f^2 + 1/3 f^3 + 1/4 f^4 + 1/5 f^5 + 1/6 f^6 + 1/7 f^7)
#             = 2*f*(1           + 1/3 f^2           + 1/5 f^4           + 1/7 f^6)
#  (aR-aL)/Log(xi) = (aR+aL)*f/(2*f*(1 + 1/3 f^2 + 1/5 f^4 + 1/7 f^6)) = (aR+aL)/(2 + 2/3 f^2 + 2/5 f^4 + 2/7 f^6)
#  (aR-aL)/Log(xi) = 0.5*(aR+aL)*(105/ (105+35 f^2+ 21 f^4 + 15 f^6)
# @muladd @inline function ln_mean(_value1, _value2)
#   value1, value2 = ifelse(_value1 > _value2, (_value1, _value2), (_value2, _value1))
@muladd @inline function ln_mean(value1, value2)
# function ln_mean(value1, value2)
  epsilon_f2 = 1.0e-4
  ratio = value2 / value1
  # f2 = f^2
  f2 = (ratio * (ratio - 2) + 1) / (ratio * (ratio + 2) + 1)
  if f2 < epsilon_f2
    return (value1 + value2) * 52.5 / (105 + f2 * (35 + f2 * (21 + f2 * 15)))
  else
    return (value2 - value1) / inline_log(ratio)
    # return (value2 - value1) / log(ratio)
  end
end

# @inline function ln_mean(a, b)
#   r = a / b
#   f = (r - 1) / (r + 1)
#   u = f * f
#   if r < 3.355e-1
#       (a + b) * evalpoly(u, (1.0, -21/13, 105/143, -35/429)) /
#                 evalpoly(u, (2.0, -100/39, 566/715, -512/15015))
#   else
#       (a + b) * f / log(r)
#   end
# end

# @inline function ln_mean(x, y)
#   mi, ma = ifelse(x > y, (y, x), (x, y))
#   #mi, ma = x, y
#   r = mi / ma
#   @fastmath ma * ifelse(r == one(r), one(r), (r - 1) / log(r))
# end

@inline function inline_log(x::Float64)
  if x > 0.0
      x == Inf && return x

      # Step 2
      if 0.9394130628134757 < x < 1.0644944589178595
          f = x-1.0
          return Base.Math.log_proc2(f)
      end

      # Step 3
      xu = reinterpret(UInt64,x)
      m = Int(xu >> 52) & 0x07ff
      if m == 0 # x is subnormal
          x *= 1.8014398509481984e16 # 0x1p54, normalise significand
          xu = reinterpret(UInt64,x)
          m = Int(xu >> 52) & 0x07ff - 54
      end
      m -= 1023
      y = reinterpret(Float64,(xu & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000)

      mf = Float64(m)
      F = (y + 3.5184372088832e13) - 3.5184372088832e13 # 0x1p-7*round(0x1p7*y)
      f = y-F
      jp = Base.Math.unsafe_trunc(Int,128.0*F)-127

      return Base.Math.log_proc1(y,mf,F,f,jp)
  elseif x == 0.0
      -Inf
  else
      NaN
  end
end
