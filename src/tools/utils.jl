apply(f::Base.Callable, X) = [f(x) for x in X]

#nmod(k, n) = mod1((k - 1), n) + 1
nmod(k, n) = mod1(k, n)