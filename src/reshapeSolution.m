function vars = reshapeSolution(sol, m, n, q, w, p, t, g)
% Helper function to reshape solution vector
% ============================================================
vars.X     = reshape(sol(1:m*n),                                                          m, n);
vars.Y     = reshape(sol(m*n+1:2*m*n),                                                    m, n);
vars.y     = -reshape(sol(2*m*n+1:2*m*n+m),                                                m, 1);
vars.z     = reshape(sol(2*m*n+m+1:2*m*n+m+n),                                            n, 1);
vars.U     = reshape(sol(2*m*n+m+n+1:2*m*n+m+n+m*q),                                      q, m);
vars.H     = reshape(sol(2*m*n+m+n+m*q+1:2*m*n+m+n+m*q+m*w),                              w, m);
vars.S     = reshape(sol(2*m*n+m+n+m*q+m*w+q+w+1:2*m*n+m+n+m*q+m*w+q+w+m*p),              p, m);
vars.beta  = reshape(sol(2*m*n+m+n+m*q+m*w+q+w+m*p+1:2*m*n+m+n+m*q+m*w+q+w+m*p+m),        1, m);
vars.gamma = reshape(sol(2*m*n+m+n+m*q+m*w+q+w+m*p+m+1:2*m*n+m+n+m*q+m*w+q+w+m*p+2*m),    1, m);
vars.F     = reshape(sol(2*m*n+m+n+m*q+m*w+q+w+m*p+2*m+1:2*m*n+m+n+m*q+m*w+q+w+m*p+2*m+m*t),         m, t);
vars.W     = reshape(sol(2*m*n+m+n+m*q+m*w+q+w+m*p+2*m+m*t+1:2*m*n+m+n+m*q+m*w+q+w+m*p+2*m+m*t+m*g), m, g);
end