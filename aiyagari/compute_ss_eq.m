function ss = compute_ss_eq(params)

r = 0.005;
w = 1;

[kpPolicy,v] = solve_decision_grid_search(w,r,params);

ss = [];
end