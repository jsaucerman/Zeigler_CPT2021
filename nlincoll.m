function dy = nlincoll(t,y,kgen,kdeg,days,fibro,kgt,kdt)

% if y(1) >= 100
%     dy = 0;
% else

cf = interp1(days,fibro,t);
kg = kgen*interp1(days,kgt,t);
kd = kdeg*interp1(days,kdt,t);
dy = kg*cf - kd*y(1);
% end

end