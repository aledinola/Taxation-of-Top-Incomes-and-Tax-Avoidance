function y = fun_entre_output(theta,nval,kval,gamma,vi)
% Output of entrepreneurial business given ability theta, labor and capital
% inputs nval, kval.
% theta_val*(k**gamma*(le+n)**(1.0d0-gamma))**vi


y = theta*(kval.^gamma.*nval.^(1-gamma)).^vi;

end