
function [Di] = calgradient(net, x)

tol=0.001;
Di= (net(x'+tol)-net(x'-tol))/tol/2;
Di=Di';
% figure(7)
% plot(x,Di);
end