
function f = objectivefcn(params,x,data,absol,phase)
iterator = numel(data.data_subj1.C.w)
f = 0
for i = 4:(iterator-1)
      f = f + abs(absol(i)*exp(1j*phase(i)*pi/180) - (params(1)*(params(2)*x(i)*1j+1)*exp(-params(3)*1j*x(i))*(params(4)^2)/((1j*x(i))^2+2*params(4)*params(5)*1j*x(i)+params(4)^2) ))^2

end

f = f/iterator


end