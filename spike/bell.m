function y = bell(x)
  x_0 = 0.5;
  width = 0.25;
  y = std_bell((x - x_0) / width);
end

function y = std_bell(x)
  epsilon = 1e-5;
  y = zeros(size(x));
  i = find(abs(x) < 1 - epsilon);
  y(i) = exp(1 - 1 ./ (1 - x(i).*x(i)));
end
