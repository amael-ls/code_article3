function [ out ] = v( t, s, s_star )
	s_star = sqrt(1 - 2*log(exp(d(t, s)*t) + exp(-0.5*s^2/ - s))) - 1;
	if (any(imag(s_star)))
		s_star = 0;
	end
	out = (2 + s_star)*exp(-s_star) ./ (1 + s); % cf Demography.h++
end
