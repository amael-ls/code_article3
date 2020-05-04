function [ out ] = d( t, s )
	% d is negative here, to get a real solution for s*. Although the density
	% eta will diverge in time (since a negative mortality becomes a source!),
	% eta and the other integrals are still bounded in size. It implies that I
	% will integrate only on a bounded time domain.
	% Not that because I know d is constant, I do not add the parameter s*
	out = -0.4; % cf Demography.h++

end
