

function Y = Greville(X)
	[m, n] = size(X);
	E = eye(n);
	a1 = X(1, :).';
	Y = a1 / (a1.' * a1);
	for i = 2:m
		a = X(i, :).';
		Z = E - Y * X(1:(i - 1), :);
		R = Y * Y.';
		if a.' * Z * a > 0
			Y = Y - (Z * a * a.' * Y) / (a.' * Z * a);
			Y(:, i) = (Z * a) / (a.' * Z * a);
		else
			Y = Y - (R * a * a.' * Y) / (1 + a.' * R * a);
			Y(:, i) = (R * a) / (1 + a.' * R * a);
		end
	end
end